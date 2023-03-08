#!/usr/local/bin/python2.7
# encoding: utf-8
'''
A stand-alone program to read HGS hydrograph files, resample to daily time-steps and count occurences of low/high flow days.

The program is a command line utility that can be executed on any HGS hydrograph file and will write output to
CSV files.

@author:     Andre R. Erler

@copyright:  2018 Aquanty Inc. All rights reserved.

@license:    GPL v3

@contact:    aerler@aquanty.com
@deffield    updated: 06/03/2019
'''

import sys
import os
import numpy as np
from scipy.interpolate import interp1d
from warnings import warn

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

__all__ = []
__version__ = '0.3'
__date__ = '2018-11-09'
__updated__ = '2023-01-11'


class CLIError(Exception):
    '''Generic exception to raise and log different fatal errors.'''
    def __init__(self, msg):
        super(CLIError).__init__(type(self))
        self.msg = "E: {:s}".format(msg)

    def __str__(self):
        return self.msg

    def __unicode__(self):
        return self.msg


# interpolation function for HGS hydrographs etc.
def interpolateIrregular(old_time, data, new_time, start_date=None, lkgs=True, lcheckComplete=True,
                         interp_kind='linear', fill_value=np.NaN):
    ''' a mass-conservative function to interpolate irregular HGS timeseries output to a regular time axis
        the function works by first integrating the flux, then interpolating the cumulative flux, and subsequently
        differentiating to obtain the average instantaneous flux; this method is mass-conservative, assuming
        piece-wise linearly varying flux. (N.B.: the 'new_time' array elements bracket the averaging time periods
        and the array has to be one element longer than the desired number of interpolated time steps.)'''
    # convert monthly time series to regular array of seconds since start date
    if start_date is None:
        start_date = new_time[0]
    new_time = (new_time.astype('datetime64[s]') - start_date.astype('datetime64[s]')) / np.timedelta64(1, 's')
    # N.B.: cast as seconds, but without time units (hence division); origin at zero
    if new_time[0] != 0:
        warn('New time axis does not start at zero ({}).'.format(new_time[0]))
    # also shift old time origin to zero
    if old_time[0] != 0:
        old_time -= old_time[0]
    # do some checks
    end_time_diff = (new_time[-1] - old_time[-1]) / 86400.
    if end_time_diff > 3 and lcheckComplete:
        warn("Data record ends more than 3 days befor end of period: {} days".format(end_time_diff))
    elif end_time_diff > 5:
        if lcheckComplete:
            raise ValueError("Data record ends more than 5 days befor end of period: {} days".format(end_time_diff))
        else:
            warn("Data record ends more than 5 days befor end of period: {} days".format(end_time_diff))
    # original time deltas in seconds
    time_diff = np.zeros_like(old_time)
    time_diff[1:] = np.diff(old_time)  # time period between time steps
    # find and remove leading zero timesteps (can happen, if internal timesteps are smaller than numerical precision)
    i_trim = 0
    while time_diff[i_trim] == 0:
        i_trim += 1
    i_trim -= 1  # keep last zero
    if i_trim > 0:
        time_diff = time_diff[i_trim:]
        old_time = old_time[i_trim:]
    if not np.all(time_diff[1:] > 0):  # except first, which is always zero
        imin = time_diff[1:].argmin() + 1
        raise ValueError(time_diff.min(), imin, old_time[imin - 1:imin + 1])
    # reshape to make sure broadcasting works
    time_diff = time_diff.reshape((len(time_diff), 1))
    if data.ndim > 2:
        oldshp = data.shape[1:]
        ncols = np.prod(oldshp)  # remember for recovery
        data = data.reshape((len(old_time), ncols))  # now simple 2D array
    else:
        oldshp = None
        ncols = data.shape[1]
    if i_trim > 0:
        data = data[i_trim:, :]
    # integrate flow over time steps before resampling
    data[1:, :] -= np.diff(data, axis=0) / 2.  # get average flow between time steps
    data *= time_diff  # integrate flow in time interval by multiplying average flow with time period
    data = np.cumsum(data, axis=0)  # integrate by summing up total flow per time interval
    # interpolate integrated flow to new time axis
    # old_time = np.concatenate(([-1], old_time), axis=0)  # integrated flow at time zero must be zero...
    # data = np.concatenate(([[0] * ncols], data), axis=0)  # ... this is probably better than interpolation
    # N.B.: we are adding zeros here so we don't have to extrapolate to the left; on the right we just fill in NaN's
    flow_interp = interp1d(x=old_time, y=data, kind=interp_kind, axis=0, copy=False,
                           bounds_error=False, fill_value=fill_value, assume_sorted=True)
    data = flow_interp(new_time)  # evaluate with call
    # compute monthly flow rate from interpolated integrated flow
    data = np.diff(data, axis=0) / np.diff(new_time, axis=0).reshape((len(new_time) - 1, 1))
    if lkgs:
        data *= 1000  # convert from m^3/s to kg/s
    # return interpolated flow data
    if oldshp:
        data = data.reshape((len(new_time) - 1,) + oldshp)  # return to old shape
    return data


def main(argv=None):
    '''Command line options and program execution: read hydrograph, resample to daily, and process flow levels.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    # Setup argument parser
    parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
    # hydrograph file (only positional argument, since it is required)
    parser.add_argument('file', metavar='file', type=str, help="the HGS hydrograph/timeseries file to be used")
    # optional arguments
    parser.add_argument("--variable", dest="variable", default=None, type=str,
                        help="variable name; overrides '--var-col' and determines column based on name in file header")
    parser.add_argument("--var-col", dest="varcol", default=1, type=int,  # should be streamflow
                        help="variable column in HGS timeseries output (0 is time; streamflow should be 1) [default: %(default)s]")
    parser.add_argument("--daily-output", dest="daily_output", default=None, type=str,
                        help="save daily resampled timeseries; the argument is the destination file name [default: %(default)s]")
    parser.add_argument("--low-flow", dest="lowflow", default=None, type=float,
                        help="count occurences and duration (number of days) below this threshold")
    parser.add_argument("--high-flow", dest="hiflow", default=None, type=float,
                        help="count occurences and duration (number of days) above this threshold")
    parser.add_argument("--min-duration", dest="min_days", default=None, type=int,
                        help="minimum duration (in days) of high/low flow conditions that will be recorded")
    parser.add_argument("--duration-out", dest="duration_out", default='duration_output.dat', type=str,
                        help="minimum duration of high/low flow conditions to be recorded [default: %(default)s]")
    # misc options
    parser.add_argument('-V', '--version', action='version', version=program_version_message)
    parser.add_argument('-d', "--debug", dest="debug", action="store_true", help="print debug output [default: %(default)s]")
    parser.add_argument('-v', "--verbose", dest="verbose", action="store_true", help="print status output [default: %(default)s]")

    # Process arguments
    args = parser.parse_args()
    filepath     = args.file
    variable     = args.variable
    varcol       = args.varcol
    daily_output = args.daily_output
    lowflow      = args.lowflow
    hiflow       = args.hiflow
    mindays      = args.min_days
    duration_out = args.duration_out
    ldebug       = args.debug
    lverbose     = args.verbose

    if not filepath:
        raise CLIError("No input hydrograph file specified!")
    elif not os.path.exists(filepath):
        raise IOError("Input file '{:s}' not found!".format(filepath))

    # determine column (based on variable name in header)
    # parse header
    with open(filepath, 'r') as f:
        f.readline(); line = f.readline(); lline = line.lower()  # skip to 2nd line and lower case
        # parse variables and determine columns
        if "variables" not in lline:
            raise IOError(line)
        variable_list = [v for v in line[line.find('=') + 1:].strip().split(',') if len(v) > 0]
        # clean up a little and remove some problematic characters
        variable_list = [v.strip().strip('"').strip().lower() for v in variable_list]
        if variable:
            # assign variable column
            try:
                varcol = variable_list.index(variable)
            except ValueError:
                raise CLIError("Variable '{}' not found in file '{}'.".format(variable, filepath))
        else:
            if varcol >= len(variable_list):
                raise CLIError("Invalid colum index '{:d}': only {:d} columns present".format(varcol, len(variable_list)))
            variable = variable_list[varcol]
        if lverbose:
            print(("\nFound variable '{:s}' in column {:d}".format(variable, varcol)))

    # read hydrograph timeseries
    if lverbose:
        print(("\nLoading hydrograph/timeseries data from file:\n '{:s}'".format(filepath)))
    data = np.genfromtxt(filepath, dtype=np.float64, delimiter=None, skip_header=3, usecols=(0, varcol))
    assert data.shape[1] == 2, data.shape
    time = data[:, 0]
    flow = data[:, 1:]
    assert flow.shape == (len(time), 1), data.shape

    # construct regular daily timeseries
    time_start = time[0]; time_end = time[-1]
    te = np.ceil((time_end - time_start) / 86400.)
    time_resampled = np.arange(0, te + 1) * 86400
    if ldebug:
        print(("\nResampling timeseries to daily intervals (in seconds):\n " + repr(time_resampled)))
        # N.B.: the array elements represent the boundaries of averaging periods,
        #       i.e. start and end of each day
    # call function to interpolate irregular HGS timeseries to regular daily timseries
    flow = interpolateIrregular(old_time=time, lkgs=False, data=flow, new_time=time_resampled,
                                start_date=None, interp_kind='linear',
                                lcheckComplete=False, fill_value=np.NaN).squeeze()
    assert len(flow) == len(time_resampled) - 1, (flow.shape, len(time_resampled - 1))

    # output hydrograph timeseries resampled to daily output
    if daily_output:
        header = "file='{:s}'\nvariable='{:s}',column={:d}\nresampled to daily averages".format(filepath, variable, varcol)
        if lverbose:
            print(("\nSaving resampled daily timeseries to:\n '{:s}'".format(daily_output)))
        np.savetxt(daily_output, np.stack((time_resampled[:-1], flow), axis=1), delimiter=',', header=header, comments='#')

    ## count and output low-flow occurences and durations
    if lowflow:
        # loop over timeseries and count low-flow occurences
        occurence = []; duration = []; c = 0  # hi/low flow counter
        for i, f in enumerate(flow):
            if f < lowflow:
                c += 1  # simply increment - we are in a low flow period
            elif c > 0:
                # a low flow period is terminating
                if c >= mindays:
                    # if of sufficient length, store duration and start day
                    duration.append(c)
                    occurence.append(i - c + 1)  # index i is zero-based, but we count one-based
                # reset counter
                c = 0
        if ldebug:
            print(("Low flow occurences:\n " + repr(occurence)))
            print(("Low flow durations:\n " + repr(duration)))
        # header for low flow
        header = "Occurence and duration of flow below {}".format(lowflow)
    # N.B.: the code below is duplicated, to avoid unnecessary branching in the loop
    if hiflow:
        # loop over timeseries and count hi-flow occurences
        occurence = []; duration = []
        c = 0  # hi/low flow counter
        for i, f in enumerate(flow):
            if f > hiflow:
                c += 1  # simply increment - we are in a low flow period
            elif c > 0:
                # a high flow period is terminating
                if c >= mindays:
                    # if of sufficient length, store duration and start day
                    duration.append(c)
                    occurence.append(i - c + 1)  # index i is zero-based, but we count one-based
                # reset counter
                c = 0
        if ldebug:
            print(("High flow occurences:\n " + repr(occurence)))
            print(("High flow durations:\n " + repr(duration)))
        # header for low flow
        header = "Occurence and duration of high above {}".format(hiflow)

    if hiflow or lowflow:
        # save hi/low flow data to CSV file
        header += "\noccurence (day): line 1; duration (days): line 2"
        header += "\nsource='{:s}',variable='{:s}',column={:d},resampled to daily averages".format(filepath, variable, varcol)
        if lverbose:
            print(("\nSaving hi/low flow occurences and durations to:\n '{:s}'".format(duration_out)))
        np.savetxt(duration_out, np.vstack((occurence, duration)), fmt='%d',
                   delimiter=',', header=header, comments='#')

    if lverbose:
        print('')  # output looks cleaner with an empty line at the end
    return 0


if __name__ == "__main__":

    # some meta data
    program_name = os.path.basename(sys.argv[0])
    program_version = "v{:s}".format(__version__)
    program_build_date = str(__updated__)
    program_version_message = '%(prog)s {:s} ({:s})'.format(program_version, program_build_date)
    program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    program_license = '''{:s}

  Created by Andre R. Erler on {:s}.
  Copyright 2018 Aquanty Inc. All rights reserved.

  Licensed under the GNU General Public License 3.0
  https://www.gnu.org/licenses/gpl.html

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.


USAGE

python -u Path/to/HGS-Tools/Python/hgs/hydrograph.py --verbose --daily_output daily_hydrograph.dat --low-flow=45 \\
                                                     --min-duration=30 --duration-out=low_flow_duration.dat \\
                                                     prefixo.hydrograph.some_station.dat

OPTIONS
'''.format(program_shortdesc, str(__date__))

    DEBUG = 1
    TEST  = 1

    # debug and testing settings
    if DEBUG:
        sys.argv.append("-d")
        sys.argv.append("-v")
    if TEST:
        sys.argv.append('1d_column_testo.observation_well_flow.100cm.dat')
        # test resampling
        resample_file = 'daily_hydrograph.csv'
        sys.argv.append('--daily-output=' + resample_file)
        # test low flow output
        duration_file = 'durations.dat'
        sys.argv.append('--duration-out=' + duration_file)
        sys.argv.append('--low-flow=45')
    # execute program
    if DEBUG:
        main()
    else:
        try:
            sys.exit(main())
        except KeyboardInterrupt:
            ### handle keyboard interrupt ###
            sys.exit(0)
        except Exception as e:
            indent = len(program_name) * " "
            sys.stderr.write(program_name + ": " + repr(e) + "\n")
            sys.stderr.write(indent + "  for help use --help\n")
            sys.exit(2)

    if TEST:

        print("\nResampled timeseries data:")
        data = np.loadtxt(resample_file, delimiter=',', comments='#')
        print(data)

        print("\nExceedance and durations:")
        occurence, duration = np.loadtxt(duration_file, dtype=np.int64, delimiter=',', comments='#')
        # N.B.: you can also use one target variable, which will be a 2D array: data = np.loadtxt(...)
        print(occurence)
        print(duration)
