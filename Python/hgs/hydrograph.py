#!/usr/local/bin/python2.7
# encoding: utf-8
'''
A program to read HGS hydrograph files, resample to daily time-steps and count occurences of low/high flow days.

The program is a command line utility that can be executed on any HGS hydrograph file and will write output to 
CSV files.

@author:     Andre R. Erler

@copyright:  2018 Aquanty Inc. All rights reserved.

@license:    GPL v3

@contact:    aerler@aquanty.com
@deffield    updated: 09/11/2018
'''

import sys, os
import numpy as np

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

# internal imports
from misc import interpolateIrregular, ArgumentError

__all__ = []
__version__ = 0.2
__date__ = '2018-11-09'
__updated__ = '2018-11-10'

DEBUG = 0
TEST  = 0

class CLIError(Exception):
    '''Generic exception to raise and log different fatal errors.'''
    def __init__(self, msg):
        super(CLIError).__init__(type(self))
        self.msg = "E: {:s}".format(msg)
    def __str__(self):
        return self.msg
    def __unicode__(self):
        return self.msg

def main(argv=None): 
    '''Command line options and program execution: read hydrograph, resample to daily, and process flow levels.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    # Setup argument parser
    parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
    # hydrograph file (only positional argument, since it is required)
    parser.add_argument('file', metavar='file', type=str, help="the hydrograph file to be used")
    # optional arguments
    parser.add_argument("--variable", dest="variable", default=None, type=str,
                        help="variable name; overrides '--var-col' and determines column based on name in file header")
    parser.add_argument("--var-col", dest="varcol", default=1, type=int, # should be streamflow
                        help="variable column in hydrograph file (zero is time; streamflow should be 1) [default: %(default)s]")
    parser.add_argument("--daily-output", dest="daily_output", default=False, type=str, 
                        help="save resampled hydrograph timeseries to this file [default: %(default)s]")
    parser.add_argument("--low-flow", dest="lowflow", default=None, type=float, 
                        help="count occurences and duration of flow below this threshold")                
    parser.add_argument("--high-flow", dest="hiflow", default=None, type=float, 
                        help="count occurences and duration of flow above this threshold")                
    parser.add_argument("--min_duration", dest="mindays", default=None, type=int, 
                        help="minimum duration of high/low flow conditions to be recorded")                
    parser.add_argument("--flow-output", dest="flow_output", default='flow_duration.dat', type=str, 
                        help="minimum duration of high/low flow conditions to be recorded [default: %(default)s]")                
    # misc options
    parser.add_argument('-V', '--version', action='version', version=program_version_message)
    parser.add_argument("--debug", dest="debug", action="store_true", help="print debug output [default: %(default)s]")
    parser.add_argument("--verbose",'-v', dest="verbose", action="store_true", help="print status output [default: %(default)s]")

    # Process arguments
    args = parser.parse_args()
    filepath     = args.file
    variable     = args.variable
    varcol       = args.varcol
    daily_output = args.daily_output
    lowflow      = args.lowflow
    hiflow       = args.hiflow
    mindays      = args.mindays
    flow_output  = args.flow_output
    ldebug       = args.debug
    lverbose     = args.verbose
    
    if not filepath:
        raise CLIError("No input hydrograph file specified!")
    elif not os.path.exists(filepath):
        raise IOError("Input file '{:s}' not found!".format(filepath))

    # determine column (based on variable name in header)
    # parse header
    with open(filepath, 'r') as f:
        f.readline(); line = f.readline(); lline = line.lower() # skip to 2nd line and lower case
        # parse variables and determine columns
        if not "variables" in lline: raise IOError(line)
        variable_list = [v for v in line[line.find('=')+1:].strip().split(',') if len(v) > 0]
        # clean up a little and remove some problematic characters
        variable_list = [v.strip().strip('"').strip().lower() for v in variable_list]
        if variable:
            # assign variable column
            try:
                varcol = variable_list.index(variable)
            except ValueError:
                raise ArgumentError("Variable '{}' not found in file '{}'.".format(variable,filepath))
        else:
            if varcol >= len(variable_list):
                raise ArgumentError("Invalid colum index '{:d}': only {:d} columns present".format(varcol,len(variable_list)))
            variable = variable_list[varcol]
        if lverbose: print("\nFound variable '{:s}' in column {:d}".format(variable,varcol))

    # read hydrograph timeseries
    if lverbose: print("\nLoading hydrograph data from file:\n '{:s}'".format(filepath))
    data = np.genfromtxt(filepath, dtype=np.float64, delimiter=None, skip_header=3, usecols = (0,varcol))
    assert data.shape[1] == 2, data.shape
    time = data[:,0]; flow = data[:,1:]
    assert flow.shape == (len(time),1), data.shape
    
    # construct regular daily timeseries    
    time_start = time[0]; time_end = time[-1]
    te = np.ceil( ( time_end - time_start ) / 86400. )
    time_resampled = np.arange(0,te+1)*86400
    if ldebug:
        print("\nResampling timeseries to daily intervals (in seconds):\n "+repr(time_resampled))
        # N.B.: the array elements represent the boundaries of averaging periods, 
        #       i.e. start and end of each day
    # call function to interpolate irregular HGS timeseries to regular daily timseries  
    flow = interpolateIrregular(old_time=time, lkgs=False, data=flow, new_time=time_resampled, 
                                start_date=None, interp_kind='linear', 
                                lcheckComplete=False, usecols=1, fill_value=np.NaN).squeeze()
    assert len(flow) == len(time_resampled)-1, (flow.shape,len(time_resampled-1))
    
    # output hydrograph timeseries resampled to daily output
    if daily_output:
        header = "file='{:s}'\nvariable='{:s}',column={:d}\nresampled to daily averages".format(filepath,variable,varcol)
        if lverbose: print("\nSaving resampled timeseries to:\n '{:s}'".format(daily_output))
        np.savetxt(daily_output, flow, delimiter=',', header=header, comments='#')

    ## count and output low-flow occurences and durations
    if lowflow:
        # loop over timeseries and count low-flow occurences
        occurence = []; duration = []; c = 0 # hi/low flow counter
        for i,f in enumerate(flow):
            if f < lowflow: 
                c += 1 # simply increment - we are in a low flow period
            elif c > 0: 
                # a low flow period is terminating
                if c >= mindays:
                    # if of sufficient length, store duration and start day
                    duration.append(c)
                    occurence.append(i-c+1) # index i is zero-based, but we count one-based
                # reset counter
                c = 0
        if ldebug:
            print("Low flow occurences:\n "+repr(occurence))
            print("Low flow durations:\n "+repr(duration))
        # header for low flow
        header = "Occurence and duration of flow below {}".format(lowflow)
    # N.B.: the code below is duplicated, to avoid unnecessary branching in the loop
    if hiflow:
        # loop over timeseries and count hi-flow occurences
        occurence = []; duration = []; c = 0 # hi/low flow counter
        for i,f in enumerate(flow):
            if f > hiflow: 
                c += 1 # simply increment - we are in a low flow period
            elif c > 0: 
                # a high flow period is terminating
                if c >= mindays:
                    # if of sufficient length, store duration and start day
                    duration.append(c)
                    occurence.append(i-c+1) # index i is zero-based, but we count one-based
                # reset counter
                c = 0
        if ldebug:
            print("High flow occurences:\n "+repr(occurence))
            print("High flow durations:\n "+repr(duration))
        # header for low flow
        header = "Occurence and duration of high above {}".format(hiflow)
        
    if hiflow or lowflow:
        # save hi/low flow data to CSV file
        header += "\noccurence (day): line 1; duration (days): line 2"
        header += "\nsource='{:s}',variable='{:s}',column={:d},resampled to daily averages".format(filepath,variable,varcol)
        if lverbose: print("\nSaving hi/low flow occurences and durations to:\n '{:s}'".format(flow_output))
        np.savetxt(flow_output, np.vstack((occurence,duration)), fmt='%d',
                   delimiter=',', header=header, comments='#')
                

    return 0


if __name__ == "__main__":
  
    # some meta data
    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    program_license = '''{:s}

  Created by Andre R. Erler on {:s}.
  Copyright 2018 Aquanty Inc. All rights reserved.

  Licensed under the GNU General Public License 3.0
  https://www.gnu.org/licenses/gpl.html

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE

python -u ~/path/to/HGS\ Tools/Python/hgs/hydrograph.py --verbose --resample_out daily_hydrograph.dat \\
                                                        --low-flow=45 --min_duration=30 --flow_output=low_flow_duration.dat \\
                                                        prefixo.hydrograph.some_station.dat

N.B.: the program depends on the misc.py module in the same package, so make sure the package is in your PYTHONPATH
      or the program is executed from the package directory (as shown above)

OPTIONS
'''.format(program_shortdesc, str(__date__))

    # debug and testing settings
    if DEBUG:
        sys.argv.append("--debug")
        sys.argv.append("-v")
    if TEST:
        sys.argv.append('D:/Data/HGS/GRW/grw2/NRCan/timeseries/hgs_run_v3/grw_omafrao.hydrograph.2GB001.dat')
        # test resampling
        resample_file = '{:s}/daily_hydrograph.dat'.format(os.getcwd())
        sys.argv.append('--resample_output='+resample_file)
        # test low flow output
        duration_file = '{:s}/durations.dat'.format(os.getcwd())
        sys.argv.append('--flow_output='+duration_file)
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
        except Exception, e:
            indent = len(program_name) * " "
            sys.stderr.write(program_name + ": " + repr(e) + "\n")
            sys.stderr.write(indent + "  for help use --help\n")
            sys.exit(2)
            
    if TEST:
        
        print("\nResampled flow data:")
        data = np.loadtxt(resample_file, delimiter=',', comments='#')
        print(data)
        
        print("\nOccurence and durations:")
        occurence,duration = np.loadtxt(duration_file, dtype=np.int64, delimiter=',', comments='#')
        # N.B.: you can also use one target variable, which will be a 2D array: data = np.loadtxt(...)
        print(occurence)
        print(duration)
        