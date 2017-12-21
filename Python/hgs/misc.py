'''
Created on Jul 17, 2017

A module contianing various utility functions for HGS data

@author: Andre R. Erler, GPL v3
'''

# external imports
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from warnings import warn
# internal imports
from geodata.misc import DataError, ArgumentError
from datasets.WSC import GageStationError


# a function to resolve a data
def convertDate(date):
  ''' convert a data into a common format '''
  if isinstance(date, (tuple,list)):
    year = date[0] # need these for time axis (day not really...)
    month = date[1] if len(date) > 1 else 1
    day = date[2] if len(date) > 2 else 1
  elif isinstance(date, dict):
    year = date['year']; month = date.get('month',1); day = date.get('day',1)
  else:
    year=int(date); month=1; day=1
  return year,month,day

# a function to conver to month since a certain year
def monthSince(year, month, start_year=1979, start_month=1):
  ''' compute the number of month since the start date '''
  return 12*(year - start_year) + month - start_month
  

# a date parser function for the pandas ascii reader
def date_parser(seconds, year=1979, month=1, day=1, **kwargs):
    ''' convert seconds into a datetime object '''
    return pd.Index(pd.datetime(year=year, month=month, day=day, **kwargs) + seconds * pd.offsets.Second())


# interpolation function for HGS hydrographs etc.
def interpolateIrregular(old_time, data, new_time, start_date=None, lkgs=True, lcheckComplete=True,  
                         usecols=None, interp_kind='linear', fill_value=np.NaN):
    ''' a mass-conservative function to interpolate irregular HGS timeseries output to a regular time axis 
        the function works by first integrating the flux, then interpolating the cumulative flux, and subsequently
        differentiating to obtain the average instantaneous flux; this method is mass-conservative, assuming 
        piece-wise linearly varying flux '''
    # convert monthly time series to regular array of seconds since start date
    if start_date is None: start_date = new_time[0] 
    new_time = ( new_time.astype('datetime64[s]') - start_date.astype('datetime64[s]') ) / np.timedelta64(1,'s') 
    # N.B.: cast as seconds, but without time units (hence division); center at zero
    if new_time[0] != 0: warn('New time axis does not start at zero ({}).'.format(new_time[0]))
    # original time deltas in seconds
    time_diff = old_time.copy(); time_diff[1:] = np.diff(old_time) # time period between time steps
    if not np.all( time_diff > 0 ):
      imin = time_diff.argmin()
      raise ValueError(time_diff.min(),imin,old_time[imin-1:imin+1])
    # reshape to make sure broadcasting works
    time_diff = time_diff.reshape((len(time_diff),1))
    if data.ndim > 2:
        oldshp = data.shape[1:]; ncols = np.prod(oldshp) # remember for recovery
        data = data.reshape((len(old_time),ncols)) # now simple 2D array
    else: 
      oldshp = None; ncols = data.shape[1]
    # integrate flow over time steps before resampling
    data[1:,:] -= np.diff(data, axis=0)/2. # get average flow between time steps
    data *= time_diff # integrate flow in time interval by multiplying average flow with time period
    data = np.cumsum(data, axis=0) # integrate by summing up total flow per time interval
    # interpolate integrated flow to new time axis
    #data = np.interp(new_time, xp=old_time[:,0], fp=data[:,0],).reshape((len(new_time),1))
    old_time = np.concatenate(([0],old_time), axis=0) # integrated flow at time zero must be zero...
    data = np.concatenate(([[0,]*ncols],data), axis=0) # ... this is probably better than interpolation
    # N.B.: we are adding zeros here so we don't have to extrapolate to the left; on the right we just fill in NaN's
    if ( new_time[-1] - old_time[-1] ) > 3*86400. and lcheckComplete: 
        warn("Data record ends more than 3 days befor end of period: {} days".format((new_time[-1]-old_time[-1])/86400.))
    elif (new_time[-1]-old_time[-1]) > 5*86400.: 
        if lcheckComplete: 
          raise DataError("Data record ends more than 5 days befor end of period: {} days".format((new_time[-1]-old_time[-1])/86400.))
        else:
          warn("Data record ends more than 5 days befor end of period: {} days".format((new_time[-1]-old_time[-1])/86400.))
    flow_interp = interp1d(x=old_time, y=data, kind=interp_kind, axis=0, copy=False, 
                           bounds_error=False, fill_value=fill_value, assume_sorted=True) 
    data = flow_interp(new_time) # evaluate with call
    # compute monthly flow rate from interpolated integrated flow
    data = np.diff(data, axis=0) / np.diff(new_time, axis=0).reshape((len(new_time)-1,1))
    if lkgs: data *= 1000 # convert from m^3/s to kg/s
    # return interpolated flow data
    if oldshp:
        data = data.reshape((len(new_time)-1,)+oldshp) # return to old shape
    return data
  
  
    # function to parse observation well output
def parseObsWells(filepath, variables=None, constants=None, layers=None, lskipNaN=True):
    ''' a function to parse observation well output and return a three dimansional array, similar to
        genfromtxt, except with an additional depth/layer dimension at the end (time,variable,layer) '''
    # load file contents (all at once)
    with open(filepath,'r') as f:
        lines = f.readlines()
    ln = len(lines)
    # validate header
    line = lines[0].lower()
    if "title" not in line: raise GageStationError((line,filepath))
    line = lines[1].lower()
    if "variables" not in line: raise GageStationError((line,filepath))
    varlist = [v for v in line[line.find('=')+1:].strip().split(',') if len(v) > 0]
    lv = len(varlist)
    # validate constants (variables with no time-dependence))
    if constants is None: 
        ce = 0 # assume no constants
        constants = []
    elif isinstance(constants, (list,tuple)):
        ce = len(constants)
    elif isinstance(constants, (int,np.integer)):
        constants = (constants,)
        ce = 1
    else: TypeError(constants)
    if lv < len(constants): raise ArgumentError((constants,varlist))
    # validate variables (variables with time-dependence)
    if variables is None: 
        ve = lv # process all columns
        variables = range(lv)
    elif isinstance(variables, (list,tuple)):
        ve = len(variables)
    elif isinstance(variables, (int,np.integer)):
        variables = (variables,)
        ve = 1
    else: TypeError(variables)
    if lv < len(variables)+ce: raise ArgumentError((variables,varlist))
    line = lines[2].lower()
    if "zone" != line[:4]: raise GageStationError((line[:4],filepath))
    if "solutiontime" != line[22:34] : raise GageStationError((line,filepath))
    # determine number of layers
    i = 3
    while 'zone' not in lines[i]: i += 1
    line = lines[i].lower()
    if "zone" != line[:4]: raise GageStationError(line,filepath)
    if "solutiontime" != line[22:34] : raise GageStationError((line,filepath))
    nlay = i - 3; te = (ln-2)/(nlay+1)
    assert (nlay+1)*te == ln-2, (ln,nlay,te)
    if layers is None: 
        #layers = tuple(range(1,nlay+1))
        le = nlay # process all layers/rows
    elif isinstance(layers, (list,tuple)):
        if nlay < len(layers): raise ArgumentError((layers,nlay))
        le = len(layers)
    elif isinstance(layers, (int,np.integer)):
        layers = (layers,)
        le = 1
    else: raise TypeError(layers)
    assert le <= nlay, (le,nlay)
    # allocate array
    time = np.zeros((te,)) # time coordinate
    data = np.zeros((te,ve,le)) # time-dependent variables
    const = np.zeros((ce,le)) # variables that are not time-dependent
    # walk through lines and pick values
    line_iter = lines.__iter__() # return an iterator over all lines
    line_iter.next(); line_iter.next() # discard first two lines
    # initialize first time-step 
    n = -1; k = nlay; te1 = te-1
    if layers: 
        layers.sort() # make sure the list is sorted
        layers.append(-1) # this is a dummy to match nothing at the end
        llayers = True
    else: llayers = False
    while n < te1 or k < nlay:
      line = line_iter.next() # should work flawlessly, if I got the indices right...
      if k == nlay:
          assert line[:4] == 'zone', line[:4]
          assert "SOLUTIONTIME" == line[22:34], line[22:34]
          # prepare next time step
          k = 0 # initialize a new time-step
          if llayers: 
              l = 0; klay = layers[0]
          n += 1
          time[n] = np.float64(line[36:])          
      else:
          if not llayers or k == klay:
              elements = line.split() # split fields, but only process those that we want
              values = [np.float64(elements[i]) for i in variables]
              # insert values into array
              if not llayers: l = k
              data[n,:,l] = values
              # process constants during first time-step
              if n == 0:
                  values = [np.float64(elements[i]) for i in constants]
                  # insert values into array
                  const[:,l] = values
              # advance layer counter
              if llayers:
                  l += 1; klay = layers[l]
                  
          # else just skip 
          k += 1
    # make sure we read all lines/time steps
    try: 
      print(line_iter.next()) # this should raise a StopIteration exception
      raise ValueError("It appears the file was not read completely...")
    except StopIteration: pass # exception is right here
    # return array 
    return time,data,const
  

if __name__ == '__main__':
    pass