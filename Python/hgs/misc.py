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
from geodata.misc import DataError


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
def interpolateIrregular(old_time, flow_data, new_time, start_date=None, lkgs=True, lcheckComplete=True,  
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
    assert np.all( time_diff > 0 ), time_diff.min()
    time_diff = time_diff.reshape((len(time_diff),1)) # reshape to make sure broadcasting works
    # integrate flow over time steps before resampling
    flow_data[1:,:] -= np.diff(flow_data, axis=0)/2. # get average flow between time steps
    flow_data *= time_diff # integrate flow in time interval by multiplying average flow with time period
    flow_data = np.cumsum(flow_data, axis=0) # integrate by summing up total flow per time interval
    # interpolate integrated flow to new time axis
    #flow_data = np.interp(new_time, xp=old_time[:,0], fp=flow_data[:,0],).reshape((len(new_time),1))
    old_time = np.concatenate(([0],old_time), axis=0) # integrated flow at time zero must be zero...
    flow_data = np.concatenate(([[0,]*len(usecols)],flow_data), axis=0) # ... this is probably better than interpolation
    # N.B.: we are adding zeros here so we don't have to extrapolate to the left; on the right we just fill in NaN's
    if ( new_time[-1] - old_time[-1] ) > 3*86400. and lcheckComplete: 
        warn("Data record ends more than 3 days befor end of period: {} days".format((new_time[-1]-old_time[-1])/86400.))
    elif (new_time[-1]-old_time[-1]) > 5*86400.: 
        if lcheckComplete: 
          raise DataError("Data record ends more than 5 days befor end of period: {} days".format((new_time[-1]-old_time[-1])/86400.))
        else:
          warn("Data record ends more than 5 days befor end of period: {} days".format((new_time[-1]-old_time[-1])/86400.))
    flow_interp = interp1d(x=old_time, y=flow_data, kind=interp_kind, axis=0, copy=False, 
                           bounds_error=False, fill_value=fill_value, assume_sorted=True) 
    flow_data = flow_interp(new_time) # evaluate with call
    # compute monthly flow rate from interpolated integrated flow
    flow_data = np.diff(flow_data, axis=0) / np.diff(new_time, axis=0).reshape((len(new_time)-1,1))
    if lkgs: flow_data *= 1000 # convert from m^3/s to kg/s
    # return interpolated flow data
    return flow_data


if __name__ == '__main__':
    pass