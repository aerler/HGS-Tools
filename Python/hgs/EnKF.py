'''
Created on 2018-01-04

This module contains meta data and access functions for EnKF output files.

@author: Andre R. Erler, GPL v3
'''

# external imports
import datetime as dt
import numpy as np
import pandas as pd
import os
# GeoPy imports
from datasets.common import getRootFolder
from geodata.misc import ArgumentError, VariableError, DataError, isNumber, DatasetError,\
  AxisError
from geodata.base import Dataset, Variable, Axis
from datasets.WSC import getGageStation, GageStationError, loadWSC_StnTS, updateScalefactor
# EnKF imports
from hgs.misc import convertDate
from enkf_utils.enkf_input import readKister
from enkf_utils.enkf_output import loadHydro, loadObs 

## EnKF Meta-vardata

dataset_name = 'EnKF'
root_folder = getRootFolder(dataset_name='HGS', fallback_name='HGS') # get dataset root folder based on environment variables
prefix_file = 'batch.pfx' # text file that contians the HGS problem prefix (also HGS convention)

variable_attributes = dict()
# hydrograph variables
class Hydro(object):
    name = 'hydro'
    atts = dict(discharge_mms = dict(name='discharge', units='m^3/s', atts=dict(long_name='Surface Flow Rate')), 
                discharge_kgs = dict(name='discharge', units='kg/s', atts=dict(long_name='Surface Flow Rate')),)
variable_attributes.update(Hydro.atts)
# observation wells
class Obs(object):
    name = 'obs'
    atts = dict(head = dict(name='head', units='m', atts=dict(long_name='Pressure Head at Well')),)
variable_attributes.update(Obs.atts)
# axes
class Axes(object):    
    name = 'axes'
    atts = dict(time = dict(name='time', units='n/a', atts=dict(long_name='Index Time')),
                # variables that are not time-dependent
                node = dict(name='node', units='#', dtype=np.int64, atts=dict(long_name='Node Number')),
                ensemble = dict(name='ensemble', units='#', dtype=np.int64, atts=dict(long_name='Ensemble Member')),
                observation = dict(name='observation', units='#', dtype=np.int64, atts=dict(long_name='Observation Number')),
                )
variable_attributes.update(Axes.atts)

# helper functions

def timeAxis(start_date=None, end_date=None, sampling=None, date_range=None, time_axis=None, 
             llastIncl=True, ntime=None, varatts=None):
    ''' figure out type and dimensions of time axis '''
    # check time input
    if date_range: start_date,end_date,sampling = date_range
    if start_date and end_date and sampling:
        start_year,start_month,start_day = convertDate(start_date)
        start_datetime = np.datetime64(dt.datetime(year=start_year, month=start_month, day=start_day), sampling)
        end_year,end_month,end_day = convertDate(end_date)
        end_datetime = np.datetime64(dt.datetime(year=end_year, month=end_month, day=end_day), sampling)
        if llastIncl: end_datetime += np.timedelta64(1, sampling)
        date_range = np.arange(start_datetime, end_datetime, dtype='datetime64[{}]'.format(sampling))
        assert date_range[0] == start_datetime, date_range[0]
        if ntime:
            if ntime > len(date_range): 
                raise ArgumentError(date_range)
            else:
                # trim 
                date_range = date_range[0:ntime]
        else: ntime = len(date_range)
    elif time_axis == 'datetime':
        raise ArgumentError('Insufficient time axis information!')
    # construct time axis
    atts = varatts['time']
    if time_axis.lower() == 'simple':
        time = Axis(atts=atts, coord=np.arange(1, ntime+1))
    elif time_axis.lower() == 'datetime':
        if   sampling.lower() == 'y' or sampling.lower() == '1y': units = 'year'
        elif sampling.lower() == 'm' or sampling.lower() == '1m': units = 'month'
        elif sampling.lower() == 'd' or sampling.lower() == '1d': units = 'day'
        elif sampling.lower() == 'h' or sampling.lower() == '1h': units = 'hour'
        else: units = sampling
        long_name = '{}s since {}'.format(units.title(),str(date_range[0])) # hope this makes sense...
        atts.update(long_name=long_name, units=units)
        time = Axis(atts=atts, coord=date_range)
    else:
        raise ArgumentError(time_axis)
    # return time axis
    return time


## function to load HGS station timeseries
def loadEnKF_StnTS(folder=None, varlist='default', varatts=None, name='enkf', title='EnKF', basin=None,  
                   start_date=None, end_date=None, sampling=None, period=None, date_range=None,  
                   llastIncl=True, WSC_station=None, basin_list=None, filenames=None, prefix=None, 
                   time_axis='datetime', scalefactors=None, metadata=None, lkgs=False, out_dir='out/', 
                   nreal=None, ntime=None, **kwargs):
    ''' load EnKF ensemble data as formatted GeoPy Dataset '''
    out_folder = os.path.join(folder,'out/') # default output folder
    if not os.path.exists(out_folder): raise IOError(out_folder)
    # default values
    if isinstance(varlist,basestring) and varlist == 'default': 
        varlist = ['discharge']
    if varatts is None: varatts = variable_attributes.copy()
    # figure out time axis
    time = timeAxis(start_date=start_date, end_date=end_date, sampling=sampling, date_range=date_range, 
                    time_axis=time_axis, llastIncl=llastIncl, ntime=ntime, varatts=varatts)
    ntime = len(time)
    # load WSC station meta data
    pass 
    # initialize Dataset
    dataset = Dataset(name=name, title=title if title else name.title(), atts=metadata)
    ensemble = None
    # load observation/innovation data
    pass
    # load discharge/hydrograph data
    if 'discharge' in varlist:
        data = loadHydro(folder=out_folder, nreal=nreal, ntime=ntime)
        ntime,nreal = data.shape
        if ensemble is None:
            # construct ensemble axis
            ensemble = Axis(atts=varatts['ensemble'], coord=np.arange(1, nreal+1))
        elif len(ensemble) != nreal:
            raise AxisError(ensemble)
        if lkgs: data *= 1000.
        dataset += Variable(atts=varatts['discharge_kgs' if lkgs else 'discharge_mms'],
                            data=data, axes=(time,ensemble))  
    
    # return formatted Dataset 
    if scalefactors is not None and scalefactors != 1: raise NotImplementedError
    return dataset

## function to load HGS station timeseries
def loadKister_StnTS(station=None, well=None, folder=None, varlist='default', varatts=None, 
                     name='observations', title=None, 
                     basin=None, start_date=None, end_date=None, sampling=None, period=None, date_range=None,  
                     llastIncl=True, WSC_station=None, basin_list=None, filenames=None, 
                     time_axis='datetime', scalefactors=None, metadata=None, lkgs=False, ntime=None, **kwargs):
    ''' load EnKF ensemble data as formatted GeoPy Dataset '''
    if folder and not os.path.exists(folder): raise IOError(folder)
    # default values
    if isinstance(varlist,basestring) and varlist == 'default':
        varlist = [] 
        if station: varlist += ['discharge']
        if well: varlist += ['head']
    if varatts is None: varatts = variable_attributes.copy()
    # figure out time axis
    if date_range: start_date,end_date,sampling = date_range
    time = timeAxis(start_date=start_date, end_date=end_date, sampling=sampling, date_range=date_range, 
                    time_axis=time_axis, llastIncl=llastIncl, ntime=ntime, varatts=varatts)
    ntime = len(time)
    # load WSC station meta data
    pass 
    # initialize Dataset
    dataset = Dataset(name=name, title=title if title else name.title(), atts=metadata)
    # load well data
    pass
    # load discharge/hydrograph data
    if 'discharge' in varlist:
        if not station: raise ArgumentError
        if folder: filepath = os.path.join(folder,station) # default output folder
        else: filepath = station
        data = readKister(filepath=filepath, period=(start_date,end_date), resample=sampling, lvalues=True)
        ntime = len(data)
        if lkgs: data *= 1000.
        dataset += Variable(atts=varatts['discharge_kgs' if lkgs else 'discharge_mms'],
                            data=data, axes=(time,))            
    # return formatted Dataset 
    if scalefactors is not None and scalefactors != 1: raise NotImplementedError
    return dataset
  

if __name__ == '__main__':
    
    ## settings
#     folder = 'D:/Data/HGS/SNW/EnKF/TWC/enkf_test_open_dec/' # experiment folder
#     folder = 'D:/Data/HGS/SNW/EnKF/TWC/enkf_test_closed_dec/' # experiment folder
#     date_range = ('2017-12-01', '2017-12-31', 'D'); ntime = None

    folder = 'D:/Data/HGS/SNW/EnKF/TWC/enkf_test_closed_may/' # experiment folder
    date_range = ('2017-05-01', '2017-12-31', 'D'); ntime = 166
    
    # check sample folder
    if not os.path.exists(folder): raise IOError(folder)

    # select test mode
    test_mode = 'test_dataset'
#     test_mode = 'test_kister'

    if test_mode == 'test_dataset':
      
      # load single dataset
      ds = loadEnKF_StnTS(folder=folder, date_range=date_range, ntime=ntime)
      print(ds)
  
    elif test_mode == 'test_kister':
      
      # test file
      csv_file = 'D:/Data/HGS/SNW/EnKF/Kister/W268-1.csv'
      
      # load single dataset
      ds = loadKister_StnTS(station=csv_file, date_range=date_range)
      print(ds)
