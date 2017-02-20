'''
Created on Aug 27, 2016

A module to gauge station data from HGS and associate the data with station meta data from the Water 
Survey of Canada; the data is stored in human-readable text files and tables.

@author: Andre R. Erler, GPL v3
'''

# external imports
import numpy as np
import pandas as pd
import scipy.interpolate as si
import os
# internal imports
from datasets.common import data_root, BatchLoad
from geodata.base import Dataset, Variable, Axis
from geodata.misc import ArgumentError, VariableError, DataError
from datasets.WSC import getGageStation, GageStationError
import datetime as dt

## WSC (Water Survey Canada) Meta-data

dataset_name = 'HGS'
root_folder = '{:s}/{:s}/'.format(data_root,dataset_name) # the dataset root folder
station_file = '{PREFIX:s}o.hydrograph.Station_{STATION:s}.dat' # general HGS naming convention
prefix_file = 'batch.pfx' # text file that contians the HGS problem prefix (also HGS convention)

# variable attributes and name
variable_attributes = dict(sfroff = dict(name='sfroff', units='kg/m^2/s', atts=dict(long_name='Surface Runoff')),    # surface flow rate over area
                           ugroff = dict(name='ugroff', units='kg/m^2/s', atts=dict(long_name='Subsurface Runoff')), # subsurface flow rate over area
                           runoff = dict(name='runoff', units='kg/m^2/s', atts=dict(long_name='Total Runoff')),      # total flow rate over area
                           discharge = dict(name='discharge', units='kg/s', atts=dict(long_name='Surface Flow Rate')),      # surface flow rate
                           seepage   = dict(name='seepage'  , units='kg/s', atts=dict(long_name='Subsurface Flow Rate')),   # subsurface flow rate
                           flow      = dict(name='flow'     , units='kg/s', atts=dict(long_name='Total Flow Rate')),      ) # total flow rate
# list of variables to load
variable_list = variable_attributes.keys() # also includes coordinate fields    
flow_to_flux = dict(discharge='sfroff', seepage='ugroff', flow='runoff') # relationship between flux and flow variables
# N.B.: computing surface flux rates from gage flows also requires the drainage area
hgs_varlist = ['time','discharge','seepage','flow'] # internal use only; needs to have 'time'
hgs_variables = {'surface':'discharge','porous media':'seepage','total':'flow'} # mapping of HGS variable names GeoPy conventions

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

## function to load HGS station timeseries
def loadHGS_StnTS(station=None, varlist=None, varatts=None, folder=None, name=None, title=None,
                  start_date=1979, end_date=None, period=15, date_parser=date_parser, 
                  basin=None, WSC_station=None, basin_list=None, 
                  filename=station_file, prefix=None, resampling='1M', **kwargs):
  ''' Get a properly formatted WRF dataset with monthly time-series at station locations; as in
      the hgsrun module, the capitalized kwargs can be used to construct folders and/or names '''
  if folder is None or ( filename is None and station is None ): raise ArgumentError
  # prepare name expansion arguments (all capitalized)
  expargs = dict(ROOT_FOLDER=root_folder, STATION=station, NAME=name, TITLE=title,
                 PREFIX=prefix, BASIN=basin, WSC_STATION=WSC_station) 
  # exparg preset keys will get overwritten if capitalized versions are defined
  for key,value in kwargs.items():
    KEY = key.upper() # we only use capitalized keywords, and non-capitalized keywords are only used/converted
    if KEY == key or KEY not in kwargs: expargs[KEY] = value # if no capitalized version is defined
  # read folder and infer prefix, if necessary
  folder = folder.format(**expargs)
  if not os.path.exists(folder): raise IOError(folder)
  if expargs['PREFIX'] is None:
    with open('{}/{}'.format(folder,prefix_file), 'r') as pfx:
      expargs['PREFIX'] = prefix = ''.join(pfx.readlines()).strip()      
  # now assemble file name for station timeseries
  filename = filename.format(**expargs)
  filepath = '{}/{}'.format(folder,filename)
  if not os.path.exists(filepath): IOError(filepath)
  # try to find meta data for gage station from WSC
  if basin is not None and basin_list is not None:
    station_name = station
    station = getGageStation(basin=basin, station=station if WSC_station is None else WSC_station, 
                             basin_list=basin_list) # only works with registered basins
    if station_name is None: station_name = station.name # backup, in case we don't have a HGS station name
    metadata = station.getMetaData() # load station meta data
    if metadata is None: raise GageStationError(name)
  else: 
    metadata = dict(); station = None
    station_name = filename[filename.index('hydrograph.Station_')+1:-4] if station is None else station
  # set meta data (and allow keyword expansion of name and title)
  metadata['problem'] = prefix
  metadata['station_name'] = metadata['long_name']
  if name is not None: name = name.format(**expargs) # name expansion with capitalized keyword arguments
  else: name = 'HGS_{:s}'.format(station_name) if name is None else name
  metadata['name'] = name
  if title is not None: title = title.format(**expargs) # name expansion with capitalized keyword arguments
  else: title = '{long_name:s} (HGS, {problem:s})'.format(**metadata) if title is None else title
  metadata['long_name'] = metadata['title'] = title
  # now determine start data for date_parser
  start_year,start_month,start_day = convertDate(start_date)
  if end_date is not None: end_year,end_month,end_day = convertDate(end_date)
  elif period is not None: end_year = start_year + period; end_month = end_day = 1
  else: raise ArgumentError("Need to specify either 'end_date' or 'period'.")
  if start_day != 1 or end_day != 1: 
    raise NotImplementedError('Currently only monthly data is supported.')
#   import functools
#   date_parser = functools.partial(date_parser, year=start_year, month=start_month, day=start_day)
#   # now load data using pandas ascii reader
#   data_frame = pd.read_table(filepath, sep='\s+', header=2, dtype=np.float64, index_col=['time'], 
#                              date_parser=date_parser, names=ascii_varlist)
#   # resample to monthly data
#   data_frame = data_frame.resample(resampling).agg(np.mean)
#       data = data_frame[flowvar].values
  # parse header
  if varlist is None: varlist = variable_list.copy() # default list 
  with open(filepath, 'r') as f:
      line = f.readline(); lline = line.lower() # 1st line
      if not "hydrograph" in lline: raise GageStationError(line)
      # parse variables and determine columns
      line = f.readline(); lline = line.lower() # 2nd line
      if not "variables" in lline: raise GageStationError(line)
      variable_order = [v.strip('"').lower() for v in line[line.find('"'):].strip().split(',')]
  # figure out varlist and data columns
  if variable_order[0] == 'time': del variable_order[0] # only keep variables
  else: raise GageStationError(variable_order)
  variable_order = [hgs_variables[v] for v in variable_order] # replace HGS names with GeoPy names
  vardict = {v:i for i,v in enumerate(variable_order)} # column mapping
  variable_order = [v for v in variable_order if v in varlist or flow_to_flux[v] in varlist]
  usecols = tuple(vardict[v] for v in variable_order) # variable columns that need to loaded (except time, which is col 0)
  # load data as tab separated values
  data = np.genfromtxt(filepath, dtype=np.float64, delimiter=None, skip_header=3, usecols = (0,)+usecols)
  assert data.shape[1] == len(usecols)+1, data.shape
  time_series = data[:,0]; flow_data = data[:,1:]
  # original time deltas in seconds
  time_diff = time_series.copy(); time_diff[1:] = np.diff(time_series) # time period between time steps
  assert np.all( time_diff > 0 ), time_diff
  time_diff = time_diff.reshape((len(time_diff),len(variable_order))) # reshape to make sure broadcasting works
  # integrate flow over time steps before resampling
  flow_data[1:] -= np.diff(flow_data, axis=0)/2. # get average flow between time steps
  flow_data *= time_diff # integrate flow in time interval by multiplying average flow with time period
  flow_data = np.cumsum(flow_data, axis=0) # integrate by summing up total flow per time interval
  # generate regular monthly time steps
  start_datetime = np.datetime64(dt.datetime(year=start_year, month=start_month, day=start_day), 'M')
  end_datetime = np.datetime64(dt.datetime(year=end_year, month=end_month, day=end_day), 'M')
  time_monthly = np.arange(start_datetime, end_datetime+np.timedelta64(1, 'M'), dtype='datetime64[M]')
  assert time_monthly[0] == start_datetime, time_monthly[0]
  assert time_monthly[-1] == end_datetime, time_monthly[-1] 
  # convert monthly time series to regular array of seconds since start date
  time_monthly = ( time_monthly.astype('datetime64[s]') - start_datetime.astype('datetime64[s]') ) / np.timedelta64(1,'s')
  assert time_monthly[0] == 0, time_monthly[0]
  # interpolate integrated flow to new time axis
  time_series = np.concatenate(([0],time_series), axis=0) # integrated flow at time zero must be zero...
  flow_data = np.concatenate(([[0,]*len(usecols)],flow_data), axis=0) # ... this is probably better than interpolation
  if (time_monthly[-1]-time_series[-1]) > 4*86400.: 
      raise DataError("Data record ends more than 4 days befor end of period: {} days".format((time_monthly[-1]-time_series[-1])/86400.))
  flow_interp = si.interp1d(x=time_series, y=flow_data, kind='linear', axis=0, copy=False, 
                            bounds_error=False, fill_value='extrapolate', assume_sorted=True) 
  flow_data = flow_interp(time_monthly) # evaluate with call
  # compute monthly flow rate from interpolated integrated flow
  flow_data = np.diff(flow_data, axis=0) / np.diff(time_monthly).reshape((len(time_monthly)-1,1))
  flow_data *= 1000 # convert from m^3/s to kg/s
  # construct time axis
  start_time = 12*(start_year - 1979) + start_month -1
  end_time = 12*(end_year - 1979) + end_month -1
  time = Axis(name='time', units='month', atts=dict(long_name='Month since 1979-01'), 
              coord=np.arange(start_time, end_time)) # not including the last, e.g. 1979-01 to 1980-01 is 12 month
  assert len(time_monthly) == end_time-start_time+1
  assert flow_data.shape == (len(time),len(variable_order)), (flow_data.shape,len(time),len(variable_order))
  # construct dataset
  dataset = Dataset(atts=metadata)
  dataset.station = station # add gage station object, if available (else None)
  for i,flowvar in enumerate(variable_order):
      data = flow_data[:,i]
      fluxvar = flow_to_flux[flowvar]
      if flowvar in varlist:
        flowatts = variable_attributes[flowvar]
        # convert variables and put into dataset (monthly time series)
        if flowatts['units'] != 'kg/s': 
          raise VariableError("Hydrograph data is read as kg/s; flow variable does not match.\n{}".format(flowatts))
        dataset += Variable(data=data, axes=(time,), **flowatts)
      if fluxvar in varlist and 'shp_area' in metadata:
        # compute surface flux variable based on drainage area
        fluxatts = variable_attributes[fluxvar]
        if fluxatts['units'] == 'kg/s' and fluxatts['units'] != 'kg/m^2/s': raise VariableError(fluxatts)
        data = data / metadata['shp_area'] # need to make a copy
        dataset += Variable(data=data, axes=(time,), **fluxatts)
  # return completed dataset
  return dataset


# wrapper to load HGS ensembles, otherwise the same
@BatchLoad
def loadHGS_StnEns(station=None, varlist=None, varatts=None, folder=None, name=None, title=None,
                   start_date=1979, end_date=None, period=15, date_parser=date_parser, 
                   basin=None, WSC_station=None, basin_list=None, 
                   filename=station_file, prefix=None, resampling='1M', **kwargs):
  ''' a wrapper for loadHGS_StnTS that supports full inner and outer product expansion of arguments  '''
  return loadHGS_StnTS(station=station, varlist=varlist, varatts=varatts, folder=folder, name=name, 
                       title=name, start_date=start_date, end_date=end_date, period=period, 
                       date_parser=date_parser, basin=basin, WSC_station=WSC_station, prefix=prefix,
                       basin_list=basin_list, filename=filename, resampling=resampling, **kwargs)
  

## abuse for testing
if __name__ == '__main__':

  from projects.WSC_basins import basin_list
#   basin_list = dict(GRW=BasinSet(name='GRW', long_name='Grand River Watershed', rivers=['Grand River'], 
#                                  data_source='Aquanty', stations={'Grand River':['Brantford']}, 
#                                  subbasins=['WholeGRW','UpperGRW','LowerGRW','NorthernGRW','SouthernGRW','WesternGRW']))

  # settings
  basin_name = 'GRW'
  name = 'HGS-{BASIN:s}' # will be expanded
  WSC_station = 'Grand River_Brantford'
#   hgs_folder = '{ROOT_FOLDER:s}/GRW/grw2/erai-g3_d01/clim_15/hgs_run'
  hgs_folder = '{ROOT_FOLDER:s}/GRW/grw2/g-ensemble-2050_d01/annual_15/hgs_run'
  hgs_station = 'GR_Brantford'
  
  # load dataset
  dataset = loadHGS_StnTS(name=name, station=hgs_station, folder=hgs_folder, start_date=1979, 
                          basin=basin_name, WSC_station=WSC_station, basin_list=basin_list, 
                          varlist=['discharge'])
  # and print
  print(dataset)
  print('')
  print(dataset.name)
  print(dataset.prettyPrint(short=True))
  
  # some common operations
  print('')
  clim = dataset.climMean()
  print(clim)
  print(clim.discharge[:]*86400)
#   print(clim.sfroff[:]*86400)