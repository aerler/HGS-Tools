'''
Created on Aug 27, 2016

A module to gauge station data from HGS and associate the data with station meta data from the Water 
Survey of Canada; the data is stored in human-readable text files and tables.

@author: Andre R. Erler, GPL v3
'''

# external imports
import numpy as np
import os
# internal imports
from geodata.misc import ArgumentError, VariableError, DataError, isNumber, DatasetError
from datasets.common import BatchLoad, getRootFolder
from geodata.base import Dataset, Variable, Axis, concatDatasets
from datasets.WSC import getGageStation, GageStationError, loadWSC_StnTS, updateScalefactor
import datetime as dt
# local imports
from hgs.misc import interpolateIrregular, convertDate

## HGS Meta-data

dataset_name = 'HGS'
root_folder = getRootFolder(dataset_name=dataset_name, fallback_name='WRF') # get dataset root folder based on environment variables
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


## function to load HGS station timeseries
def loadHGS_StnTS(station=None, varlist=None, varatts=None, folder=None, name=None, title=None, lcheckComplete=True,
                  start_date=None, end_date=None, run_period=None, period=None, lskipNaN=False, basin=None,  
                  WSC_station=None, basin_list=None, filename=None, prefix=None, scalefactors=None, **kwargs):
  ''' Get a properly formatted WRF dataset with monthly time-series at station locations; as in
      the hgsrun module, the capitalized kwargs can be used to construct folders and/or names '''
  if folder is None or ( filename is None and station is None ): raise ArgumentError
  # try to find meta data for gage station from WSC
  HGS_station = station
  if basin is not None and basin_list is not None:
    station_name = station
    station = getGageStation(basin=basin, station=station if WSC_station is None else WSC_station, 
                             basin_list=basin_list) # only works with registered basins
    if station_name is None: station_name = station.name # backup, in case we don't have a HGS station name
    metadata = station.getMetaData() # load station meta data
    if metadata is None: raise GageStationError(name)
  else: 
    metadata = dict(); station = None; station_name =  None    
  # prepare name expansion arguments (all capitalized)
  expargs = dict(ROOT_FOLDER=root_folder, STATION=HGS_station, NAME=name, TITLE=title,
                 PREFIX=prefix, BASIN=basin, WSC_STATION=WSC_station)
  for key,value in metadata.items():
      if isinstance(value,basestring):
          expargs['WSC_'+key.upper()] = value # in particular, this includes WSC_ID
  if 'WSC_ID' in expargs: 
      if expargs['WSC_ID'][0] == '0': expargs['WSC_ID0'] = expargs['WSC_ID'][1:]
      else: raise DatasetError('Expected leading zero in WSC station ID: {}'.format(expargs['WSC_ID']))
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
  if station_name is None: 
      station_name = filename[filename.index('hydrograph.')+1:-4] if station is None else station
  # set meta data (and allow keyword expansion of name and title)
  metadata['problem'] = prefix
  metadata['station_name'] = metadata.get('long_name', station_name)
  if name is not None: name = name.format(**expargs) # name expansion with capitalized keyword arguments
  else: name = 'HGS_{:s}'.format(station_name)
  metadata['name'] = name
  if title is None: 
    title = ' (HGS, {problem:s})'.format(**metadata)
    if name == name.lower(): title = name.title() + title  # capitalize
    else: title = name + title # assume already correctly capitalized
  else: title = title.format(**expargs) # name expansion with capitalized keyword arguments
  metadata['long_name'] = metadata['title'] = title
  # now determine start data for date_parser
  if end_date is None: 
      if start_date and run_period: end_date = start_date + run_period 
      elif period: end_date = period[1]
      else: raise ArgumentError("Need to specify either 'start_date' & 'run_period' or 'period' to infer 'end_date'.")
  end_year,end_month,end_day = convertDate(end_date)
  if start_date is None: 
      if end_date and run_period: start_date = end_date - run_period 
      elif period: start_date = period[0]
      else: raise ArgumentError("Need to specify either 'end_date' & 'run_period' or 'period' to infer 'start_date'.")
  start_year,start_month,start_day = convertDate(start_date)
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
  if varlist is None: varlist = variable_list[:] # default list 
  with open(filepath, 'r') as f:
      line = f.readline(); lline = line.lower() # 1st line
      if not "hydrograph" in lline: raise GageStationError(line,filepath)
      # parse variables and determine columns
      line = f.readline(); lline = line.lower() # 2nd line
      if not "variables" in lline: raise GageStationError(line)
      variable_order = [v.strip('"').lower() for v in line[line.find('"'):].strip().split(',')]
  # figure out varlist and data columns
  if variable_order[0] == 'time': del variable_order[0] # only keep variables
  else: raise GageStationError(variable_order)
  variable_order = [hgs_variables[v] for v in variable_order] # replace HGS names with GeoPy names
  vardict = {v:i+1 for i,v in enumerate(variable_order)} # column mapping; +1 because time was removed
  variable_order = [v for v in variable_order if v in varlist or flow_to_flux[v] in varlist]
  usecols = tuple(vardict[v] for v in variable_order) # variable columns that need to be loaded (except time, which is col 0)
  assert 0 not in usecols, usecols
  
  # generate regular monthly time steps
  start_datetime = np.datetime64(dt.datetime(year=start_year, month=start_month, day=start_day), 'M')
  end_datetime = np.datetime64(dt.datetime(year=end_year, month=end_month, day=end_day), 'M')
  time_monthly = np.arange(start_datetime, end_datetime+np.timedelta64(1, 'M'), dtype='datetime64[M]')
  assert time_monthly[0] == start_datetime, time_monthly[0]
  assert time_monthly[-1] == end_datetime, time_monthly[-1] 
  
  # load data as tab separated values
  data = np.genfromtxt(filepath, dtype=np.float64, delimiter=None, skip_header=3, usecols = (0,)+usecols)
  assert data.shape[1] == len(usecols)+1, data.shape
  if lskipNaN:
      data = data[np.isnan(data).sum(axis=1)==0,:]
  elif np.any( np.isnan(data) ):
      raise DataError("Missing values (NaN) encountered in hydrograph file; use 'lskipNaN' to ignore.\n('{:s}')".format(filepath))    
  time_series = data[:,0]; flow_data = data[:,1:]
  assert flow_data.shape == (len(time_series),len(usecols)), flow_data.shape
  
  # call function to interpolate irregular HGS timeseries to regular monthly timseries  
  flow_data = interpolateIrregular(old_time=time_series, flow_data=flow_data, new_time=time_monthly, start_date=start_datetime, 
                                   lkgs=True, lcheckComplete=lcheckComplete, usecols=usecols, interp_kind='linear', fill_value=np.NaN)
  
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
  # apply analysis period
  if period is not None:
      dataset = dataset(years=period)
  # adjust scalefactors, if necessary
  if scalefactors:
      if isinstance(scalefactors,dict):
          dataset = updateScalefactor(dataset, varlist=scalefactors, scalefactor=None)
      elif isNumber(scalefactors):
          scalelist = ('discharge','seepage','flow')
          dataset = updateScalefactor(dataset, varlist=scalelist, scalefactor=scalefactors)
      else: 
          raise TypeError(scalefactors) 
        
  # return completed dataset
  return dataset


# an enhanced ensemble loader that supports argument expansion and construction of ensemble datasets
@BatchLoad
def loadHGS_StnEns(ensemble=None, station=None, varlist=None, varatts=None, name=None, title=None, 
                   period=None, run_period=15, folder=None, obs_period=None,  
                   ensemble_list=None, ensemble_args=None, observation_list=None, # ensemble and obs lists for project
                   loadHGS_StnTS=loadHGS_StnTS, loadWSC_StnTS=loadWSC_StnTS, # these can also be overloaded
                   prefix=None, WSC_station=None, basin=None, basin_list=None, **kwargs):
  ''' a wrapper for the regular HGS loader that can also load gage stations and assemble ensembles '''
  if observation_list is None: observation_list = ('obs','observations')
  if ensemble_list is None: ensemble_list = dict() # empty, i.e. no ensembles
  elif not isinstance(ensemble_list, dict): raise TypeError(ensemble_list)
  if ensemble is None: raise ArgumentError("Mandatory argument 'ensemble' is not defined!")
  # decide what to do, based on inputs
  if ensemble.lower() in observation_list:
      # translate parameters
      station = station if WSC_station is None else WSC_station
      period = period if obs_period is None else obs_period
      filetype = 'monthly'
      # load gage station with slightly altered parameters
      dataset = loadWSC_StnTS(station=station, name=name, title=title, basin=basin, basin_list=basin_list, 
                              varlist=varlist, varatts=varatts, period=period, filetype=filetype)
  elif ensemble.lower() in ensemble_list:
      if ensemble_args is None: ensemble_args = dict()
      # loop over list of experiments in ensemble
      ens = []
      for exp in ensemble_list[ensemble]:
          # load individual HGS simulation
          ds = loadHGS_StnTS(station=station, varlist=varlist, varatts=varatts, name=name, title=title, 
                             period=period, ENSEMBLE=exp, run_period=run_period, folder=folder, prefix=prefix, 
                             WSC_station=WSC_station, basin=basin, basin_list=basin_list, **kwargs)
          ens.append(ds)
      # construct ensemble by concatenating time-series
      ensemble_args.setdefault('name',ds.name.replace(exp,ensemble).replace(exp.title(),ensemble.title()))
      ensemble_args.setdefault('title',ds.title.replace(exp,ensemble).replace(exp.title(),ensemble.title())) 
      # N.B.: the ensemble name is constructed by replacing the experiment name in specific dataset names with the ensemble name
      ensemble_args.setdefault('axis','time')
      dataset = concatDatasets(ens, **ensemble_args)
  else:
      # load HGS simulation
      dataset = loadHGS_StnTS(station=station, varlist=varlist, varatts=varatts, name=name, title=title, period=period, 
                              ENSEMBLE=ensemble, run_period=run_period, folder=folder, prefix=prefix, 
                              WSC_station=WSC_station, basin=basin, basin_list=basin_list, **kwargs)
  return dataset


## abuse for testing
if __name__ == '__main__':

#   from projects.WSC_basins import basin_list
  from datasets.WSC import BasinSet
  basin_list = dict(GRW=BasinSet(name='GRW', long_name='Grand River Watershed', rivers=['Grand River'], 
                                 data_source='Aquanty', stations={'Grand River':['Brantford']}, 
                                 subbasins=['WholeGRW','UpperGRW','LowerGRW','NorthernGRW','SouthernGRW','WesternGRW']))

  # settings
  basin_name = 'GRW'
  WSC_station = 'Grand River_Brantford'
  hgs_station = 'GR_Brantford'
  filename = '{PREFIX:s}o.hydrograph.Station_{STATION:s}.dat'

#   test_mode = 'gage_station'
  test_mode = 'dataset'
#   test_mode = 'ensemble'


  if test_mode == 'gage_station':
    
    # load single dataset
    ds = loadWSC_StnTS(station=WSC_station, basin=basin_name, period=(1974,2004), 
                            basin_list=basin_list, filetype='monthly', scalefactors=1e-4)
    print(ds)
    
  elif test_mode == 'dataset':

    hgs_name = 'HGS-{BASIN:s}' # will be expanded
#       hgs_folder = '{ROOT_FOLDER:s}/GRW/grw2/erai-g3_d01/clim_15/hgs_run'
#       hgs_folder = '{ROOT_FOLDER:s}/GRW/grw2/NRCan/annual_15/hgs_run'
    hgs_folder = '{ROOT_FOLDER:s}/GRW/grw2/{EXP:s}{PRD:s}_d{DOM:02d}/{BC:s}{CLIM:s}/hgs_run'

    # load dataset
    dataset = loadHGS_StnTS(station=hgs_station, folder=hgs_folder, filename=filename, 
                            start_date=1979, run_period=10, PRD='', DOM=2, CLIM='clim_15', BC='AABC_', 
                            basin=basin_name, WSC_station=WSC_station, basin_list=basin_list, 
                            lskipNaN=True, lcheckComplete=True, varlist=None, scalefactors=1e-4,
                            EXP='erai-g', name='{EXP:s} ({BASIN:s})')
    # N.B.: there is not record of actual calendar time in HGS, so periods are anchored through start_date/run_period
    # and print
    print(dataset)
    print('')
    print(dataset.name)
    print(dataset.prettyPrint(short=True))
    
    # some common operations
    print('')
    clim = dataset.climMean()
    print(clim)
    # test exact results
    if dataset.name == 'erai-g (GRW)':
        print(clim.discharge[:])
        assert np.allclose(clim.discharge[:], 
                        np.asarray([24793.523584608138, 25172.635322536684, 39248.71087752686, 73361.80217956303, 64505.67974315114, 
                                    32456.80709658126, 18431.93890164255, 15018.095766333918, 16045.543845416256, 17636.665822798554,
                                    18529.952477226405, 22288.711837028015]))
#     print(clim.sfroff[:]*86400)

    
  elif test_mode == 'ensemble':
    
    ens_name = '{ENSEMBLE:s}{PRDSTR:s}'
    ens_folder = '{ROOT_FOLDER:s}/GRW/grw2/{ENSEMBLE:s}{PRDSTR:s}_d{DOM:02d}/{CLIM:s}/hgs_run'
    # actual ensemble definition
    ensemble_list = {'g-mean':('g-ctrl','g-ens-A','g-ens-B','g-ens-C'),
                     't-mean':('t-ctrl','t-ens-A','t-ens-B','t-ens-C')}

    # load an esemble of datasets
    ens = loadHGS_StnEns(ensemble=['g-mean','t-mean'], DOM=1, CLIM='clim_15', run_period=15,
                         period=[(1984,1994),(2050,2060),(2090,2100)], PRDSTR=['','-2050','-2100'],
                         station=hgs_station, name=ens_name, basin=basin_name, filename=filename,
                         WSC_station=WSC_station, basin_list=basin_list, folder=ens_folder,
                         lskipNaN=True, lcheckComplete=True, varlist=None, obs_period=(1974,2004),
                         ens_name='HGS Ensemble', ens_title='HGS Ensemble based on WRF Ensemble',
                         ensemble_list=ensemble_list,
                         outer_list=['ensemble',('period','PRDSTR')], lensemble=True)
    # N.B.: all need to have unique names... whihc is a problem with obs...
    print(ens)
    print('\n')
    print(ens[0])