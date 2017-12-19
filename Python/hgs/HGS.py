'''
Created on Aug 27, 2016

A module to gauge station vardata from HGS and associate the vardata with station meta vardata from the Water 
Survey of Canada; the vardata is stored in human-readable text files and tables.

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
from hgs.misc import interpolateIrregular, convertDate, parseObsWells
# import filename patterns
from hgsrun.misc import hydro_files, well_files, newton_file, water_file

## HGS Meta-vardata

dataset_name = 'HGS'
root_folder = getRootFolder(dataset_name=dataset_name, fallback_name='WRF') # get dataset root folder based on environment variables
prefix_file = 'batch.pfx' # text file that contians the HGS problem prefix (also HGS convention)

# variable attributes and name
variable_attributes_mms = dict(# hydrograph variables
                               surface      = dict(name='discharge', units='m^3/s', atts=dict(long_name='Surface Flow Rate')),      # surface flow rate
                               porous_media = dict(name='seepage'  , units='m^3/s', atts=dict(long_name='Subsurface Flow Rate')),   # subsurface flow rate
                               total        = dict(name='flow'     , units='m^3/s', atts=dict(long_name='Total Flow Rate')),        # total flow rate
                               sfroff = dict(name='sfroff', units='mm/s', atts=dict(long_name='Surface Runoff')),    # surface flow rate over area
                               ugroff = dict(name='ugroff', units='mm/s', atts=dict(long_name='Subsurface Runoff')), # subsurface flow rate over area
                               runoff = dict(name='runoff', units='mm/s', atts=dict(long_name='Total Runoff')),      # total flow rate over area
                               # Newton iteration diagnostics
                               absolute_error = dict(name='error', units='n/a', atts=dict(long_name='Absolute Error')), # I don't know the units...
                               residual_error = dict(name='residual', units='n/a', atts=dict(long_name='Residual Error')),
                               number_of_iterations = dict(name='niter', units='#', atts=dict(long_name='Number of Iterations')),
                               number_of_solver_iterations = dict(name='nsiter', units='#', atts=dict(long_name='Number of Solver Iterations')),
                               time_step = dict(name='delta_t', units='s', atts=dict(long_name='Time Step')),
                               # water balance variables
                               outerboundary = dict(name='outflow', units='m^3/s', atts=dict(long_name='Outer Boundary Flow')),
                               rainfall = dict(name='tot_precip', units='m^3/s', atts=dict(long_name='Basin-integrated Precipitation')),
                               pet = dict(name='tot_pet', units='m^3/s', atts=dict(long_name='Basin-integrated Potential ET')),
                               tot_et = dict(name='tot_et', units='m^3/s', atts=dict(long_name='Basin-integrated Actual ET')),
                               infilt = dict(name='infil', units='m^3/s', atts=dict(long_name='Basin-integrated Infiltration')),
                               exfilt = dict(name='exfil', units='m^3/s', atts=dict(long_name='Basin-integrated Exfiltration')),
                               delta_stor_int = dict(name='delta_storage', units='m^3/s', atts=dict(long_name='Basin-integrated Storage Change')),
                               #TODO: add remaining 11 water balance variables...
                               # observation wells
                               h    = dict(name='head', units='m', atts=dict(long_name='Head at Well')),
                               s    = dict(name='sat', units='', atts=dict(long_name='Relative Saturation')),
                               x    = dict(name='x', units='m', atts=dict(long_name='X Coord.')),
                               y    = dict(name='y', units='m', atts=dict(long_name='Y Coord.')),
                               z    = dict(name='z', units='m', atts=dict(long_name='Elevation (M.S.L.)')),
                               node = dict(name='node', units='', atts=dict(long_name='Node Number')),
                               )
# optional unit conversion
mms_to_kgs = {'mm/s':'kg/m^2/s', 'm^3/s':'kg/s'}
variable_attributes_kgs = dict()
for varname,varatts in variable_attributes_mms.items():
    varatts = varatts.copy()
    varatts['units'] = mms_to_kgs.get(varatts['units'],varatts['units'])
    variable_attributes_kgs[varname] = varatts
# list of variables to load
variable_list = variable_attributes_mms.keys() # also includes coordinate fields    
flow_to_flux = dict(discharge='sfroff', seepage='ugroff', flow='runoff') # relationship between flux and flow variables
# N.B.: computing surface flux rates from gage flows also requires the drainage area
hgs_varmap = {value['name']:key for key,value in variable_attributes_mms.items()}

## function to load HGS station timeseries
def loadHGS_StnTS(station=None, well=None, varlist='default', layers=None, varatts=None, folder=None, name=None, title=None, lcheckComplete=True,
                  start_date=None, end_date=None, run_period=None, period=None, lskipNaN=False, basin=None, lkgs=True,
                  WSC_station=None, basin_list=None, filename=None, prefix=None, scalefactors=None, metadata=None, **kwargs):
  ''' Get a properly formatted WRF dataset with monthly time-series at station locations; as in
      the hgsrun module, the capitalized kwargs can be used to construct folders and/or names '''
  if folder is None or ( filename is None and station is None and well is None ): raise ArgumentError
  # distinguish between special timeseries files, hydrographs, and observation wells
  if metadata is None: metadata = dict()
  if station == 'water_balance' or well == 'water_balance': 
    filename = water_file; well = None; zone = None
    name_tag = 'water_balance'; long_name = 'Integrated Water Balance'
    file_title = 'transient water balance summary'
  elif station == 'newton_info' or well == 'newton_info': 
    filename = newton_file; well = None; zone = None;  lkgs = False # don't change units... too messy
    name_tag = 'newton_info'; long_name = 'Newton Iteration Information'
    file_title = 'transient newton iteration summary'
  elif station and well: raise ArgumentError
  elif well is not None:
    filename = well_files; zone = well # zone is used for file verification later on
    name_tag = well; long_name = well
    file_title = 'flow data at observation well:'
  elif station is not None:
    filename = hydro_files; zone = station
    name_tag = station; long_name = station
    file_title = station+' Hydrograph'
    # try to find meta vardata for gage station from WSC
    if basin is not None and basin_list is not None:
      station = getGageStation(basin=basin, station=station if WSC_station is None else WSC_station, 
                               basin_list=basin_list) # only works with registered basins
      if long_name is None: long_name = station.name # backup, in case we don't have a HGS station name
      metadata = station.getMetaData() # load station meta vardata
      #if metadata is None: raise GageStationError(name)
  else: 
    if filename is None: raise ArgumentError
    long_name = None; name_tag = None; file_title = None
  # prepare name expansion arguments (all capitalized)
  expargs = dict(ROOT_FOLDER=root_folder, STATION=name_tag, WELL=name_tag, NAME=name, TITLE=title,
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
    with open(os.path.join(folder,prefix_file), 'r') as pfx:
      expargs['PREFIX'] = prefix = ''.join(pfx.readlines()).strip()  
  # some more argument expansions
  zone = prefix if zone is None else zone.format(**expargs) # this only applies to newton and water balance files      
  if file_title: file_title = file_title.format(**expargs) # used to validate file header below
  name_tag = name_tag.format(**expargs) # inserted filename as TAG

  # set meta vardata (and allow keyword expansion of name and title)
  metadata['problem'] = prefix
  metadata['station_name'] = metadata.get('long_name', long_name)
  metadata['basin'] = basin if basin else 'n/a'
  if name is not None: name = name.format(**expargs) # name expansion with capitalized keyword arguments
  else: name = 'HGS_{:s}'.format(long_name)
  metadata['name'] = name
  if title is None: 
    title = ' (HGS, {problem:s})'.format(**metadata)
    if name == name.lower(): title = name.title() + title  # capitalize
    else: title = name + title # assume already correctly capitalized
  else: title = title.format(**expargs) # name expansion with capitalized keyword arguments
  metadata['long_name'] = metadata['title'] = title

  # now determine start vardata for date_parser
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
    raise NotImplementedError('Currently only monthly vardata is supported.')
  # generate regular monthly time steps
  start_datetime = np.datetime64(dt.datetime(year=start_year, month=start_month, day=start_day), 'M')
  end_datetime = np.datetime64(dt.datetime(year=end_year, month=end_month, day=end_day), 'M')
  time_monthly = np.arange(start_datetime, end_datetime+np.timedelta64(1, 'M'), dtype='datetime64[M]')
  assert time_monthly[0] == start_datetime, time_monthly[0]
  assert time_monthly[-1] == end_datetime, time_monthly[-1] 
  # construct time axis
  start_time = 12*(start_year - 1979) + start_month -1
  end_time = 12*(end_year - 1979) + end_month -1
  time = Axis(name='time', units='month', atts=dict(long_name='Month since 1979-01'), 
              coord=np.arange(start_time, end_time)) # not including the last, e.g. 1979-01 to 1980-01 is 12 month
  assert len(time_monthly) == end_time-start_time+1

  ## load vardata
  # now assemble file name for station timeseries
  filename = filename.format(PROBLEM=prefix,TAG=name_tag)
  filepath = os.path.join(folder,filename)
  if not os.path.exists(filepath): IOError(filepath)
  # parse header
  with open(filepath, 'r') as f:
      line = f.readline(); lline = line.lower() # 1st line
      if "title" not in lline : raise GageStationError(line,filepath)
      if file_title and file_title.lower() not in lline : 
          raise GageStationError((file_title, line, filepath))
      # parse variables and determine columns
      line = f.readline(); lline = line.lower() # 2nd line
      if not "variables" in lline: raise GageStationError(line)
      variable_order = [v for v in line[line.find('=')+1:].strip().split(',') if len(v) > 0]
      # clean up a little and remove some problematic characters
      variable_order = [v.strip().strip('"').strip().lower() for v in variable_order]
      for c,r in zip((' ','-','(',')'),('_','_','','')):
          variable_order = [v.replace(c,r) for v in variable_order]
      # this line is just for verification
      line = f.readline(); lline = line.lower() # 3rd line
      if "zone" not in lline: raise GageStationError(line,filepath)
      if zone.lower() not in lline: raise GageStationError(line,filepath)
  # figure out varlist and vardata columns
  if variable_order[0].lower() == 'time':
      offset = 1 
      del variable_order[0] # only keep variables
  elif well is not None: 
      offset = 0 # observation wells have different time stamps
  else: raise GageStationError(variable_order)
  if varlist is None: 
      varlist = variable_order[:] # load all in the file
  elif varlist.lower() == 'default': 
      varlist = [varname for varname in variable_attributes_mms.keys() if varname in variable_order] # load all that are known
  else:
      varlist = [hgs_varmap[varname] for varname in varlist] # translate back to internal HGS names
  vardict = {v:i+offset for i,v in enumerate(variable_order)} # column mapping; +1 because time was removed
  variable_order = [v for v in variable_order if v in varlist or flow_to_flux.get(v,False) in varlist]
  usecols = tuple(vardict[v] for v in variable_order) # variable columns that need to be loaded (except time, which is col 0)
  assert offset-1 not in usecols, usecols
  
  # load vardata as tab separated values
  if well:
      # well files need special treatment
      lsqueeze = isinstance(layers,(int,np.integer))
      time_series,data = parseObsWells(filepath, variables=usecols, layers=layers, lskipNaN=lskipNaN)
      assert data.shape[1] == len(usecols), data.shape
      nlay = data.shape[2] # number of layers
      if lsqueeze:
          assert nlay == 1, data.shape 
          data = data.squeeze()
          layer = None
          assert data.shape == (len(time_series),len(usecols),), data.shape
      else:
          layer = Axis(name='i_lay', units='', coord=np.arange(1,nlay+1))
          assert data.shape == (len(time_series),len(usecols),len(layer)), data.shape
  else:
      # all other files follow the same format
      data = np.genfromtxt(filepath, dtype=np.float64, delimiter=None, skip_header=3, usecols = (0,)+usecols)
      assert data.shape[1] == len(usecols)+1, data.shape
      time_series = data[:,0]; data = data[:,1:]
      assert data.shape == (len(time_series),len(usecols)), data.shape
      layer = None # no layer axis
      if lskipNaN:
          data = data[np.isnan(data).sum(axis=1)==0,:]
      elif np.any( np.isnan(data) ):
          raise DataError("Missing values (NaN) encountered in timeseries file; use 'lskipNaN' to ignore.\n('{:s}')".format(filepath))    
  
  # call function to interpolate irregular HGS timeseries to regular monthly timseries  
  data = interpolateIrregular(old_time=time_series, data=data, new_time=time_monthly, start_date=start_datetime, 
                                   lkgs=lkgs, lcheckComplete=lcheckComplete, usecols=usecols, interp_kind='linear', fill_value=np.NaN)
  assert data.shape[0] == len(time), (data.shape,len(time),len(variable_order))
  
  
  ## construct dataset
  dataset = Dataset(atts=metadata)
  dataset.station = station # add gage station object, if available (else None)
  
  # unit options: cubic meters or kg
  if name_tag == 'newton_info' or well is not None:
    flow_units = None; flux_units = None; den = None
    variable_attributes = variable_attributes_mms    
  elif lkgs:
    flow_units = 'kg/s'; flux_units = 'kg/m^2/s'
    variable_attributes = variable_attributes_kgs
    den = metadata.get('shp_area',None) 
  else:
    flow_units = 'm^3/s'; flux_units = 'mm/s'
    variable_attributes = variable_attributes_mms
    den = metadata['shp_area'] / 1000. if 'shp_area' in metadata else None
    
  # create variables
  for i,varname in enumerate(variable_order):
      if layer: 
        vardata = data[:,i,:]
        axes = (time,layer)
      else: 
        vardata = data[:,i]
        axes = (time,)
      # process variable as is first 
      # N.B.: we need to check again, because sometimes we only want the flux variable
      if varname in varlist:
        if varname in variable_attributes:
            varatts = variable_attributes[varname]
        else:
            varatts = dict(name=varname, units=flow_units)
        # convert variables and put into dataset (monthly time series)
        if flow_units and varatts['units'] != flow_units: 
          raise VariableError("Hydrograph vardata is read as kg/s; flow variable does not match.\n{}".format(varatts))
        dataset += Variable(data=vardata, axes=axes, **varatts)
      # process possible flux variable
      fluxvar = flow_to_flux.get(varname,None)      
      if ( fluxvar and fluxvar in varlist ) and ( den and den > 0 ):
        # compute surface flux variable based on drainage area
        if fluxvar in variable_attributes:
            fluxatts = variable_attributes[fluxvar]
        else:
            fluxatts = dict(name=fluxvar, units=flux_units)
        if flux_units and fluxatts['units'] != flux_units: 
          raise VariableError("Hydrograph vardata is read as kg/s; flux variable does not match.\n{}".format(fluxatts))
        vardata = vardata / den # need to make a copy
        dataset += Variable(vardata=vardata, axes=axes, **fluxatts)
        
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
def loadHGS_StnEns(ensemble=None, station=None, well=None, varlist='default', layers=None, varatts=None, 
                   name=None, title=None, period=None, run_period=15, folder=None, obs_period=None,  
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
          ds = loadHGS_StnTS(station=station, well=well, varlist=varlist, layers=layers, varatts=varatts, 
                             name=name, title=title, 
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
      dataset = loadHGS_StnTS(station=station, well=well, varlist=varlist, layers=layers, varatts=varatts, 
                              name=name, title=title, period=period, 
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
  hgs_well = hgs_station = WSC_station= None
  # V1 GRW model
#   hgs_folder = '{ROOT_FOLDER:s}/GRW/grw2/{EXP:s}{PRD:s}_d{DOM:02d}/{BC:s}{CLIM:s}/hgs_run'
#   hgs_station = 'Station_GR_Brantford'; WSC_station = 'Grand River_Brantford'
  # V3 GRW model
  hgs_folder = '{ROOT_FOLDER:s}/GRW/grw2/{EXP:s}{PRD:s}_d{DOM:02d}/{BC:s}{CLIM:s}/hgs_run_v3_wrfpet'
#   hgs_station = '{WSC_ID0:s}'; WSC_station = 'Grand River_Brantford'
#   hgs_station = 'water_balance'; WSC_station = None
#   hgs_station = 'newton_info'; WSC_station = None
  hgs_well = 'W0000347_3'



#   test_mode = 'gage_station'
  test_mode = 'dataset'
#   test_mode = 'ensemble'


  if test_mode == 'gage_station':
    
    # load single dataset
    ds = loadWSC_StnTS(station=WSC_station, basin=basin_name, period=(1974,2004), 
                            basin_list=basin_list, filetype='monthly', scalefactors=1e-3)
    print(ds)
    
  elif test_mode == 'dataset':

    hgs_name = 'HGS-{BASIN:s}' # will be expanded
#       hgs_folder = '{ROOT_FOLDER:s}/GRW/grw2/erai-g3_d01/clim_15/hgs_run'
#       hgs_folder = '{ROOT_FOLDER:s}/GRW/grw2/NRCan/annual_15/hgs_run'
    

    # load dataset
    lkgs = True
    dataset = loadHGS_StnTS(station=hgs_station, well=hgs_well, folder=hgs_folder, layers=None, #[16,17,18], 
                            start_date=1979, run_period=10, PRD='', DOM=2, CLIM='clim_15', BC='AABC_', 
                            basin=basin_name, WSC_station=WSC_station, basin_list=basin_list, lkgs=lkgs,
                            lskipNaN=True, lcheckComplete=True, varlist='default', scalefactors=1e-4,
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
    
    if hgs_station == 'Station_GR_Brantford':
        test_results = np.asarray([24793.523584608138, 25172.635322536684, 39248.71087752686, 73361.80217956303, 64505.67974315114, 
                                   32456.80709658126, 18431.93890164255, 15018.095766333918, 16045.543845416256, 17636.665822798554,
                                   18529.952477226405, 22288.711837028015])
        if not lkgs: test_results /= 1000.
        # test exact results
        if dataset.name == 'erai-g (GRW)':
            print(clim.discharge[:])
            assert np.allclose(clim.discharge[:], test_results)
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