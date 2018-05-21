'''
Created on Aug 27, 2016

A module to gauge station vardata from HGS and associate the vardata with station meta vardata from the Water 
Survey of Canada; the vardata is stored in human-readable text files and tables.

@author: Andre R. Erler, GPL v3
'''

# external imports
import datetime as dt
import numpy as np
import os.path as osp
import os, glob
# internal imports
from geodata.misc import ArgumentError, VariableError, DataError, isNumber, DatasetError, translateSeasons
from datasets.common import BatchLoad, getRootFolder
from geodata.base import Dataset, Variable, Axis, concatDatasets
from datasets.WSC import getGageStation, GageStationError, loadWSC_StnTS, updateScalefactor
from geodata.gdal import loadPickledGridDef, addGDALtoDataset, GridDefinition
from geodata.gdal import grid_folder as common_grid_folder
# Graham's package
from hgs_output import binary
# local imports
from hgs.misc import interpolateIrregular, convertDate, parseObsWells
from hgs.PGMN import loadMetadata, loadPGMN_TS
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
                               outerboundary = dict(name='outflow', units='m^3/s', atts=dict(long_name='Outer Boundary Flow', flip_sign=True)),
                               outer_edge = dict(name='outflow', units='m^3/s', atts=dict(long_name='Outer Boundary Flow', flip_sign=True)),
                               rainfall = dict(name='precip_tot', units='m^3/s', atts=dict(long_name='Basin-integrated Precipitation')),
                               pet_pet = dict(name='pet_pet', units='m^3/s', atts=dict(long_name='Basin-integrated Potential ET')),
                               pet = dict(name='pet_tot', units='m^3/s', atts=dict(long_name='Basin-integrated Potential ET')),
                               tot_et = dict(name='et_tot', units='m^3/s', atts=dict(long_name='Basin-integrated Actual ET', flip_sign=True)),
                               canopy_evap = dict(name='et_can', units='m^3/s', atts=dict(long_name='Basin-integrated Canopy Evaporation', flip_sign=True)),
                               surf_evap = dict(name='et_sfc', units='m^3/s', atts=dict(long_name='Basin-integrated Surface Evaporation', flip_sign=True)),
                               pm_evap = dict(name='evap_pm', units='m^3/s', atts=dict(long_name='Basin-integrated PM Evaporation', flip_sign=True)),
                               pm_trans = dict(name='trans_pm', units='m^3/s', atts=dict(long_name='Basin-integrated PM Transpiration', flip_sign=True)),
                               infilt = dict(name='infil', units='m^3/s', atts=dict(long_name='Basin-integrated Infiltration')),
                               exfilt = dict(name='exfil', units='m^3/s', atts=dict(long_name='Basin-integrated Exfiltration', flip_sign=True)),
                               overland = dict(name='d_olf', units='m^3/s', atts=dict(long_name='Total Change in Surface Flow',)),
                               pm = dict(name='d_pm', units='m^3/s', atts=dict(long_name='Total Subsurface Storage Change',)),
                               #delta_stor_int = dict(name='delta_storage', units='m^3/s', atts=dict(long_name='Basin-integrated Storage Change')),
                               #TODO: add remaining 11 water balance variables...
                               # observation wells
                               h    = dict(name='head', units='m', atts=dict(long_name='Total Head at Well')),
                               s    = dict(name='sat', units='', atts=dict(long_name='Relative Saturation')),
                               )
binary_attributes_mms = dict(# 3D porous medium variables (scalar)
                             p_pm  = dict(name='p_pm', units='m', atts=dict(long_name='Pressure Head (PM)', function='calculate_pressure_head',
                                                                            dependencies=['z_pm','head_pm'], pm=True),),
                             head_pm  = dict(name='head_pm', units='m', atts=dict(long_name='Total Head (PM)', pm=True),),
                             sat_pm   = dict(name='sat', units='', atts=dict(long_name='Relative Saturation', pm=True),),
                             # 3D porous medium variables (vector)
                             v_pm  = dict(name='flow_pm', units='m/s', atts=dict(long_name='Flow Velocity (PM)', pm=True, elemental=True, vector=True)),
                             q_pm  = dict(name='dflx', units='m/s', atts=dict(long_name='Darcy Flux', pm=True, elemental=True, vector=True)),
                             # Overland Flow (2D) variables (scalar)
                             p_olf        = dict(name='p_olf', units='m', atts=dict(long_name='Pressure Head (OLF)', function='calculate_pressure_head',
                                                                                    dependencies=['z','head_olf'])),
                             head_olf     = dict(name='head_olf', units='m', atts=dict(long_name='Total Head (OLF)'),),
                             ExchFlux_olf = dict(name='exflx', units='m/s', atts=dict(long_name='Exchange Flux'),),
                             ETEvap_olf   = dict(name='evap', units='m/s', atts=dict(long_name='Surface Evaporation'),),
                             ETPmEvap_olf = dict(name='evap_pm', units='m/s', atts=dict(long_name='Porous Media Evaporation'),),
                             ETTotal_olf  = dict(name='ET', units='m/s', atts=dict(long_name='Total Evapo-Transpiration'),),
                             ETPmTranspire_olf = dict(name='trans', units='m/s', atts=dict(long_name='Plant Transpiration'),),
                             # Overland Flow (2D) variables (vector)
                             v_olf  = dict(name='flow_olf', units='m/s', atts=dict(long_name='Flow Velocity (OLF)', elemental=True, vector=True)),
                             # derived variables
                             depth2gw  = dict(name='d_gw', units='m', atts=dict(long_name='Depth to Groundwater Table', function='calculate_depth2gw', 
                                                                                dependencies=['coordinates_pm','coordinates_olf','head_pm']),),
                             olf_depth = dict(name='d_olf', units='m', atts=dict(long_name='Overland Flow Depth', function='calculate_olf_depth', 
                                                                                 dependencies=['coordinates_olf','head_olf']),),
                             )
constant_attributes = dict(# variables that are not time-dependent (mostly coordinate variables)
                           vector = dict(name='vector', units='', dtype=np.int64, atts=dict(
                                         long_name='Vector Component', order='0:x, 1:y, 2:z', vector=True)),
                           tensor = dict(name='tensor', units='', dtype=np.int64, atts=dict(long_name='Vector Component', 
                                         order='0:xx, 1:yy, 2:zz, 0:xy, 1:yz, 2:zx', tensor=True)),
                           # Nodal OLF
                           x    = dict(name='x', units='m', atts=dict(long_name='X Coord.')),
                           y    = dict(name='y', units='m', atts=dict(long_name='Y Coord.')),
                           z    = dict(name='zs', units='m', atts=dict(long_name='Surface Elevation')),
                           node = dict(name='node', units='', dtype=np.int64, atts=dict(long_name='Node Number')),
                           # Nodal PM
                           x_pm = dict(name='x_pm', units='m', atts=dict(long_name='X Coord. (PM)', pm=True)),
                           y_pm = dict(name='y_pm', units='m', atts=dict(long_name='Y Coord. (PM)', pm=True)),
                           z_pm = dict(name='z', units='m', atts=dict(long_name='Z Coord.', pm=True)),
                           sheet    = dict(name='sheet', units='', dtype=np.int64, atts=dict(long_name='Sheet Number', pm=True)),
                           nodes_pm = dict(name='node_pm', units='', dtype=np.int64, atts=dict(long_name='3D (PM) Node Number', pm=True)),
                           # Elemental (all)
                           x_elm = dict(name='x_elm', units='m', atts=dict(long_name='X Coord. (Elemental)', elemental=True)),
                           y_elm = dict(name='y_elm', units='m', atts=dict(long_name='Y Coord. (Elemental)', elemental=True)),
                           z_elm = dict(name='zs_elm', units='m', atts=dict(long_name='Surface Elevation (Elemental)', elemental=True)),
                           z_pmelm = dict(name='z_elm', units='m', atts=dict(long_name='Z Coord. (Elemental)', elemental=True, pm=True)),
                           layer   = dict(name='layer', units='', dtype=np.int64, atts=dict(long_name='Layer Number', elemental=True, pm=True)),
                           element = dict(name='element', units='', dtype=np.int64, atts=dict(long_name='2D Element Number', elemental=True, pm=True)),
                           elements_pm = dict(name='elemements_pm', units='', dtype=np.int64, atts=dict(long_name='3D Element Number', elemental=True, pm=True)),
                           elem_k = dict(name='K', units='m/s', atts=dict(long_name='Elemental k', elemental=True, tensor=True, pm=True)),
                           # Time
                           time_ts    = dict(name='time', units='month', dtype=np.int64, atts=dict(long_name='Time since 1979-01')),
                           time_clim  = dict(name='time', units='month', dtype=np.int64, atts=dict(long_name='Time of the Year')),
                           model_time = dict(name='model_time', units='s', atts=dict(long_name='Time since Simulation Start')),
                            )
# optional unit conversion
mms_to_kgs = {'mm/s':'kg/m^2/s', 'm^3/s':'kg/s', 'm/s':'kg/m^2/s'}
variable_attributes_kgs = dict() 
for varname,varatts in variable_attributes_mms.items():
    varatts = varatts.copy()
    varatts['units'] = mms_to_kgs.get(varatts['units'],varatts['units'])
    variable_attributes_kgs[varname] = varatts
binary_attributes_kgs = dict() 
for varname,varatts in binary_attributes_mms.items():
    varatts = varatts.copy()
    if 'flow' not in varatts['name']: # flow stays m/s regardless
        varatts['units'] = mms_to_kgs.get(varatts['units'],varatts['units'])
    binary_attributes_kgs[varname] = varatts
# list of variables to load
variable_list = variable_attributes_mms.keys()
binary_list = binary_attributes_mms.keys()    
flow_to_flux = dict(discharge='sfroff', seepage='ugroff', flow='runoff') # relationship between flux and flow variables
# N.B.: computing surface flux rates from gage flows also requires the drainage area
hgs_varmap = {value['name']:key for key,value in variable_attributes_mms.items()}
bin_varmap = {value['name']:key for key,value in binary_attributes_mms.items()}
const_varmap = {value['name']:key for key,value in constant_attributes.items()}

## function to load HGS station timeseries
def loadHGS_StnTS(station=None, well=None, varlist='default', layers=None, z_layers=None, varatts=None, 
                  folder=None, name=None, title=None, lcheckComplete=True, start_date=None, end_date=None, 
                  run_period=None, period=None, lskipNaN=False, basin=None, lkgs=False, z_axis='z', 
                  time_axis='simple', resample='M', llastIncl=False, WSC_station=None, Obs_well=None, 
                  basin_list=None, filename=None, scalefactors=None, metadata=None, 
                  z_aggregation=None, correct_z=20., conservation_authority=None, **kwargs):
  ''' Get a properly formatted WRF dataset with monthly time-series at station locations; as in
      the hgsrun module, the capitalized kwargs can be used to construct folders and/or names '''
  if folder is None or ( filename is None and station is None and well is None ): raise ArgumentError
  if metadata is None: metadata = dict()
  meta_pfx = None # prefix for station/well attibutes from observations
  # distinguish between special timeseries files, hydrographs, and observation wells
  if station == 'water_balance' or well == 'water_balance': 
    filename = water_file; well = None; zone = None
    name_tag = 'water_balance'; long_name = 'Integrated Water Balance'
    file_title = 'transient water balance summary'
  elif station == 'newton_info' or well == 'newton_info': 
    filename = newton_file; well = None; zone = None;  lkgs = False # don't change units... too messy
    name_tag = 'newton_info'; long_name = 'Newton Iteration Information'
    file_title = 'transient newton iteration summary'
  elif station and well: raise ArgumentError
  elif well:
    filename = well_files; zone = well # zone is used for file verification later on
    name_tag = well; long_name = well; lkgs = False # don't change units!
    file_title = 'flow data at observation well:'
    if conservation_authority:
        meta_pfx = 'PGMN_'
        metadata = loadMetadata(Obs_well if Obs_well else well, conservation_authority=conservation_authority)
  elif station:
    filename = hydro_files; zone = station
    name_tag = station; long_name = station
    file_title = station+' Hydrograph'
    # try to find meta vardata for gage station from WSC
    if basin is not None and basin_list is not None:
        meta_pfx = 'WSC_'
        station = getGageStation(basin=basin, station=station if WSC_station is None else WSC_station, 
                                 basin_list=basin_list) # only works with registered basins
        if long_name is None: long_name = station.name # backup, in case we don't have a HGS station name
        metadata = station.getMetaData() # load station meta vardata
        #if metadata is None: raise GageStationError(name)
  else: 
    if filename is None: raise ArgumentError
    long_name = None; name_tag = None; file_title = None; zone = None
  # prepare name expansion arguments (all capitalized)
  expargs = dict(ROOT_FOLDER=root_folder, STATION=name_tag, WELL=name_tag, NAME=name, TITLE=title,
                 BASIN=basin, WSC_STATION=WSC_station)
  for key,value in metadata.items():
      if isinstance(value,basestring):
          expargs[meta_pfx+key.upper()] = value # in particular, this includes WSC_ID
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
  if expargs.get('PREFIX',None) is None:
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
      if start_date and run_period: 
          start_year,start_month,start_day = convertDate(start_date); del start_month,start_day
          end_date = start_year + run_period 
      elif period: end_date = period[1]
      else: raise ArgumentError("Need to specify either 'start_date' & 'run_period' or 'period' to infer 'end_date'.")
  end_year,end_month,end_day = convertDate(end_date)
  if start_date is None: 
      if end_date and run_period: start_date = end_date - run_period 
      elif period: start_date = period[0]
      else: raise ArgumentError("Need to specify either 'end_date' & 'run_period' or 'period' to infer 'start_date'.")
  start_year,start_month,start_day = convertDate(start_date)
  # generate regular monthly time steps
  start_datetime = np.datetime64(dt.datetime(year=start_year, month=start_month, day=start_day), resample)
  end_datetime = np.datetime64(dt.datetime(year=end_year, month=end_month, day=end_day), resample)
  if llastIncl: end_datetime += np.timedelta64(1, resample)
  time_resampled = np.arange(start_datetime, end_datetime+np.timedelta64(1, resample), dtype='datetime64[{}]'.format(resample))
  assert time_resampled[0] == start_datetime, time_resampled[0]
  assert time_resampled[-1] == end_datetime, time_resampled[-1] 
  # construct time axis
  if time_axis.lower() == 'simple':
      start_time = 12*(start_year - 1979) + start_month -1
      end_time = 12*(end_year - 1979) + end_month -1
      time = Axis(name='time', units='month', atts=dict(long_name='Month since 1979-01'), 
                  coord=np.arange(start_time, end_time)) # not including the last, e.g. 1979-01 to 1980-01 is 12 month
      assert len(time_resampled) == end_time-start_time+1
  elif time_axis.lower() == 'datetime':
      if resample.lower() == 'y': units = 'year'
      elif resample.lower() == 'm': units = 'month'
      elif resample.lower() == 'd': units = 'day'
      elif resample.lower() == 'h': units = 'hour'
      else: units = resample
      long_name = '{}s since {}'.format(units.title(),str(time_resampled[0])) # hope this makes sense...
      time = Axis(name='time', units=units, atts=dict(long_name=long_name), coord=time_resampled[:-1])
  else:
      raise ArgumentError(time_axis)

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
      for c,r in {' ':'_','-':'_','(':'',')':''}.items():
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
  if isinstance(varlist,basestring) and varlist.lower() == 'all': 
      varlist = variable_order[:] # load all in the file
  elif varlist is None or ( isinstance(varlist,basestring) and varlist.lower() == 'default' ): 
      varlist = [varname for varname in variable_attributes_mms.keys() if varname in variable_order] # load all that are known
      varlist += [varname for varname in constant_attributes.keys() if varname in variable_order] # load all constants
  else:
      varlist = [hgs_varmap.get(varname,varname) for varname in varlist] # translate back to internal HGS names
      if z_axis.lower() == 'z' and 'z' not in varlist: varlist.append('z') # needed for vertical coordinate
  vardict = {v:i+offset for i,v in enumerate(variable_order)} # column mapping; +1 because time was removed
  variable_order = [v for v in variable_order if v in varlist or flow_to_flux.get(v,False) in varlist]
  constant_order = [v for v in variable_order if v in constant_attributes]
  variable_order = [v for v in variable_order if v not in constant_attributes]
  varcols = tuple(vardict[v] for v in variable_order) # variable columns that need to be loaded (except time, which is col 0)
  constcols = tuple(vardict[v] for v in constant_order) # constants columns that need to be loaded
  assert offset-1 not in varcols, varcols
  
  # load vardata as tab separated values
  if well:
      # resolve layers
      if isinstance(layers,basestring) and z_layers is None: z_layers = layers
      if isinstance(z_layers,basestring):
          if z_layers.lower() == 'screen':
              if 'z_t' in metadata and 'z_b' in metadata:
                  if 'Screen' not in metadata and 'screen' not in metadata and 'SCREEN' not in metadata: 
                      raise ValueError("Values for 'z_b' and 'z_t' likely do not refer to a screen interval")
                  z_layers = (metadata['z_b'],metadata['z_t'])
              else:
                  raise ValueError("Unable to infer screen interval from well metadata.") 
          else: raise ArgumentError(z_layers)
      # well files need special treatment
      lsqueeze = isinstance(layers,(int,np.integer))
      time_series,data,const,z_s = parseObsWells(filepath, variables=varcols, constants=constcols, 
                                                  layers=layers, z_layers=z_layers, lskipNaN=lskipNaN,
                                                  lelev=True)
      assert data.shape[1] == len(varcols), data.shape
      nlay = data.shape[2] # number of layers
      if isinstance(correct_z, (int,np.integer,float,np.inexact)) and not isinstance(correct_z, (bool,np.bool_)):
          if 'screen_depth' not in metadata: 
              raise AttributeError('Need screen_depth attribute to correct based on elevation and depth!')
          correct_z =  ( metadata['screen_depth'] < correct_z )
      if correct_z:
          i_h = variable_order.index('h')
          bias = metadata['ELVA_GROUN'] - z_s
          if data.ndim == 2: data[:,i_h] += bias
          else: data[:,i_h,:] += bias
      if lsqueeze:
          assert nlay == 1, data.shape 
          data = data.squeeze()
          sheet = None
          assert data.shape == (len(time_series),len(varcols),), data.shape
      else:
          # create vertical axis
          if z_axis.lower() == 'z':
              iz = constant_order.index('z',)
              sheet = Axis(coord=const[iz,:], **constant_attributes['z_pm'])
          elif z_axis.lower() == 'i': 
              sheet = Axis(name='i_lay', units='', coord=np.arange(1,nlay+1))
          else:
              raise ArgumentError(z_axis)
          assert data.shape == (len(time_series),len(varcols),len(sheet)), data.shape
  else:
      # all other files follow the same format
      assert len(constcols) == 0, constcols
      data = np.genfromtxt(filepath, dtype=np.float64, delimiter=None, skip_header=3, usecols = (0,)+varcols)
      assert data.shape[1] == len(varcols)+1, data.shape
      if lskipNaN:
          data = data[np.isnan(data).sum(axis=1)==0,:]
      elif np.any( np.isnan(data) ):
          raise DataError("Missing values (NaN) encountered in timeseries file; use 'lskipNaN' to ignore.\n('{:s}')".format(filepath))    
      time_series = data[:,0]; data = data[:,1:]
      assert data.shape == (len(time_series),len(varcols)), data.shape
      sheet = None # no sheet axis
  
  # call function to interpolate irregular HGS timeseries to regular monthly timseries  
  data = interpolateIrregular(old_time=time_series, lkgs=lkgs, data=data, new_time=time_resampled, 
                              start_date=start_datetime, interp_kind='linear', 
                              lcheckComplete=lcheckComplete, usecols=varcols, fill_value=np.NaN)
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
        
  # add constants to dataset (only for wells at the moment)
  for i,varname in enumerate(constant_order):
    if varname != sheet.name: # skip z-axis
      if sheet: 
          vardata = const[i,:]
          # check if there is a sheet-dependence
          if np.all(vardata == vardata[0]): vardata = vardata[0]
      else: 
          vardata = const[i] # an actual constant...
      #dataset.atts[varname] = np.asscalar(vardata) # just a constant... 
      axes = () if vardata.size == 1 else (sheet,)
      if varname in constant_attributes: varatts = constant_attributes[varname]
      else: varatts = dict(name=varname, units=flow_units)
      dataset += Variable(data=vardata, axes=axes, plotatts_dict={}, **varatts) # add variable
  
  # create variables
  for i,varname in enumerate(variable_order):
      if sheet: 
        vardata = data[:,i,:]
        axes = (time,sheet)
      else: 
        vardata = data[:,i]
        axes = (time,)
      # process variable as is first 
      # N.B.: we need to check again, because sometimes we only want the flux variable
      if varname in varlist:
        if varname in variable_attributes: varatts = variable_attributes[varname]
        else: varatts = dict(name=varname, units=flow_units)
        # convert variables and put into dataset (monthly time series)
        if flow_units and varatts['units'] != flow_units: 
          raise VariableError("Hydrograph vardata is read as kg/s; flow variable does not match.\n{}".format(varatts))
        # flip sign for some variables
        if varatts['atts'].get('flip_sign',False): vardata *= -1
        dataset += Variable(data=vardata, axes=axes, plotatts_dict={}, **varatts)
      # process possible flux variable
      fluxvar = flow_to_flux.get(varname,None)      
      if ( fluxvar and fluxvar in varlist ) and ( den and den > 0 ):
        # compute surface flux variable based on drainage area
        if fluxvar in variable_attributes: fluxatts = variable_attributes[fluxvar]
        else: fluxatts = dict(name=fluxvar, units=flux_units)
        if flux_units and fluxatts['units'] != flux_units: 
          raise VariableError("Hydrograph vardata is read as kg/s; flux variable does not match.\n{}".format(fluxatts))
        if varatts['atts'].get('flip_sign',False) and not fluxatts['atts'].get('flip_sign',False):
            raise VariableError("If the variable sign has been flipped, the sign of the flux variable has to be flipped, too.")
        vardata = vardata / den # need to make a copy
        dataset += Variable(vardata=vardata, axes=axes, plotatts_dict={}, **fluxatts)
        
  # apply analysis period
  if period is not None:
      dataset = dataset(years=period)
      
  # aggregate vertical axis
  if z_aggregation and sheet: 
      if not hasattr(np, z_aggregation):
          raise ArgumentError("'z_aggregation' has to be a valid numpy function.")
      dataset = getattr(dataset, z_aggregation)(axis=sheet.name)      
      
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
                   ensemble_list=None, ensemble_args=None, observation_list=None, conservation_authority=None,# ensemble and obs lists for project
                   loadHGS_StnTS=loadHGS_StnTS, loadWSC_StnTS=loadWSC_StnTS, # these can also be overloaded
                   WSC_station=None, Obs_well=None, basin=None, basin_list=None, **kwargs):
  ''' a wrapper for the regular HGS loader that can also load gage stations and assemble ensembles '''
  if observation_list is None: observation_list = ('obs','observations')
  if ensemble_list is None: ensemble_list = dict() # empty, i.e. no ensembles
  elif not isinstance(ensemble_list, dict): raise TypeError(ensemble_list)
  if ensemble is None: raise ArgumentError("Mandatory argument 'ensemble' is not defined!")
  # decide what to do, based on inputs
  if ensemble.lower() in observation_list:
      # translate parameters
      station = station if WSC_station is None else WSC_station
      well = well if Obs_well is None else Obs_well
      obs_period = obs_period or period
      filetype = 'monthly'
      # load gage station with slightly altered parameters
      if station and well: raise ArgumentError()
      elif station:
          dataset = loadWSC_StnTS(station=station, name=name, title=title, basin=basin, basin_list=basin_list, 
                                  varlist=varlist, varatts=varatts, period=obs_period, filetype=filetype)
      elif well:
          dataset = loadPGMN_TS(well=well, name=name, title=title, varlist=varlist, varatts=varatts,
                                conservation_authority=conservation_authority,)
          if obs_period: dataset = dataset(years=obs_period)
      else:
          raise ArgumentError(ensemble)
  elif ensemble.lower() in ensemble_list:
      if ensemble_args is None: ensemble_args = dict()
      # loop over list of experiments in ensemble
      ens = []
      for exp in ensemble_list[ensemble]:
          # load individual HGS simulation
          ds = loadHGS_StnTS(station=station, well=well, varlist=varlist, layers=layers, varatts=varatts, 
                             name=name, title=title, conservation_authority=conservation_authority,
                             period=period, experiment=exp, run_period=run_period, folder=folder, 
                             WSC_station=WSC_station, Obs_well=Obs_well, basin=basin, basin_list=basin_list, 
                             **kwargs)
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
                              name=name, title=title, period=period, conservation_authority=conservation_authority, 
                              experiment=ensemble, run_period=run_period, folder=folder,
                              WSC_station=WSC_station, Obs_well=Obs_well, basin=basin, basin_list=basin_list,
                               **kwargs)
  return dataset


## load functions for binary data


## function to interpolate nodal/elemental datasets to a regular grid
def gridDataset(dataset, griddef=None, basin=None, subbasin=None, shape_file=None,  
                basin_list=None, grid_folder=None):
  ''' interpolate nodal/elemental datasets to a regular grid, add GDAL, and mask to basin outlines '''
  if isinstance(griddef,basestring): 
      if grid_folder is None: grid_folder = common_grid_folder  
      griddef = loadPickledGridDef(grid=griddef, folder=grid_folder)
  elif not isinstance(griddef,GridDefinition):
      raise ArgumentError(griddef)
  # identify shape file
  if basin and basin_list and shape_file is None:
    basininfo = basin_list[basin]
    shape_file = basininfo.shapefiles[subbasin if subbasin else basininfo.outline]
    if not os.path.exists(shape_file): 
        raise IOError(shape_file)
  # interpolate (regular nodal values first)
  if 'x' in dataset and 'y' in dataset: 
      dataset = dataset.gridDataset(grid_axes=(griddef.ylat,griddef.xlon), method='cubic')
  # also interpolate elemental values, if present
  if 'x_elm' in dataset and 'y_elm' in dataset: 
      dataset = dataset.gridDataset(grid_axes=(griddef.ylat,griddef.xlon), 
                                    coord_map=dict(x='x_elm',y='y_elm'), method='cubic')
  # add GDAL
  dataset = addGDALtoDataset(dataset=dataset, griddef=griddef, )
  # mask basin shape
  if shape_file and dataset.gdal:
      dataset.maskShape(name=basin, filename=shape_file, invert=True)
  # return gridded and masked dataset
  return dataset


## function to compute pressure head in loadHGS
def calculate_pressure_head(head_pm=None, head_olf=None, z_pm=None, z=None):
    ''' a function to calculate pressure head from total head and elevation '''
    if head_pm is not None and z_pm is not None:
        assert isinstance(z_pm,np.ndarray)
        if hasattr(head_pm, 'values'): head_pm = head_pm.values.reshape(z_pm.shape)
        p = head_pm - z_pm
    elif head_olf is not None and z is not None:
        assert isinstance(z,np.ndarray)
        if hasattr(head_olf, 'values'): head_olf = head_olf.values.reshape(z.shape)
        p = head_olf - z
    else:
        raise ArgumentError()
    return p
  
## function to load HGS binary data
@BatchLoad
def loadHGS(varlist=None, folder=None, name=None, title=None, basin=None, season=None, 
            lgrid=False, griddef=None, subbasin=None, shape_file=None, t_list=None, lflipdgw=False, 
            mode='climatology', file_mode='last_12', file_pattern='{PREFIX}o.head_olf.????',  
            lkgs=False, varatts=None, constatts=None, lstrip=True, lxyt=True, grid_folder=None, 
            basin_list=None, metadata=None, conservation_authority=None, 
            override_k_option='Anisotropic Elemental K', lallelem=False, **kwargs):
  ''' Get a properly formatted WRF dataset with monthly time-series at station locations; as in
      the hgsrun module, the capitalized kwargs can be used to construct folders and/or names '''
  if folder is None: raise ArgumentError
  if metadata is None: metadata = dict()
  # unit options: cubic meters or kg  
  varatts = ( varatts or ( binary_attributes_kgs if lkgs else binary_attributes_mms ) ).copy()
  constatts = constatts or constant_attributes
  if varlist is None: varlist = binary_list

  # prepare name expansion arguments (all capitalized)
  expargs = dict(ROOT_FOLDER=root_folder, NAME=name, TITLE=title, BASIN=basin,)
  for key,value in metadata.items():
      if isinstance(value,basestring): expargs[key.upper()] = value 
  # exparg preset keys will get overwritten if capitalized versions are defined
  for key,value in kwargs.items():
    KEY = key.upper() # we only use capitalized keywords, and non-capitalized keywords are only used/converted
    if KEY == key or KEY not in kwargs: expargs[KEY] = value # if no capitalized version is defined
  # read folder and infer prefix, if necessary
  folder = folder.format(**expargs)
  if not os.path.exists(folder): raise IOError(folder)
  if expargs.get('PREFIX',None) is None:
    with open(os.path.join(folder,prefix_file), 'r') as pfx:
      expargs['PREFIX'] = prefix = ''.join(pfx.readlines()).strip()  
  # set meta vardata (and allow keyword expansion of name and title)
  metadata['problem'] = prefix
  metadata['basin'] = basin if basin else 'n/a'
  if name is not None: name = name.format(**expargs) # name expansion with capitalized keyword arguments
  else: name = 'HGS Binary Fields'
  metadata['name'] = name
  if title is None: 
    title = ' (HGS, {problem:s})'.format(**metadata)
    if name == name.lower(): title = name.title() + title  # capitalize
    else: title = name + title # assume already correctly capitalized
  else: title = title.format(**expargs) # name expansion with capitalized keyword arguments
  metadata['long_name'] = metadata['title'] = title

  # find files/time-steps to load
  if not t_list:
      glob_folder = osp.join(folder.format(**expargs),file_pattern.format(**expargs))
      file_list = glob.glob(glob_folder)
      if len(file_list) == 0: 
          raise DataError("No binary output files found:\n '{}'".format(glob_folder))
      t_list = [int(f[-4:]) for f in file_list]
      t_list.sort()      
  if mode.lower()[:4] == 'clim':
      if file_mode.lower() == 'last_12':
          if len(t_list) < 12: 
            raise ValueError("Need at least 12 time-steps to assemble Monthly climatology; {} given".format(len(t_list)))
          t_list = t_list[-12:]            
      if len(t_list) != 12: 
            raise ValueError("Need at exactly 12 time-steps to assemble Monthly climatology; {} given".format(len(t_list)))
      # construct time Axis
      time = Axis(coord=np.arange(1, 13), **constatts['time_clim'])
      # extract a season
      if season: 
          season = translateSeasons(season)    
          # slice time axis
          t_list = [t_list[i] for i in season]
          time = time(time=season, lidx=True)          
  elif mode.lower()[:4] == 'time':
      if file_mode.lower() == 'last_12':
          if len(t_list) < 12: 
            raise ValueError("Need at least 12 time-steps to select last year; {} given".format(len(t_list)))
          t_list = t_list[-12:]
      # construct time axis
      time = Axis(coord=np.arange(1, len(t_list)+1), **constatts['time_ts']) 
  else: raise NotImplementedError(mode)
  te = len(time)
  assert len(t_list) == te, (len(t_list),te)
  
  # save some more metadata
  metadata['prefix'] = prefix
  metadata['HGS_folder'] = folder
  metadata['time_list'] = t_list

  # different varlist for different purposes
  if lstrip:
      # the variables we want at the end (to clean up)
      final_varlist = set(['x','y','x_elm','y_elm','model_time']) if lxyt else set() 
      for var in varlist:
          if var not in bin_varmap and var not in const_varmap:
              raise ArgumentError("When using variable stripping, only GeoPy names can be used," +
                                  " HGS names, such as '{}', do not work!".format(var)) 
          final_varlist.add(var)
  varlist = [bin_varmap.get(var,var) for var in varlist] # translate to HGS var names
  varlist = [const_varmap.get(var,var) for var in varlist] # translate to HGS var names
  # make sure variable dependencies work (derived variables last)
  all_deps = dict()
  for var in varlist[:]: # use copy to avoid interference
      if var in varatts:
          deplist = varatts[var]['atts'].get('dependencies',None)
          if deplist:
              varlist.remove(var)
              for depvar in deplist:
                  if ( depvar not in varlist and depvar not in constant_attributes and
                       depvar not in ('coordinates_pm','coordinates_olf') ): 
                      varlist.append(depvar)
                  all_deps[depvar] = None
              varlist.append(var)
      else:
          if var not in constatts: 
              raise VariableError("Variable '{}' not found in variable attributes.".format(var))
          
  ## construct dataset
  dataset = Dataset(atts=metadata)
  dataset += time
        
  ## load vardata using Graham's hgs_output package
  # load first time step to create coordinate arrays etc.
  prefixo = prefix+'o'
  reader = binary.IO(prefixo,folder,t_list[0])
  coords_pm = reader.read_coordinates_pm()
  coords_olf = reader.read_coordinates_olf(coords_pm)
  # create mesh axes
  ne = len(coords_olf)
  se = coords_pm['sheet'].max()
  # load necessary nodal coordinates
  if not lallelem or any([var in varlist for var in ('x','y','z','z_pm','nodes_pm')]):
      node_ax = Axis(coord=np.arange(1,ne+1), **constatts['node'])
      assert ne == len(node_ax)
      dataset += node_ax
  if not lallelem or any([var in varlist for var in ('z_pm','nodes_pm')]):
      sheet_ax = Axis(coord=np.arange(coords_pm['sheet'].min(),se+1), **constatts['sheet'])
      assert se == len(sheet_ax)
      dataset += sheet_ax
  # add coordinate fields (surface)
  for var in ('x','y','z'):
      if not lallelem or var in varlist:
          dataset += Variable(data=coords_olf[var].values, axes=(node_ax,), **constatts[var])
  # extract coordinate fields (porous medium)
  if not lallelem or 'z_pm' in varlist:
      dataset += Variable(data=coords_pm['z'].values.reshape((se,ne)), axes=(sheet_ax,node_ax), 
                          **constatts['z_pm'])
  if not lallelem or 'nodes_pm' in varlist:      
      dataset += Variable(data=coords_pm.index.values.reshape((se,ne)), axes=(sheet_ax,node_ax), 
                          **constatts['nodes_pm'])
  nne = len(coords_pm)
  assert ne*se == nne
  vector_ax = Axis(coord=np.arange(3), **constatts['vector'])

  # load elemental coordinate, if necessary
  lelem = False; lelem3D = False
  for hgsvar in varlist:
      atts = constatts[hgsvar]['atts'] if hgsvar in constatts else varatts[hgsvar]['atts']
      if atts.get('elemental',False) or lallelem: 
          lelem = True
          lelem3D = lelem3D or atts.get('pm',False)  
          if lallelem:
              atts['interp_elem'] = ( not atts.get('elemental',False) )
#   lelem3D = lelem3D or 'z_elm' in final_varlist or 'K' in final_varlist
#   lelem = lelem or 'zs_elm' in final_varlist or lelem3D
  
  if lelem:
      elem_olf = reader.read_elements(domain='olf')
      nelem = len(elem_olf)
      elem_ax = Axis(coord=np.arange(1,nelem+1), **constatts['element'])
      dataset += elem_ax
      elem_coords_olf = reader.compute_element_coordinates(elements=elem_olf, coords_pm=coords_pm, 
                                                           coord_list=('x','y','z'), lpd=False)
      if lallelem: 
          elem_olf_offset = elem_olf - (se-1)*ne
          # N.B.: in principle it would be possible to look up the corresponding PM and OLF nodes,
          #       and replace PM indices with OLF indices, but that is extremely inefficient, and
          #       based on some testing it appears save to assume that the order of nodes is the 
          #       same in each sheet (and the OLF domain), so that simple subtraction should work.
          #       Nevertheless, it is still saver to test this, at least a little...
          assert elem_olf_offset.values[:,1:3].min() == 1
          assert elem_olf_offset.values[:,1:3].max() == ne
      # add surface element coordinate fields (x, y, and surface elevation zs)
      for var in ('x','y','z'):
          dataset += Variable(data=elem_coords_olf[var], axes=(elem_ax,), **constatts[var+'_elm'])
          # N.B.: 'z' is actually 'zs' but this is already taken care of in constant_attributes
      # add 3D element coordinates
      if lelem3D:
          elem_pm  = reader.read_elements(domain='pm')
          nlay = len(elem_pm)/nelem
          assert nelem*nlay == len(elem_pm) 
          assert nlay+1 == se # there is one extra sheet 
          layer_ax = Axis(coord=np.arange(1,nlay+1), **constatts['layer'])
          dataset += layer_ax
          # add 3D element numbers
          dataset += Variable(data=elem_pm.index.values.reshape((nlay,nelem)), axes=(layer_ax,elem_ax), 
                              **constatts['elements_pm'])
          # add 3D element elevation (z coordinate), if requested
          if 'z_elm' in final_varlist:
              elem_coords_pm = reader.compute_element_coordinates(elements=elem_pm, coords_pm=coords_pm, 
                                                                  coord_list=('z',), lpd=False)
              dataset += Variable(data=elem_coords_pm['z'].reshape((nlay,nelem)), axes=(layer_ax,elem_ax,), 
                                  **constatts['z_pmelm']) # 'z' is already used for the 'zs' variable
              assert np.all(np.diff(dataset['z_elm'][:], axis=0) > 0) 
          # add elemental K
          if 'K' in final_varlist:
              elem_k = reader.read_k(override_k_option=override_k_option)
              nten = len(elem_k.columns)
              lk = len(elem_k)
              if lk != nlay*nelem: 
                  raise ValueError("Number of K values ({}) does not match number of elements ({}); consider overriding the k_option.".format(lk,nlay*nelem))
              if nten == 1:
                  dataset += Variable(data=elem_k.values.reshape((nlay,nelem)), axes=(layer_ax,elem_ax,), 
                                      **constatts['elem_k']) # scalar value
              else:
                  tensor_ax = Axis(coord=np.arange(1,nten+1), **constatts['tensor'])
                  dataset += Variable(data=elem_k.values.reshape((nlay,nelem,nten)), 
                                      axes=(layer_ax,elem_ax,tensor_ax), **constatts['elem_k']) # tensor value
          # N.B.: currently elements are assumed to be organized in vertically symmetric layers
          
          
  # remove constant variables from varlist (already loaded)
  varlist = [var for var in varlist if var not in constatts]  
     
  # initialize variables
  load_varlist = []
  for hgsvar in varlist:
      atts = varatts[hgsvar]
      aa = atts['atts']
      # elemental variables are currently not supported
      aa['HGS_name'] = hgsvar # needed later
      if aa.get('elemental',False) or aa.get('interp_elem',False):
          axes = (elem_ax,)
          if aa.get('pm',False): axes = (layer_ax,)+axes
      else:
          axes = (node_ax,)
          if aa.get('pm',False): axes = (sheet_ax,)+axes
      if aa.get('vector',False): axes = axes+(vector_ax,)
      axes = (time,)+axes
      shape = tuple([len(ax) for ax in axes]) 
      # save name and variable
      dataset += Variable(data=np.zeros(shape),axes=axes, **atts)
      load_varlist.append(atts['name'])
  # add simulation time variable
  dataset += Variable(data=np.zeros((te,)), axes=(time,), **constatts['model_time'])
  
  # now fill in the remaining data
  for i,t in enumerate(t_list):
      if i > 0: # reuse old reader for first step 
          reader = binary.IO(prefixo,folder,t)
      # reset dependencies
      for depvar in all_deps.keys(): all_deps[depvar] = None
      all_deps['coordinates_pm'] = coords_pm
      all_deps['coordinates_olf'] = coords_olf
      for depvar in ['z_pm','z']:
          if depvar in all_deps:
              all_deps[depvar] = dataset[constatts[depvar]['name']][:] # these are just arrays
      sim_time = reader.read_timestamp()
      # loop over variables
      for var in load_varlist:
          variable = dataset[var]; aa = variable.atts; hgsvar = aa['HGS_name']
          linterp = aa.get('interp_elem',False); l3d = variable.atts.get('pm',False)
          shp3d = (nlay,nelem) if lelem3D and ( aa.get('elemental',False) or linterp ) else (se,ne)
          # load data
          if aa.get('function',False):
              deplist = {depvar:all_deps[depvar] for depvar in aa.get('dependencies',[])}
              fct_name = aa['function']; _locals = globals()
              if fct_name in _locals:
                  data = _locals[fct_name](**deplist)
              elif hasattr(reader, fct_name):
                  df = getattr(reader,fct_name)(**deplist)                  
                  if l3d: 
                      if linterp:
                          data = reader.interpolate_node2element(df, elements=elem_pm, lpd=False)
                      else: data = df.values
                      data = data.reshape(shp3d)
                  else: 
                    if linterp:
                        data = reader.interpolate_node2element(df, elements=elem_olf_offset, lpd=False)
                    else: data = df.values
                    data = data.squeeze()
                    if data.size == variable.shape[1]+1: 
                        print("Warning: Trimming first element of {} array.".format(var))
                        data = data[1:] 
              else:
                  raise NotImplementedError(fct_name)
              variable.data_array[i,:] = data       
          elif aa.get('vector',False):
              df = reader.read_vec(hgsvar)
              if linterp:
                  if l3d: data = reader.interpolate_node2element(df, elements=elem_pm, lpd=False)
                  else: data = reader.interpolate_node2element(df, elements=elem_olf_offset, lpd=False)
              else: data = df.values
              if l3d: data = data.reshape(shp3d+(3,))   
              variable.data_array[i,:] = data
#               else:        
#                   for j in range(3): # transposing vectors
#                       variable.data_array[i,j,:] = data[:,j]
          else:
              # check simulation time
              st = reader.read_timestamp(var=hgsvar)
              if sim_time != st:
                  raise ValueError("Timestamps in output files are not consistent between variables: {} != {} ({})".format(sim_time,st,hgsvar))
              # read actual binary 2D or 3D data
              if l3d:
                  df = reader.read_var(hgsvar, nne)
                  if linterp:
                      data = reader.interpolate_node2element(df, elements=elem_pm, lpd=False)
                  else: data = df.values
                  variable.data_array[i,:] = data.reshape(shp3d)
              else:
                  df = reader.read_var(hgsvar, ne)
                  if linterp:
                      data = reader.interpolate_node2element(df, elements=elem_olf_offset, lpd=False)
                  else: data = df.values
                  variable.data_array[i,:] = data.squeeze()
          # save dependencies
          if var in all_deps: all_deps[var] = df # dataframes are not interpolated to elements
      # save timestamp
      dataset['model_time'].data_array[i] = sim_time              
    
  # now remove all unwanted variables...
  if lstrip:
      # clean up variables
      for var in dataset.variables.keys():
          if var not in final_varlist: dataset.removeVariable(var)
      # clean up axes
      for ax in dataset.axes.keys():
          if ax not in final_varlist: dataset.removeAxis(ax, force=False) 
          # N.B.: force=False means only remove unused axes
  
  # do some multi-purpose slicing
  tensor_idx = dict(x=0,y=1,z=2,xx=0,yy=1,zz=2,xy=3,yz=4,zx=5)
  slc_axes = dict()
  for axname,axval in kwargs.items():
      if axname in dataset.axes and axval is not None:
          # some special values...
          if axname in ('vector','tensor'):
              axval = tensor_idx.get(axval,axval)
          if axval < 0: 
              axval = dataset.axes[axname].max() + axval +1
          slc_axes[axname] = axval
  dataset = dataset(**slc_axes)
      
  # flip sign of depth to groundwater variable
  if lflipdgw:
      d_gw = binary_attributes_mms['depth2gw']['name']
      if d_gw in dataset:
          dataset[d_gw] *= -1

  # interpolate to regular grid      
  if lgrid:
      dataset = gridDataset(dataset, griddef=griddef, basin=basin, subbasin=subbasin, 
                            shape_file=shape_file, basin_list=basin_list, grid_folder=grid_folder) 
  
  # return completed dataset
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
#   test_mode = 'dataset_regrid'
  test_mode = 'binary_dataset'
#   test_mode = 'time_axis'
#   test_mode = 'station_dataset'
#   test_mode = 'station_ensemble'


  if test_mode == 'gage_station':
    
    # load single dataset
    if WSC_station:
        ds = loadWSC_StnTS(station=WSC_station, basin=basin_name, period=(1974,2004), 
                                basin_list=basin_list, filetype='monthly', scalefactors=1e-3)
    elif hgs_well:
        ds = loadPGMN_TS(well=hgs_well, conservation_authority='GRCA',)
    print(ds)
  
  
  elif test_mode == 'dataset_regrid':

    # load dataset
    dataset = loadHGS(folder=hgs_folder, varlist=['dflx','zs'], conservation_authority='GRCA', 
                      sheet=-2, layer=-3, vector='z', lallelem=True,
                      PRD='', DOM=2, CLIM='clim_15', BC='AABC_', basin=basin_name, basin_list=basin_list, 
                      lkgs=False, EXP='erai-g', name='{EXP:s} ({BASIN:s})',
                      lgrid=True, griddef='grw3',)
    # load griddefition
    # interpolate to regular x,y-grid
#     x = Axis(name='x', units='m', coord=np.linspace(dataset.x.min(),dataset.x.max(),200))
#     y = Axis(name='y', units='m', coord=np.linspace(dataset.y.min(),dataset.y.max(),200))
#     dataset = gridDataset(dataset, griddef='grw2', basin=basin_name, basin_list=basin_list, 
#                           grid_folder=grid_folder)
    # and print
    print('')
    print(dataset)
    print('')
    print(dataset.dflx)
    print(dataset.dflx[2,15,:])
    
  
  elif test_mode == 'binary_dataset':

    # load dataset
    vecvar = 'dflx'
    #hgs_folder = '{ROOT_FOLDER:s}/GRW/grw2/{EXP:s}{PRD:s}_d{DOM:02d}/{BC:s}{CLIM:s}/hgs_run_deep'
    dataset = loadHGS(varlist=['d_gw','zs_elm',], EXP='g-ensemble', name='{EXP:s} ({BASIN:s})', 
                      lallelem=True, sheet=-2, layer=(12,17), season='MAM', tensor='z', vector='z', 
                      folder=hgs_folder, conservation_authority='GRCA', 
                      PRD='', DOM=2, CLIM='clim_15', BC='AABC_', 
                      basin=basin_name, basin_list=basin_list, lkgs=False, )
    # N.B.: there is no record of actual calendar time in HGS, so periods are anchored through start_date/run_period
    # and print
    print('')
    print(dataset)
    print(dataset.atts.HGS_folder)
    print('')
    print(dataset.model_time)
    print(dataset.model_time[:])
    # inspect layers and depth
    if dataset.hasAxis('layer'):
        print('')
        print(dataset.layer)
        print(dataset.layer[:])
    if 'z_elm' in dataset and 'zs_elm' in dataset:
        d0 = dataset.zs_elm - dataset.z_elm(layer=dataset.layer.max())
        print(d0)
        print(d0.min(),d0.mean(),d0.max())
    if vecvar in dataset:
        print('')
        vec = dataset[vecvar].mean(axis=('element','time'))
        print(vec)
        print(vec[:])
    
  
  elif test_mode == 'time_axis':

    # load dataset
    dataset = loadHGS_StnTS(station=hgs_station, well=hgs_well, folder=hgs_folder, layers=None, #[16,17,18], 
                            start_date='1979-01-01', time_axis='datetime', resample='D', 
#                             end_date='1988-12-31', llastIncl=True,
                            end_date='1989-01-01', llastIncl=False,
                            basin=basin_name, WSC_station=WSC_station, basin_list=basin_list, lkgs=False,
                            conservation_authority='GRCA',
                            lskipNaN=True, lcheckComplete=True, varlist='default', scalefactors=1e-4,
                            PRD='', DOM=2, CLIM='clim_15', BC='AABC_', EXP='erai-g', name='{EXP:s} ({BASIN:s})')
    # N.B.: there is no record of actual calendar time in HGS, so periods are anchored through start_date/run_period
    # and print
    print(dataset)
    print('')
    print(dataset.name)
    print(dataset.prettyPrint(short=True))
    
    # view time axis
    print('')
    print(dataset.time)
    print(dataset.time[:])

    # some variable
    if 'discharge' in dataset:
        print('')
        print(dataset.discharge)
        print(dataset.discharge.plot)

#     # test climatology... currently only works with month
#     print('')
#     clim = dataset.climMean()
#     print(clim)


  elif test_mode == 'station_dataset':

    # load dataset
    lkgs = True
    dataset = loadHGS_StnTS(station=hgs_station, conservation_authority='GRCA', well=hgs_well, folder=hgs_folder, 
                            layers=None, z_layers='screen', z_axis='z', z_aggregation='max', correct_z=True,
                            start_date=1979, run_period=10, PRD='', DOM=2, CLIM='clim_15', BC='AABC_', 
                            basin=basin_name, WSC_station=WSC_station, basin_list=basin_list, lkgs=lkgs,
                            lskipNaN=True, lcheckComplete=True, varlist='default', scalefactors=1e-4,
                            EXP='erai-g', name='{EXP:s} ({BASIN:s})')
    # N.B.: there is no record of actual calendar time in HGS, so periods are anchored through start_date/run_period
    # and print
    print(dataset)
    print('')
    print(dataset.name)
    print(dataset.prettyPrint(short=True))
    
    # some common operations
    print('')
    clim = dataset.climMean()
    print(clim)
    if dataset.hasAxis('z'):
        print('')
        #print(dataset.x[:],dataset.y[:],)
        print(dataset.z[:])
        assert dataset.z.units == 'm', dataset.z
        print("Screen Interval: {}m - {}m".format(dataset.atts.z_b,dataset.atts.z_t))
    
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

    
  elif test_mode == 'station_ensemble':
    
    ens_name = '{EXPERIMENT:s}{PRDSTR:s}'
    ens_folder = '{ROOT_FOLDER:s}/GRW/grw2/{EXPERIMENT:s}{PRDSTR:s}_d{DOM:02d}/{CLIM:s}/hgs_run_v3_wrfpet'
    # actual ensemble definition
    ensemble_list = {'g-mean':('g-ctrl','g-ens-A','g-ens-B','g-ens-C'),
                     't-mean':('t-ctrl','t-ens-A','t-ens-B','t-ens-C')}

    # load an esemble of datasets
    ens = loadHGS_StnEns(ensemble=['g-mean','t-mean'], DOM=2, CLIM='AABC_clim_15', run_period=15,
                         period=[(1984,1994),(2050,2060),(2090,2100)], PRDSTR=['','-2050','-2100'],
                         station=hgs_station, well=hgs_well, conservation_authority='GRCA', 
                         name=ens_name, title=ens_name, basin=basin_name, layers=[1,2], # z_layers='screen',
                         WSC_station=WSC_station, basin_list=basin_list, folder=ens_folder,
                         lskipNaN=True, lcheckComplete=True, ensemble_list=ensemble_list, 
                         ens_name='HGS Ensemble', ens_title='HGS Ensemble based on WRF Ensemble',
                         outer_list=['ensemble',('period','PRDSTR')], lensemble=True)
    # load an esemble of datasets
    obs = loadHGS_StnEns(ensemble='Observations', basin=basin_name,
                         station=hgs_station, well=hgs_well, conservation_authority='GRCA', 
                         WSC_station=WSC_station, basin_list=basin_list, folder=None,
                         lskipNaN=True, lcheckComplete=True, varlist=None, obs_period=(1974,2015),
                         outer_list=None, lensemble=False)
    ens.insertMember(0,obs)
    # N.B.: all need to have unique names... whihc is a problem with obs...
    print(ens)
    print('\n')
    print(ens[-1])