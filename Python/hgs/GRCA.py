'''
Created on 2017-12-21

This module contains meta data and access functions for GRCA well data from xls spread sheets and
database files.

@author: Andre R. Erler, GPL v3
'''

# external imports
import os
import numpy as np
import numpy.ma as ma
# internal imports
from datasets.common import days_per_month, name_of_month, getRootFolder, loadObservations
from warnings import warn

## GRCA Meta-data

dataset_name = 'GRCA'
root_folder = getRootFolder(dataset_name=dataset_name) # get dataset root folder based on environment variables

# variable attributes and name (basically no alterations necessary...)
varatts = dict(h = dict(name='head', units='m', atts=dict(long_name='Pressure Head at Well')), # head in observation well
               # constants
               well_depth = dict(name='d', units='m', atts=dict(long_name='Well Depth')), # well depth from surface
               wel_piezom = dict(name='d_piezo', units='m', atts=dict(long_name='Depth of Piezometer')), # where the head is measured
               elva_groun = dict(name='zs', units='m', atts=dict(long_name='Elevation (M.S.L.)')), # surface elevation
               longitude  = dict(name='lon', units='deg E', atts=dict(long_name='Longitude')), # geographic longitude field
               latitude   = dict(name='lat', units='deg N', atts=dict(long_name='Latitude')), # geographic latitude field
               # axis 
               z = dict(name='z', units='m', atts=dict(long_name='Sampling Elevation (M.S.L.)')), # elevation where the measurement is taken
               time = dict(name='time', units='month', atts=dict(long_name='Month of the Year')),) # time coordinate
# N.B.: the time-series time offset is chose such that 1979 begins with the origin (time=0)
# list of variables to load
varlist = varatts.keys() # also includes coordinate fields    


## Functions that provide access to well-formatted PRISM NetCDF files

# pre-processed climatology files (varatts etc. should not be necessary)
avgfile = dataset_name.lower()+'{0:s}_clim{1:s}.nc' # formatted NetCDF file
avgfolder = root_folder + dataset_name.lower()+'avg/' # prefix
# function to load these files...
def loadGRCA(name=dataset_name, period=None, grid=None, resolution=None, varlist=None, 
              varatts=None, folder=avgfolder, filelist=None, lautoregrid=True):
  ''' Get the pre-processed monthly PRISM climatology as a DatasetNetCDF. '''
  # only the climatology is available
  if period is not None: 
    warn('Only the full climatology is currently available: setting \'period\' to None.')
    period = None
  # load standardized climatology dataset with PRISM-specific parameters  
  dataset = loadObservations(name=name, folder=folder, period=period, grid=grid, station=None, 
                             varlist=varlist, varatts=varatts, filepattern=avgfile, filelist=filelist, 
                             lautoregrid=lautoregrid, mode='climatology')
  # return formatted dataset
  return dataset

def loadGRCA_Stn(name=dataset_name, period=None, station=None, resolution=None, varlist=None, 
              varatts=None, folder=avgfolder, filelist=None):
  ''' Get the pre-processed monthly PRISM climatology at station locations as a DatasetNetCDF. '''
  # only the climatology is available
  if period is not None: 
    warn('Only the full climatology is currently available: setting \'period\' to None.')
    period = None
  # load standardized climatology dataset with PRISM-specific parameters  
  dataset = loadObservations(name=name, folder=folder, grid=None, station=station, shape=None, 
                             varlist=varlist, varatts=varatts, filepattern=avgfile, filelist=filelist, 
                             lautoregrid=False, period=period, mode='climatology')
  # return formatted dataset
  return dataset


## Dataset API

dataset_name # dataset name
root_folder # root folder of the dataset
ts_file_pattern = None
clim_file_pattern = avgfile # filename pattern
data_folder = avgfolder # folder for user data
grid_def = None # no special name, since there is only one...
LTM_grids = [] # grids that have long-term mean data 
TS_grids = [] # grids that have time-series data
grid_res = {} # approximate resolution in degrees at 45 degrees latitude
default_grid = None
# functions to access specific datasets
loadLongTermMean = None # climatology provided by publisher
loadTimeSeries = None # time-series data
loadClimatology = loadGRCA # pre-processed, standardized climatology
loadStationClimatology = loadGRCA_Stn
loadShapeClimatology = None

## Functions that handle access to GRCA XLS/DBF files

grca_folder = os.path.join(root_folder,'GRCA Transient Waterlevels')

def getWellName(name):
  ''' helper function to break down well names into ID and number '''
  well_no = 1 # default, unless specified
  if isinstance(name,basestring):
      well = name.upper()
      if '-' in well: 
          well_id,well_no = well.split('-')
          well_no = int(well_no)
      else:
          well_id = well
      if 'W' == well_id[0]:
          well_id = int(well_id[1:])
      else: 
          well_id = int(well_id)
  else: 
      well_id = int(well)
  # return components
  return well_id, well_no

# loads data from original XLS files and returns Pandas data_frame
def loadXLS(well, filename='', folder=grca_folder):
  # local imports
  from numpy.ma import zeros
  # return array
  return data

# function to meta data from database file
def loadMetadata(well, filename='metadata.dbf', format='W{WELL_ID:07d}-{WELL_NO:1d}', folder=grca_folder):
  # clean up well name
  well_id, well_no = getWellName(well)
  well = format.format(WELL_ID=well_id, WELL_NO=well_no)
  # open database and get relevant entry
  #from simpledbf import Dbf5
  from dbfread import DBF
  table = DBF(os.path.join(grca_folder,filename))
  for record in table:
      if record['PGMN_WELL'] == well: 
          meta = record.copy()
  # parse screen information
  screen_type,screen_depth= meta['SCREEN_HOL'].split(':')
  meta['Screen'] = screen_type.title()
  screen_hilo = []
  for hilo in screen_depth.split('-'):
      if hilo[-1] != 'M': raise ValueError(hilo)
      screen_hilo.append(float(hilo[:-1]))
  assert len(screen_hilo) == 2, screen_hilo
  meta['Screen_top'] = screen_hilo[0]
  meta['Screen_bottom'] = screen_hilo[1]
  meta['Screen_depth'] = ( screen_hilo[0] + screen_hilo[1] ) / 2.
  meta['z'] = meta['ELVA_GROUN'] - meta['Screen_depth']
  meta['z_t'] = meta['ELVA_GROUN'] - meta['Screen_top']
  meta['z_b'] = meta['ELVA_GROUN'] - meta['Screen_bottom']
  # return coordinate arrays (in degree)
  return meta

if __name__ == '__main__':
    
#   mode = 'test_climatology'
#   mode = 'test_timeseries'
#   mode = 'convert_XLS'
  mode = 'test_load_XLS'
  
  # do some tests
  if mode == 'test_climatology':
    
    # load point climatology
    print('')
    dataset = loadGRCA_Stn(station='')
    print(dataset)
    print('')
    print(dataset.time)
    print(dataset.time.coord)


  ## test load function for XLS files
  elif mode == 'test_load_XLS': 

    meta = loadMetadata(well='W178', )
    # inspect dictionary
    for item in meta.items(): print(item)

  ## convert XLS files to NetCDF
  elif mode == 'convert_XLS': 
    
    raise NotImplementedError