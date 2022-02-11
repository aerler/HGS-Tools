'''
Created on Feb. 26, 2019

Utility functions to create NetCDF datasets and interoperate with GDAL/rasterio and xarray/dask.
This module is an extended Python 3 port of the GeoPy utils.nctools module, and should be Python 2
and Python 3 compatible.

@author: Andre R. Erler, GPL v3
'''

# external imports
import os
import os.path as osp
import netCDF4 as nc # netCDF4-python module: Dataset is probably all we need
# N.B.: most of the functions simply take a Dataset object as the first argument
from six import string_types # for testing string in Python 2 and 3
import collections as col
import numpy as np
import numpy.ma as ma
from warnings import warn

## definitions

geospatial_netcdf_version = 1.0

# NC4 compression options
zlib_default = dict(zlib=True, complevel=1, shuffle=True) # my own default compression settings

# data error class
class NCDataError(Exception):
  ''' Exceptions related to data passed to NetCDF datasets. '''
  pass

# axis error class
class NCAxisError(Exception):
  ''' Exceptions related to axis/dimensions in NetCDF datasets. '''
  pass

# prevent failure if xarray is not installed
try:
    from xarray import DataArray
except:
    warn("Import of 'xarray' failed: cannot handle 'DataArray'.")
    DataArray = NotImplemented.__class__ # set to __class__, so that isinstance works (and usually fails)


## helper functions

# def autoChunk(shape, min_tiles=8, tile_size=256, total_size=524288):
# def autoChunk(shape, min_tiles=8, tile_size=128, total_size=65536):
def autoChunk(shape, min_tiles=16, tile_size=64, total_size=32768):
    ''' function to find good chunks based on simple heuristics
        This could be revised using the script linked here:
        https://www.unidata.ucar.edu/blogs/developer/en/entry/chunking_data_choosing_shapes
    '''
    ndim = len(shape)
    if ndim == 2: shp_2d = shape
    elif ndim == 3: shp_2d = shape[1:] # assuming time is the outermost...
    else:
        raise NotImplementedError(shape)
    chunks = (1,1); cs = tile_size*2
    # find smaller decomposition, so we have at least 8 tiles
    while np.prod(chunks) < min_tiles:
        cs /= 2
        chunks = [int(np.ceil(s/np.ceil(s/cs))) for s in shp_2d]
    if ndim == 3:
        chunks = [int(total_size/np.prod(chunks))] + chunks # assign remaining size to time dimension
    elif ndim > 3:
        raise NotImplementedError(shape)
    return chunks


def checkFillValue(fillValue, dtype):
    ''' check a fill value and either return a valid value or raise an exception '''

    if dtype.kind != 'S' and isinstance(fillValue, string_types): # string variable
        fillValue = None # invalid value
        warn("Removed malformed fill_value '{:s}'.")

    if fillValue is not None:
        if isinstance(fillValue, np.ndarray): fillValue = np.asscalar(fillValue)
        fillValue = dtype.type(fillValue) # transform into appropriate numpy scalar
        if not np.issubdtype(fillValue,dtype):
          raise TypeError(fillValue) # convert to Numpy type

    # return valid fill value
    return fillValue


def coerceAtts(atts):
    ''' Convert an attribute dictionary to a NetCDF compatible format (creates a new copy). '''

    if not isinstance(atts,dict): raise TypeError(atts)
    ncatts = dict()

    # loop over items
    for key,value in atts.items():
        if key in ('missing_value','fillValue','_FillValue'): pass
        # N.B.: these are special attributes that the NetCDF module will try to read
        elif isinstance(key,string_types) and key[0] == '_' : pass # skip (internal attributes)
        elif value is None: pass # skip (invalid value / not assigned)
        elif not isinstance(value,(string_types,np.ndarray,np.inexact,float,np.integer,int)):
            if 'name' in dir(value):
                ncatts[key] = value.name # mostly for datasets and variables
            elif isinstance(value,col.Iterable):
                if len(value) == 0: ncatts[key] = '' # empty attribute
                elif all([isinstance(val,(int,np.integer,float,np.inexact)) for val in value]):
                    # N.B.: int & float are not part of their numpy equivalents
                    ncatts[key] = np.array(value)
                else:
                    l = '(' # fake list representation
                    for elt in value[0:-1]: l += '{0:s}, '.format(str(elt))
                    l += '{0:s})'.format(str(value[-1]))
                    ncatts[key] = l
            else: ncatts[key] = str(value)
        else: ncatts[key] = value

    # return reformatted attribute dict
    return ncatts


def getNCAtts(src):
    ''' return NetCDF attributes as dictionary '''
    atts = dict()
    # loop over attribute list
    for att in src.ncattrs():
        # extract each value individually (no other option)
        atts[att] = src.getncattr(att)
    return atts

def setNCAtts(tgt, atts, lcoerce=True):
    ''' set NetCDF attributes based on dictionary, with some general checking '''
    # coerce type of attributes
    if lcoerce:
        atts = coerceAtts(atts)
    # loop over dictionary
    for key,value in atts.items():
        # set attributed individually (no other option)
        if isinstance(value,string_types): tgt.setncattr_string(key,value)
        else: tgt.setncattr(key,value)
    return atts


## generic functions to create netCDF variables

def add_coord(dst, name, data=None, length=None, atts=None, dtype=None, zlib=True, fillValue=None, **kwargs):
    ''' Function to add a Coordinate Variable to a NetCDF Dataset; returns the Variable reference. '''

    # check input
    if length is None:
        # that means unlimited record dimensions (can only have one of these)
        dst.createDimension(name, size=None) # need to create here, or add_var will assign a size
    if isinstance(length,(int,np.integer)):
        length=(length,)
    if data is not None and data.ndim == 1:
        if length is None: length = data.shape
    else:
        raise NCDataError(data)

    # basically a simplified interface for add_var
    coord = add_var(dst, name, (name,), data=data, shape=length, atts=atts, dtype=dtype,
                    zlib=zlib, fillValue=fillValue, **kwargs)

    # return handle to coordinate variable
    return coord


def add_var(dst, name, dims, data=None, shape=None, atts=None, dtype=None, zlib=True, smart_chunks=True, fillValue=None, **kwargs):
    ''' function to add a Variable to a NetCDF Dataset; returns the Variable reference;
        all remaining kwargs are passed on to dst.createVariable() '''

    # use data array to infer dimensions and data type
    if data is not None:
        if not isinstance(data,np.ndarray): raise TypeError(data)
        if len(dims) != data.ndim:
            raise NCDataError("Number of dimensions in '{:s}' does not match data array.".format(name,))
        if shape:
            if shape != data.shape:
                raise NCDataError("Shape of '{:s}' does not match data array.".format(name,))
        else: shape = data.shape
        # get dtype
        if dtype:
            if dtype != data.dtype: data = data.astype(dtype)
        else: dtype = data.dtype
    if dtype is None: raise NCDataError("Cannot construct a NetCDF Variable without a data array or an abstract data type.")
    dtype = np.dtype(dtype) # use numpy types
    if dtype is np.dtype('bool_'): dtype = np.dtype('i1') # cast numpy bools as 8-bit integers
    if dtype.kind == 'S': dtype = np.dtype(str)
    # N.B.: the dtype 'S' causes a TypeError, because NetCDF4 only allows 'S1', but Python str can be used instead

    # check/create dimensions
    if shape is None: shape = [None,]*len(dims)
    else: shape = list(shape)
    if len(shape) != len(dims):
        raise NCAxisError(shape)
    for i,dim in zip(range(len(dims)),dims):
        if dim in dst.dimensions:
          if shape[i] is None:
              shape[i] = len(dst.dimensions[dim])
          else:
              if shape[i] != len(dst.dimensions[dim]) and len(dst.dimensions[dim]) > 0:
                  raise NCAxisError('Size of dimension {:s} does not match records! {:d} != {:d}'.format(dim,shape[i],len(dst.dimensions[dim])))
        else:
            dst.createDimension(dim, size=shape[i])
    dims = tuple(dims); shape = tuple(shape)

    # figure out parameters for variable
    varargs = dict() # arguments to be passed to createVariable
    if isinstance(zlib,dict): varargs.update(zlib)
    elif zlib: varargs.update(zlib_default)
    varargs.update(kwargs)
    if fillValue is None:
        if atts and '_FillValue' in atts: fillValue = atts['_FillValue'] # will be removed later
        elif atts and 'missing_value' in atts: fillValue = atts['missing_value']
        elif data is not None and isinstance(data,ma.MaskedArray): # defaults values for numpy masked arrays
          fillValue = ma.default_fill_value(dtype)
        else: pass # if it is not a masked array and no missing value information was passed, don't assign fillValue
    else:
        if data is not None and isinstance(data,ma.MaskedArray): data._fill_value = fillValue
    # make sure fillValue is OK (there have been problems...)
    fillValue = checkFillValue(fillValue, dtype)

    # set chunks based on reasonable division of domain
    if smart_chunks and len(dims) > 1 and 'chunksizes' not in varargs:
        varargs['chunksizes'] = autoChunk(shape)

    # create netcdf variable
    var = dst.createVariable(name, dtype, dims, fill_value=fillValue, **varargs)
    # add attributes
    if atts:
        var.setncatts(coerceAtts(atts))
        if fillValue is not None:
            var.setncattr('missing_value',fillValue) # coerceAtts removes this one... but I prefer it to _FillValue
    # assign coordinate data if given
    if data is not None: var[:] = data

    # return var reference
    return var


## copy functions

def copy_ncatts(dst, src, prefix='', incl_=True):
    ''' copy attributes from a variable or dataset to another '''

    # loop over attributes
    for att in src.ncattrs():
      if att in ('missing_value','fillValue','_FillValue'): pass
      elif att[0] != '_' or incl_: # these seem to cause problems
          value = src.getncattr(att)
          if isinstance(value,string_types): dst.setncattr_string(prefix+att,value)
          else: dst.setncattr(prefix+att,value)

    # return modified dataset
    return dst


def copy_vars(dst, src, varlist=None, namemap=None, dimmap=None, remove_dims=None, copy_data=True,
              copy_atts=True, zlib=True, prefix='', incl_=True, fillValue=None, **kwargs):
    ''' copy variables from one dataset to another '''

    # prefix is passed to copy_ncatts, the remaining kwargs are passed to dst.createVariable()
    if varlist is None: varlist = src.variables.keys() # just copy all
    #if dimmap: midmap = dict(zip(dimmap.values(),dimmap.keys())) # reverse mapping
    varargs = dict() # arguments to be passed to createVariable
    if zlib: varargs.update(zlib_default)
    varargs.update(kwargs)
    dtype = varargs.pop('dtype', None)

    # loop over variable list
    for name in varlist:

        if namemap and (name in namemap.keys()): rav = src.variables[namemap[name]] # apply name mapping
        else: rav = src.variables[name]
        dims = [] # figure out dimension list
        for dim in rav.dimensions:
            #if dimmap and dim in midmap: dim = midmap[dim] # apply name mapping (in reverse)
            if dimmap and dim in dimmap: dim = dimmap[dim] # apply name mapping
            if not (remove_dims and dim in remove_dims): dims.append(dim)

        # create new variable
        dtype = dtype or rav.dtype
        if '_FillValue' in rav.ncattrs(): fillValue = rav.getncattr('_FillValue')
        var = dst.createVariable(name, dtype, dims, fill_value=fillValue, **varargs)
        if copy_atts: copy_ncatts(var, rav, prefix=prefix, incl_=incl_) # copy attributes, if desired (default)
        if copy_data: var[:] = rav[:] # copy actual data, if desired (default)

    # return modified dataset
    return dst


def copy_dims(dst, src, dimlist=None, namemap=None, copy_coords=True, **kwargs):
    ''' copy dimensions and coordinate variables from one dataset to another '''

    # all remaining kwargs are passed on to dst.createVariable()
    if dimlist is None: dimlist = src.dimensions.keys() # just copy all
    if namemap is None: namemap = dict() # a dummy - assigning pointers in argument list is dangerous!

    # loop over dimensions
    for name in dimlist:
        mid = src.dimensions[namemap.get(name,name)]
        # create actual dimensions
        dst.createDimension(name, size=len(mid))

    # create coordinate variable
    if copy_coords:
        remove_dims = [dim for dim in src.dimensions.keys() if dim not in dimlist] # remove_dims=remove_dims
        dimlist = [dim for dim in dimlist if dim in src.variables] # only the ones that have coordinates
        copy_vars(dst, src, varlist=dimlist, namemap=namemap, dimmap=namemap, remove_dims=remove_dims, **kwargs)

    # return modified dataset
    return dst


## function to create a geo-referenced NetCDF dataset

def interprete_time_units(units):
    ''' homogenize various time formats... '''
    # identify unit names and default formats
    if units == 'M' or units.lower() in ('month','months'):
        unit_name= 'Months'; units = 'month'; dt_code = 'M'; date_fmt = '%Y-%m'
    elif units.lower() in ('days','day','d'):
        unit_name= 'Days'; units = 'day'; dt_code = 'D'; date_fmt = '%Y-%m-%d'
    elif units.lower() in ('hour','hours','h'):
        unit_name= 'Hours'; units = 'h'; dt_code = 'h'; date_fmt = '%Y-%m-%d_%H'
    elif units.lower() in ('min','minute','minutes','m'):
        unit_name= 'Minutes'; units = 'min'; dt_code = 'm'; date_fmt = '%Y-%m-%d_%H:%M:%S'
    elif units.lower() in ('sec','second','seconds','s'):
        unit_name= 'Seconds'; units = 's'; dt_code = 's'; date_fmt = '%Y-%m-%d_%H:%M:%S'
    else:
        raise ValueError(units)
    # return various translations
    return unit_name,units,dt_code,date_fmt

def addTimeStamps(dataset, time='time', units=None, atts=None):
    ''' add human-readable time-stamps to NetCDF dataset based on datetime64 axis'''
    # settings
    if atts is None: atts = dict(name=time+'_stamp', units='', long_name='Time Stamp')
    # load existing time coordinate
    nc_unit_str = dataset.variables[time].getncattr('units')
    dt_coord = nc.num2date(dataset.variables[time][:], units=nc_unit_str) # load as datetime64
    # figure out target format
    if units is None: units = nc_unit_str.split()[0]
    unit_name,units,ts_dt_code,date_fmt = interprete_time_units(units); del unit_name, date_fmt
    tsdata = np.stack([str(t) for t in dt_coord.astype('datetime64[{}]'.format(ts_dt_code))], axis=0)
    # create NC variable
    tsnc = add_var(dataset, atts['name'], dims=('time',), data=tsdata, atts=atts, zlib=True, ) # daily time-stamp
    # add time-stamps
    return tsnc

def add_time_coord(dst, data, name=None, units='D', atts=None, ts_atts=None, date_fmt=None,
                   dtype=None,ltimestamp=True, zlib=True, **kwargs):
    ''' add a time coordinate and a human-readable time-stamp variables based on a datetime64 array '''

    # make time stamps
    unit_name,units,dt_code,date_format = interprete_time_units(units)
    if date_fmt is None: date_fmt = date_format
    if ltimestamp:
        if date_fmt:
            import pandas as pd
            ts_data = np.asanyarray([pd.to_datetime(date).strftime(date_fmt) for date in data])
        else:
            ts_data = np.datetime_as_string(data, unit=dt_code) # very convenient and seems to do what I want
        start_date = ts_data[0]
    else:
        if date_fmt:
            import pandas as pd
            start_date = pd.to_datetime(data[0]).strftime(date_fmt)
        else:
            start_date = np.datetime_as_string(data[0], unit=dt_code)

    # time is the outer-most (record) dimension
    if name is None and atts is None: name = 'time'
    elif name is None: name = atts.get('name','time')
    if atts is None: atts = dict(name=name)
    atts['units'] = '{} since {}'.format(unit_name.lower(), start_date)
    atts['geospatial_netcdf_version'] = geospatial_netcdf_version
    coord = (data - data[0]) / np.timedelta64(1, dt_code)
    if dtype is None: dtype = np.dtype('i4') # 32 bit int
    add_coord(dst, name, data=coord, length=None, atts=atts, dtype=dtype, zlib=zlib, **kwargs) # record dim

    # also add a time-stamp variable
    if ltimestamp:
        if ts_atts is None:
            ts_atts = atts.copy() # build based on time atts
            ts_atts['name'] = name+'_stamp'
            ts_atts['long_name'] = "Human-readable Time Stamp"
            ts_atts['units'] = ''
            ts_atts['geospatial_netcdf_version'] = geospatial_netcdf_version
        add_var(dst, name+'_stamp', dims=(name,), data=ts_data, shape=(len(data),),
                atts=ts_atts, dtype=np.dtype('U'), zlib=zlib, **kwargs) # time-stamp

    # return dataset handle
    return dst


# default attributes for month meta data variables
default_mon_name_atts = dict(name='name_of_month', units='', long_name='Name of the Month',
                             geospatial_netcdf_versio=geospatial_netcdf_version)
default_mon_len_atts = dict(name='length_of_month', units='days',long_name='Length of Month',
                             geospatial_netcdf_versio=geospatial_netcdf_version)

def addNameLengthMonth(ds, time_dim='time', mon_name_atts=default_mon_name_atts, mon_len_atts=default_mon_len_atts):
    ''' add name and length of month to monthly normals '''
    from datasets.common import name_of_month, days_per_month
    assert len(ds.dimensions[time_dim]) == 12, ds.dimensions[time_dim]
    # name of month
    varnc = add_var(ds, mon_name_atts['name'], dims=('time',), data=None, shape=(None,),
                    atts=mon_name_atts, dtype=str, zlib=True, fillValue=None) # string variable
    varnc[:] = np.stack(name_of_month, axis=0)
    # length of month
    varnc = add_var(ds, mon_len_atts['name'], dims=('time',), data=None, shape=(None,),
                    atts=mon_len_atts, dtype=np.dtype('float32'), zlib=True, fillValue=None,) # float variable
    varnc[:] = np.stack(days_per_month, axis=0)
    # sync and return dataset with new variables
    ds.sync()
    return ds


def createGeoReference(ds, crs=None, geotrans=None, size=None, xlon=None, ylat=None, griddef=None,
                       varatts=None, zlib=True, lcoords=True, loverwrite=False):
    ''' function to add georeferencing information based on rasterio/GDAL info to existing netCDF4 Dataset;
        this can include adding new axes, but does not modify existing variables '''
    from geospatial.rasterio_tools import genCRS, constructCoords

    # interpret GridDefinition object (defined in geodata.gdal, but not imported to avoid dependencies)
    if griddef:
        crs = griddef.projection.ExportToProj4() # for some reason, rasterio does not understand WKT...
        geotrans = griddef.geotransform
        size = griddef.size
        xlon = griddef.xlon.name
        ylat = griddef.ylat.name

    # determine type of coordinate reference system
    if crs is None:
        raise NotImplementedError
    crs = genCRS(crs)
    lproj = crs.is_projected
    # construct coordinate axes
    if not xlon and not ylat:
        xlon,ylat = ('x','y') if lproj else ('lon','lat')
    assert xlon and ylat, (xlon,ylat)

    # construct coordinate variables
    if lcoords:
      xlon_coord,ylat_coord = constructCoords(geotrans, size, dtype=nc_coord_dtype)
      # add coordinate variables
      if varatts is None: varatts = default_varatts
      # latitude (intermediate/regular dimension)
      yatts = varatts[ylat]; yname = yatts['name']
      if yname not in ds.variables or loverwrite:
          dtype = yatts.get('dtype',nc_coord_dtype)
          add_coord(ds, yname, data=ylat_coord, length=size[1], atts=yatts, dtype=dtype, zlib=zlib)
      # longitude is typically the inner-most dimension (continuous)
      xatts = varatts[xlon]; xname = xatts['name']
      if xname not in ds.variables or loverwrite:
          dtype = xatts.get('dtype',nc_coord_dtype)
          add_coord(ds, xname, data=xlon_coord, length=size[0], atts=xatts, dtype=dtype, zlib=zlib)

    ## save attributes, including georeferencing information
    if isinstance(geotrans,tuple) and hasattr(geotrans,'to_gdal'):
        # probably rio.transform.Affine
        dx = abs(geotrans.a); dy = abs(geotrans.e) # I'm never sure if I got this right...
        gdal_geotrans = geotrans.to_gdal()
        assert (abs(gdal_geotrans[1]),abs(gdal_geotrans[5])) == (dx,dy), gdal_geotrans
    else:
        assert len(geotrans)==6, geotrans
        dx,dy = abs(geotrans[1]), abs(geotrans[5])
    ds.setncattr('dx',dx); ds.setncattr('dy',dy)
    ds.setncattr_string('xlon',xlon)
    ds.setncattr_string('ylat',ylat)
    ds.setncattr_string('proj4',crs.to_proj4())
    wkt = crs.to_wkt() if hasattr(crs,'to_wkt') else 'n/a'
    ds.setncattr_string('WKT',wkt if wkt else 'n/a')
    epsg = crs.to_epsg() if hasattr(crs,'to_epsg') else 'n/a'
    ds.setncattr('EPSG',epsg if epsg is not None else 'n/a')
    ds.setncattr('is_projected',int(lproj))
    ds.setncattr('geospatial_netcdf_version',geospatial_netcdf_version)
    # return modified dataset
    return ds


nc_coord_dtype  = np.dtype('<f8') # little-endian 64-bit float, otherwise errors in geotransform accumulate
netcdf_dtype    = np.dtype('<f4') # little-endian 32-bit float, to save space...
default_varatts = dict( # attributed for coordinate variables
                       x  = dict(name='x', units='m', long_name='Northing'), # geographic northing
                       y  = dict(name='y', units='m', long_name='Easting'), # geographic easting
                       lon  = dict(name='lon', units='deg E', long_name='Longitude'), # geographic longitude
                       lat  = dict(name='lat', units='deg N', long_name='Latitude'), # geographic latitude
                       time = dict(name='time', units='day', long_name='Days'), # time coordinate
                       time_stamp = dict(name='time_stamp', units='', long_name='Human-readable Time Stamp'),) # readable time stamp (string)

def createGeoNetCDF(filename, atts=None, folder=None, xlon=None, ylat=None, time=None,
                    time_args=None, varatts=None, crs=None, geotrans=None, size=None, griddef=None,
                    nc_format='NETCDF4', zlib=True, loverwrite=True):
    ''' create a NetCDF-4 file, create dimensions and geographic coordinates, and set attributes '''

    # create Dataset/file
    if folder:
        if not osp.exists(folder):
            os.makedirs(folder)
        filepath = osp.join(folder,filename)
    else:
        filepath = filename
    if osp.exists(filepath) and not loverwrite:
        ds = nc.Dataset(filepath, mode='a', format=nc_format, clobber=False)
    else:
        ds = nc.Dataset(filepath, mode='w', format=nc_format, clobber=True)
    # add attributes
    if atts is None: atts = dict()
    setNCAtts(ds, atts, lcoerce=True)
    # default varatts
    if varatts is None: varatts = default_varatts

    ## construct time dimension and time stamp variable from xarray/datetime64
    if time is not None:

        # infer datetime64 array
        if isinstance(time,DataArray):
            time_coord = time.data
            time_atts = time.attrs
            time_axis = time.name
        elif isinstance(time,np.ndarray):
            if np.issubdtype(time.dtype, np.datetime64):
                time_coord = time
                time_axis = 'time'
                time_atts = varatts[time_axis]
            else:
                raise NotImplementedError(time)
        else:
            raise NotImplementedError(time)
        # add time axis and time stamps based on datetime64
        kwargs = dict(name=time_axis, atts=time_atts, ltimestamp=True, zlib=zlib)
        if time_args: kwargs.update(time_args)
        if kwargs['name'] not in ds.variables or loverwrite:
            add_time_coord(ds, time_coord, **kwargs)

    ## construct geographic coordinate axes from xarray or GDAL/rasterio georeference
    ds = createGeoReference(ds, crs=crs, geotrans=geotrans, size=size, xlon=xlon, ylat=ylat, griddef=griddef,
                            varatts=varatts, zlib=zlib, lcoords=True, loverwrite=loverwrite)

    # return dataset object
    return ds


# abuse for testing
if __name__ == '__main__':

    # define a grid (son1 for testing)
    crs = "+proj=utm +zone=17 +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
    size = (118,82) # lower resolution 5 km grid
    geotrans = (320920.,5.e3,0,4624073.,0,5.e3) # 5 km

    # create georeferenced NetCDF file for testing
    filepath = 'D:/Data/test/geo_test.nc'
    print(filepath)
    atts = dict(note='a geo-referenced test dataset')
    time = np.arange('2014-09-28','2014-10-21',dtype='datetime64[D]')
    print(time)
    ds = createGeoNetCDF(filepath, atts=atts, time=time, crs=crs, geotrans=geotrans, size=size,
                         time_args=dict(date_fmt='%Y%m%dT%H'), loverwrite=True)
    ds.close()

    # open file with xarray and print
    import xarray as xr

    print('')
    assert osp.exists(filepath), filepath
    ds = xr.open_dataset(filepath, decode_cf=True, mask_and_scale=True, decode_times=True, autoclose=True,
                         decode_coords=True, engine='netcdf4', )
    print(ds)

    print('')
    print(ds.variables['time'])

    print('')
    print(ds.variables['time_stamp'])
