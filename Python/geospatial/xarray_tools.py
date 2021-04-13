'''
Created on Feb. 23, 2019

Utility functions to extract data from xarray Dataset or DataArray classes.

@author: Andre R. Erler, GPL v3
'''

from warnings import warn
from datetime import datetime
import os
import numpy as np
import xarray as xr
import netCDF4 as nc
from dask.diagnostics import ProgressBar

# internal imports
from geospatial.netcdf_tools import getNCAtts, geospatial_netcdf_version, zlib_default # this import should be fine

## an important option I am relying on!
xr.set_options(keep_attrs=True)

# names of valid geographic/projected coordinates
default_x_coords = dict(geo=('lon','long','longitude',), proj=('x','easting','west_east') )
default_y_coords = dict(geo=('lat','latitude',),         proj=('y','northing','south_north'))
default_lon_coords = default_x_coords['geo']; default_lat_coords = default_y_coords['geo']



## helper functions

def getAtts(xvar, lraise=True):
    ''' return dictionary of attributed from netCDF4 or xarray '''
    if isinstance(xvar,(xr.DataArray,xr.Variable,xr.Dataset)):
        atts = xvar.attrs.copy()
    elif isinstance(xvar,(nc.Variable,nc.Dataset)):
        atts = getNCAtts(xvar)
    elif lraise: 
        raise TypeError(xvar)
    return atts

## functions to interface with rasterio

def getGeoDims(xvar, x_coords=None, y_coords=None, lraise=True):
    ''' helper function to identify geographic/projected dimensions by name ''' 
    if x_coords is None: x_coords = default_x_coords
    if y_coords is None: y_coords = default_y_coords
    
    xlon,ylat = None,None # return None, if nothing is found    

    if isinstance(xvar,(xr.DataArray,xr.Dataset)):
        # test geographic grid and projected grids separately
        for coord_type in x_coords.keys():
            for name in xvar.dims.keys() if isinstance(xvar,xr.Dataset) else xvar.dims:
                if name.lower() in x_coords[coord_type]: 
                    xlon = name; break
            for name in xvar.dims.keys() if isinstance(xvar,xr.Dataset) else xvar.dims:
                if name.lower() in y_coords[coord_type]: 
                    ylat = name; break
            if xlon is not None and ylat is not None: break
            else: xlon,ylat = None,None
    elif isinstance(xvar,(nc.Dataset,nc.Variable)):
        # test geographic grid and projected grids separately
        for coord_type in x_coords.keys():
            for name in xvar.dimensions:
                if name.lower() in x_coords[coord_type]: 
                    xlon = name; break
            for name in xvar.dimensions:
                if name.lower() in y_coords[coord_type]: 
                    ylat = name; break
            if xlon is not None and ylat is not None: break
            else: xlon,ylat = None,None      
    elif lraise: # optionally check input
        raise TypeError("Can only infer coordinates from xarray or netCDF4 - not from {}".format(xvar.__class__))
    else:
        pass # return None,None
    
    return xlon,ylat

def getGeoCoords(xvar, x_coords=None, y_coords=None, lraise=True, lvars=True):
    '''  helper function to extract geographic/projected coordinates from xarray'''
    
    # find dim names
    xlon_dim,ylat_dim = getGeoDims(xvar, x_coords=x_coords, y_coords=y_coords, lraise=lraise)
    
    # find coordinates
    if isinstance(xvar,(xr.DataArray,xr.Dataset)):
        if xlon_dim in xvar.coords:
            xlon = xvar.coords[xlon_dim] if lvars else xlon_dim
        else: xlon = None
        if ylat_dim in xvar.coords:
            ylat = xvar.coords[ylat_dim] if lvars else ylat_dim 
        else: ylat = None
    elif isinstance(xvar,nc.Variable) and lraise:
        raise TypeError("Cannot infer coordinates from netCDF4 Variable - only Dataset!")
    elif isinstance(xvar,nc.Dataset):
        if xlon_dim in xvar.variables:
            xlon = xvar.variables[xlon_dim] if lvars else xlon_dim
        else: xlon = None
        if ylat_dim in xvar.variables:
            ylat = xvar.variables[ylat_dim] if lvars else ylat_dim 
        else: ylat = None
        
    # optionally raise error if no coordinates are found, otherwise just return None
    if lraise and (xlon is None or ylat is None):
        raise ValueError("No valid pair of geographic coodinates found:\n {}".format(xvar.dims))
      
    # return a valid pair of geographic or projected coordinate axis
    return xlon,ylat

  
def isGeoVar(xvar, x_coords=None, y_coords=None, lraise=True):
    ''' helper function to identify variables that have geospatial coordinates (geographic or 
        projected), based on xarray or netCDF4 dimension names '''
    
    if x_coords is None: x_coords = default_x_coords
    if y_coords is None: y_coords = default_y_coords

    if isinstance(xvar,(xr.DataArray,xr.Dataset)):
        dims = xvar.coords.keys()
    elif isinstance(xvar,(nc.Dataset,nc.Variable)):
        dims = xvar.dimensions
    elif lraise:
        raise TypeError("Can only infer coordinate system from xarray or netCDF4 - not from {}".format(xvar.__class__))
    else:
        return None # evaluates as False, but allows checking
        
    # test geographic grid and projected grids separately
    for coord_type in x_coords.keys():
        xlon,ylat = False,False
        for name in dims:
            if name.lower() in x_coords[coord_type]: 
                xlon = True; break
        for name in dims:
            if name.lower() in y_coords[coord_type]: 
                ylat = True; break
        if xlon and ylat: break
    
    # if it has a valid pair of geographic or projected coordinate axis
    return ( xlon and ylat )

  
def isGeoCRS(xvar, lat_coords=None, lon_coords=None, lraise=True):
    ''' helper function to determine if we have a simple geographic lat/lon CRS (based on xarray dimension names) '''
    lat,lon = False,False
    
    if lon_coords is None: lon_coords = default_x_coords['geo']
    if lat_coords is None: lat_coords = default_y_coords['geo']
    
    if isinstance(xvar,(xr.DataArray,xr.Dataset)):
        dims = xvar.coords.keys()
    elif isinstance(xvar,(nc.Dataset,nc.Variable)):
        dims = xvar.dimensions
    elif lraise:
        raise TypeError("Can only infer coordinate system from xarray or netCDF4- not from {}".format(xvar.__class__))
    else:
        return None # evaluates as False, but allows checking
      
    # check dimension names
    for name in dims:
        if name.lower() in lon_coords: 
            lon = True; break
    for name in dims:
        if name.lower() in lat_coords: 
            lat = True; break
    
    # it is a geographic coordinate system if both, lat & lon are present
    return ( lat and lon )


def getTransform(xvar=None, x=None, y=None, lcheck=True):
    ''' generate an affine transformation from xarray coordinate axes '''
    from rasterio.transform import Affine # to generate Affine transform
    
    if isinstance(xvar,(xr.DataArray,xr.Dataset,nc.Dataset)):
        x,y = getGeoCoords(xvar, lraise=True)
    elif xvar is None and isinstance(x,(xr.DataArray,nc.Variable)) and isinstance(y,(xr.DataArray,nc.Variable)):
        pass # x and y axes are supplied directly
    elif xvar:
        raise TypeError('Can only infer GeoTransform from xarray Dataset or DataArray or netCDF4 Dataset\n - not from {}.'.format(xvar))
    
    # check X-axis
    if isinstance(x,xr.DataArray): x = x.data
    elif isinstance(x,nc.Variable): x = x[:]
    if not isinstance(x,np.ndarray): 
        raise TypeError(x)
    diff_x = np.diff(x); dx = diff_x.min()
    if lcheck and not np.isclose(dx, diff_x.max(), rtol=1.e-2): 
        raise ValueError("X-axis is not regular: {} - {}".format(dx, diff_x.max()))
    
    # check Y-axis
    if isinstance(y,xr.DataArray): y = y.data
    elif isinstance(y,nc.Variable): y = y[:]
    if not isinstance(y,np.ndarray): 
        raise TypeError(y)
    diff_y = np.diff(y); dy = diff_y.min()
    if lcheck and not np.isclose(dy, diff_y.max(), rtol=1.e-2): 
        raise ValueError("Y-axis is not regular. {} - {}".format(dy, diff_y.max()))
    
    # generate transform
    return Affine.from_gdal(x[0]-dx/2.,dx,0.,y[0]-dy/2.,0.,dy), (len(x),len(y))


def readCFCRS(xds, grid_mapping=None, lraise=True, lproj4=False):
    ''' function to generate CRS from CF-Convention grid mapping variable; only works with Datasets '''
    # read CF convention string
    if not isinstance(xds,(nc.Dataset,xr.Dataset)):
        raise TypeError("Only xarray of netCDF4 Datasets are supported.")
    atts = getAtts(xds) # works for xarray or netCDF4
    if 'Conventions' in atts:
        cf_str = atts['Conventions']
        if cf_str[:3] != 'CF-' or float(cf_str[3:]) < 1:
            raise ValueError("Only CF convection version 1 or later is supported; not '{}'.".format(cf_str))
    elif lraise:
        raise ValueError("No CF convention attribute found; this Dataset may not adhere to CF conventions.")
    else:
        return None # return without CRS
    # find grid mapping variable
    if grid_mapping:
        if grid_mapping in xds.variables:
            grid_type = grid_mapping
            grid_atts = getAtts(xds.variables[grid_mapping])
        else: 
            raise ValueError("Grid mapping '{}' not found in dataset.".format(grid_mapping))
    else:
        grid_type = None
        grid_varlist = ['Lambert_Conformal']
        for grid_var in grid_varlist:
            if grid_var in xds.variables:
                if grid_type is None:
                    grid_type = grid_var
                    grid_atts = getAtts(xds.variables[grid_var])
                else:
                    raise ValueError("Multiple grid_mapping variables detected:",grid_type,grid_var)
    if grid_type is None:
        if lraise:
            raise NotImplementedError("No supported grid_mapping variable detected:\n",grid_varlist)
        else:
            return None # return without CRS
    elif grid_type == 'Lambert_Conformal':
        assert grid_atts['grid_mapping_name'] == "lambert_conformal_conic", grid_atts
        proj4 = ('+proj=lcc +lat_1={lat_1} +lat_2={lat_1} '.format(lat_1=grid_atts['standard_parallel'])
                 + '+lat_0={lat_0} +lon_0={lon_0} '.format(lat_0=grid_atts['latitude_of_projection_origin'],
                                                           lon_0=grid_atts['longitude_of_central_meridian'])
                 + '+x_0=0 +y_0=0 +a=6371229 +b=6371229 +units=m +no_defs' )
    else:
        raise NotImplementedError("The grid_mapping '{}' is currently not implemented/supported.".format(grid_type))
    import rasterio as rio
    # return either string or CRS object
    if lproj4: crs = proj4
    else: crs = rio.crs.CRS.from_string(proj4) # initialize from Proj4 string
    return crs

def getCRS(xvar, lraise=True):
    ''' infer projection from a xarray Dataset or DataArray; this function assumes that either a proj4 string or
        an EPSG designation is stored in the attributes of the dataset/variable. '''
    from geospatial.rasterio_tools import genCRS # used to generate CRS object
    
    if isinstance(xvar,(xr.DataArray,xr.Dataset)):
        atts = xvar.attrs
    elif isinstance(xvar,(nc.Variable,nc.Dataset)):        
        atts = getAtts(xvar)
    elif lraise:
        raise TypeError("Can only infer coordinate system from xarray or netCDF4 - not from {}".format(xvar.__class__))
    else: 
        return None # no projection
          
    crs = None
    # check CF convention
    if isinstance(xvar,(xr.Dataset,nc.Dataset)):
        crs = readCFCRS(xvar, lraise=False, lproj4=False)        
    # search for EPSG number
    if crs is None:
        for key,value in atts.items():
            if key.upper() == 'EPSG' and value != 'n/a': crs = genCRS(value); break    
    # search for Proj4 string
    if crs is None:
        for key,value in atts.items():
            if key.lower() == 'proj4' and value != 'n/a': crs = genCRS(value); break    
    # check for simple geographic lat/lon system
    if crs is None:
        if isGeoCRS(xvar, lraise=False): # error will be raised below (if desired)
            crs = genCRS() # no arguments for default lat/lon    
    # return values
    if lraise and crs is None:
        raise ValueError("No projection information found in attributes.")
    
    # return a GDAL/rasterio CRS instance
    return crs


def inferGeoInfo(xvar, varname=None, crs=None, transform=None, size=None, lraise=True, lcheck=True):
    ''' infere geo-reference information from xarray DataArray or Dataset and netCDF4 Dataset '''
    
    # CRS
    _crs = getCRS(xvar, lraise=lraise)
    if crs is None: crs = _crs
    elif crs != _crs: 
        from geospatial.rasterio_tools import genCRS # used to generate CRS object
        crs = genCRS(crs)
        if crs != _crs:
            raise ValueError("Prescribed CRS and inferred CRS are incompatible:\n{}\n{}".format(crs,_crs))
    crs = _crs # for some reason EPSG ints also pass the equality test...
    
    # geotransform & grid size
    xlon,ylat = getGeoCoords(xvar, lraise=True, lvars=False)
    _transform, _size = getTransform(xvar, lcheck=lraise)
    if transform is None: transform = _transform
    elif not transform is _transform:
        raise ValueError("Prescribed and inferred Geotransform are incompatible:\n{}\n{}".format(transform,_transform))
    if size is None: size = _size
    elif not size is _size:
        raise ValueError("Prescribed and inferred grid sizes are incompatible:\n{}\n{}".format(size,_size))
    
    # do some checks
    if lcheck:
        if crs.is_projected and isGeoCRS(xvar):
            raise ValueError(crs,xvar) # simple check
        if isinstance(xvar,xr.Dataset) and varname: 
            xvar = xvar[varname]
        shape = None; dims = None
        if isinstance(xvar,xr.DataArray): 
            shape = xvar.data.shape; dims = xvar.dims
            if xvar.attrs.get('dim_order',None) is False:
                raise NotImplementedError("The x/lon and y/lat axes of this xarray have to be swapped:\n {}".format(xvar))
        elif isinstance(xvar,nc.Dataset) and varname:
            xvar = xvar.variables[varname]
            shape = xvar.shape; dims = xvar.dimensions
        if shape:
            if shape[-2:] != (size[1],size[0]):
                raise ValueError(xvar)
        if dims:
            if dims[-2] != ylat or dims[-1] != xlon:
                raise ValueError(xvar)
          
    # return verified georef info
    return crs, transform, size    


## functions that modify a dataset

def _inferVarmap(varmap=None, varatts=None, linvert=False):
    ''' simple function that infers a varmap using varatts, if necessary '''
    if varmap is None:
        varmap = dict()
        if varatts is not None:
            for varname,atts in varatts.items():
                if 'name' in atts: varmap[varname] = atts['name']  
    elif not isinstance(varmap,dict): 
        raise TypeError(varmap)
    if linvert:
        varmap = {value:key for key,value in varmap.items()}
    # return varmap (guaranteed to be a dict)
    return varmap

def updateVariableAttrs(xds, varatts=None, varmap=None, varlist=None, **kwargs):
    ''' a helper function to update variable attributes, rename variables, and apply scaling factors '''
    # update varatts
    if varatts is None: 
        varatts = dict()
    elif isinstance(varatts,dict): 
        varatts = varatts.copy()
    else: 
        raise TypeError(varatts)
    varatts.update(kwargs) # add kwargs
    # generate varmap
    varmap = _inferVarmap(varmap=varmap, varatts=varatts, linvert=False)
    # drop variables
    if varlist is not None:
        drop_list = []
        for varname in xds.data_vars.keys():
            name = varmap.get(varname,varname)
            if name not in varlist: drop_list.append(varname)
        xds = xds.drop_vars(drop_list)                    
    # update attributes (using old names)
    date_str = datetime.today().strftime('%Y%m%d')
    for varname,atts in varatts.items():
        if varname in xds.variables:
            if varname == 'time':
                warn("The 'time' coordinate is handled automatically by xarray using numpy datetime64; "
                      + "changing attributes can break this functionality when the dataset is saved to file. ")
            var = xds.variables[varname]
            attrs = var.attrs.copy()
            if 'updated' not in attrs:
                if 'units' in atts:
                  if 'units' not in attrs or attrs['units'] != atts['units']:
                    if 'scalefactor' in atts and atts['scalefactor'] != 1:
                        var *= atts['scalefactor'] # this should execute lazily...
                    if 'offset' in atts and atts['offset'] != 0:
                        var += atts['offset'] # this should execute lazily...
                # transfer attributes
                for key,value in atts.items():
                    if key not in ('scalefactor','offset'):
                        if key in attrs: attrs['old_'+key] = attrs[key]
                        attrs[key] = value
                attrs['updated'] = date_str # indicate we have updated with date string
                var.attrs = attrs
    # actually rename (but only vars that are present and need to be renamed...)
    xds = xds.rename({key:val for key,val in varmap.items() if key in xds.variables and key != val})
    xds = xds.rename_dims({key:val for key,val in varmap.items() if key in xds.dims and key != val})
    xds.attrs['updated'] = date_str
    return xds


def addGeoReference(xds, proj4_string=None, x_coords=None, y_coords=None, lcreate=False, xlon_coord=None, ylat_coord=None):
    ''' helper function to add GDAL/rasterio-style georeferencing information to an xarray dataset;
        note that this only refers to attributed, not axes, but also includes variables '''
    xlon,ylat = getGeoCoords(xds, x_coords=x_coords, y_coords=y_coords, lvars=lcreate, lraise=not lcreate)
    if lcreate:
        if (xlon is None and ylat is None):
            assert xlon_coord is not None and ylat_coord is not None
            # need to find names again...
            xlon_dim,ylat_dim = getGeoDims(xds, x_coords=x_coords, y_coords=y_coords, lraise=True)
            # create new xlon/ylat coordinates, based on coordinates passed down
            coords = {xlon_dim:xlon_coord, ylat_dim:ylat_coord}
            xds = xds.assign_coords(**coords)
        elif (xlon is not None) and (ylat is not None):
            xlon = xlon.name; ylat = ylat.name # from here on only need names
        else:
            raise ValueError("No valid pair of geographic coodinates found:\n {}".format(xds.dims))
    xds.attrs['xlon'] = xlon
    xds.attrs['ylat'] = ylat
    if proj4_string is None:
        if isGeoVar(xds, x_coords, y_coords, lraise=True):
            proj4_string = '+proj=longlat +lon_0=0 +lat_0=0 +ellps=WGS84 +datum=WGS84' # default geographic, also EPSG 4326
        else:
            raise ValueError("Cannot infer projection - need to provide proj4 string!")
    elif isinstance(proj4_string,str):
        xds.attrs['proj4'] = proj4_string
    else:
        raise TypeError("Cannot infer projection - need to provide proj4 string!")
    for xvar in list(xds.data_vars.values()): 
        if isGeoVar(xvar):
            xvar.attrs['proj4'] = proj4_string
            xvar.attrs['xlon'] = xlon
            xvar.attrs['ylat'] = ylat
            xvar.attrs['dim_order'] = int( xvar.dims[-2:] == (ylat, xlon) )
            # N.B.: the NetCDF-4 backend does not like Python bools
    return xds


def rechunkTo2Dslices(xvar, **other_chunks):
    ''' convenience function to rechunk an xarray so that the horizontal dimensions are contiguous (not chunked)
        N.B.: rechunking in a way that does not simply combine existing chunks seems to cause all chunks/data
              to be loaded into memory (we want to avoid that); also, chunks are defined by their size, not by 
              their number, i.e. the definition for one large 2D chunk is (len(y),len(x)) and *not* (1,1) '''
    if not isinstance(xvar,(xr.DataArray,xr.Dataset)):
        raise TypeError(xvar)
    # old chunk sizes
    if 'chunksizes' in xvar.encoding:
        chunks = {dim:cs for dim,cs in zip(xvar.sizes,xvar.encoding['chunksizes'])}
    else: chunks = dict()
    chunks.update(other_chunks)
    # find horizontal/map dimensions
    xlon = xvar.attrs['xlon']; ylat = xvar.attrs['ylat'] 
    chunks[xlon] = xvar.sizes[xlon]; chunks[ylat] = xvar.sizes[ylat]
    return xvar.chunk(chunks=chunks) # rechunk x/lon and y/lat
  
def autoChunkXArray(xds, chunks=None, dims=None, **kwargs):
    ''' apply auto-chunking to an xarray object, like a Dataset or DataArray (chunks kw arg can override) '''
    from geospatial.netcdf_tools import autoChunk
    if dims is None:
        xlon,ylat = getGeoCoords(xds)
        dims = ('time', ylat.name, xlon.name)
    dims = [dim for dim in dims if dim in xds.sizes]
    shape = [xds.sizes[dim] for dim in dims]
    cks = autoChunk(shape, **kwargs)
    cks = {dim:c for dim,c in zip(dims,cks)}
    if chunks: cks.update(chunks) # manually/explicitly specified chunks override 
    return xds.chunk(chunks=cks)

def getCommonChunks(xds, method='min'):
    ''' get smallest/largest/mean common denominator for chunks in dataset '''
    chunk_list = dict()
    # collect chunks
    if isinstance(xds,xr.Dataset):
        for xvar in xds.data_vars.values():
            if 'chunksizes' in xvar.encoding:
                for dim,cks in zip(xvar.dims,xvar.encoding['chunksizes']):
                    if dim in chunk_list: chunk_list[dim].append(cks)
                    else: chunk_list[dim] = [cks]
    elif isinstance(xds,nc.Dataset):
        for ncvar in xds.variables.values():
            if ncvar.chunking():
                for dim,cks in zip(ncvar.dimensions,ncvar.chunking()):
                    if dim in chunk_list: chunk_list[dim].append(cks)
                    else: chunk_list[dim] = [cks]
    else:
        raise TypeError(xds)
    # reduce chunks
    chunks = dict()
    for dim,cks in list(chunk_list.items()):
        chunks[dim] = getattr(np,method)(cks)
    # return dict with chunksize for each dimension
    return chunks        


def computeNormals(xds, aggregation='month', time_stamp='time_stamp', lresample=False, time_name='time'):
    ''' function invoking lazy groupby() call and replacing the resulting time axis with a new time axis '''  
    
    lts = time_stamp and time_stamp in xds
    # time stamp variable for meta data
    if lts:
        import pandas as pd
        ts_var = xds[time_stamp].load()
        period = (pd.to_datetime(ts_var.data[0]).year, (pd.to_datetime(ts_var.data[-1])+pd.Timedelta(31, unit='D')).year)
        prdstr = '{:04d}-{:04d}'.format(*period)
    
    # resample data to aggregation interval
    if lresample:
        if aggregation.lower() == 'month': rsi = 'MS'
        else:
            raise NotImplementedError(aggregation)
        xds = xds.resample(time=rsi,skipna=True,).mean()
    # N.B.: I am not sure to which extent resampling is necessary
    
    # compute monthly normals
    cds = xds.groupby('time.'+aggregation).mean('time')
    assert len(cds['month']) == 12, cds
    
    # convert time axis
    cds = cds.rename({aggregation:time_name}) # the new time axis is named 'month'
    tm = cds.coords[time_name]
    tm.attrs['name']       = time_name
    tm.attrs['long_name']  = 'Calendar '+aggregation.title()
    tm.attrs['units']      = aggregation
    # add period info for quick identification
    if lts:        
        tm.attrs['start_date'] = str(ts_var.data[0])
        tm.attrs['end_date']   = str(ts_var.data[-1])
        tm.attrs['period']     = prdstr
        # add attributes to dataset
        cds.attrs['start_date'] = str(ts_var.data[0])
        cds.attrs['end_date']   = str(ts_var.data[-1])
        cds.attrs['period']     = prdstr    
        
    # return formatted climatology dataset      
    return cds         
         
         
## function to load a dataset         

def _multichunkPresets(multi_chunks):
    ''' translate string identifiers into valid multichunk dicts, based on presets '''
    if isinstance(multi_chunks,str):
        if multi_chunks.lower() == 'regular': # 256 MB
            multi_chunks = {dim:16 for dim in ('lat','lon','latitude','longitude','x','y',)}
            multi_chunks['time'] = 8
        elif multi_chunks.lower() == 'small': # 64 MB
            multi_chunks = {dim:8 for dim in ('lat','lon','latitude','longitude','x','y','time')}
        elif multi_chunks.lower() == 'time': # 184 MB
            multi_chunks = {dim:4 for dim in ('lat','lon','latitude','longitude','x','y')}
            multi_chunks['time'] = 92 # for reductions along time, we can use a higher value (8 days * 92 ~ 2 years)
        else:
            raise NotImplementedError(multi_chunks)
    elif ( multi_chunks is not None ) and not isinstance(multi_chunks, dict):
        raise TypeError(multi_chunks)
    # return valid multi_chunks (dict)
    return multi_chunks
         
def loadXArray(varname=None, varlist=None, folder=None, varatts=None, filename_pattern=None, filelist=None, default_varlist=None, 
               varmap=None, mask_and_scale=True, grid=None, lgeoref=True, geoargs=None, chunks=True, multi_chunks=None, 
               ldropAtts=False, lskip=False, filetypes=None, 
               compat='override', join='inner', fill_value=np.NaN, combine_attrs='no_conflicts', **kwargs):
    ''' function to open a dataset in one of two modes: 1) variables are stored in separate files, but in the same folder (this mainly 
        applies to high-resolution, high-frequency (daily) observations, e.g. SnoDAS) or 2) multiple variables are stored in different
        filetypes and each is opened and then merged (usually model output); datasets are opened using xarray '''
    # load variables
    if filetypes is None: 
        lopt1 = True
        # option 1: one variable per file
        if varname and varlist: 
            raise ValueError(varname,varlist)
        elif varname:
            varlist = [varname] # load a single variable
        elif varlist is None:
            varlist = default_varlist
        # add variable filetypes
        # if there is a (implied) varmap, we need to apply that to variable-filetypes
        ravmap = _inferVarmap(varmap=varmap, varatts=varatts, linvert=True)
        filetypes = [ravmap.get(varname,varname) for varname in varlist]
        # now also transform varatts and varmap
        varmap_single = None if varmap is None else varmap.copy()
        varatts_single = None if varatts is None else varatts.copy()
        varatts = {filetype:varatts_single for filetype in filetypes}
        varmap = {filetype:varmap_single for filetype in filetypes}
    else:
        lopt1 = False # just to remember when using option 2
    ## now use option 2: multiple variables per file
    # expand varmap to filetypes
    if varmap is None: 
        varmap = {filetype:None for filetype in filetypes} # no varmap
    elif isinstance(varmap,dict):
        filetypes_set = set(filetypes); varmap_set = set(varmap.keys())
        if varmap_set.issubset(filetypes_set) or filetypes_set.issubset(varmap_set): # expand to filetypes using None
            for filetype in filetypes:
                if filetype in varmap_set:
                    if not isinstance(varmap[filetype],dict) and varmap[filetype] is not None:
                        raise TypeError(filetype,varmap[filetype])
                else:
                    varmap[filetype] = None
        elif any([key in filetypes for key in varmap.keys()]):
            raise ValueError("It is unclear if varmap is a dict containing varmap dicts for each filetype or just one varmap dict.",varmap.keys())

        if all([key in filetypes for key in varmap.keys()]): # one varmap per filetype
            if not all([isinstance(value,dict) or value is None for value in varmap.values()]):
                raise TypeError(varmap)
        elif any([key in filetypes for key in varmap.keys()]):
            raise ValueError(varmap.keys())
        else: 
            varmap = {filetype:varmap for filetype in filetypes} # same varmap for all
    else:
        raise TypeError(varmap)
    # expand varatts to filetypes
    if varatts is None: 
        varatts = {filetype:None for filetype in filetypes} # no varatts
    elif isinstance(varatts,dict):
        filetypes_set = set(filetypes); varatts_set = set(varatts.keys())
        if varatts_set.issubset(filetypes_set) or filetypes_set.issubset(varatts_set): # expand to filetypes using None
            for filetype in filetypes:
                if filetype in varatts_set:
                    if not isinstance(varatts[filetype],dict) and varatts[filetype] is not None:
                        raise TypeError(filetype,varatts[filetype])
                else:
                    varatts[filetype] = None
        elif any([key in filetypes for key in varatts.keys()]):
            raise ValueError("It is unclear if varatts is a dict containing varatts dicts for each filetype or just one varatts dict.",varatts.keys())
        else: 
            varatts = {filetype:varatts for filetype in filetypes} # same varatts for all
    else:
        raise TypeError(varatts)
    # expand filename/pattern to filetypes
    if filename_pattern and not filelist: 
        filelist = filename_pattern
    if isinstance(filelist, dict):
        if len(filelist) != len(filetypes):
            raise ValueError(filelist)
    elif isinstance(filelist, str):
        filelist = {filetype:filelist for filetype in filetypes}
    else:
        raise ValueError(filelist)
    # just some default settings that will produce chunks larger than 100 MB on 8*64*64 float chunks
    multi_chunks = _multichunkPresets(multi_chunks)
    orig_chunks = chunks.copy() if isinstance(chunks, dict) else chunks # deep copy or True or None
    # construct dataset
    ds_list = []
    for filetype in filetypes:
        filename = filelist[filetype].lower().format(var=filetype.lower(), type=filetype.lower()) # all lower case
        filepath = '{}/{}'.format(folder,filename)
        chunks = orig_chunks # reset
        # apply varmap in reverse to varlist
        if os.path.exists(filepath):
            # load dataset
            if chunks is True:
                # infer chunks from NetCDF-4 file (not sure why xarray doesn't do this automatically...)
                with nc.Dataset(filepath, 'r') as ncds : # open in read-only using NetCDF4 module
                    chunks = dict()
                    for varname,ncvar in ncds.variables.items():
                        for dim,size in zip(ncvar.dimensions,ncvar.chunking()):
                            chunks[dim] = size # this just selects the last value... not necessarily always the same
                            if dim in chunks and chunks[dim] != size:
                                print("WARNING: Chunks for dimension '{}' not coherent in file:\n '{}'".format(dim, filepath))
            if multi_chunks: # enlarge chunks with multiplier
                chunks = {dim:(val*multi_chunks.get(dim,1)) for dim,val in chunks.items()}                  
            # open dataset with xarray
            #print(varname,chunks)
            ds = xr.open_dataset(filepath, chunks=chunks, mask_and_scale=mask_and_scale, **kwargs)
            # N.B.: the use of open_mfdataset is problematic, because it does not play nicely with chunking - 
            #       by default it loads everything as one chunk, and it only respects chunking, if chunks are 
            #       specified explicitly at the initial load time (later chunking seems to have no effect!)
            #       That being said, I don't know if this is still the case...
            # rename, prune/drop vars and apply attributes
            if ldropAtts: ds.attrs = dict() # drop original attributes from NC file (still add georef etc.)
            if varatts or varmap:
                ds = updateVariableAttrs(ds, varatts=varatts[filetype], varmap=varmap[filetype], 
                                         varlist=None if lopt1 else varlist)
            ds_list.append(ds)
        else:
            if lskip:
                print("Skipping missing dataset file '{}' ('{}')".format(filename,folder))
            else:
                raise IOError("The dataset file '{}' was not found in folder:\n '{}'".format(filename,folder))
    # merge into new dataset
    if len(ds_list) == 0:
        raise ValueError("Dataset is empty - aborting! Folder: \n '{}'".format(folder))
    # resolve a very common conflict caused by NCO logging
    if np.sum(['history' in ds.attrs for ds in ds_list]) > 1:
        for ds in ds_list: 
            if 'history' in ds.attrs: ds.attrs['history'] = 'conflicting sources'
    xds = xr.merge(ds_list, compat=compat, join=join, fill_value=fill_value, combine_attrs=combine_attrs)
    # add projection info
    if lgeoref:
        if geoargs is not None:
            # check
            if 'proj4' in xds.attrs and 'proj4_string' in geoargs:
                if xds.attrs['proj4'] != geoargs['proj4_string']:
                    raise ValueError(xds.attrs['proj4'])
            # custom options 
            xds = addGeoReference(xds, **geoargs)
        # default options            
        elif 'proj4' in xds.attrs: 
            # read projection string
            xds = addGeoReference(xds, proj4_string=xds.attrs['proj4'])
        elif grid:
            # load griddef from pickle
            from geodata.gdal import loadPickledGridDef
            griddef = loadPickledGridDef(grid=grid)
            xds = addGeoReference(xds, proj4_string=griddef.projection.ExportToProj4(),) 
        else: 
            # use default lat/lon, if it works...
            xds = addGeoReference(xds,) 
    return xds
  
    
def saveXArray(xds, filename=None, folder=None, mode='overwrite', varlist=None, chunks=None, encoding=None, laddTime=None, 
               time_dim='time', time_agg=None, ltmpfile=True, lcompute=True, lprogress=True, lfeedback=True, **kwargs):
    ''' function to save a xarray dataset to disk, with options to add/overwrite variables, choose smart encoding, 
        add timstamps, use a temp file, and handle dask functionality '''
    from geospatial.netcdf_tools import addTimeStamps, addNameLengthMonth
    # file path and tmp file
    if folder: 
        filepath = '{}/{}'.format(folder,filename)
    # if file exists, get varlist and chunks
    if not os.path.exists(filepath) or mode.lower() in ('overwrite','write'):
        # create a new file
        nc_mode = 'w'
        if lfeedback: print("\nExporting to new NetCDF-4 file:")
    else:
        # if file exists and we are appending...
        nc_mode = 'a' # in most cases
        ltmpfile = not lcompute # only works with new file (or dummy...)
        if mode.lower() in ('add_new',):
            if lfeedback: print("\nAppending to existing NetCDF-4 file (only adding new variables):")
        elif mode.lower() in ('add_all',):
            if lfeedback: print("\nAppending to existing NetCDF-4 file (overwriting existing variables):")
        else:
            raise ValueError(mode)
    # determine tmp file
    if ltmpfile: 
        tmp_filepath = filepath + ( '.tmp' if lcompute else '.test' ) # use temporary file during creation
    else:
        tmp_filepath = filepath
    if lfeedback: print(" '{}'".format(tmp_filepath))
    ## handle varlist and existing variables in file
    # user-supplied varlist
    if varlist:
        drop_vars = [xvar for xvar in xds.data_vars.keys() if xvar not in varlist]
        xds = xds.drop_vars(drop_vars) # returns a shallow copy with vars removed
    # handle existing 
    if nc_mode == 'a':
        # open existing file and get encoding
        with nc.Dataset(filepath, 'r') as ncds:
            if chunks is None: chunks = getCommonChunks(ncds)
            if mode.lower() == 'add_new':
                nc_varlist = [var for var in ncds.variables.keys() if var not in ncds.dimensions]
                drop_vars = [xvar for xvar in xds.data_vars.keys() if xvar in nc_varlist]
                xds = xds.drop_vars(drop_vars) # returns a shallow copy with vars removed
        # adding all variables and overwriting existing ones, requires no changes except nc_mode = 'a'
    # setup encoding
    if encoding is None: 
        encoding = dict(); default = None
    else:
        default = encoding.pop('DEFAULT',None)
    for varname,xvar in xds.data_vars.items():
        tmp = zlib_default.copy()
        cks = tuple(1 if dim == 'time' else chunks[dim] for dim in xvar.dims)
        tmp['chunksizes'] = cks # depends on variable
        # N.B.: use chunk size 1 for time and as before for space; monthly chunks make sense, since
        #       otherwise normals will be expensive to compute (access patterns are not sequential)            
        if isinstance(xvar.dtype,np.inexact): encoding[varname]['_FillValue'] = np.NaN
        if default: tmp.update(default)
        if varname not in encoding: 
            encoding[varname] = tmp
        else:
            tmp.update(encoding[varname])
            encoding[varname] = tmp
        #print(varname,cks,rvar.encoding)
    # write to NetCDF
    ## write to file (with progress)
    
    # write results to file (actually just create file)
    task = xds.to_netcdf(tmp_filepath, mode=nc_mode, format='NETCDF4', unlimited_dims=['time'], 
                         engine='netcdf4', encoding=encoding, compute=False)
    if lcompute:
        # execute with or without progress bar
        if lprogress:
            with ProgressBar():
                task.compute()
        else:
            task.compute()

        ## add extras
        with nc.Dataset(tmp_filepath, mode='a') as ncds:
            if laddTime:
                time_coord = ncds.variables[time_dim]
                tatts = getNCAtts(time_coord)
                tname = tatts.get('long_name','')
                if tname.lower().startswith('calendar '):
                    # info on month for climatology
                    from geospatial.netcdf_tools import default_mon_name_atts
                    if default_mon_name_atts['name'] in ncds.variables:
                        if lfeedback: print("\nName of months variable alrady exists.")
                    else:
                        if lfeedback: print("\nAdding name and length of months.")
                        assert tatts.get('units','').lower().startswith('month'), tatts # this assumes monthly time aggregation
                        assert not time_agg or time_agg.lower().startswith('month')
                        addNameLengthMonth(ncds, time_dim=time_dim)
                else:
                    # time stamps for transient
                    if time_dim+'_stamp' in ncds.variables:
                        if lfeedback: print("\nTime-stamp variable ('{}_stamp') already exists.".format(time_dim))
                    else:
                        time_agg = time_agg.lower()
                        if time_agg.endswith('ly'): time_agg = time_agg[:-2]
                        if lfeedback: print("\nAdding human-readable time-stamp variable ('_stamp').".format(time_dim))
                        addTimeStamps(ncds, units=time_agg) # add time-stamps
            ## make sure the spatial units are present!!! xarray seems to loose the spatial coordinate units
            lgeo = isGeoCRS(ncds, lraise=True)
            for coord in getGeoCoords(ncds, lvars=True, lraise=True):
                if 'units' not in coord.ncattrs():
                    coord.setncattr('units','deg' if lgeo else 'm')
            # store geospatial code version
            ncds.setncattr('geospatial_netcdf_version',geospatial_netcdf_version) 
        
        # replace original file
        if ltmpfile:
            if lfeedback: print("\nMoving file to final destination (overwriting old file):\n '{}'".format(filepath))
            if os.path.exists(filepath): os.remove(filepath)
            os.rename(tmp_filepath, filepath)            
    else:
        # just show some info and save task graph
        if lfeedback:
            print("\nEncoding info:") 
            print(encoding)
            print(task)
            print("\nSaving task graph to:\n '{}.svg'".format(filepath))
        task.visualize(filename=filepath+'.svg')  # This file is never produced
    
    # return file path
    return filepath


if __name__ == '__main__':
    pass