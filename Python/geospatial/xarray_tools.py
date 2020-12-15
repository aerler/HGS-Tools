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

# internal imports
from geospatial.netcdf_tools import getNCAtts # this import shpuld be fine

## an important option I am relying on!
xr.set_options(keep_attrs=True)

# names of valid geographic/projected coordinates
default_x_coords = dict(geo=('lon','long','longitude',), proj=('x','easting') )
default_y_coords = dict(geo=('lat','latitude',),         proj=('y','northing'))
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

def getGeoCoords(xvar, x_coords=None, y_coords=None, lraise=True, lvars=True):
    '''  helper function to extract geographic/projected coordinates from xarray'''
    
    if x_coords is None: x_coords = default_x_coords
    if y_coords is None: y_coords = default_y_coords
    
    xlon,ylat = None,None # return None, if nothing is found
    
    if isinstance(xvar,(xr.DataArray,xr.Dataset)):
        # test geographic grid and projected grids separately
        for coord_type in x_coords.keys():
            for name,coord in xvar.coords.items():
                if name.lower() in x_coords[coord_type]: 
                    xlon = coord if lvars else name; break
            for name,coord in xvar.coords.items():
                if name.lower() in y_coords[coord_type]: 
                    ylat = coord if lvars else name; break
            if xlon is not None and ylat is not None: break
            else: xlon,ylat = None,None
    elif isinstance(xvar,nc.Variable) and lraise:
        raise TypeError("Cannot infer coordinates from netCDF4 Variable - only Dataset!")
    elif isinstance(xvar,nc.Dataset):
        # test geographic grid and projected grids separately
        for coord_type in x_coords.keys():
            for name in xvar.dimensions:
                if name.lower() in x_coords[coord_type]: 
                    if name in xvar.variables:
                        xlon = xvar.variables[name] if lvars else name; break
            for name in xvar.dimensions:
                if name.lower() in y_coords[coord_type]: 
                    if name in xvar.variables:
                        ylat = xvar.variables[name] if lvars else name; break
            if xlon is not None and ylat is not None: break
            else: xlon,ylat = None,None      
    elif lraise: # optionally check input
        raise TypeError("Can only infer coordinates from xarray or netCDF4 - not from {}".format(xvar.__class__))
    else:
        pass # return None,None
        
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

def updateVariableAttrs(xds, varatts=None, varmap=None, varlist=None, **kwargs):
    ''' a helper function to update variable attributes, rename variables, and apply scaling factors '''
    if varatts is None: 
        varatts = dict()
    elif isinstance(varatts,dict): 
        varatts = varatts.copy()
    else: 
        raise TypeError(varatts)
    varatts.update(kwargs) # add kwargs
    # generate var map
    if varmap is None:
        varmap = dict()
        for varname,atts in varatts.items():
            if varname in xds and 'name' in atts: varmap[varname] = atts['name']  
    elif not isinstance(varmap,dict): 
        raise TypeError(varmap)
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
    xds = xds.rename({key:val for key,val in varmap.items() if key in xds and key != val})
    xds.attrs['updated'] = date_str
    return xds


def addGeoReference(xds, proj4_string=None, x_coords=None, y_coords=None):
    ''' helper function to add GDAL/rasterio-style georeferencing information to an xarray dataset;
        note that this only refers to attributed, not axes, but also includes variables '''
    xlon,ylat = getGeoCoords(xds, x_coords=x_coords, y_coords=y_coords, lvars=False)
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
         
def loadXArray(varname=None, varlist=None, folder=None, varatts=None, filename_pattern=None, default_varlist=None, varmap=None, 
               mask_and_scale=True, grid=None, lgeoref=True, geoargs=None, chunks=True, lautoChunk=False, lskip=False, **kwargs):
    ''' function to open a dataset where variables are stored in separate files (but in the same folder); this mainly applies 
        to high-resolution, high-frequency (daily) observations (e.g. SnoDAS); datasets are opened using xarray '''
    # load variables
    if varname and varlist: 
        raise ValueError(varname,varlist)
    elif varname:
        varlist = [varname] # load a single variable
    elif varlist is None:
        varlist = default_varlist
    # apply varmap in reverse to varlist
    if varmap is None and varatts is not None:
        varmap = {name:atts.get('name',name) for name,atts in varatts.items()}
    if varmap is not None:
        ravmap = {value:key for key,value in varmap.items()}
        varlist = [ravmap.get(varname,varname) for varname in varlist]
    # construct dataset
    xds = None
    for varname in varlist:
        filename = filename_pattern.lower().format(var=varname.lower()) # all lower case
        filepath = '{}/{}'.format(folder,filename)
        if os.path.exists(filepath):
            # load dataset
            cks = {} if chunks is True else chunks
            ds = xr.open_dataset(filepath, chunks=cks, mask_and_scale=mask_and_scale, **kwargs)
            # N.B.: the use of open_mfdataset is problematic, because it does not play nicely with chunking - 
            #       by default it loads everything as one chunk, and it only respects chunking, if chunks are 
            #       specified explicitly at the initial load time (later chunking seems to have no effect!)
            if chunks is True:
                # apply chunks from NetCDF-4 file (not sure why xarray doesn't do this automatically...)
                tmpvar = ds[varname]
                warn("Post-hoc chunking does not seem to work properly!")
                ds[varname] = tmpvar.chunk(chunks=tmpvar.encoding['chunksizes'])
            # merge into new dataset
            if xds is None: xds = ds
            else: xds.update(ds)
        else:
            if lskip:
                print("Skipping missing dataset file '{}' ('{}')".format(filename,folder))
            else:
                raise IOError("The dataset file '{}' was not found in folder:\n '{}'".format(filename,folder))
    assert xds is not None, "Dataset is empty - aborting"
    # rewrite chunking, if desired
    if lautoChunk:
        cks = None if not isinstance(chunks, dict) else chunks
        xds = autoChunkXArray(xds, chunks=cks)
    # rename and apply attributes
    if varatts or varmap:
        xds = updateVariableAttrs(xds, varatts=varatts, varmap=varmap)
    # add projection info
    if lgeoref:
        if geoargs is not None:
            # check
            if 'proj4' in xds.attrs and 'proj4' in geoargs:
                if xds.attrs['proj4'] != geoargs['proj4']:
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


if __name__ == '__main__':
    pass