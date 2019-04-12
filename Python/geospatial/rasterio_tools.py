'''
Created on Feb. 23, 2019

Utility tools to work with rasterio or that use rasterio as their backend.

@author: Andre R. Erler, GPL v3
'''

import functools
import os
import os.path as osp
import numpy as np
from six import string_types # for testing string in Python 2 and 3
import pandas as pd
import rasterio as rio

# prevent failue if xarray is not installed
try: 
    from xarray import DataArray
    from geospatial.xarray_tools import getTransform, isGeoCRS, getProj
except: 
    DataArray = NotImplemented
    getTransform  = NotImplemented
    isGeoCRS = NotImplemented
    getProj = NotImplemented


## functions to interface with xarray

def genProj(*args,**kwargs):
    ''' generate a rasterio CRS object, based on Proj4/pyproj convention '''
 
    if args:
        if len(args) > 1: raise ValueError(args)
        arg = args[0]
        if isinstance(arg,rio.crs.CRS):
            crs = arg # nothing to do
        elif isinstance(arg,string_types):
            if kwargs:
                for key,value in kwargs.items():
                    arg += ' +{:s}={}'.format(key,value)
            crs = rio.crs.CRS.from_string(arg) # initialize from Proj4 string
        elif isinstance(arg,(int,np.integer)):
            crs = rio.crs.CRS.from_epsg(arg) # initialize from EPSG reference number
            if kwargs: raise ValueError("kwargs don't work with EPSG.")
        else: raise TypeError(arg)
    else:

        # initialize from dictionary, with some defaults
        if not kwargs:
            kwargs = dict(proj='longlat',lon_0=0,lat_0=0,x_0=0,y_0=0) # default geographic, wrapping at 180 E
        if 'proj' not in kwargs: kwargs['proj'] = 'longlat'
        elif kwargs['proj'].lower() in ('lonlat','latlon','latlong'):kwargs['proj'] = 'longlat'
        if 'ellps' not in kwargs: kwargs['ellps'] = 'WGS84'
        if 'datum' not in kwargs: kwargs['datum'] = 'WGS84'
        # create CRS object
        crs = rio.crs.CRS(**kwargs)

    # return a rasterio CRS instance (should be based on GDAL CRS)
    return crs


def constructCoords(geotrans, size, dtype='float32'):
    ''' construct coordinate arrays from geotransform and size; coordinates are centered w.r.t. pixels '''
    
    # use GDAL conventions
    if isinstance(geotrans,rio.transform.Affine):
        geotrans = geotrans.to_gdal()
    # can't handle tilted axes
    if geotrans[2] != 0 or geotrans[4] != 0:
        raise NotImplementedError(geotrans)
    
    # compute x/lon
    xlon = np.arange(size[0], dtype=np.dtype(dtype))
    xlon *= geotrans[1]; xlon += geotrans[0] + geotrans[1]/2.
    # compute y/lat
    ylat = np.arange(size[1], dtype=np.dtype(dtype)) 
    ylat *= geotrans[5]; ylat += geotrans[3] + geotrans[5]/2.
    
    # return numpy arrays     
    return xlon,ylat


## functions for Dask execution

def regrid_array(data, tgt_crs=None, tgt_transform=None, tgt_size=None, resampling='bilinear', 
                 src_crs=None, src_transform=None, src_size=None, lxarray=False, **kwargs):
    ''' a function to regrid/reproject a data array to a new grid; src attributes can be inferred from xarray '''
    
    # infer source attributes
    if isinstance(data,DataArray):
        if src_transform is None:
            src_transform, size = getTransform(data,)
        else: size = (data.shape[-1],data.shape[-2])
        if src_size and src_size != size: raise ValueError(src_size,size)
        else: src_size = size
        if data.attrs.get('dim_order',None) is False:
            raise NotImplementedError("This the x/lon and y/lat axes of this xarray have to be swapped:\n {}".format(data))
        if src_crs is None:
            if isGeoCRS(data): src_crs = getProj()
            else: raise ValueError("No source projection 'src_crs' supplied and can't infer projection from source data.")
        data = data.data
    
    # N.B.: GDAL convention for data arrays is band,y/lat,x/lon,
    #       but the GDAL size tuple is stored as (x/lon,y/lat)
    if src_size: 
        if data.shape[-2:] != src_size[::-1]:
            raise ValueError(data.shape,src_size[::-1])
        src_shp = data.shape[:-2] + src_size[::-1] 
    else:
        src_shp = data.shape[-2:]
    if data.ndim > 3:
        raise NotImplementedError 
    
    # prepare data
    tgt_shp = src_shp[:-2] + tgt_size[::-1] # GDAL convention for data arrays is band,y/lat,x/lon
    # reshape data
    if len(src_shp) > 3:
        bnds = np.prod(src_shp[:-2])
        src_data = data.reshape((bnds,)+src_shp[-2:])
        tgt_data = np.zeros((bnds,)+tgt_shp[-2:])
    else:
        src_data = data
        tgt_data = np.zeros(tgt_shp)
    
    # prepare reprojection
    if isinstance(resampling,string_types):
        resampling = getattr(rio.warp.Resampling,resampling)
    # do GDAL reprojection
    rio.warp.reproject(src_data, tgt_data, src_transform=src_transform, src_crs=src_crs,
                       dst_transform=tgt_transform, dst_crs=tgt_crs, resampling=resampling, **kwargs)
    
    # restore shape
    if len(tgt_shp) > 3:
        tgt_data = tgt_data.reshape(tgt_shp)
    # cast as xarray
    if lxarray:
        raise NotImplementedError
        # we need to assume new geospatial dimensions and construct new coordinate axes,
        # then construct the DataArray and copy the attributes...
    
    # return regridded data
    return tgt_data


def write_raster(filename, data, crs=None, transform=None, driver='AAIGrid', missing_value=-9999, 
                 missing_flag=None, lmask_invalid=True, lecho=True, loverwrite=True):
    ''' write an array to a raster '''
    if osp.exists(filename) and not loverwrite: 
        if lecho: print("Skipping existing file: {}".format(filename))
    else:
        if lecho: print(filename)
        data = data.squeeze()
        if data.ndim == 2 :
            count = 1
            height,width = data.shape
        else:
            if driver == 'AAIGrid':
                raise ValueError(driver)
            height,width = data.shape[-2:]
            count = np.prod(data.shape[:-2])
        
        # prep data
        if isinstance(data,DataArray): 
            data = data.data
            # optionally infer grid from xarray and defaults
            if crs is None:
                crs = getProj(data, lraise=True) # default lat/lon
            if transform is None:
                transform, size = getTransform(data, lcheck=False) 
                assert data.shape == size, data.shape
        elif crs is None or transform is None:
            raise ValueError("Can only infer grid/projection from xarray DataArray.")
        if isinstance(data,np.ma.MaskedArray): data = data.filled(missing_value)
        if lmask_invalid: data[~np.isfinite(data)] = missing_value
        
        # fix transform
        if transform.e > 0:
            geotrans = list(transform.to_gdal())
            # move origin to upper left corner
            geotrans[3] += geotrans[5]*height
            geotrans[5] *= -1
            transform = rio.transform.Affine.from_gdal(*geotrans)
            data = np.flip(data, axis=-2) # and now flip data along y-axis!
        
        # write data
        if missing_flag is None: missing_flag = missing_value
        with rio.open(filename, mode='w', driver=driver, crs=crs, transform=transform,
                      width=width, height=height, count=count, dtype=str(data.dtype), nodata=missing_flag) as dst:
            if count == 1:
                dst.write(data,1) # GDAL/rasterio bands are one-based
            else:
                # loop over bands (if there are bands)
                data = data.reshape((count,height,width)) # merge non-geospatial dimensions
                for i in range(1,count+1):
                    dst.write(data[i,:,:],i) # GDAL/rasterio bands are one-based
    # done...  
    return None
  
    
def write_time_rasters(filepath_pattern, xvar, time_coord=None, crs=None, transform=None, driver='AAIGrid', 
                       missing_value=-9999, missing_flag=None, lecho=True, loverwrite=True):
    ''' a wrapper that iterates over write_raster in order to write multi-dimensional arrays to 2D rasters; 
        note, however, that this loads all data into memory!
        also note that this implementation assumes iteration over a datetime axis and uses the date
        to construct the file name for each file '''
   
    # determine time coordinate
    if time_coord is None:
        if isinstance(xvar,DataArray):
            time_coord = xvar.coords['time'].data
        else: TypeError(xvar)
    elif not np.issubdtype(time_coord.dtype,np.datetime64): 
        raise TypeError(time_coord)
    
    # loop over time axis
    for i,date in enumerate(time_coord):
    
        # use date to construct file name
        filepath = filepath_pattern.format(pd.to_datetime(date).strftime('%Y%m%d'))
        # command to write raster
        write_raster(filepath, xvar[i,:,:],crs=crs, transform=transform, driver=driver, 
                     missing_value=missing_value, missing_flag=missing_flag, lecho=lecho, loverwrite=loverwrite)    
    # done...
    return None


dummy_array = np.zeros((1,1,1), dtype=np.int8)
def regrid_and_export(data, block_id=None, time_chunks=None, src_crs=None, src_geotrans=None, src_size=None, 
                      tgt_crs=None, tgt_geotrans=None, tgt_size=None, mode='raster2D', resampling='bilinear',
                      filepath=None, time_coord=None, driver='AAIGrid',missing_value=0., missing_flag=-9999.,
                      ncvar=None, bias_correction=None, bc_varname=None,
                      lecho=True, loverwrite=True, return_dummy=dummy_array):
    ''' a function for use with Dask lazy execution that regrids an array and exports/writes the results 
        to either NetCDF or any GDAL/rasterio raster format; the array has to be chunked in a way that  
        full 2D horizontal surfaces are processed, otherwise regridding does not work '''
    
    # regrid array
    if tgt_crs != src_crs or tgt_geotrans != src_geotrans or tgt_size != src_size:
        # reproject data
        data = regrid_array(data, tgt_crs=tgt_crs, tgt_transform=tgt_geotrans, tgt_size=tgt_size, 
                            src_crs=src_crs, src_transform=src_geotrans, src_size=src_size, 
                            resampling=resampling, lxarray=False)

    # figure out time coordinates/dates
    ts = time_chunks[block_id[0]]; te = ts + data.shape[0] 
    
    # bias correction
    if bias_correction:
        time_chunk = time_coord[ts:te] # chunk of time axis that corresponds to this chunk
        # N.B.: bias_correction is assumed to be an object that is similar to the BiasCorrection class,
        #       which is defined in the bc_methods module and exposes a 'correctArray' method; the bc_varname
        #       variabel indicates, which correction array from bias_correction should be used
        data = bias_correction.correctArray(data, varname=bc_varname, time=time_chunk, ldt=True, time_idx=0)
        
    if mode.lower() == 'raster2d':    
        # write raster files, one per time step
        if not bias_correction:
            time_chunk = time_coord[ts:te] # chunk of time axis that corresponds to this chunk 
        #print(time_chunk)
        write_time_rasters(filepath, data, time_coord=time_chunk, driver=driver, crs=tgt_crs, 
                           transform=tgt_geotrans, missing_value=missing_value, missing_flag=missing_flag, 
                           lecho=lecho, loverwrite=loverwrite)
    elif mode.upper() == 'NETCDF':
        # append to existing NetCDF variable
        ncvar[ts:te,:,:] = data
    else:
        raise NotImplementedError(mode)
    
    # return a small dummy array, because we need to return something array-like
    return return_dummy


def generate_regrid_and_export(xvar, mode='raster2D', time_coord='time', folder=None, filename=None,
                               tgt_crs=None, tgt_geotrans=None, tgt_size=None, resampling='bilinear', 
                               driver='AAIGrid', missing_value=0., missing_flag=-9999., 
                               bias_correction=None, bc_varname=None, 
                               lecho=True, loverwrite=True,):
    ''' a function that returns another function, which is suitable for regridding and direct export to disk 
        using Dask lazy execution; some parameters can be inferred from xarray attributes '''
    
    # infer xarray meta data
    if not isinstance(xvar,DataArray): 
        raise NotImplementedError("Can only infer grid/projection from xarray DataArray.")
    if len(xvar.dims) != 3:
        raise NotImplementedError(xvar)
    
    # infer grid from xarray and defaults
    src_crs = getProj(xvar, lraise=True) # default lat/lon
    src_geotrans, src_size = getTransform(xvar, lcheck=True) 
    
    # time axis
    if isinstance(time_coord,str):
        assert xvar.dims[0] == time_coord, xvar
        time_coord = xvar.coords[time_coord].data
    elif isinstance(time_coord,DataArray):
        if len(xvar.coords[time_coord.name].data) != len(time_coord):
            raise ValueError(time_coord)
    elif isinstance(time_coord,np.ndarray):
        if xvar.shape[0] != len(time_coord):
            raise ValueError(time_coord)
    else:
        raise TypeError(time_coord)        
        
    # ensure correct chunking 
    if not all([len(c)==1 for c in xvar.chunks[-2:]]): 
        raise ValueError("Regridding can only work, if horizontal dimensions are contiguous;" + 
                         "\n however, the horizontal dimensions appear to be chunked:\n {}".format(xvar.chunks) +
                         "\n Consider using the 'rechunkTo2Dslices' function to fix chunking.")
    
    # create folder
    filepath = osp.join(folder,filename)
    if not os.path.exists(folder):
        os.makedirs(folder)
        
    # also create the NetCDF file/dataset/variable, if necessary
    if mode.upper() == 'NETCDF':
        
        from geospatial.netcdf_tools import createGeoNetCDF, add_var
        ncds = createGeoNetCDF(filepath, time=time_coord, crs=tgt_crs, geotrans=tgt_geotrans, size=tgt_size, 
                               zlib=True, loverwrite=loverwrite)
        
        # create variable, if necessary
        if xvar.name in ncds.variables:
            ncvar = ncds.variables[xvar.name]
        else:
            dims = xvar.dims[:-2]+(ncds.ylat,ncds.xlon) # geographic dimensions need to be replaced
            shp = xvar.shape[:-2]+tgt_size[::-1]
            ncvar = add_var(ncds, name=xvar.name, dims=dims, data=None, shape=shp, 
                            atts=xvar.attrs.copy(), dtype=xvar.dtype, zlib=True,)
        dataset = ncds 
    else:
        ncvar = None
        dataset = folder

    # define options for dask execution
    time_chunks = np.concatenate([[0],np.cumsum(xvar.chunks[0][:-1], dtype=np.int)])
    dummy_array = np.zeros((1,1,1), dtype=np.int8)
    dask_fct = functools.partial(regrid_and_export, time_chunks=time_chunks, resampling=resampling,   
                                 src_crs=src_crs, src_geotrans=src_geotrans, src_size=src_size, 
                                 tgt_crs=tgt_crs, tgt_geotrans=tgt_geotrans, tgt_size=tgt_size, 
                                 mode=mode, filepath=filepath, time_coord=time_coord, ncvar=ncvar,
                                 driver=driver,missing_value=missing_value, missing_flag=missing_flag, 
                                 bias_correction=bias_correction, bc_varname=bc_varname,                                
                                 lecho=lecho, loverwrite=loverwrite, return_dummy=dummy_array)
    
    # return function with dummy array (for measure)
    return dask_fct,dummy_array,dataset


if __name__ == '__main__':
    pass