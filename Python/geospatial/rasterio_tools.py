'''
Created on Feb. 23, 2019

Utility tools to work with rasterio or that use rasterio as their backend.

@author: Andre R. Erler, GPL v3
'''

import os.path as osp
import numpy as np
import pandas as pd
import rasterio as rio

# prevent failue if xarray is not installed
try: 
    from xarray import DataArray
    from xarray_tools import getTransform, isGeoCRS, getProj
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
        if isinstance(arg,basestring):
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
    if isinstance(resampling,basestring):
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
                      width=width, height=height, count=count, dtype=str(data.dtype), nodata=missing_value) as dst:
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


if __name__ == '__main__':
    pass