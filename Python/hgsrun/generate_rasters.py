'''
Created on Dec. 3, 2018

A module to export ASCII raster data for use with HGS from NetCDF-4 and related datasets. Also includes 
functionality for clipping and reprojection/regridding.

This module is designed to make use of rasterio georeferencing/gdal functionality and xarray/dask lazy 
execution and chunking.

@author: Andre R. Erler, GPL v3
'''

# external imports
import os
import os.path as osp
import numpy as np
import pandas as pd
import rasterio as rio
import xarray as xr
from importlib import import_module


## functions to interface rasterio and xarray

def getGeoTransform(xvar=None, x=None, y=None, lcheck=True):
    ''' generate an affine transformation from coordinate axes '''
    if isinstance(xvar,xr.DataArray):
        for name in ('lon','long','longitude','x'):
            if name in xvar.coords: 
                x = xvar.coords[name]; break
        if x is None: raise ValueError('No x/lon coordinate found!')
        for name in ('lat','latitude','y'):
            if name in xvar.coords: 
                y = xvar.coords[name]; break
        if y is None: raise ValueError('No y/lat coordinate found!')
    elif isinstance(xvar,xr.Dataset):
        raise NotImplementedError(xvar)
    # check X-axis
    if isinstance(x,xr.DataArray): x = x.data
    elif not isinstance(x,np.ndarray): 
        raise TypeError(x)
    diff_x = np.diff(x); dx = diff_x.min()
    if lcheck and not np.isclose(dx, diff_x.max(), rtol=1.e-3): 
        raise ValueError("X-axis is not regular: {} - {}".format(dx, diff_x.max()))
    # check Y-axis
    if isinstance(y,xr.DataArray): y = y.data
    elif not isinstance(y,np.ndarray): 
        raise TypeError(y)
    diff_y = np.diff(y); dy = diff_y.min()
    if lcheck and not np.isclose(dy, diff_y.max(), rtol=1.e-3): 
        raise ValueError("Y-axis is not regular. {} - {}".format(dy, diff_y.max()))
    # generate transform
    return rio.transform.Affine.from_gdal(x[0],dx,0.,y[0],0.,dy)

def getProjection(*args,**kwargs):
    ''' generate a rasterio CRS object, based on Proj4/pyproj convention '''
    if args:
        if len(args) > 1: raise ValueError(args)
        if kwargs: raise ValueError("Can only parse one set of options.")
        arg = args[0]
        if isinstance(arg,basestring):
            crs = rio.crs.CRS(arg) # initialize from Proj4 string
        elif isinstance(arg,(int,np.integer)):
            crs = rio.crs.CRS.from_epsg(arg) # initialize from EPSG reference number
        else: raise TypeError(arg)
    else:
        # initialize from dictionary, with some defaults
        if not kwargs:
            kwargs = dict(proj='longlat',lon_0=0,lat_0=0,x_0=0,y_0=0) # default geographic, wrapping at 180 E
        if 'proj' not in kwargs: kwargs['proj'] = 'longlat'
        elif kwargs['proj'].lower() in ('lonlat','latlon','latlong'):kwargs['proj'] = 'longlat'
        if 'ellps' not in kwargs: kwargs['ellps'] = 'WGS84'
        if 'datum' not in kwargs: kwargs['datum'] = 'WGS84'
        crs = rio.crs.CRS(**kwargs)
    return crs


## functions for Dask execution

def write_raster(filename, data, crs=None, transform=None, driver='AAIGrid', missing_value=-9999):
    ''' write an array to a raster '''
    print(filename)
    data = data.squeeze()
    if data.ndim == 2 :
        count = 1
        height,width = data.shape
    else:
        if driver == 'AAIGrid':
            raise ValueError(driver)
        count,height,width = data.shape
    # optionally infer grid from xarray and defaults
    if crs is None:
        crs = getProjection() # default lat/lon
    if transform is None:
        transform = getGeoTransform(data) 
    # write data
    if isinstance(data,xr.DataArray): data = data.data
    with rio.open(filename, mode='w', driver=driver, crs=crs, transform=transform,
                  width=width, height=height, count=count, dtype=str(data.dtype), nodata=missing_value) as dst:
        if count == 1:
            dst.write(data,1) # GDAL/rasterio bands are one-based
        else:
            # loop over bands (if there are bands)
            for i in range(1,count+1):
                dst.write(data,i) # GDAL/rasterio bands are one-based
    # done...  
    return
  
def batch_write_rasters(filepath_pattern, xvar, crs=None, transform=None, driver='AAIGrid', missing_value=-9999):
    ''' a wrapper that loops over write_raster; note, however, that this loads all data into memory! '''
   
    # print time coordinate
    print(xvar.coords['time']) 
    
    # loop over time axis
    for i,date in enumerate(xvar.coords['time'].data):
    
        # use date to construct file name
        filepath = filepath_pattern.format(pd.to_datetime(date).strftime('%Y%m%d'))
        #print(filepath)
        # command to write raster
        write_raster(filepath, xvar[i,:,:],crs=crs, transform=transform, driver=driver, missing_value=missing_value)


## execute raster export
if __name__ == '__main__':

    import dask
    from time import time
    
    start = time()

    ## define target data/projection
    start_date = '2011-01-01'; end_date = '2012-01-16'
    root_folder = '{:s}/SON/test/'.format(os.getenv('HGS_ROOT'))
    raster_name = '{:s}_{:s}_{:s}.asc'
    raster_format = 'AAIGrid'
    # preliminary bounds
    lat_min, lat_max = 35,45; lon_min,lon_max = -95,-75
    
    ## define source data
    # SnoDAS
    dataset = 'SnoDAS'
    varname = 'snwmlt'
    filename_pattern = raster_name.format(dataset.lower(),varname.lower(),'{:s}')
    print("\n***   Exporting '{}' from '{}' to raster format {}   ***\n".format(varname,dataset,raster_format))
    
    ## get dataset
    ds_mod = import_module('datasets.{0:s}'.format(dataset))
    xds = ds_mod.loadDataset_Daily(varname=varname, time_chunks=8)
    # clip source data
    xvar = xds[varname].loc[start_date:end_date, lat_min:lat_max, lon_min:lon_max]
    print(xvar)
    # get georeference
    src_crs = getProjection(name=dataset, **ds_mod.projdict)
    src_geotrans = getGeoTransform(xvar)
    
    
    ## reproject data to target grid
    tgt_crs = src_crs
    tgt_geotrans = src_geotrans
    
    
    ## write rasters
    subfolder = '{:s}/transient/'.format(dataset,)
    target_folder = root_folder + subfolder
    if not osp.exists(target_folder): os.mkdir(target_folder)
    filepath_pattern = target_folder + filename_pattern
    
    # generate workload for lazy execution
    print("\n***   Constructing Workload from {} to {}.   ***\n".format(start_date,end_date))
    print("Output folder: '{}'\nRaster pattern: '{}'\n".format(target_folder,filename_pattern)) 
    # loop over time coordinate
    work_load = []
    for i,date in enumerate(xvar.coords['time'].data):
     
        # use date to construct file name
        filepath = filepath_pattern.format(pd.to_datetime(date).strftime('%Y%m%d'))
        #print(filepath)
        # command to write raster
        work_load.append( dask.delayed( write_raster )(filepath, xvar[i,:,:], 
                                                       crs=tgt_crs, transform=tgt_geotrans, driver=raster_format) )
    
#     # try batch-executing everything at once in a single workload
#     work_load = [ dask.delayed( batch_write_rasters )(filepath_pattern, xvar, 
#                                                       crs=tgt_crs, transform=tgt_geotrans, driver=raster_format) ]
#     ## N.B.: a single workload leads to everything being loaded at once...
    
    # execute delayed computation
    print("\n***   Executing Workload using Dask   ***\n")
    dask.compute(*work_load)
    
    end = time()
    print("\n***   Completed in {:.2f} seconds   ***\n".format(end-start))
    