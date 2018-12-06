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
from rasterio.warp import Resampling, reproject, calculate_default_transform
import xarray as xr
from importlib import import_module

## functions to interface rasterio and xarray

def getGeoTransform(xvar=None, x=None, y=None, lcheck=True):
    ''' generate an affine transformation from coordinate axes '''
    if isinstance(xvar,(xr.DataArray,xr.Dataset)):
        for name in ('lon','long','longitude','x'):
            if name in xvar.coords: 
                x = xvar.coords[name]; break
        if x is None: raise ValueError('No x/lon coordinate found!')
        for name in ('lat','latitude','y'):
            if name in xvar.coords: 
                y = xvar.coords[name]; break
        if y is None: raise ValueError('No y/lat coordinate found!')
#     elif isinstance(xvar,xr.Dataset):
#         raise NotImplementedError(xvar)
    # check X-axis
    if isinstance(x,xr.DataArray): x = x.data
    elif not isinstance(x,np.ndarray): 
        raise TypeError(x)
    diff_x = np.diff(x); dx = diff_x.min()
    if lcheck and not np.isclose(dx, diff_x.max(), rtol=1.e-2): 
        raise ValueError("X-axis is not regular: {} - {}".format(dx, diff_x.max()))
    # check Y-axis
    if isinstance(y,xr.DataArray): y = y.data
    elif not isinstance(y,np.ndarray): 
        raise TypeError(y)
    diff_y = np.diff(y); dy = diff_y.min()
    if lcheck and not np.isclose(dy, diff_y.max(), rtol=1.e-2): 
        raise ValueError("Y-axis is not regular. {} - {}".format(dy, diff_y.max()))
    # generate transform
    return rio.transform.Affine.from_gdal(x[0],dx,0.,y[0],0.,dy), (len(x),len(y))

def getProjection(*args,**kwargs):
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
        crs = rio.crs.CRS(**kwargs)
    return crs


## functions for Dask execution

def regrid_array(data, tgt_crs=None, tgt_transform=None, tgt_size=None, resampling='bilinear', 
                 src_crs=None, src_transform=None, src_size=None, **kwargs):
    ''' a function to regrid/reproject a data array to a new grid; src attributes can be inferred from xarray '''
    # infer source attributes
    if isinstance(data,xr.DataArray):
        src_transform, src_size = getGeoTransform(data,)
        data = data.data
    if src_size: 
        src_shp = src_size[::-1] # GDAL convention is band,y/lat,x/lon
        if data.shape[-2:] != src_shp:
            raise ValueError(data.shape,src_shp)
    else:
        src_shp = data.shape[-2:]
    if data.ndim > 2:
        raise NotImplementedError 
    # prepare
    tgt_shp = tgt_size[::-1] # GDAL convention is band,y/lat,x/lon
    tgt_data = np.zeros(tgt_shp)
    if isinstance(resampling,basestring):
        resampling = getattr(Resampling,resampling)
    # do GDAL reprojection
    reproject(data, tgt_data, src_transform=src_transform, src_crs=src_crs,
              dst_transform=tgt_transform, dst_crs=tgt_crs, resampling=resampling, **kwargs)
    # return regridded data
    return tgt_data

def write_raster(filename, data, crs=None, transform=None, driver='AAIGrid', missing_value=-9999, 
                 lmask_invalid=True):
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
    # prep data
    if isinstance(data,xr.DataArray): data = data.data
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
    with rio.open(filename, mode='w', driver=driver, crs=crs, transform=transform,
                  width=width, height=height, count=count, dtype=str(data.dtype), nodata=missing_value) as dst:
        if count == 1:
            dst.write(data,1) # GDAL/rasterio bands are one-based
        else:
            # loop over bands (if there are bands)
            for i in range(1,count+1):
                dst.write(data[i,:,:],i) # GDAL/rasterio bands are one-based
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
    start_date = '2011-01-01'; end_date = '2011-02-01'
    root_folder = '{:s}/SON/test/'.format(os.getenv('HGS_ROOT'))
    raster_name = '{:s}_{:s}_{:s}.asc'
    raster_format = 'AAIGrid'
    # southern Ontario grid (UTM 17)
    tgt_crs = getProjection("+proj=utm +zone=17 +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs", name='son1',)
#     tgt_geotrans = rio.transform.Affine.from_gdal(320920.,5.e3,0,4624073.,0,5.e3); tgt_size = (118,82) # 5 km
    tgt_geotrans = rio.transform.Affine.from_gdal(320920.,1.e3,0,4624073.,0,1.e3); tgt_size = (590,410) # 1 km    
    
    ## define source data
    # SnoDAS
    dataset = 'SnoDAS'
    varname = 'snwmlt'
    filename_pattern = raster_name.format(dataset.lower(),varname.lower(),'{:s}')
    print("\n***   Exporting '{}' from '{}' to raster format {}   ***\n".format(varname,dataset,raster_format))
    
    ## get dataset
    ds_mod = import_module('datasets.{0:s}'.format(dataset))
    xds = ds_mod.loadDataset_Daily(varname=varname, time_chunks=8)
    # get georeference
    src_crs = getProjection(name=dataset, **ds_mod.projdict)
    src_geotrans,src_size = getGeoTransform(xds)
    # figure out bounds for clipping
    left,bottom = tgt_geotrans*(0,0); right,top = tgt_geotrans*tgt_size
    trans, w,h = calculate_default_transform(src_crs=tgt_crs, dst_crs=src_crs, width=tgt_size[0], height=tgt_size[1], 
                                             left=left, bottom=bottom, right=right, top=top,)
    left,top = trans*(0,0); rigth,bottom = trans*(w,h)
    # clip source data
    xvar = xds[varname].loc[start_date:end_date, bottom:top, left:right]
    print(xvar)
    
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
        data = xvar[i,:,:] # chunk to work on
        # command to regrid data
        data = dask.delayed( regrid_array )(data, tgt_crs=tgt_crs, tgt_transform=tgt_geotrans, tgt_size=tgt_size,
                                            src_crs=src_crs, src_transform=src_geotrans, src_size=src_size, 
                                            resampling='average')
        # command to write raster
        work_load.append( dask.delayed( write_raster )(filepath, data, driver=raster_format, 
                                                       crs=tgt_crs, transform=tgt_geotrans) )
    
#     # try batch-executing everything at once in a single workload
#     work_load = [ dask.delayed( batch_write_rasters )(filepath_pattern, xvar, 
#                                                       crs=tgt_crs, transform=tgt_geotrans, driver=raster_format) ]
#     ## N.B.: a single workload leads to everything being loaded at once...
    
    # execute delayed computation
    print("\n***   Executing Workload using Dask   ***\n")

#     with dask.set_options(scheduler='processes'):  
    from multiprocessing.pool import ThreadPool
    with dask.set_options(pool=ThreadPool(4)):    
        dask.compute(*work_load)
    
    end = time()
    print("\n***   Completed in {:.2f} seconds   ***\n".format(end-start))
    