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

# WindowsError is not defined on Linux - need a dummy
try: 
    lWin = True
    WindowsError
except NameError:
    lWin = False
    WindowsError = None


## functions to interface rasterio and xarray

# valid geographic/projected coordinates
x_coords = (('lon','long','longitude',), ('x','easting') )
y_coords = (('lat','latitude',),         ('y','northing'))


def getGeoCoords(xvar, x_coords=x_coords, y_coords=y_coords, lraise=True):
    '''  helper function to extract geographic/projected coordinates from xarray'''
    if not isinstance(xvar,(xr.DataArray,xr.Dataset)):
        if lraise: # optionally check input
            raise TypeError("Can only infer coordinates from xarray - not from {}".format(xvar.__class__))
    else:
        # test geographic grid and projected grids separately
        xlon,ylat = None,None # return None, if nothing is found
        for i in range(len(x_coords)):
            for name,coord in xvar.coords.items():
                if name.lower() in x_coords[i]: 
                    xlon = coord; break
            for name,coord in xvar.coords.items():
                if name.lower() in y_coords[i]: 
                    ylat = coord; break
            if xlon is not None and ylat is not None: break
            else: xlon,ylat = None,None
        # optionally raise error if no coordinates are found, otherwise just return None
        if lraise and (xlon is None or ylat is None):
            raise ValueError("No valid pair of geographic coodinates found:\n {}".format(xvar.dims))
    # return a valid pair of geographic or projected coordinate axis
    return xlon,ylat

  
def isGeoVar(xvar, x_coords=x_coords, y_coords=y_coords, lraise=True):
    ''' helper function to identify variables that have geographic or projected coordintes '''
    if not isinstance(xvar,(xr.DataArray,xr.Dataset)):
        if lraise: # optionally check input
            raise TypeError("Can only infer coordinate system from xarray - not from {}".format(xvar.__class__))
        else: 
            return None # evaluates as False, but allows checking
    else:
        # test geographic grid and projected grids separately
        for i in range(len(x_coords)):
            xlon,ylat = False,False
            for name in xvar.coords.keys():
                if name.lower() in x_coords[i]: 
                    xlon = True; break
            for name in xvar.coords.keys():
                if name.lower() in y_coords[i]: 
                    ylat = True; break
            if xlon and ylat: break
    # if it has a valid pair of geographic or projected coordinate axis
    return ( xlon and ylat )

  
def isGeoCRS(xvar, lraise=True):
    ''' helper function to determine if we have a simple geographic lat/lon CRS '''
    lat,lon = False,False
    if not isinstance(xvar,(xr.DataArray,xr.Dataset)):
        if lraise:
            raise TypeError("Can only infer coordinate system from xarray - not from {}".format(xvar.__class__))
        else: 
            return None # evaluates as False, but allows checking
    else:
        for name in xvar.coords.keys():
            if name.lower() in ('lon','long','longitude',): 
                lon = True; break
        for name in xvar.coords.keys():
            if name.lower() in ('lat','latitude',): 
                lat = True; break
    # it is a geographic coordinate system if both, lat & lon are present
    return ( lat and lon )


def getTransform(xvar=None, x=None, y=None, lcheck=True):
    ''' generate an affine transformation from coordinate axes '''
    if isinstance(xvar,(xr.DataArray,xr.Dataset)):
        x,y = getGeoCoords(xvar, lraise=True)
    elif xvar:
        raise TypeError('Can only infer GeoTransform from xarray Dataset or DataArray - not from {}.'.format(xvar))
    # check X-axis
    if isinstance(x,xr.DataArray): x = x.data
    if not isinstance(x,np.ndarray): 
        raise TypeError(x)
    diff_x = np.diff(x); dx = diff_x.min()
    if lcheck and not np.isclose(dx, diff_x.max(), rtol=1.e-2): 
        raise ValueError("X-axis is not regular: {} - {}".format(dx, diff_x.max()))
    # check Y-axis
    if isinstance(y,xr.DataArray): y = y.data
    if not isinstance(y,np.ndarray): 
        raise TypeError(y)
    diff_y = np.diff(y); dy = diff_y.min()
    if lcheck and not np.isclose(dy, diff_y.max(), rtol=1.e-2): 
        raise ValueError("Y-axis is not regular. {} - {}".format(dy, diff_y.max()))
    # generate transform
    return rio.transform.Affine.from_gdal(x[0],dx,0.,y[0],0.,dy), (len(x),len(y))


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
        crs = rio.crs.CRS(**kwargs)
    return crs


def getProj(xvar, lraise=True):
    ''' infer projection from a xarray Dataset of DataArray '''
    if not isinstance(xvar,(xr.DataArray,xr.Dataset)):
        if lraise:
            raise TypeError("Can only infer coordinate system from xarray - not from {}".format(xvar.__class__))
        else: 
            return None # no projection
    proj = None
    # search for Proj4 string
    for key,value in xvar.attr.items():
        if key.lower() == 'proj4': proj = genProj(value); break
    # search for EPSG number
    if proj is None:
        for key,value in xvar.attr.items():
            if key.upper() == 'EPSG': proj = genProj(value); break
    # check for simple geographic lat/lon system
    if proj is None:
        if isGeoCRS(xvar, lraise=False): # error will be raised below (if desired)
            proj = genProj() # no arguments for default lat/lon
    # return values
    if lraise and proj is None:
        raise ValueError("No projection information found in attributes.")
    return proj


## functions for Dask execution

def regrid_array(data, tgt_crs=None, tgt_transform=None, tgt_size=None, resampling='bilinear', 
                 src_crs=None, src_transform=None, src_size=None, lxarray=False, **kwargs):
    ''' a function to regrid/reproject a data array to a new grid; src attributes can be inferred from xarray '''
    # infer source attributes
    if isinstance(data,xr.DataArray):
        if src_transform is None:
            src_transform, size = getTransform(data,)
        if src_size and src_size != size: raise ValueError(src_size,size)
        else: src_size = size
        if data.attrs.get('dim_order',None) is False:
            raise NotImplementedError("This the x/lon and y/lat axes of this xarray have to be swapped:\n {}".format(data))
        if src_crs is None:
            if isGeoCRS(data): src_crs = getProj()
            else: raise ValueError("No source projection 'src_crs' supplied and can't infer projection from source data.")
        data = data.data
    if src_size: 
        if data.shape[-2:] != src_size[::-1]:
            raise ValueError(data.shape,src_size[::-1])
        src_shp = data.shape[:-2] + src_size[::-1] # GDAL convention is band,y/lat,x/lon
    else:
        src_shp = data.shape[-2:]
    if data.ndim > 3:
        raise NotImplementedError 
    # prepare data
    tgt_shp = src_shp[:-2] + tgt_size[::-1] # GDAL convention is band,y/lat,x/lon
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
        resampling = getattr(Resampling,resampling)
    # do GDAL reprojection
    reproject(src_data, tgt_data, src_transform=src_transform, src_crs=src_crs,
              dst_transform=tgt_transform, dst_crs=tgt_crs, resampling=resampling, **kwargs)
    # restore shape
    if len(tgt_shp) > 3:
        tgt_data = tgt_data.reshape(tgt_shp)
    # cast as xarray
    if lxarray:
        raise NotImplementedError
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
        crs = getProj(data, lraise=True) # default lat/lon
    if transform is None:
        transform, size = getTransform(data, lcheck=False) 
        assert data.shape == size, data.shape
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

    ## fast test config
#     loverwrite = True
#     start_date = '2011-01-01'; end_date = '2011-02-01'
#     grid_name  = 'son1'
    ## operational test config
    loverwrite = False
    start_date = '2011-06-01'; end_date = '2011-07-01'
    grid_name  = 'son2'
    ## operational config
#     loverwrite = False
#     start_date = None; end_date = None
#     grid_name  = 'son2'

    ## define target data/projection
    root_folder = '{:s}/SON/{:s}/'.format(os.getenv('HGS_ROOT'),grid_name)
    raster_name = '{dataset:s}_{variable:s}_{grid:s}_{date:s}.asc'
    raster_format = 'AAIGrid'
    # southern Ontario grid (UTM 17)
    if grid_name == 'son1':
        tgt_size = (118,82) # lower resolution 5 km grid
        tgt_geotrans = rio.transform.Affine.from_gdal(320920.,5.e3,0,4624073.,0,5.e3) # 5 km
    elif grid_name == 'son2':
        tgt_size = (590,410) # higher resolution 1 km grid (~ 1 MB per day)
        tgt_geotrans = rio.transform.Affine.from_gdal(320920.,1.e3,0,4624073.,0,1.e3) # 1 km    
    # southern Ontario projection
    tgt_crs = genProj("+proj=utm +zone=17 +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs", name=grid_name)
    
    ## define source data
    # SnoDAS
    dataset = 'SnoDAS'
    varname = 'liqwatflx'
    filename_pattern = raster_name.format(dataset=dataset.lower(), variable=varname.lower(),
                                          grid=grid_name.lower(), date='{:s}') # no date for now...
    print("\n***   Exporting '{}' from '{}' to raster format {}   ***\n".format(varname,dataset,raster_format))
    
    ## get dataset
    ds_mod = import_module('datasets.{0:s}'.format(dataset))
    xds = ds_mod.loadDataset_Daily(varname=varname, time_chunks=8)
    # get georeference
    src_crs = getProj(xds)
    # figure out bounds for clipping
    left,bottom = tgt_geotrans*(0,0); right,top = tgt_geotrans*tgt_size
    trans, w,h = calculate_default_transform(src_crs=tgt_crs, dst_crs=src_crs, width=tgt_size[0], height=tgt_size[1], 
                                             left=left, bottom=bottom, right=right, top=top,)
    left,top = trans*(0,0); rigth,bottom = trans*(w,h)
    # clip source data
    xvar = xds[varname].loc[:, bottom:top, left:right]
    if start_date and end_date:
        xvar = xds[varname].loc[start_date:end_date,:,:]
    print(xvar)
    src_geotrans,src_size = getTransform(xvar)
    
    ## write rasters
    subfolder = '{:s}/transient/'.format(dataset,)
    target_folder = root_folder + subfolder
    # make sure path exists
    try:
        if not osp.exists(target_folder): os.mkdir(target_folder)
    except (WindowsError,OSError):
        os.makedirs(target_folder)
    filepath_pattern = target_folder + filename_pattern
    
    # generate workload for lazy execution
    print("\n***   Constructing Workload from {} to {}.   ***\n".format(start_date,end_date))
    print("Output folder: '{}'\nRaster pattern: '{}'\n".format(target_folder,filename_pattern)) 
    # loop over time coordinate
    work_load = []
    for i,date in enumerate(xvar.coords['time'].data):
      
        # use date to construct file name
        filename = filename_pattern.format(pd.to_datetime(date).strftime('%Y%m%d'))
        filepath = target_folder + filename
        if osp.exists(filepath) and not loverwrite: 
            print("Skipping existing file: {}".format(filename))
        else:
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
    