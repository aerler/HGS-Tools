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
    for key,value in xvar.attrs.items():
        if key.lower() == 'proj4': proj = genProj(value); break
    # search for EPSG number
    if proj is None:
        for key,value in xvar.attrs.items():
            if key.upper() == 'EPSG': proj = genProj(value); break
    # check for simple geographic lat/lon system
    if proj is None:
        if isGeoCRS(xvar, lraise=False): # error will be raised below (if desired)
            proj = genProj() # no arguments for default lat/lon
    # return values
    if lraise and proj is None:
        raise ValueError("No projection information found in attributes.")
    return proj


def rechunkXlonYalt(xvar):
    ''' convenience function to rechunk an xarray so that the horizontal dimensions are contiguous (not chunked)
        N.B.: rechunking in a way that does not simply combine existing chunks seems to cause all chunks/data
              to be loaded into memory (we want to avoid that); also, chunks are defined by their size, not by 
              their number '''
    if not isinstance(xvar,(xr.DataArray,xr.Dataset)):
        raise TypeError(xvar)
    # find horizontal/map dimentions
    xlon = xvar.coords[xvar.attrs['xlon']]; ylat = xvar.coords[xvar.attrs['ylat']]
    return xvar.chunk(chunks={xlon.name:len(xlon),ylat.name:len(ylat)}) # rechunk x/lon and y/lat
         

## functions for Dask execution

def regrid_array(data, tgt_crs=None, tgt_transform=None, tgt_size=None, resampling='bilinear', 
                 src_crs=None, src_transform=None, src_size=None, lxarray=False, **kwargs):
    ''' a function to regrid/reproject a data array to a new grid; src attributes can be inferred from xarray '''
    # infer source attributes
    if isinstance(data,xr.DataArray):
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
                 missing_flag=None, lmask_invalid=True, lecho=True):
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
        if missing_flag is None: missing_flag = missing_value
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
  
def batch_write_rasters(filepath_pattern, xvar, time_coord=None, crs=None, transform=None, 
                        driver='AAIGrid', missing_value=-9999, missing_flag=None, lecho=True):
    ''' a wrapper that loops over write_raster; note, however, that this loads all data into memory! '''
   
    # determine time coordinate
    if time_coord is None:
        if isinstance(xvar,xr.DataArray):
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
                     missing_value=missing_value, missing_flag=missing_flag, lecho=lecho)


## execute raster export
if __name__ == '__main__':

    import dask
    from multiprocessing.pool import ThreadPool
    from time import time
    
    start = time()

    ## fast test config
#     project = 'SON'
#     loverwrite = True
#     start_date = '2011-01-01'; end_date = '2011-02-01'
#     grid_name  = 'son1'
    ## operational test config
#     loverwrite = True
#     start_date = '2010-11-01'; end_date = '2011-01-01'
#     grid_name  = 'son2'
    ## operational config for SON2
#     loverwrite = True
#     start_date = '2011-01-01'; end_date = None
#     grid_name  = 'son2'
    ## operational config for ASB2
    project = 'ASB'
    loverwrite = True
    start_date = '2009-12-13'; end_date = None
    grid_name  = 'asb2'
    # HGS include file
    lexec = True
    inc_file = 'precip.inc'

    ## define target data/projection
    root_folder = '{:s}/{:s}/{:s}/'.format(os.getenv('HGS_ROOT'),project,grid_name)
    raster_name = '{dataset:s}_{variable:s}_{grid:s}_{date:s}.asc'
    raster_format = 'AAIGrid'
    # projection/UTM zone
    if project.upper() == 'SON':
        # southern Ontario projection
        tgt_crs = genProj("+proj=utm +zone=17 +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs", name=grid_name)
    elif project.upper() == 'ASB':
        # Assiniboin projection
        tgt_crs = genProj("+proj=utm +zone=14 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs", name=grid_name)
    # UTM grid definition
    if grid_name == 'son1':
        tgt_size = (118,82) # lower resolution 5 km grid
        tgt_geotrans = rio.transform.Affine.from_gdal(320920.,5.e3,0,4624073.,0,5.e3) # 5 km
    elif grid_name == 'son2':
        tgt_size = (590,410) # higher resolution 1 km grid (~ 1 MB per day)
        tgt_geotrans = rio.transform.Affine.from_gdal(320920.,1.e3,0,4624073.,0,1.e3) # 1 km 
    elif grid_name == 'asb2':
        tgt_size = (955,675) # higher resolution 1 km grid (> 1 MB per day)
        tgt_geotrans = rio.transform.Affine.from_gdal(-159.e3, 1.e3, 0., 5202.e3, 0., 1.e3) # 1 km 
    
    ## define source data
    # SnoDAS
    dataset = 'SnoDAS'
    varname = 'liqwatflx'
    scalefactor = 1000. # divide by this (convert kg/m^2 to m^3/m^2, SI units to HGS internal)
#     varname = 'snow'
    filename_pattern = raster_name.format(dataset=dataset.lower(), variable=varname.lower(),
                                          grid=grid_name.lower(), date='{:s}') # no date for now...
    print("\n***   Exporting '{}' from '{}' to raster format {}   ***\n".format(varname,dataset,raster_format))
    
    ## get dataset
    ds_mod = import_module('datasets.{0:s}'.format(dataset))
    xds = ds_mod.loadDailyTimeSeries(varname=varname, time_chunks=2)
    # get georeference
    src_crs = getProj(xds)
    # figure out bounds for clipping
    left,bottom = tgt_geotrans*(0,0); right,top = tgt_geotrans*tgt_size
    trans, w,h = calculate_default_transform(src_crs=tgt_crs, dst_crs=src_crs, width=tgt_size[0], height=tgt_size[1], 
                                             left=left, bottom=bottom, right=right, top=top,)
    left,top = trans*(0,0); rigth,bottom = trans*(w,h)
    # clip source data
    xvar = xds[varname].loc[:, bottom:top, left:right]
    if start_date or end_date:
        xvar = xds[varname].loc[start_date:end_date,:,:]
    print(xvar)
    src_geotrans,src_size = getTransform(xvar)
    time_coord = xvar.coords['time'].data
    
    
    ## write rasters
    subfolder = '{:s}/snow/'.format(dataset,)
    target_folder = root_folder + subfolder
    # make sure path exists
    try:
        if not osp.exists(target_folder): os.mkdir(target_folder)
    except (WindowsError,OSError):
        os.makedirs(target_folder)
    filepath_pattern = target_folder + filename_pattern
    
    
    ## generate inc file
    start_inc = time()
    inc_filepath = target_folder+inc_file
    print("\nWriting HGS include file:\n '{:s}'\n".format(inc_filepath))
    with open(inc_filepath, mode='w') as incf:
        
        # first line
        filename = filename_pattern.format(pd.to_datetime(time_coord[0]).strftime('%Y%m%d'))
        line = '{:15d}     {:s}\n'.format(0,filename)
        incf.write(line)
        
        # loop over time coordinate
        for i,date in enumerate(time_coord):
             
            # use date to construct file name
            filename = filename_pattern.format(pd.to_datetime(date).strftime('%Y%m%d'))
            line = '{:15d}     {:s}\n'.format(43200+86400*i,filename)
            incf.write(line)
            
        # last line
        filename = filename_pattern.format(pd.to_datetime(time_coord[-1]).strftime('%Y%m%d'))
        line = '{:15d}     {:s}\n'.format(86400*len(time_coord),filename)
        incf.write(line)
        
        ## N.B.: are the first and last lines really necessary???
    end_inc = time()
    print("Timing to write include file: {} seconds\n".format(end_inc-start_inc))
        
    if not lexec: exit()
           
           
    # generate workload for lazy execution
    start_load = time()
    print("\n***   Constructing Workload from {} to {}.   ***\n".format(start_date,end_date))
    print("Output folder: '{}'\nRaster pattern: '{}'\n".format(target_folder,filename_pattern)) 
       
    # explicitly determine chunking to get complete 2D lat/lon slices
    xvar = rechunkXlonYalt(xvar)
    assert all([len(c)==1 for c in xvar.chunks[-2:]]), xvar.chunks
    time_chunks = np.concatenate([[0],np.cumsum(xvar.chunks[0][:-1], dtype=np.int)])
    
    # apply a scaling factor
    xvar /= scalefactor
    # N.B.: apply scalefactor 'in-place' so that xarray variable attributes 
    #       are preserved (it will still execute delayed); applying the scale-
    #       factor after regridding is slightly faster, but this is cleaner

    # define function to apply to blocks
    dummy = np.zeros((1,1,1), dtype=np.int8)
    def regrid_and_write(data, block_id=None):
        ''' function to apply regridding and subsequent export to raster on blocks '''
        # figure out time coordinates/dates
        ts = time_chunks[block_id[0]]; te = ts + data.shape[0] 
        time_chunk = time_coord[ts:te] # chunk of time axis that corresponds to this chunk 
        #print(time_chunk)
        # regrid array
        data = regrid_array(data, tgt_crs=tgt_crs, tgt_transform=tgt_geotrans, tgt_size=tgt_size, 
                            src_crs=src_crs, src_transform=src_geotrans, src_size=src_size, 
                            resampling='bilinear', lxarray=False)
        # write raster files
        batch_write_rasters(filepath_pattern, data, time_coord=time_chunk, crs=tgt_crs, lecho=True,
                            transform=tgt_geotrans, driver='AAIGrid', missing_value=0., missing_flag=-9999.)
        # return data
        return data
    # now map regridding operation to blocks
    n_loads = len(xvar.chunks[0])
    dummy_output = xvar.data.map_blocks(regrid_and_write, chunks=dummy.shape, dtype=dummy.dtype)
    work_load = [dummy_output]
    
    end_load = time()
    print("Timing to construct workload: {} seconds\n".format(end_load-start_load))
    
    # execute delayed computation
    print("\n***   Executing {} Workloads using Dask   ***".format(n_loads))
    print("Chunks (time only): {}\n".format(xvar.chunks[0]))

#     with dask.set_options(scheduler='processes'):      
    with dask.config.set(pool=ThreadPool(4)):    
        dask.compute(*work_load)
        
    print("\nDummy output:")
    print(dummy_output)
    print("Size in memory: {} MB\n".format(dummy_output.nbytes/1024./1024.))

    
    end = time()
    print("\n***   Completed in {:.2f} seconds   ***\n".format(end-start))
    