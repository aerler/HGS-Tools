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
from rasterio.warp import calculate_default_transform
from importlib import import_module
# internal imports
from geospatial.rasterio_tools import genProj,regrid_array, write_time_rasters
from geospatial.xarray_tools import getProj, getTransform, rechunkTo2Dslices

# WindowsError is not defined on Linux - need a dummy
try: 
    lWin = True
    WindowsError
except NameError:
    lWin = False
    WindowsError = None


## execute raster export
if __name__ == '__main__':

    import dask
    from multiprocessing.pool import ThreadPool
    from time import time
    
    start = time()

    ## generate a full SnoDAS raster
#     project = 'SnoDAS'
#     loverwrite = True
#     start_date = '2014-01-01'; end_date = '2014-01-05'
#     grid_name  = 'snodas'
    ## fast test config
    project = 'SON'
    loverwrite = True
    start_date = '2011-01-01'; end_date = '2011-02-01'
    grid_name  = 'son1'
    ## operational test config
#     loverwrite = True
#     start_date = '2010-11-01'; end_date = '2011-01-01'
#     grid_name  = 'son2'
    ## operational config for SON2
#     loverwrite = True
#     start_date = '2011-01-01'; end_date = None
#     grid_name  = 'son2'
    ## operational config for ASB2
#     project = 'ASB'
#     loverwrite = True
#     start_date = '2010-01-01'; end_date = None
#     grid_name  = 'asb2'
    # HGS include file
    lexec = True
    inc_file = 'precip.inc'

    ## define target data/projection
    root_folder = '{:s}/{:s}/{:s}/'.format(os.getenv('HGS_ROOT'),project,grid_name)
    raster_name = '{dataset:s}_{variable:s}_{grid:s}_{date:s}.asc'
    raster_format = 'AAIGrid'
    # projection/UTM zone
    if project == 'SnoDAS':
        tgt_crs = None # native grid
    elif project.upper() == 'SON':
        # southern Ontario projection
        tgt_crs = genProj("+proj=utm +zone=17 +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs", name=grid_name)
    elif project.upper() == 'ASB':
        # Assiniboin projection
        tgt_crs = genProj("+proj=utm +zone=14 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs", name=grid_name)
    # UTM grid definition
    if grid_name == 'son1':
        tgt_size = (118,82) # lower resolution 5 km grid
        tgt_geotrans = (320920.,5.e3,0,4624073.,0,5.e3) # 5 km
    elif grid_name == 'son2':
        tgt_size = (590,410) # higher resolution 1 km grid (~ 1 MB per day)
        tgt_geotrans = (320920.,1.e3,0,4624073.,0,1.e3) # 1 km 
    elif grid_name == 'asb2':
        tgt_size = (955,675) # higher resolution 1 km grid (> 1 MB per day)
        tgt_geotrans = (-159.e3, 1.e3, 0., 5202.e3, 0., 1.e3) # 1 km 
    elif grid_name == 'snodas': # original SnoDAS grid
        tgt_size = None; tgt_geotrans = None # native grid
    # cast geotransform into an Affine object
    if tgt_geotrans is not None:
        tgt_geotrans = rio.transform.Affine.from_gdal(*tgt_geotrans)
    
    ## define source data
    # SnoDAS
    dataset = 'SnoDAS'
#     varname = 'liqwatflx'; scalefactor = 1000. # divide by this (convert kg/m^2 to m^3/m^2, SI units to HGS internal)
    varname = 'snow'; scalefactor = 1.
    filename_pattern = raster_name.format(dataset=dataset.lower(), variable=varname.lower(),
                                          grid=grid_name.lower(), date='{:s}') # no date for now...
    print("\n***   Exporting '{}' from '{}' to raster format {}   ***\n".format(varname,dataset,raster_format))
    
    ## get dataset
    ds_mod = import_module('datasets.{0:s}'.format(dataset))
    xds = ds_mod.loadDailyTimeSeries(varname=varname, time_chunks=2)
    # get georeference
    src_crs = getProj(xds)
    if tgt_crs is None: tgt_crs = src_crs
    # figure out bounds for clipping
    if tgt_geotrans is not None and tgt_size is not None:
        # get outline
        left,bottom = tgt_geotrans*(0,0); right,top = tgt_geotrans*tgt_size
        if tgt_crs != src_crs:
            # reproject outline, if necessary
            trans, w,h = calculate_default_transform(src_crs=tgt_crs, dst_crs=src_crs, 
                                                     width=tgt_size[0], height=tgt_size[1], 
                                                     left=left, bottom=bottom, right=right, top=top,)
            left,top = trans*(0,0); rigth,bottom = trans*(w,h)
        # clip source data
        xvar = xds[varname].loc[:, bottom:top, left:right]
    # clip time axis as well
    if start_date or end_date:
        xvar = xds[varname].loc[start_date:end_date,:,:]
    print(xvar)
    src_geotrans,src_size = getTransform(xvar)
    if tgt_geotrans is None: tgt_geotrans = src_geotrans
    if tgt_size is None: tgt_size = src_size
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
    xvar = rechunkTo2Dslices(xvar)
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
        if tgt_crs != src_crs:
            # reproject data
            data = regrid_array(data, tgt_crs=tgt_crs, tgt_transform=tgt_geotrans, tgt_size=tgt_size, 
                                src_crs=src_crs, src_transform=src_geotrans, src_size=src_size, 
                                resampling='bilinear', lxarray=False)

        # write raster files
        write_time_rasters(filepath_pattern, data, time_coord=time_chunk, driver='AAIGrid', crs=tgt_crs,
                            transform=tgt_geotrans, missing_value=0., missing_flag=-9999.,
                            lecho=True, loverwrite=loverwrite)
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
    