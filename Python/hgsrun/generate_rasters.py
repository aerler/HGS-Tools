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
from rasterio.warp import calculate_default_transform # need to import separately...
from importlib import import_module
# internal imports
from geospatial.rasterio_tools import genProj, generate_regrid_and_export
from geospatial.xarray_tools import getProj, getTransform, rechunkTo2Dslices

# WindowsError is not defined on Linux - need a dummy
try: 
    lWin = True
    WindowsError
except NameError:
    lWin = False
    WindowsError = None


# write HGS include file
def writeIncFile(filepath, time_coord, filename_pattern, time_fmt='{:15.0f}', date_fmt='%Y%m%d'):
    ''' write an HGS include file based on a time coordinate and a filename pattern;
        the date is used as part of the filename '''
    line_fmt = time_fmt+'     {:s}\n'
    # open file and write entries
    with open(filepath, mode='w') as incf:
        
        # first line
        filename = filename_pattern.format(pd.to_datetime(time_coord[0]).strftime(date_fmt))
        line = line_fmt.format(0,filename)
        incf.write(line)
        
        # construct array of seconds from simulation start
        sim_time = (time_coord - time_coord[0])/np.timedelta64(1, 's')
        sim_time[:-1] += ( np.diff(sim_time) / 2. ) # recenter between time-steps
        sim_time[-1] += ( sim_time[-1] - sim_time[-2] ) # extrapolate offset for last element
        
        # loop over time coordinate
        for stime,date in zip(sim_time,time_coord):
             
            # use date to construct file name
            filename = filename_pattern.format(pd.to_datetime(date).strftime(date_fmt))
            line = line_fmt.format(stime,filename)
            incf.write(line)
            
        # last line
        last_time = sim_time[-1] + (sim_time[-1] - sim_time[-2])/2. # the official end of validity
        filename = filename_pattern.format(pd.to_datetime(time_coord[-1]).strftime(date_fmt))
        line = line_fmt.format(last_time,filename)
        incf.write(line)
        
    ## N.B.: are the first and last lines really necessary???
    # done...
    return None


## execute raster export
if __name__ == '__main__':

    import dask
    from multiprocessing.pool import ThreadPool
    from time import time
    
    start = time()

    loverwrite = True
    time_chunks = 8 # typically not much speed-up beyond 8
    mode = 'raster2d'
    ## WRF grids
#     project = 'WRF'
# #     start_date = '2014-01-01'; end_date = '2015-01-01'
#     start_date = None; end_date = None
#     grid_name  = 'wc2_d02'    
#     project = 'CMB'
# #     start_date = '2014-01-01'; end_date = '2014-02-01'
#     start_date = None; end_date = None
#     grid_name  = 'cmb1'    
#     mode = 'NetCDF'
    ## generate a full SnoDAS raster
#     project = 'SnoDAS'
#     start_date = '2014-01-01'; end_date = '2014-01-05'
#     grid_name  = 'snodas'
    ## fast test config
    project = 'SON'
    start_date = '2013-01-01'; end_date = '2013-01-31'
    grid_name  = 'son1'
    ## operational test config
#     project = 'SON'
#     start_date = '2010-11-01'; end_date = '2011-01-01'
#     grid_name  = 'son2'
    ## operational config for SON2
#     project = 'SON'
#     start_date = '2011-01-01'; end_date = None
#     grid_name  = 'son2'
    ## operational config for ASB2
#     project = 'ASB'
#     start_date = '2010-01-01'; end_date = None
#     grid_name  = 'asb2'

    ## define target grid/projection
    resampling = 'cubic'
    # projection/UTM zone
    tgt_size = None; tgt_geotrans = None # valid for native grid
    if project == 'WRF':
        # load pickled GridDef
        from geodata.gdal import loadPickledGridDef
        griddef = loadPickledGridDef(grid=grid_name, encoding='latin1')
        print(griddef)
        tgt_crs = genProj(griddef.projection.ExportToProj4(), name=grid_name)
        tgt_geotrans = griddef.geotransform; tgt_size = griddef.size
    elif project == 'SnoDAS':
        tgt_crs = None # native grid
    elif project.upper() == 'SON':
        # southern Ontario projection
        tgt_crs = genProj("+proj=utm +zone=17 +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs", name=grid_name)
    elif project.upper() == 'CMB':
        # southern Ontario projection
        tgt_crs = genProj("+proj=utm +zone=11 +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs", name=grid_name)
    elif project.upper() == 'ASB':
        # Assiniboin projection
        tgt_crs = genProj("+proj=utm +zone=14 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs", name=grid_name)
    # grid definition (mostly UTM grids for HGS)
    if tgt_geotrans is not None and tgt_size is not None:
        pass # already assigned above
    elif grid_name == 'son1':
        tgt_size = (118,82) # lower resolution 5 km grid
        tgt_geotrans = (320920.,5.e3,0,4624073.,0,5.e3) # 5 km
    elif grid_name == 'cmb1':
        tgt_size = (640,826) # higher resolution 500 m grid
        tgt_geotrans = (292557.,500,0,5872251.,0,-500.) # 500 m 
    elif grid_name == 'son2':
        tgt_size = (590,410) # higher resolution 1 km grid (~ 1 MB per day)
        tgt_geotrans = (320920.,1.e3,0,4624073.,0,1.e3) # 1 km 
    elif grid_name == 'asb2':
        tgt_size = (955,675) # higher resolution 1 km grid (> 1 MB per day)
        tgt_geotrans = (-159.e3, 1.e3, 0., 5202.e3, 0., 1.e3) # 1 km 
    elif grid_name == 'snodas': # original SnoDAS grid
        time_chunks = 1 # this is pretty big!
        tgt_size = None; tgt_geotrans = None # native grid
    # cast geotransform into an Affine object
    if tgt_geotrans is not None:
        tgt_geotrans = rio.transform.Affine.from_gdal(*tgt_geotrans)
    
    
    ## define source data
    # SnoDAS
    dataset = 'SnoDAS'
    ds_mod = import_module('datasets.{0:s}'.format(dataset))
    
    ## define export parameters
    raster_format = None; scalefactor = 1.; lexec = True
    # modes
    if mode.lower() == 'raster2d':
        # raster output using rasterio
        #varname = 'liqwatflx'; scalefactor = 1000. # divide by this (convert kg/m^2 to m^3/m^2, SI units to HGS internal)
        varlist = ['liqwatflx','precip']
        target_folder = '{root:s}/{proj:s}/{grid:s}/{name:s}/'.format(root=os.getenv('HGS_ROOT'),
                                                                      proj=project,grid=grid_name,name=dataset)
        raster_format = 'AAIGrid'
        raster_name = '{dataset:s}_{variable:s}_{grid:s}_{date:s}.asc'
        filename_novar = raster_name.format(dataset=dataset.lower(), variable='{var:s}',
                                            grid=grid_name.lower(), date='{date:s}') # no date for now...
        print(("\n***   Exporting '{}' to raster format {}   ***\n".format(dataset,raster_format)))
        # HGS include file only?
        lexec = True        
    elif mode.upper() == 'NETCDF':
        # NetCDF output using netCDF4
        #varlist = ['snow','liqwatflx']
        varlist = ds_mod.netcdf_varlist # all primary and secondary variables, excluding coordinate variables
        target_folder = osp.join(ds_mod.daily_folder,grid_name,resampling)
        filename_novar = ds_mod.netcdf_filename.format('{:s}_{:s}'.format('{var:s}',grid_name))
        print(("\n***   Regridding '{}' to '{}' (NetCDF format)   ***".format(dataset,grid_name)))
    else:
        raise NotImplementedError
    print(("   Variable list: {}\n".format(str(varlist))))
    
    
    # lazily load dataset (assuming xarray)
    xds = ds_mod.loadDailyTimeSeries(varlist=varlist, time_chunks=time_chunks)
    
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
            left,top = trans*(0,0); right,bottom = trans*(w,h)
        # clip source data
        space_slice = {xds.attrs['xlon']:slice(left,right), xds.attrs['ylat']:slice(bottom,top)}
        xds = xds.loc[space_slice]
    # clip time axis as well
    if start_date or end_date:
        time_slice = {'time':slice(start_date,end_date)}
        xds = xds.loc[time_slice]
    print(xds)
    src_geotrans,src_size = getTransform(xds)
    if tgt_geotrans is None: tgt_geotrans = src_geotrans
    if tgt_size is None: tgt_size = src_size
    time_coord = xds.coords['time'].data
    if start_date is None: start_date = pd.to_datetime(time_coord[0]).strftime('%Y-%m-%d')
    if end_date is None: end_date = pd.to_datetime(time_coord[0]).strftime('%Y-%m-%d')
    
    
    # make sure target path exists
    try:
        if not osp.exists(target_folder): os.mkdir(target_folder)
    except (WindowsError,OSError):
        os.makedirs(target_folder)
   
    
    ## loop over variables (for rasters, that will typically not be necessary)
    
    for varname in varlist:  
        
        start_var = time()
        
        # select variable
        xvar = xds[varname]
        filename = filename_novar.format(var=varname.lower(),date='{:s}')
        
        print(("\n\n###   Processing Variable '{:s}'   ###".format(varname,start_date,end_date)))
                
        ## generate inc file
        if mode.lower() == 'raster2d':
            start_inc = time()
            inc_filepath = target_folder+varname.lower()+'.inc'
            print(("\nWriting HGS include file:\n '{:s}'".format(inc_filepath)))
            writeIncFile(filepath=inc_filepath, time_coord=time_coord, 
                          filename_pattern=filename, date_fmt='%Y%m%d')
            end_inc = time()
            #print("\nTiming to write include file: {} seconds".format(end_inc-start_inc))
            
        if not lexec: exit()
               
        ## generate workload for lazy execution
        start_load = time()
        print(("\nConstructing Workload for '{:s}' from {:s} to {:s}.   ***".format(varname,start_date,end_date)))
        if mode.lower() == 'raster2d':
            print(("Output folder: '{:s}'\nRaster pattern: '{:s}'".format(target_folder,filename)))
        elif mode.upper() == 'NETCDF':
            print(("NetCDF file: '{:s}'".format(osp.join(target_folder,filename))))        
         
        # explicitly determine chunking to get complete 2D lat/lon slices
        xvar = rechunkTo2Dslices(xvar)    
          
        # apply a scaling factor
        xvar /= scalefactor
        # N.B.: apply scalefactor 'in-place' so that xarray variable attributes 
        #       are preserved (it will still execute delayed); applying the scale-
        #       factor after regridding is slightly faster, but this is cleaner
    
       
        # generate dask execution function
        dask_fct,dummy,dataset = generate_regrid_and_export(xvar, time_coord=time_coord,
                                                    tgt_crs=tgt_crs, tgt_geotrans=tgt_geotrans, tgt_size=tgt_size, 
                                                    mode=mode, resampling=resampling, 
                                                    folder=target_folder, filename=filename, driver=raster_format,
                                                    lecho=True, loverwrite=loverwrite,)
                  
        # now map regridding operation to blocks
        n_loads = len(xvar.chunks[0])
        dummy_output = xvar.data.map_blocks(dask_fct, chunks=dummy.shape, dtype=dummy.dtype)
        work_load = [dummy_output]
        
        end_load = time()
        #print("\nTiming to construct workload: {:.2f} seconds".format(end_load-start_load))
        
        
        # execute delayed computation
        print(("\n***   Executing {:d} Workloads for '{:s}' using Dask   ***".format(n_loads,varname)))
        print(("Chunks (time only): {}".format(xvar.chunks[0])))
    
    #     with dask.set_options(scheduler='processes'):      
        with dask.config.set(pool=ThreadPool(2)):    
            dask.compute(*work_load)
            
        #print("\nDummy output:")
        #print(dummy_output)
        print(("\nDummy Size in memory: {:f} MB".format(dummy_output.nbytes/1024./1024.)))
    
        if mode.upper() == 'NETCDF':
            dataset.close()
        
        end_var = time()
        print(("\n\n***   Completed '{:s}' in {:.2f} seconds   ***\n".format(varname,end_var-start_var)))
        
    end = time()
    print(("\n***   Overall Timing: {:.2f} seconds   ***\n".format(end-start)))

        