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
from geospatial.rasterio_tools import genCRS, generate_regrid_and_export
from geospatial.xarray_tools import getCRS, getTransform, rechunkTo2Dslices

# WindowsError is not defined on Linux - need a dummy
try: 
    lWin = True
    WindowsError
except NameError:
    lWin = False
    class WindowsError(BaseException): pass


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
    # N.B.: it seems dask multi-processing is prone to memory leaks...
    #from multiprocessing.pool import ThreadPool
    from time import time
    
    start = time()

    loverwrite = True
    lhourly = False
    time_chunks = 1 # typically not much speed-up beyond 8
    resampling = 'bilinear'
    lexec = True # actually write rasters or just include file
    # WRF grids
#     project = 'WRF'
#     #grid_name  = 'wc2_d01'    
#     #project = 'CMB'
#     project = 'ARB'
#     grid_name  = 'arb3'    
    ## Fraser's Ontario domain
#     project = 'WRF' # load grid from pickle
# #     grid_name = 'glb1_d02'    
#     start_date = None; end_date = None
#     grid_name = 'glb1_d01'    
#     project = 'Geo'
#     grid_name = 'on1'
#     start_date = None; end_date = None
#     grid_name  = 'cmb1'    
#     mode = 'NetCDF'
    ## generate a full SnoDAS raster
#     project = 'native'
#     grid_name  = 'native'
    ## fast test config
#     project = 'SON'
#     grid_name  = 'son1'
    ## config for Hugo's domain in Quebec
#     project = 'Hugo'
#     grid_name = 'hd1'
#     mode = 'NetCDF'
    ## operational config for GRW
#     project = 'GRW'
#     grid_name  = 'grw1'
    ## test config for GRW
#     project = 'GRW'
#     grid_name  = 'grw2'; resampling = 'nearest'; #source_grid = 'grw1'
    ## operational config for SON2
    project = 'SON'
    grid_name  = 'son2'
    ## 
#     project = 'SNW'
#     grid_name  = 'snw2'
    ## operational config for ASB2
#     project = 'ASB'
#     grid_name  = 'asb1'
#     #grid_name  = 'asb2'
    ## CA12 NRCan grid
#     project = 'Geo'
#     grid_name = 'ca12'
    ## Queensland (Australia) grid
#     project = 'QEL'
#     grid_name = 'qel1' # 10 km Queensland grid

    ## define target grid/projection
    # projection/UTM zone
    tgt_size = None; tgt_geotrans = None # valid for native grid
    if project == 'WRF':
        # load pickled GridDef
        from geodata.gdal import loadPickledGridDef
        griddef = loadPickledGridDef(grid=grid_name, encoding='latin1')
        print(griddef)
        tgt_crs = genCRS(griddef.projection.ExportToProj4(), name=grid_name)
        tgt_geotrans = griddef.geotransform; tgt_size = griddef.size
    elif project == 'Geo':
        # generic geographic lat/lon
        tgt_crs = genCRS(name=grid_name)
    elif project == 'ARB':
        # Projection for ARB model
        tgt_crs = genCRS('+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +ellps=sphere +units=m +no_defs', name=grid_name)
    elif project == 'Hugo':
        # Hugo's projection for Quebec
        tgt_crs = genCRS("+proj=lcc +ellps=NAD83 +datum=NAD83 +lat_0=44.0 +lat_1=46.0 +lat_2=60.0 +lon_0=-68.5  +x_0=0 +y_0=0 +units=m +no_defs", name=grid_name)
        # N.B.: something seems to be wrong with this... rasterio can't read it...
    elif project.upper() in ('SON','GRW'):
        # southern Ontario projection
        tgt_crs = genCRS("+proj=utm +zone=17 +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs", name=grid_name)
    elif project.upper() in ('SNW'):
        # South Nation projection
        tgt_crs = genCRS("+proj=utm +zone=18 +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs", name=grid_name)
    elif project.upper() == 'CMB':
        # Columbia Mass Balance projection
        tgt_crs = genCRS("+proj=utm +zone=11 +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs", name=grid_name)
    elif project.upper() == 'ASB':
        # Assiniboin projection
        tgt_crs = genCRS("+proj=utm +zone=14 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs", name=grid_name)
    elif project.upper() in ('QEL'):
        # Queensland (Australia) projection
        tgt_crs = genCRS("+proj=utm +zone=55 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs", name=grid_name)
    else:
        tgt_crs = None # native grid
    # grid definition (mostly UTM grids for HGS)
    if tgt_geotrans is not None and tgt_size is not None:
        pass # already assigned above
    elif grid_name == 'ca12': # the NRCan 1/12t deg Canada grid
        tgt_geotrans = (-141.0, 1./12., 0.0, 41.0, 0.0, 1./12.); tgt_size = (1068, 510) # (x,y) map size of NRCan grid
        resampling = 'cubic_spline'
    elif grid_name == 'on1':
        tgt_geotrans = [-87.87916564941406,0.008331298828125,0.0,41.995832443237305,0.0,0.008335113525390625,]
        resampling = 'cubic_spline'
    elif grid_name == 'arb2':
        tgt_geotrans = [-1460500,5e3,0,810500,0,5e3]; tgt_size = (284,258)
        resampling = 'average' # it's a fairly coarse grid...
    elif grid_name == 'arb3':
        tgt_geotrans = [-1280e3,5e3,0,900e3,0,5e3]; tgt_size = (172,144)
        resampling = 'average' # it's a fairly coarse grid...
    elif grid_name == 'hd1':
        tgt_size = (70,49) # lower resolution 5 km grid
        tgt_geotrans = (-479184.769227,5.e3,0,68508.4877898,0,5.e3) # 5 km
        resampling = 'average' # it's a fairly coarse grid...
    elif grid_name == 'son1':
        tgt_size = (118,82) # lower resolution 5 km grid
        tgt_geotrans = (320920.,5.e3,0,4624073.,0,5.e3) # 5 km
        resampling = 'average' # it's a fairly coarse grid...
    elif grid_name == 'son2':
        tgt_size = (590,410) # higher resolution 1 km grid (~ 1 MB per day)
        tgt_geotrans = (320920.,1.e3,0,4624073.,0,1.e3) # 1 km 
    elif grid_name == 'grw1':
        tgt_size = (132,162) # smaller, higher resolution 1 km grid for GRW
        tgt_geotrans = (500.e3,1.e3,0,4740.e3,0,1.e3) # 1 km 
    elif grid_name == 'grw2':
        tgt_size = (27,33) # smaller, lower resolution 5 km grid for GRW
        tgt_geotrans = (500.e3,5.e3,0,4740.e3,0,5.e3) # 5 km 
        resampling = 'average' # it's a fairly coarse grid...
    elif grid_name == 'snw2':
        tgt_size = (44,55) # 2 km resolution SNW grid for CaLDAS/CaPA data
        tgt_geotrans = (438.e3,2.e3,0,4940.e3,0,2.e3) # 2 km 
        resampling = 'average' # it's a fairly coarse grid...
    elif grid_name == 'cmb1':
        tgt_size = (640,826) # higher resolution 500 m grid
        tgt_geotrans = (292557.,500,0,5872251.,0,-500.) # 500 m 
    elif grid_name == 'asb1':
        tgt_size = (191,135) # lower resolution 5 km grid
        tgt_geotrans = (-159.e3, 5.e3, 0., 5202.e3, 0., 5.e3) # 5 km
    elif grid_name == 'asb2':
        tgt_size = (955,675) # higher resolution 1 km grid (> 1 MB per day)
        tgt_geotrans = (-159.e3, 1.e3, 0., 5202.e3, 0., 1.e3) # 1 km
    elif grid_name == 'qel1':
        tgt_size = (121,191) # 10 km grid, similar to ERA5
#         tgt_geotrans = (-26770., 10400, 0., 6902510, 0., 10970) # ~10 km
        tgt_geotrans = (-26772., 10398.9, 0., 6902511, 0., 10968.7) # ~10 km  
        # N.B.: this grid attempts to achieve direct grid point correspondence to the ERA5-Land lat/lon grid
        #       at 148E/23S
    elif grid_name == 'native': # original grid
        time_chunks = 1 # this can be pretty big!
        tgt_size = None; tgt_geotrans = None # native grid
    else:
        raise NotImplementedError(grid_name)
    # cast geotransform into an Affine object
    if tgt_geotrans is not None:
        tgt_geotrans = rio.transform.Affine.from_gdal(*tgt_geotrans)
    
    
    ## define source data
    lexec = True
    loverwrite = True  
    lwarp = True # set False to suppress reprojection    
    # some defaults for most datasets
    time_chunks = 8 # SnoDAS & CaSPAr only! typically not much speed-up beyond 8
    dataset_kwargs = dict(); subdataset = None 
    bias_correction = None; bc_varmap = dict(); obs_name = None; bc_method = None
    target_folder_ascii = '{root:s}/{proj:s}/{grid:s}/{name:s}/{bc:s}transient_{int:s}/'
    raster_name = '{dataset:s}_{variable:s}_{grid:s}_{date:s}.asc'
    target_folder_netcdf = '{daily:s}/{grid:s}/{smpl:s}/'
    
    # dummy variables for WRF
    exp_name = None; domain = None; WRF_exps = None; filetype = None
    
    ## SnoDAS        
#     dataset = 'SnoDAS' # default...
# #     bias_correction = 'SMBC'; obs_name = 'NRCan' 
#     bc_varmap = dict(liqprec='liqwatflx') # just for testing...
#     varlist = ['snow']; bc_method = 'rfbc'
#     dataset_kwargs = dict(grid='on1', bias_correction=bc_method)
#     resampling = 'bilinear'
#     end_date = '2011-02-01'
    ## CaSPAr
    #dataset = 'CaSPAr'; lhourly = True; dataset_kwargs = dict(grid='lcc_snw')
    ## MergedForcing
#     varlist = []
#     dataset = 'MergedForcing'
#     subdataset = dataset; varlist = ['liqwatflx_ne5']
# #     subdataset = dataset; varlist = ['liqwatflx_ne5','pet_har','pet_hog',]
# # #     subdataset = dataset; varlist = ['pet_pts',] # PET based on Priestley-Taylor with solar radiation only
# # #     subdataset = dataset; varlist = ['liqwatflx','pet_pts',] # assorted forcing
# # #     subdataset = 'NRCan'; varlist += ['Tmin',] # base variables
# # #     subdataset = 'NRCan'; varlist += ['precip','Tmin','Tmax','T2',] # base variables
# # #     subdataset = 'NRCan'; varlist += ['precip_adj',] # adjusted precip data (up to 2016)
# #     subdataset = 'NRCan'; varlist += ['pet_har',] # PET based on Hargreaves' method
# # #     subdataset = 'NRCan'; varlist += ['pet_haa',] # PET based on Hargreaves' with Allen's correction
# # #     subdataset = 'NRCan'; varlist += ['pet_th',] # PET based on Thornthwaite method
# #     subdataset = 'NRCan'; varlist += ['pet_hog',] # PET based on simple Hogg method
# #     dataset_kwargs = dict(dataset=subdataset)
#     dataset_kwargs['resolution'] = 'CA12'; resampling = 'cubic_spline'; dataset_kwargs['chunks'] = dict(time=8, lon=63, lat=64)
# #     dataset_kwargs['resolution'] = 'SON60'; resampling = 'bilinear'
#     dataset_kwargs['dataset_index'] = dict(liqwatflx='MergedForcing',liqwatflx_ne5='MergedForcing',pet_har='NRCan',pet_hog='NRCan')
# #     dataset_kwargs['grid'] = 'son2'; resampling = None; lwarp = False
# #     dataset_kwargs['grid'] = 'snw2'; resampling = None; lwarp = False
    ## ERA5
    dataset = 'ERA5'; subdataset = 'ERA5L'
#     varlist = ['snow','dswe',]
#     varlist = ['precip','pet_era5','liqwatflx','snow','dswe',]
    varlist = ['pet_era5','liqwatflx','snow']
    dataset_kwargs = dict(filetype=subdataset)
#     dataset_kwargs['resolution'] = 'NA10'; chunks = dict(time=8,latitude=61,longitude=62)
    dataset_kwargs['resolution'] = 'NA10'; dataset_kwargs['grid'] = 'son2'; chunks = dict(time=9,y=59,x=59)
#     dataset_kwargs['resolution'] = 'AU10'; chunks = dict(time=8, latitude=59, longitude=62)  
    resampling = 'cubic_spline'; dataset_kwargs['chunks'] = chunks # apparently we need to pre-chunk or there is a memory leak..   
    

#     start_date = '1997-01-01'; end_date = '2017-12-31' # SON/SNW full period
#     start_date = '1981-01-01'; end_date = '2017-12-31' # SON/SNW full period
    start_date = '1981-01-01'; end_date = '2020-08-31' # full ERA5-Land period
#     start_date = '2000-01-01'; end_date = '2018-01-01'
#     start_date = '2011-01-01'; end_date = '2017-12-31' # combined NRCan-SnoDAS period
#     start_date = '2016-01-01'; end_date = '2017-12-31'
#     start_date = '2016-01-01'; end_date = '2016-01-31' # for testing NRCan
#     start_date = '1997-01-01'; end_date = '1997-02-01'; resampling = 'nearest' # for testing...

    ## output type: ASCII raster or NetCDF-4
#     mode = 'NetCDF'
    mode = 'raster2d'

    
    ## WRF requires special treatment
#     dataset = 'WRF';  lhourly = False; bias_correction = None; resampling = 'bilinear'
#     if project in ('ARB','CMB','ASB'): from projects.WesternCanada import WRF_exps
#     else: from projects.GreatLakes import WRF_exps
#     exp_name = os.getenv('WRFEXP')
# #     exp_name = 'max-ctrl'
# #     exp_name = 'ctrl-1'  
#     domain = 2; filetype = 'hydro'
#     dataset_kwargs = dict(experiment=exp_name, domain=domain, filetypes=filetype, exps=WRF_exps)
# #     start_date = '1979-01-01'; end_date = '1979-12-31'
#     start_date = None; end_date = None    
#     target_folder_ascii = '{root:s}/{proj:s}/{grid:s}/{exp_name:s}_d{dom:0=2d}/{bc:s}transient_{int:s}/climate_forcing/'
#     target_folder_netcdf = '{exp_folder:s}/{grid:s}/{smpl:s}/'  
# #     bias_correction = 'MyBC'; bc_varmap = dict(liqwatflx=None); obs_name = 'CRU'
#     bias_correction = 'AABC'; bc_varmap = dict(liqwatflx='precip'); obs_name = 'CRU'
#     varlist = ['liqwatflx','pet',]
    
    # import dataset module
    ds_mod = import_module('datasets.{0:s}'.format(dataset))
    
    # get some experiment/dataset info
    exp_folder = None; exp = None; bc_folder = None
    if dataset == 'WRF':
        exp_folder,exp,exp_name,domain = ds_mod.getFolderNameDomain(experiment=exp_name, domains=domain, exps=WRF_exps, lreduce=True)
        print('{exp_name:s}_d{dom:0=2d}'.format(exp_name=exp_name,dom=domain))
        netcdf_name = ds_mod.fileclasses[filetype].dailyfile.format(domain, '{grid:s}')
        bc_folder = exp_folder
        daily_folder = ds_mod.daily_folder
    elif dataset == 'MergedForcing':
        daily_folder,netcdf_name =ds_mod.getFolderFileName(varname='{var_str:s}', mode='daily', aggregation='daily', **dataset_kwargs)
        print(daily_folder,netcdf_name)
    else:
        if hasattr(ds_mod,'getFolderFileName',):
            daily_folder,netcdf_name =ds_mod.getFolderFileName(varname='{var_str:s}', dataset=dataset, mode='daily', aggregation='daily', **dataset_kwargs)
        else: 
            netcdf_name = ds_mod.netcdf_filename.format('{var_str:s}', VAR='{var_str:s}')
            daily_folder = ds_mod.daily_folder
        bc_folder = ds_mod.avgfolder
        
    ## bias correction
    if bias_correction:
        assert bc_method is None, "'bc_method' can only be set manually as a varname extension,  if no explicit bias_correction is applied"
        # load pickle from file
        from processing.bc_methods import loadBCpickle
        bc_method = bias_correction
        bias_correction = loadBCpickle(method=bias_correction, obs_name=obs_name, gridstr=grid_name, domain=domain, folder=bc_folder)
        
    ## define export parameters
    driver_args = dict(); scalefactor = 1.; raster_format = None
    # modes
    if mode.lower() == 'raster2d':
        # raster output using rasterio
        scalefactor = 1000. # divide by this (convert kg/m^2 to m^3/m^2, SI units to HGS internal)
        gridstr = dataset.lower() if grid_name.lower() == 'native' else grid_name.lower()
        bc_str = bc_method+'_' if bc_method else ''
        time_interval = 'hourly' if lhourly else 'daily'
        hgs_root = os.getenv('HGS_ROOT', os.getenv('DATA_ROOT')+'HGS/')
        target_folder = target_folder_ascii.format(root=hgs_root, proj=project, grid=gridstr, name=dataset, int=time_interval, 
                                                   bc=bc_str, exp_name=exp_name, dom=domain, exp_folder=exp_folder)          
        raster_format = 'AAIGrid'        
        filename_novar = raster_name.format(dataset=dataset.lower(), variable='{var:s}',
                                            grid=grid_name.lower(), date='{date:s}') # no date for now...
        driver_args = dict(significant_digits=4, fill_value=0., nodata_flag=-9999.)
        print(("\n***   Exporting '{}' to raster format {}   ***\n".format(dataset,raster_format)))
    elif mode.upper() == 'NETCDF':
        # NetCDF output using netCDF4
        gridstr = '' if grid_name.lower() == 'native' else '_'+grid_name.lower()
        bc_str1 = bc_method+'_' if bc_method else ''
        bc_str2 = '_'+bc_method if bc_method else ''
        target_folder = target_folder_netcdf.format(daily=daily_folder,grid=grid_name,smpl=resampling,
                                                    exp_name=exp_name, dom=domain, exp_folder=exp_folder)        
        filename_novar = netcdf_name.format(var_str=bc_str1+'{var:s}'+gridstr, grid=bc_str2.lower()+gridstr,) # usually either var or grid
        #driver_args = dict(least_significant_digit=4)
        print(("\n***   Regridding '{}' to '{}' (NetCDF format)   ***".format(dataset,grid_name)))
    else:
        raise NotImplementedError
    print(("   Variable list: {}\n".format(str(varlist))))
    
    
    # lazily load dataset (assuming xarray)
#     if dataset == 'MergedForcing':
#         xds = getattr(ds_mod,'load{}_Daily'.format(subdataset))(varlist=varlist, **dataset_kwargs)
    if dataset == 'SnoDAS':
        xds = ds_mod.loadDailyTimeSeries(varlist=varlist, time_chunks=time_chunks, **dataset_kwargs)
    elif lhourly:
        xds = ds_mod.loadHourlyTimeSeries(varlist=varlist, time_chunks=time_chunks, **dataset_kwargs)
    else:
        xds = ds_mod.loadDailyTimeSeries(varlist=varlist, **dataset_kwargs)
    
    # get georeference
    src_crs = getCRS(xds)
    src_geotrans,src_size = getTransform(xds)
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
            left,top = trans*(-1,-1); right,bottom = trans*(w+1,h+1) # need some padding to avoid margins 
        # clip source data
        if src_geotrans.e < 0: bottom,top = top,bottom
        space_slice = {xds.attrs['xlon']:slice(left,right), xds.attrs['ylat']:slice(bottom,top)}
        print(space_slice)
        xds = xds.loc[space_slice]
    # clip time axis as well
    if start_date or end_date:
        time_slice = {'time':slice(start_date,end_date)}
        xds = xds.loc[time_slice]
    print(xds)
    src_geotrans,src_size = getTransform(xds) # recalculate after clipping
    if tgt_geotrans is None: tgt_geotrans = src_geotrans
    if tgt_size is None: tgt_size = src_size
    time_coord = xds.coords['time'].data
    # time format
    if lhourly: date_fmt = '%Y%m%d%H'; date_longfmt = '%Y-%m-%dT%H'
    else: date_fmt = '%Y%m%d'; date_longfmt = '%Y-%m-%d'
    if start_date is None: start_date = pd.to_datetime(time_coord[0]).strftime(date_longfmt)
    if end_date is None: end_date = pd.to_datetime(time_coord[-1]).strftime(date_longfmt)
    
    
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
                          filename_pattern=filename, date_fmt=date_fmt)
            end_inc = time()
            #print("\nTiming to write include file: {} seconds".format(end_inc-start_inc))
            
        if not lexec: 
            print("\nNot executing workload --- set 'lexec' to 'True' to execute workload")
            exit()
               
        ## generate workload for lazy execution
        start_load = time()
        print(("\nConstructing Workload for '{:s}' from {:s} to {:s}.   ***".format(varname,start_date,end_date)))
        if mode.lower() == 'raster2d':
            print(("Output folder: '{:s}'\nRaster pattern: '{:s}'".format(target_folder,filename)))
        elif mode.upper() == 'NETCDF':
            print(("NetCDF file: '{:s}'".format(osp.join(target_folder,filename))))        
        
        
        # explicitly determine chunking to get complete 2D lat/lon slices
        xvar = rechunkTo2Dslices(xvar, time=time_chunks)    
          
        # apply a scaling factor
        xvar /= scalefactor
        # N.B.: apply scalefactor 'in-place' so that xarray variable attributes 
        #       are preserved (it will still execute delayed); applying the scale-
        #       factor after regridding is slightly faster, but this is cleaner
        
        # Bias-correction parameters
        if bias_correction:
            bc_varname = bc_varmap.get(varname,varname)
            bc_object = None if bc_varname is None else bias_correction
        else:
            bc_object = None; bc_varname = None 
        
        
        
        # generate dask execution function
        dask_fct,dummy,dataset = generate_regrid_and_export(xvar, time_coord=time_coord, lwarp=lwarp,
                                                    tgt_crs=tgt_crs, tgt_geotrans=tgt_geotrans, tgt_size=tgt_size, 
                                                    mode=mode, resampling=resampling, time_fmt=date_fmt,
                                                    folder=target_folder, filename=filename, driver=raster_format,
                                                    bias_correction=bc_object, bc_varname=bc_varname,
                                                    lecho=True, loverwrite=loverwrite, **driver_args)
        # N.B.: the dataset returned here is a NetCDF dataset, not a xarray dataset!
        
        # switch of overwrite/deletion if filename is not variable-specific
        if varname.lower() not in filename: 
            loverwrite = False
            assert mode.upper() == 'NETCDF' 
        # N.B.: this allows writing of all variables in the varlist to a single file
                  
        # now map regridding operation to blocks
        n_loads = len(xvar.chunks[0])
        dummy_output = xvar.data.map_blocks(dask_fct, chunks=dummy.shape, dtype=dummy.dtype)
        work_load = [dummy_output]
        
        end_load = time()
        #print("\nTiming to construct workload: {:.2f} seconds".format(end_load-start_load))
        
        
        # execute delayed computation
        print(("\n***   Executing {:d} Workloads for '{:s}' using Dask   ***".format(n_loads,varname)))
        print(("Chunks (time only): {}".format(xvar.chunks[0])))
    
        #with dask.set_options(scheduler='processes'):      
        #with dask.config.set(pool=ThreadPool(4)):    
        dask.compute(*work_load)
            
        #print("\nDummy output:")
        #print(dummy_output)
        print(("\nDummy Size in memory: {:f} MB".format(dummy_output.nbytes/1024./1024.)))
    
        if mode.upper() == 'NETCDF':
            dataset.setncattr('resampling',resampling) # netCDF4 dataset, not xarray
            xlon = dataset.getncattr('xlon')
            ylat = dataset.getncattr('ylat')
            for var in dataset.variables.values():
                if  xlon in var.dimensions and ylat in var.dimensions:
                    var.setncattr('resampling',resampling)
            dataset.close()
        
        end_var = time()
        print(("\n\n***   Completed '{:s}' in {:.2f} seconds   ***\n".format(varname,end_var-start_var)))
        
    end = time()
    print(("\n***   Overall Timing: {:.2f} seconds   ***\n".format(end-start)))

        