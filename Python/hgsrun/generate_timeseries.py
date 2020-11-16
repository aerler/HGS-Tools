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
    from multiprocessing.pool import ThreadPool
    from time import time
    
    start = time()

    loverwrite = True
    lexec = True # actually write rasters or just include file
    lhourly = False
    interpolation = 'nearest'
    
    # project parameters
    project = 'SON'
    point_set = 'Elora'
    
    # project settings
    tgt_crs = genCRS(); lgeo = True # generic geographic lat/lon    
    # target coordinates
    if point_set == 'Elora':
        tgt_coords = (43.64,-80.4) # Elora RCS


    ## define source data
    # some defaults for most datasets
    time_chunks = 8 # SnoDAS & CaSPAr only! typically not much speed-up beyond 8
    subdataset = None; src_grid = None; dataset_kwargs = dict(grid=src_grid)
    bias_correction = None; bc_varmap = dict(); obs_name = None; bc_method = None
    target_folder_ascii = '{root:s}/{proj:s}/{grid:s}/{name:s}/{bc:s}transient_{int:s}/climate_forcing/'
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
    dataset = 'MergedForcing'; subdataset = dataset; varlist = []
#     varlist = ['liqwatflx',]
#     varlist = ['pet_pts',] # PET based on Priestley-Taylor with solar radiation only
#     varlist = ['liqwatflx','pet_har',] # assorted forcing    
    subdataset = 'NRCan'
    varlist += ['precip','Tmin','Tmax','T2',] # base variables
    varlist += ['precip_adj',] # adjusted precip data (up to 2016)
    varlist += ['pet_har',] # PET based on Hargreaves' method
    varlist += ['pet_haa',] # PET based on Hargreaves' with Allen's correction
    varlist += ['pet_th',] # PET based on Thornthwaite method
    varlist += ['pet_hog',] # PET based on simple Hogg method
#     dataset_kwargs['resolution'] = 'CA12'
    dataset_kwargs['resolution'] = 'SON60'
    
    ## Time Slice
    start_date = '1997-01-01'; end_date = '2017-12-31'
#     start_date = '2000-01-01'; end_date = '2018-01-01'
#     start_date = '2011-01-01'; end_date = '2011-02-01'    
#     start_date = '2016-01-01'; end_date = '2017-12-31'
#     start_date = '2011-01-01'; end_date = '2011-02-01'
    
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
        daily_folder,netcdf_name =ds_mod.getFolderFileName(varname='{var_str:s}', dataset=subdataset, resolution=dataset_kwargs.get('resolution',None))
        print(daily_folder,netcdf_name)
    else:
        netcdf_name = ds_mod.netcdf_filename.format('{var_str:s}', VAR='{var_str:s}')
        bc_folder = ds_mod.avgfolder
        daily_folder = ds_mod.daily_folder
        
    ## bias correction
    if bias_correction:
        assert bc_method is None, "'bc_method' can only be set manually as a varname extension,  if no explicit bias_correction is applied"
        # load pickle from file
        from processing.bc_methods import loadBCpickle
        bc_method = bias_correction
        bias_correction = loadBCpickle(method=bias_correction, obs_name=obs_name, gridstr=src_grid, domain=domain, folder=bc_folder)
        
    ## define export parameters
    driver_args = dict()
    mode = 'NetCDF'
#     mode = 'ascii'
    # modes
    if mode.lower() == 'ascii':
        # raster output using rasterio
        scalefactor = 1000. # divide by this (convert kg/m^2 to m^3/m^2, SI units to HGS internal)
        gridstr = '_'+point_set.lower()
        bc_str = bc_method+'_' if bc_method else ''
        time_interval = 'hourly' if lhourly else 'daily'
        hgs_root = os.getenv('HGS_ROOT', os.getenv('DATA_ROOT')+'HGS/')
        target_folder = target_folder_ascii.format(root=hgs_root, proj=project, grid=gridstr, name=dataset, int=time_interval, 
                                                   bc=bc_str, exp_name=exp_name, dom=domain, exp_folder=exp_folder)          
        raster_format = 'AAIGrid'        
        filename_novar = raster_name.format(dataset=dataset.lower(), variable='{var:s}',
                                            grid=point_set.lower(), date='{date:s}') # no date for now...
        driver_args = dict(significant_digits=4, fill_value=0., nodata_flag=-9999.)
        print(("\n***   Exporting point set '{}' from '{}'   ***\n".format(point_set,dataset)))
    elif mode.upper() == 'NETCDF':
        # NetCDF output using netCDF4
        gridstr = '_'+point_set.lower()
        bc_str1 = bc_method+'_' if bc_method else ''
        bc_str2 = '_'+bc_method if bc_method else ''
        target_folder = target_folder_netcdf.format(daily=daily_folder,grid=point_set,smpl=interpolation,
                                                    exp_name=exp_name, dom=domain, exp_folder=exp_folder)        
        filename = netcdf_name.format(var_str=bc_str1+gridstr, grid=bc_str2.lower()+gridstr,) # usually either var or grid
        driver_args = dict(zlib=True, complevel=1,) # compression
        #driver_args = dict(least_significant_digit=4)
        print(("\n***   Extracting point set '{}' from '{}' (NetCDF format)   ***".format(point_set,dataset)))
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
        
    # clip time axis first
    if start_date or end_date:
        time_slice = {'time':slice(start_date,end_date)}
        xds = xds.loc[time_slice]
    time_coord = xds.coords['time'].data
    # time format
    if lhourly: date_fmt = '%Y%m%d%H'; date_longfmt = '%Y-%m-%dT%H'
    else: date_fmt = '%Y%m%d'; date_longfmt = '%Y-%m-%d'
    if start_date is None: start_date = pd.to_datetime(time_coord[0]).strftime(date_longfmt)
    if end_date is None: end_date = pd.to_datetime(time_coord[-1]).strftime(date_longfmt)    
    
    ## extract point set
    # get georeference
    src_crs = getCRS(xds)
    if tgt_crs is None: tgt_crs = src_crs
    # figure outpoint coordinates in source projection
    if tgt_crs != src_crs:
        # reproject coordinates to source system
        raise NotImplementedError()
    # extract points by using slicing logic
    is_projected = int(xds.attrs.get('is_projected',0))
    ylat = xds.attrs.get('ylat','y' if is_projected else 'lat')
    xlon = xds.attrs.get('xlon','x' if is_projected else 'lon')
    point_slice = {xlon:tgt_coords[1], ylat:tgt_coords[0]}
    if interpolation == 'nearest':
        dataset = xds.sel(point_slice, method='nearest')
    else:
        raise NotImplementedError(interpolation)
    
    # print dataset
    print(dataset)
    
    # make sure target path exists
    try:
        if not osp.exists(target_folder): os.mkdir(target_folder)
    except (WindowsError,OSError):
        os.makedirs(target_folder)
        
        
    ## Export point data
   
    if mode.upper() == 'NETCDF':
      
        # xarray can just write the new dataset to disk
        nc_filepath = osp.join(target_folder,filename)
        print(("Writing to NetCDF file: '{:s}'".format(nc_filepath)))
        
        # save some meta data
        for name,coord in point_slice.items():
            dataset.setncattr(name,coord) # extracted (target) coordinate values
        dataset.setncattr('interpolation',interpolation) # netCDF4 dataset, not xarray
        for var in dataset.variables.values():
            if  xlon in var.dimensions and ylat in var.dimensions:
                var.setncattr('interpolation',interpolation)
        
        # write dataset to NetCDF file
        var_enc = {varname:driver_args for varname in varlist}
        task = xds.to_netcdf(nc_filepath, mode='w', format='NETCDF4', unlimited_dims=['time'], engine='netcdf4',
                             encoding=var_enc, compute=False)
        if lexec:
            task.compute()
        else:
            print(driver_args)
            print(task)
            task.visualize(filename=target_folder+'/netcdf.svg')  # This file is never produced

    elif mode.upper() == 'ASCII':
      
        # for ASCII, loop over variables
    
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
            
            
            # switch of overwrite/deletion if filename is not variable-specific
            if varname.lower() not in filename: 
                loverwrite = False
                assert mode.upper() == 'NETCDF' 
            # N.B.: this allows writing of all variables in the varlist to a single file
                                             
            end_var = time()
            print(("\n\n***   Completed '{:s}' in {:.2f} seconds   ***\n".format(varname,end_var-start_var)))
        
    end = time()
    print(("\n***   Overall Timing: {:.2f} seconds   ***\n".format(end-start)))

        