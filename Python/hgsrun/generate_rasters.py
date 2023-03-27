"""
Created on Dec. 3, 2018

A module to export ASCII raster data for use with HGS from NetCDF-4 and related datasets. Also includes
functionality for clipping and reprojection/regridding.

This module is designed to make use of rasterio georeferencing/gdal functionality and xarray/dask lazy
execution and chunking.

@author: Andre R. Erler, GPL v3
"""

# external imports
import os
import os.path as osp
import numpy as np
import pandas as pd
import rasterio as rio
from rasterio.warp import calculate_default_transform  # need to import separately...
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

    class WindowsError(BaseException):
        pass


# write HGS include file
def writeIncFile(
    filepath,
    time_coord,
    filename_pattern,
    time_fmt="{:15.0f}",
    date_fmt="%Y%m%d",
    cycles=1,
):
    """write an HGS include file based on a time coordinate and a filename pattern;
    the date is used as part of the filename"""
    line_fmt = time_fmt + "     {:s}\n"
    # open file and write entries
    with open(filepath, mode="w") as incf:

        # first line
        filename = filename_pattern.format(
            pd.to_datetime(time_coord[0]).strftime(date_fmt)
        )
        line = line_fmt.format(0, filename)
        incf.write(line)

        # construct array of seconds from simulation start
        sim_time = (time_coord - time_coord[0]) / np.timedelta64(1, "s")
        sim_time[:-1] += np.diff(sim_time) / 2.0  # recenter between time-steps
        sim_time[-1] += (
            sim_time[-1] - sim_time[-2]
        )  # extrapolate offset for last element

        full_sim_time = sim_time[-1] + (sim_time[-1] - sim_time[-2]) / 2.0

        # loop over cycles
        for i in range(cycles):

            # loop over time coordinate
            for stime, date in zip(sim_time, time_coord):

                # use date to construct file name
                filename = filename_pattern.format(pd.to_datetime(date).strftime(date_fmt))
                line = line_fmt.format(stime + (i * full_sim_time), filename)
                incf.write(line)

        # last line
        last_time = cycles * full_sim_time
        # the official end of validity
        filename = filename_pattern.format(
            pd.to_datetime(time_coord[-1]).strftime(date_fmt)
        )
        line = line_fmt.format(last_time, filename)
        incf.write(line)

    ## N.B.: are the first and last lines really necessary???
    # done...
    return None


## execute raster export
if __name__ == "__main__":

    import dask
    from dask.diagnostics import ProgressBar

    # N.B.: it seems dask multi-processing is prone to memory leaks...
    # from multiprocessing.pool import ThreadPool
    from time import time

    start = time()

    loverwrite = True
    ltest = False  # prepend 'test' to output
    lwarp = True  # set False to suppress reprojection
    lclip = None  # None: clip rasters, but not NetCDF
    time_interval = "daily"
    time_chunks = 8  # typically not much speed-up beyond 8
    output_chunks = None
    resampling = "bilinear"
    lexec = True  # actually write rasters or just include file
    ## WRF grids
    # project = 'WRF'
    # grid_name  = 'wc2_d01'
    ## generate a full SnoDAS raster
    #     project = 'native'
    #     grid_name  = 'native'
    ## fast test config
    # project = 'SNW'
    # grid_name = 'snw1'
    ## explicitly defined grids
    project = 'C1W'
    grid_name = 'bfi1'
    # grid_name = 'c1w1'
    # project = 'GLB'
    # grid_name = 'dog2'
    # project = 'ASB'
    # grid_name = 'swan1'    
    # project = 'Geo'
    # grid_name = 'snodas'

    ## define target grid/projection
    # projection/UTM zone
    tgt_size = None
    tgt_geotrans = None  # valid for native grid
    if project == "WRF":
        # load pickled GridDef
        from geodata.gdal import loadPickledGridDef

        griddef = loadPickledGridDef(grid=grid_name, encoding="latin1")
        print(griddef)
        tgt_crs = genCRS(griddef.projection.ExportToProj4(), name=grid_name)
        tgt_geotrans = griddef.geotransform
        tgt_size = griddef.size
    elif project == "Geo":
        # generic geographic lat/lon
        tgt_crs = genCRS(name=grid_name)
    elif project == "C1W":
        # Projection for the Canada1Water models
        tgt_crs = genCRS(
            "+proj=aea +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
            name=grid_name,
        )
    elif project == "ARB":
        # Projection for ARB model
        tgt_crs = genCRS(
            "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +ellps=sphere +units=m +no_defs",
            name=grid_name,
        )
    elif project == "Hugo":
        # Hugo's projection for Quebec
        tgt_crs = genCRS(
            "+proj=lcc +ellps=NAD83 +datum=NAD83 +lat_0=44.0 +lat_1=46.0 +lat_2=60.0 +lon_0=-68.5  +x_0=0 +y_0=0 +units=m +no_defs",
            name=grid_name,
        )
        # N.B.: something seems to be wrong with this... rasterio can't read it...
    elif project.upper() in ("SNW"):
        # South Nation projection
        tgt_crs = genCRS(
            "+proj=utm +zone=18 +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs",
            name=grid_name,
        )
    elif project.upper() in ("SON", "GRW"):
        # southern Ontario projection
        tgt_crs = genCRS(
            "+proj=utm +zone=17 +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs",
            name=grid_name,
        )
    elif project.upper() == "GLB":
        # Assiniboin projection
        tgt_crs = genCRS(
            "+proj=utm +zone=16 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
            name=grid_name,
        )
    elif project.upper() == "ASB":
        # Assiniboin projection
        tgt_crs = genCRS(
            "+proj=utm +zone=14 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
            name=grid_name,
        )
    elif project.upper() == "CMB":
        # Columbia Mass Balance projection
        tgt_crs = genCRS(
            "+proj=utm +zone=11 +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs",
            name=grid_name,
        )
    elif project.upper() in ("QEL"):
        # Queensland (Australia) projection
        tgt_crs = genCRS(
            "+proj=utm +zone=55 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs",
            name=grid_name,
        )
    else:
        tgt_crs = None  # native grid
    # grid definition (mostly UTM grids for HGS)
    if tgt_geotrans is not None and tgt_size is not None:
        pass  # already assigned above
    elif grid_name == "ca12":  # the NRCan 1/12t deg Canada grid
        tgt_geotrans = (-141.0, 1.0 / 12.0, 0.0, 41.0, 0.0, 1.0 / 12.0)
        tgt_size = (1068, 510)  # (x,y) map size of NRCan grid
        resampling = "cubic_spline"
    elif grid_name == "na12":  # the NRCan 1/12t deg Canada grid
        tgt_geotrans = (-168.0, 1.0 / 12.0, 0.0, 25.0, 0.0, 1.0 / 12.0)
        tgt_size = (1392, 720)  # (x,y) map size of NRCan grid
        resampling = "cubic_spline"
        output_chunks = (8, 60, 58)  # time, ylat, xlon
    elif grid_name == 'snodas':
        tgt_geotrans = (-130.516666666667-0.00416666666666052, 0.00833333333333333, 0,
                        24.1000000000000-0.00416666666666052, 0, 0.00833333333333333)
        tgt_size = (8192, 4096) # (x,y) map size of SnoDAS grid
        resampling = "cubic_spline"
        output_chunks = (8, 256, 256)  # time, ylat, xlon
    elif grid_name == 'c1w1':
        tgt_size = (1205, 808)  # lower resolution 5 km grid (> 1 MB per raster)
        tgt_geotrans = (-2895.e3, 5.e3, 0., -8.e3, 0., 5.e3)  # 5 km
    elif grid_name == 'bfi1':  # Baffin Island
        tgt_size = (388, 297)
        tgt_geotrans = (-42863.98, 5.e3, 0., 2468861, 0., 5.e3)  # 5 km
    elif grid_name == "on1":
        tgt_geotrans = [-87.87916564, 0.008331298, 0.0, 41.995832443, 0.0, 0.008335113525]
        resampling = "cubic_spline"
    elif grid_name == "hd1":
        tgt_size = (70, 49)  # lower resolution 5 km grid
        tgt_geotrans = (-479184.769227, 5.0e3, 0, 68508.4877898, 0, 5.0e3)  # 5 km
        resampling = "average"  # it's a fairly coarse grid...
    elif grid_name == "son1":
        tgt_size = (118, 82)  # lower resolution 5 km grid
        tgt_geotrans = (320920.0, 5.0e3, 0, 4624073.0, 0, 5.0e3)  # 5 km
        resampling = "average"  # it's a fairly coarse grid...
    elif grid_name == "son2":
        tgt_size = (590, 410)  # higher resolution 1 km grid (~ 1 MB per day)
        tgt_geotrans = (320920.0, 1.0e3, 0, 4624073.0, 0, 1.0e3)  # 1 km
    elif grid_name == "grw1":
        tgt_size = (132, 162)  # smaller, higher resolution 1 km grid for GRW
        tgt_geotrans = (500.0e3, 1.0e3, 0, 4740.0e3, 0, 1.0e3)  # 1 km
    elif grid_name == "grw2":
        tgt_size = (27, 33)  # smaller, lower resolution 5 km grid for GRW
        tgt_geotrans = (500.0e3, 5.0e3, 0, 4740.0e3, 0, 5.0e3)  # 5 km
        resampling = "average"  # it's a fairly coarse grid...
    elif grid_name == "snw1":
        tgt_size = (18, 22)  # 5 km resolution SNW grid for WRF data
        tgt_geotrans = (438.0e3, 5.0e3, 0, 4940.0e3, 0, 5.0e3)  # 5 km
        resampling = "average"  # it's a fairly coarse grid...
    elif grid_name == "snw2":
        tgt_size = (44, 55)  # 2 km resolution SNW grid for CaLDAS/CaPA data
        tgt_geotrans = (438.0e3, 2.0e3, 0, 4940.0e3, 0, 2.0e3)  # 2 km
        resampling = "average"  # it's a fairly coarse grid...
    elif grid_name == "cmb1":
        tgt_size = (640, 826)  # higher resolution 500 m grid
        tgt_geotrans = (292557.0, 500, 0, 5872251.0, 0, -500.0)  # 500 m
    elif grid_name == "dog1":
        tgt_size = (12, 11)  # lower resolution 10 km grid
        tgt_geotrans = (247.e3, 10.e3, 0., 5375.e3, 0., 10.e3)
    elif grid_name == "dog2":
        tgt_size = (24, 22)  # medium resolution 5 km grid
        tgt_geotrans = (247.e3, 5.e3, 0., 5375.e3, 0., 5.e3)
    elif grid_name == 'pem1':
        tgt_size = (55,40)  # lower resolution 5 km grid
        tgt_geotrans = (388.e3, 5.e3, 0., 5336.e3, 0., 5.e3)  # 5 km
    elif grid_name == 'pem2':
        tgt_size = (275,200)  # higher resolution 1 km grid
        tgt_geotrans = (388.e3, 1.e3, 0., 5336.e3, 0., 1.e3)  # 1 km
    elif grid_name == "asb1":
        tgt_size = (191, 135)  # lower resolution 5 km grid
        tgt_geotrans = (-159.0e3, 5.0e3, 0.0, 5202.0e3, 0.0, 5.0e3)  # 5 km
    elif grid_name == "asb2":
        tgt_size = (955, 675)  # higher resolution 1 km grid (> 1 MB per day)
        tgt_geotrans = (-159.0e3, 1.0e3, 0.0, 5202.0e3, 0.0, 1.0e3)  # 1 km
    elif grid_name == "swan1":
        tgt_size = (10, 7)
        tgt_geotrans = (501.e3, 1.0e3, 0.0, 5464.0e3, 0.0, 1.0e3)        
    elif grid_name == "arb2":
        tgt_geotrans = [-1460500, 5e3, 0, 810500, 0, 5e3]
        tgt_size = (284, 258)
        resampling = "average"  # it's a fairly coarse grid...
    elif grid_name == "arb3":
        tgt_geotrans = [-1280e3, 5e3, 0, 900e3, 0, 5e3]
        tgt_size = (172, 144)
        resampling = "average"  # it's a fairly coarse grid...
    elif grid_name == "ccj1":
        tgt_geotrans = (-1618000, 5000, 0, 3058000, 0, 5000)
        tgt_size = (22, 33)  # 5 km
    elif grid_name == "qel1":
        tgt_size = (121, 191)  # approx. 10 km grid, similar to ERA5
        # tgt_geotrans = (-26770., 10400, 0., 6902510, 0., 10970)  # ~10 km
        tgt_geotrans = (-26772.0, 10398.9, 0.0, 6902511, 0.0, 10968.7)  # ~10 km
        # N.B.: this grid attempts to achieve direct grid point correspondence to the ERA5-Land lat/lon grid
        #       at 148E/23S
    elif grid_name == "qel2":
        tgt_size = (126, 211)  # exactly 10 km grid, similar to ERA5
        tgt_geotrans = (-26772.0, 10e3, 0.0, 6902511, 0.0, 10e3)  # 10 km
    elif grid_name == "native":  # original grid
        time_chunks = 1  # this can be pretty big!
        tgt_size = None
        tgt_geotrans = None  # native grid
    else:
        raise NotImplementedError(grid_name)
    # cast geotransform into an Affine object
    if tgt_geotrans is not None:
        tgt_geotrans = rio.transform.Affine.from_gdal(*tgt_geotrans)

    ## define source data

    lexec = True
    ltest = False  # prefix with 'test' - don't overwrite exiting data

    # some defaults for most datasets
    dataset_kwargs = dict()
    subdataset = None
    dataset_name = None  # the dataset name in the target folder; defaults to dataset
    bias_correction = None
    bc_varmap = dict()
    obs_name = None
    bc_method = None
    period = None
    sim_cycles = 1
    fill_masked = True
    fill_max_search = 5
    raster_name = "{dataset:s}_{variable:s}_{grid:s}_{date:s}.asc"
    target_folder_ascii = "{root:s}/{proj:s}/{grid:s}/{name:s}/{bc:s}_{int:s}/"
    # target_folder_ascii = '//aquanty-nas/share/temp_data_exchange/Erler/{proj:s}/{grid:s}/{name:s}/{bc:s}_{int:s}/'
    # target_folder_ascii = '//aquanty-nas/share/temp_data_exchange/Erler/{proj:s}/{grid:s}/{name:s}/{int:s}/'
    if ltest:
        target_folder_ascii += "test/"  # store in subfolder
    data_mode = "daily"

    # dummy variables for WRF
    exp_name = None
    domain = None
    WRF_exps = None
    filetype = None

    ## SnoDAS & NRCan
    #     dataset = 'SnoDAS' # default...
    # #     bias_correction = 'SMBC'; obs_name = 'NRCan'
    #     bc_varmap = dict(liqprec='liqwatflx') # just for testing...
    #     varlist = ['snow']; bc_method = 'rfbc'
    #     dataset_kwargs = dict(grid='on1', bias_correction=bc_method)
    #     resampling = 'bilinear'
    #     end_date = '2011-02-01'

    ## CaSPAr
    # dataset = 'CaSPAr'; dataset_kwargs = dict(grid='lcc_snw'); time_interval = 'hourly'

    ## MergedForcing Daily
    # varlist = []
    # dataset = 'MergedForcing'
    # # subdataset = 'NRCan'; varlist = ['pet_hog']
    # # subdataset = dataset; varlist = ['liqwatflx_ne5']
    # subdataset = None; varlist = ['liqwatflx_ne5', 'pet_hog']
    # # subdataset = dataset; varlist = ['pet_pts',] # PET based on Priestley-Taylor with solar radiation only
    # dataset_kwargs = dict(dataset=subdataset)
    # dataset_kwargs['resolution'] = 'NA12'  # for NRCan: SON60, NA12
    # resampling = 'bilinear'
    # dataset_kwargs['dataset_index'] = dict(liqwatflx='MergedForcing',
    #                                        liqwatflx_ne5='MergedForcing',
    #                                        pet_har='NRCan', pet_hog='NRCan')
    # dataset_kwargs['grid'] = 'son2'
    # dataset_kwargs['grid'] = 'snw2'; resampling = None; lwarp = False

    ## MergedForcing Monthly (incl. ERA5)
    dataset = 'MergedForcing'  # to load module and use in this script
    dataset_args = dict(ERA5=dict(grid='NA10'.lower(), subset='ERA5L'),
                        NRCan=dict(resolution='NA12'))
    # subdataset = 'MergedForcing';  varlist = ['liqwatflx_ne5']
    # subdataset = 'MergedForcing';  varlist = ['liqwatflx_sno',]
    # subdataset = 'NRCan';  varlist = ['pet_hog']
    subdataset = 'NRCan';  varlist = ['precip']
    # subdataset = 'SnoDAS';  varlist = ['dswe',]
    # subdataset = 'ERA5';  dataset_name = 'ERA5';  varlist = ['liqwatflx', 'pet_era5']
    src_grid = 'na12'  # can be either resolution or grid, depending on source dataset
    # src_grid = None
    # src_grid = 'snodas'
    # time_interval = "clim"; period = (1981, 2011)
    # time_interval = "clim"; period = (2000, 2020)
    time_interval = "monthly"
    # start_date = end_date = None; sim_cycles = 20  # cycles/repetitions in include file for periodic forcing
    # data_mode = "daily"  # averages computed from daily data
    # time_interval = "daily"
    dataset_kwargs = dict(period=period, grid=src_grid, ldt64=True, chunks=True,
                          mode=data_mode, aggregation=time_interval, multi_chunks=None,
                          dataset=subdataset, dataset_args=dataset_args)

    # ## ERA5 Daily
    # dataset = 'ERA5'; subdataset = 'ERA5L'
    # #time_chunks = 92 # for small grids only!
    # varlist = ['snow','dswe',]
    # # varlist = ['precip','pet_era5','liqwatflx','snow','dswe',]
    # # varlist = ['pet_era5','liqwatflx',]
    # data_mode = "daily"
    # time_interval = "daily"
    # sim_cycles = 10  # cycles/repetitions in include file for periodic forcing
    # start_date = end_date = None
    # dataset_kwargs = dict(subset=subdataset, combine_attrs='override')
    # dataset_kwargs['resolution'] = 'NA10'  # NA10, AU10, SON10
    # resampling = 'bilinear'  # apparently we need to pre-chunk or there is a memory leak..
    # fill_masked = True

    ## C1W data
    # dataset = 'C1W'; subdataset = 'C1W_Soil'
    # time_interval = 'monthly'; data_mode = 'avg'
    # varlist = ['Tsl1','Tsl2','Tsl3']
    # dataset_kwargs = dict(subset=subdataset)
    # dataset_kwargs['resolution'] = 'NA005'
    # resampling = 'cubic_spline' # apparently we need to pre-chunk or there is a memory leak..
    # fill_masked = True

    # ## WRF requires special treatment
    # dataset = "WRF"
    # lhourly = False
    # bias_correction = None
    # resampling = "bilinear"
    # if project in ('ARB', 'CMB', 'ASB'):
    #     from projects.WesternCanada import WRF_exps
    # else:
    #     from projects.GreatLakes import WRF_exps
    # exp_name = os.getenv('WRFEXP', "g-ens")
    # domain = 1
    # filetype = "aux"
    # data_mode = "avg"
    # time_interval = "clim"
    # sim_cycles = 10  # cycles/repetitions in include file for periodic forcing
    # start_date = end_date = None
    # # start_date = '1979-01-01'; end_date = '1979-12-31'
    # dataset_kwargs = dict(experiment=exp_name, domains=domain, filetypes=filetype,
    #                       exps=WRF_exps, lconst=False,)
    # # target_folder_ascii = "//aquanty-nas/share/temp_data_exchange/Erler/{proj:s}/{grid:s}/{exp_name:s}_d{dom:0=2d}/{bc:s}_{int:s}/climate_forcing/"
    # target_folder_ascii = "{root:s}/{proj:s}/{grid:s}/{exp_name:s}_d{dom:0=2d}/{bc:s}_{int:s}/climate_forcing/"
    # target_folder_netcdf = "{exp_folder:s}/{grid:s}/{smpl:s}/"
    # # bias_correction = 'MyBC'; bc_varmap = dict(liqwatflx=None); obs_name = 'CRU'
    # # bias_correction = 'AABC'; bc_varmap = dict(liqwatflx='precip'); obs_name = 'CRU'
    # varlist = ["liqwatflx", "pet"]

    # start_date = '1997-01-01'; end_date = '2017-12-31' # SON/SNW full period
    # start_date = '1981-01-01'; end_date = '2017-12-31' # SON/SNW full period
    # start_date = '1981-01-01'; end_date = '2020-09-01' # MergedForcing period
    # start_date = '1981-01-01'; end_date = '2020-12-31' # full ERA5-Land period
    # start_date = '2000-01-01'; end_date = '2018-01-01'
    start_date = '2009-12-14'; end_date = '2020-12-31' # combined NRCan-SnoDAS period
    # start_date = '2010-10-01'; end_date = '2020-12-31' # Prairie NRCan-SnoDAS period    
    # start_date = '2010-01-01'; end_date = '2020-12-31' # conservative NRCan-SnoDAS period
    # start_date = '2016-01-01'; end_date = '2017-12-31'
    # start_date = '2016-01-01'; end_date = '2016-01-31' # for testing NRCan

    if ltest:
        varlist = varlist[:1]
        resampling = "nearest"  # faster, for testing...
        if start_date and end_date:
            start_date = "2011-01-01"
            end_date = "2012-01-01" if time_interval.lower() == "monthly" else "2011-01-15"

    ## output type: ASCII raster or NetCDF-4
    # mode = 'NetCDF'
    mode = "raster2d"

    # import dataset module
    ds_mod = import_module("datasets.{0:s}".format(dataset))

    # set target folder/file names based on experiment/dataset
    exp_folder = None
    exp = None
    bc_folder = None
    target_kwargs = dataset_kwargs.copy()
    for key in ('grid', 'resampling', 'mode', 'aggregation', 'domain', 'filetypes', 'lconst'):
        target_kwargs.pop(key, None)
    if dataset == "WRF":
        exp_folder, exp, exp_name, domain = ds_mod.getFolderNameDomain(grid=grid_name,
                                                                       resampling=resampling,
                                                                       lreduce=True,
                                                                       mode=data_mode,
                                                                       **target_kwargs)
        print("{exp_name:s}_d{dom:0=2d}".format(exp_name=exp_name, dom=domain))
        netcdf_folder = exp_folder
        bc_folder = exp_folder
        if data_mode.lower() == "avg":
            netcdf_name = ds_mod.fileclasses[filetype].tsfile.format(domain, "{grid:s}")
        else:
            netcdf_name = ds_mod.fileclasses[filetype].dailyfile.format(domain, "{grid:s}")
    elif hasattr(ds_mod, "getFolderFileName"):
        netcdf_folder, netcdf_name = ds_mod.getFolderFileName(grid=grid_name,
                                                              resampling=resampling,
                                                              mode=data_mode,
                                                              aggregation=time_interval,
                                                              varname="{var_str:s}",
                                                              lcreateFolder=False,
                                                              **target_kwargs)
        bc_folder = ds_mod.avgfolder
    else:
        raise NotImplementedError('''Without a 'getFolderFileName' function, handling
 of grids and resampling in filenames has to be implemented explicitly''')

    if ltest:
        netcdf_name = "test_" + netcdf_name
        print(netcdf_folder, netcdf_name)

    ## bias correction
    if bias_correction:
        assert (bc_method is None), "'bc_method' can only be set manually as a varname extension,  if no explicit bias_correction is applied"
        # load pickle from file
        from processing.bc_methods import loadBCpickle

        bc_method = bias_correction
        bias_correction = loadBCpickle(
            method=bias_correction,
            obs_name=obs_name,
            gridstr=grid_name,
            domain=domain,
            folder=bc_folder,
        )

    ## define export parameters
    driver_args = dict()
    m3factor = None
    raster_format = None
    # modes
    if mode.lower() == "raster2d":
        # raster output using rasterio
        if lclip is None:
            lclip = True  # remove negative values
        m3factor = 1000.0  # divide by this (convert kg/m^2 to m^3/m^2, SI units to HGS internal)
        gridstr = (
            dataset.lower() if grid_name.lower() == "native" else grid_name.lower())
        bc_str = bc_method + "_" if bc_method else "noBC"
        hgs_root = os.getenv("HGS_ROOT", os.getenv("DATA_ROOT") + "HGS/")
        target_folder = target_folder_ascii.format(root=hgs_root,
                                                   proj=project,
                                                   grid=gridstr,
                                                   name=dataset_name or dataset,
                                                   int=time_interval,
                                                   bc=bc_str,
                                                   exp_name=exp_name,
                                                   dom=domain,
                                                   exp_folder=exp_folder,)
        raster_format = "AAIGrid"
        filename_novar = raster_name.format(dataset=(subdataset or dataset).lower(),
                                            variable="{var:s}",
                                            grid=grid_name.lower(),
                                            date="{date:s}",)  # no date for now...
        driver_args = dict(significant_digits=4, fill_value=0.0, nodata_flag=np.NaN)
        print("\n***   Exporting '{}' to raster format {}   ***\n".format(dataset, raster_format))

    elif mode.upper() == "NETCDF":
        # NetCDF output using netCDF4
        if lclip is None:
            lclip = False
        gridstr = "" if grid_name.lower() == "native" else "_" + grid_name.lower()
        bc_str1 = bc_method + "_" if bc_method else ""
        bc_str2 = "_" + bc_method if bc_method else ""
        target_folder = netcdf_folder
        filename_novar = netcdf_name.format(var_str=bc_str1 + "{var:s}",  # usually either var or grid
                                            grid=bc_str2.lower() + gridstr)
        # driver_args = dict(least_significant_digit=4)
        print("\n***   Regridding '{}' to '{}' (NetCDF format)   ***".format(dataset, grid_name))
    else:
        raise NotImplementedError
    print(("   Variable list: {}\n".format(str(varlist))))

    # lazily load dataset (assuming xarray)
    if time_interval.lower() == "hourly":
        xds = ds_mod.loadHourlyTimeSeries(
            varlist=varlist, time_chunks=time_chunks, **dataset_kwargs)
    elif time_interval.lower() == "daily":
        xds = ds_mod.loadDailyTimeSeries(
            varlist=varlist, **dataset_kwargs)
    elif time_interval.lower() == "monthly":
        xds = ds_mod.loadTimeSeries(varlist=varlist, lxarray=True, **dataset_kwargs)
    elif time_interval.lower().startswith("clim"):
        xds = ds_mod.loadClimatology(varlist=varlist, lxarray=True, **dataset_kwargs)
    else:
        raise ValueError(time_interval)

    # get georeference
    src_crs = getCRS(xds)
    src_geotrans, src_size = getTransform(xds)
    if tgt_crs is None:
        tgt_crs = src_crs
    # figure out bounds for clipping
    if tgt_geotrans is not None and tgt_size is not None:
        # get outline
        left, bottom = tgt_geotrans * (0, 0)
        right, top = tgt_geotrans * tgt_size
        if tgt_crs != src_crs:
            # reproject outline, if necessary
            trans, w, h = calculate_default_transform(src_crs=tgt_crs,
                                                      dst_crs=src_crs,
                                                      width=tgt_size[0],
                                                      height=tgt_size[1],
                                                      left=left,
                                                      bottom=bottom,
                                                      right=right,
                                                      top=top,)
            left, top = trans * (-1, -1)
            right, bottom = trans * (w + 1, h + 1)  # need some padding to avoid margins
        # clip source data
        if src_geotrans.e < 0:
            bottom, top = top, bottom
        space_slice = {xds.attrs["xlon"]: slice(left, right),
                       xds.attrs["ylat"]: slice(bottom, top)}
        print(space_slice)
        xds = xds.loc[space_slice]
    # clip time axis as well
    if start_date or end_date:
        time_slice = {"time": slice(start_date, end_date)}
        xds = xds.loc[time_slice]
    print(xds)
    src_geotrans, src_size = getTransform(xds)  # recalculate after clipping
    if tgt_geotrans is None:
        tgt_geotrans = src_geotrans
    if tgt_size is None:
        tgt_size = src_size
    # time format
    if time_interval.lower() == "hourly":
        date_fmt = "%Y%m%d%H"
        date_longfmt = "%Y-%m-%dT%H"
    elif time_interval.lower() == "daily":
        date_fmt = "%Y%m%d"
        date_longfmt = "%Y-%m-%d"
    elif time_interval.lower() == "monthly":
        date_fmt = "%Y%m"
        date_longfmt = "%Y-%m"
    elif time_interval.lower().startswith("clim"):
        # fix timedelta for climatologies: this will center each month in the year 2000
        xds = xds.assign_coords({'time': xds.coords["time"].data + np.datetime64('1999-12-16', 'D')})
        date_fmt = "%m"
        date_longfmt = "%m"

    time_coord = xds.coords["time"].data
    if start_date is None:
        start_date = pd.to_datetime(time_coord[0]).strftime(date_longfmt)
    if end_date is None:
        end_date = pd.to_datetime(time_coord[-1]).strftime(date_longfmt)

    # make sure target path exists
    try:
        if not osp.exists(target_folder):
            os.mkdir(target_folder)
    except (WindowsError, OSError):
        os.makedirs(target_folder)

    ## loop over variables (for rasters, that will typically not be necessary)

    for varname in varlist:

        start_var = time()

        # select variable
        xvar = xds[varname]
        filename = filename_novar.format(var=varname.lower(), date="{:s}")

        print("\n\n###   Processing Variable '{:s}'   ###".format(varname, start_date, end_date))

        ## generate inc file
        if mode.lower() == "raster2d":
            start_inc = time()
            inc_filepath = target_folder + varname.lower() + ".inc"
            print(("\nWriting HGS include file:\n '{:s}'".format(inc_filepath)))
            writeIncFile(filepath=inc_filepath,
                         time_coord=time_coord,
                         filename_pattern=filename,
                         date_fmt=date_fmt,
                         cycles=sim_cycles,
                         )
            end_inc = time()
            # print("\nTiming to write include file: {} seconds".format(end_inc-start_inc))

        if not lexec:
            print("\nNot executing workload --- set 'lexec' to 'True' to execute workload")
            exit()

        ## generate workload for lazy execution
        start_load = time()
        print("\nConstructing Workload for '{:s}' from {:s} to {:s}.   ***".format(varname, start_date, end_date))
        filepath = target_folder + filename
        if mode.lower() == "raster2d":
            print("Output folder: '{:s}'\nRaster pattern: '{:s}'".format(target_folder, filename))
        elif mode.upper() == "NETCDF":
            print(("NetCDF file: '{:s}'".format(filepath)))
            original_filepath = filepath
            filepath += ".tmp"
            filename += ".tmp"
            # determine chunks
            driver_args["smart_chunks"] = False if output_chunks else True
            driver_args["chunksizes"] = output_chunks if output_chunks else None

        # explicitly determine chunking to get complete 2D lat/lon slices
        xvar = rechunkTo2Dslices(xvar, time=time_chunks)

        # apply a scaling factor
        if (m3factor is not None) and ("kg/m^2" in xvar.units or "mm" in xvar.units):
            xvar /= m3factor
        # N.B.: apply scalefactor 'in-place' so that xarray variable attributes
        #       are preserved (it will still execute delayed); applying the scale-
        #       factor after regridding is slightly faster, but this is cleaner

        # Bias-correction parameters
        if bias_correction:
            bc_varname = bc_varmap.get(varname, varname)
            bc_object = None if bc_varname is None else bias_correction
        else:
            bc_object = None
            bc_varname = None

        # switch of overwrite/deletion if filename is not variable-specific
        if varname.lower() not in filename:
            loverwrite = False
            assert mode.upper() == "NETCDF"
        # N.B.: this allows writing of all variables in the varlist to a single file,
        #       note, however, that it does not allow appending to existing files,
        #       since the files are created with a .tmp extension and later renamed.

        # generate dask execution function
        dask_fct, dummy, dataset = generate_regrid_and_export(xvar,
                                                              time_coord=time_coord,
                                                              lwarp=lwarp,
                                                              tgt_crs=tgt_crs,
                                                              tgt_geotrans=tgt_geotrans,
                                                              tgt_size=tgt_size,
                                                              mode=mode,
                                                              resampling=resampling,
                                                              time_fmt=date_fmt,
                                                              folder=target_folder,
                                                              filename=filename,
                                                              driver=raster_format,
                                                              bias_correction=bc_object,
                                                              bc_varname=bc_varname,
                                                              lclip=lclip,
                                                              fill_masked=fill_masked,
                                                              fill_max_search=fill_max_search,
                                                              lecho=True,
                                                              loverwrite=loverwrite,
                                                              **driver_args)
        # N.B.: the dataset returned here is a NetCDF dataset, not a xarray dataset!

        # now map regridding operation to blocks
        n_loads = len(xvar.chunks[0])
        dummy_output = xvar.data.map_blocks(
            dask_fct, chunks=dummy.shape, dtype=dummy.dtype
        )
        work_load = [dummy_output]

        end_load = time()
        # print("\nTiming to construct workload: {:.2f} seconds".format(end_load-start_load))

        # execute delayed computation
        print("\n***   Executing {:d} Workloads for '{:s}' using Dask   ***".format(n_loads, varname))
        # print(("Chunks (time only): {}".format(xvar.chunks[0])))

        # Dask scheduler settings - threading can make debugging very difficult
        # if ltest:
        #     dask.config.set(scheduler='synchronous')  # single-threaded for small workload and debugging
        # else:
        #     # dask.config.set(scheduler='threading')  # default scheduler - some parallelization
        #     dask.config.set(scheduler='synchronous')  # for very large fields

        with ProgressBar():
            dask.compute(*work_load, scheduler='threads', num_workers=4)

        # print("\nDummy output:")
        # print(dummy_output)
        print("\nDummy Size in memory: {:f} MB".format(dummy_output.nbytes / 1024.0 / 1024.0))

        # finalize
        if mode.upper() == "NETCDF":
            dataset.setncattr("resampling", resampling)  # netCDF4 dataset, not xarray
            dataset.variables[varname].setncattr("resampling", resampling)
            dataset.close()

        end_var = time()
        print("\n\n***   Completed '{:s}' in {:.2f} seconds   ***\n".format(varname, end_var - start_var))

        # replace original file
        if mode.upper() == "NETCDF":
            if os.path.exists(original_filepath):
                os.remove(original_filepath)
            os.rename(filepath, original_filepath)

    end = time()
    print(("\n***   Overall Timing: {:.2f} seconds   ***\n".format(end - start)))
