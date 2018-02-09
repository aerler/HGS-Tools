'''
Created on Dec 30, 2017

A utility program to average rasters over a shape and generate a time value table 
for HGS input, based on a time raster table.

@author: Andre R. Erler, GPL v3
'''

# imports
import os
import pandas as pd
# GeoPy imports
from geodata.gdal import Shape, GridDefinition
from utils.ascii import readRasterArray


# load and average raster
def averageRasters(filelist=None, projection=None, shape=None, lfeedback=True):
    ''' load a list of rasters into an array and average over a shape; return timeseries '''
    # read rasters
    data, geotransform = readRasterArray(file_pattern='{FILENAME}', lfeedback=lfeedback, lgeotransform=True, 
                                             lskipMissing=False, lna=False, 
                                             axes=('FILENAME',), FILENAME=filelist)
    if lfeedback: print(data.shape,geotransform)
    # create GridDef to rasterize shape
    griddef = GridDefinition(name='raster', projection=projection, geotransform=geotransform, 
                             size=(data.shape[-1],data.shape[-2]), convention='proj4')
    if lfeedback: print(griddef)
    # create mask
    mask = shape.rasterize(griddef=griddef, invert=False)
    assert mask.shape == data.shape[1:], griddef
    mask = mask.reshape((1,)+mask.shape).repeat(data.shape[0],axis=0)
    assert mask.shape == data.shape, mask
    #print(mask)  
    # extract and average data
    #print(mask[1,:].sum(),)
    data.mask = mask
    #print(data[1,:])
    #print(data[1,:].sum(),data[1,:].mean())
    ts_data = data.mean(axis=-1).mean(axis=-1)
    assert len(filelist) == data.shape[0], data.shape
    # return timeseries of shape averages
    return ts_data

# read time raster table
def readTimeRaster(filepath, folder=None, lvalidate=True):
    ''' parse a HGS include file containing a time raster table and return a list of times and rasters '''
    # read file
    with open(filepath) as f:
        lines = f.readlines()
    # parse lines
    times = []; filelist = []
    for line in lines:
        line = line.strip().split()
        # time is the first element
        times.append(line.pop(0))
        # the file is the second
        if len(line) == 1: filename = line[0]
        else: filename = ' '.join(line)
        if lvalidate:
            filepath = os.path.join(folder,filename) if folder else filename
            if not os.path.exists(filepath):
                raise IOError(filepath)
        filelist.append(filename)
    # return
    return times, filelist
  
# write time raster table based on datetime and file pattern
def writeTimeRaster(inc_file, date_range=None, filepattern=None, date_format='%Y%m%d', 
                    lvalidate=True, lfeedback=False):
    ''' take a datetime and filepattern and write an include file '''
    # date settings
    begin = pd.to_datetime(date_range[0])
    end   = pd.to_datetime(date_range[1])
    freq  = date_range[2]
    datelist = pd.date_range(begin, end, freq=freq)
    model_time = ( datelist - begin ).astype('timedelta64[s]')
    # open file and write
    if lfeedback: 
        print("Writing include file '{}':".format(inc_file))
    with open(inc_file, 'w') as f:
        # loop over dates
        for delta,date in zip(model_time,datelist):
            # assemble raster filename
            date_str = date.strftime(date_format)
            raster = filepattern.format(date_str)
            if lvalidate and not os.path.exists(raster):
                raise IOError(raster)
            # assemble line
            line = '{:>20s}   {:s}\n'.format(str(delta), raster)
            f.write(line)
            if lfeedback: print(line[:-2]) # omit linebreak
    # done...


# write simple time table based on datetime
def writeTimeTable(inc_file, date_range=None, lfeedback=False):
    ''' take a datetime write an include file with time values (for output times) '''
    # date settings
    begin = pd.to_datetime(date_range[0])
    end   = pd.to_datetime(date_range[1])
    freq  = date_range[2]
    datelist = pd.date_range(begin, end, freq=freq)
    model_time = ( datelist - begin ).astype('timedelta64[s]')
    # open file and write
    if lfeedback: 
        print("Writing include file '{}':".format(inc_file))
    with open(inc_file, 'w') as f:
        # loop over dates
        for time in model_time:
            # assemble line
            line = '{:>20s}\n'.format(str(time))
            f.write(line)
            if lfeedback: print(line[:-2]) # omit linebreak
    # done...


# write time value table
def writeTimeValue(filepath=None, times=None, values=None):
    ''' write a time value file for HGS input using a list of time-stamps and values'''
    with open(filepath, 'w') as f:     
        # loop over times/values
        for time,value in zip(times,values):        
            # assemble line
            line = '{:>20s}   {:e}\n'.format(time, value)
            f.write(line)
            

if __name__ == '__main__':
  
    # patch symlink on Windows
    from hgsrun.misc import symlink_ms
    if os.name == "nt":
      os.symlink = symlink_ms # replace os symlink with this function

    ## settings
    # available range
    folder = 'D:/Data/HGS/SNW/EnKF/TWC/enkf_may/' # folder where files a written
    date_range = ('2017-05-01', '2018-01-31', 'D') # date range for files
    # just december
#     folder = 'D:/Data/HGS/SNW/EnKF/TWC/enkf_december/' # folder where files a written
#     date_range = ('2017-12-01', '2017-12-31', 'D') # date range for files
    # just november and december
#     folder = 'D:/Data/HGS/SNW/EnKF/TWC/enkf_november/' # folder where files a written
#     date_range = ('2017-11-01', '2018-12-31', 'D') # date range for files
    
    # work folder setup
    if not os.path.exists(folder): os.mkdir(folder)
    forcing_folder = 'D:/Data/HGS/SNW/EnKF/TWC/forcing/HistoricDailyTransientRaster_TWC/'
    local_forcing_folder = 'climate_forcing'
    if not os.path.exists(os.path.join(folder,local_forcing_folder)):
        os.symlink(forcing_folder,os.path.join(folder,local_forcing_folder))
    
    # task execution
    tasks = []
    tasks += ['write_time_inc'  ]
    tasks += ['write_raster_inc']
    tasks += ['raster_average'  ]

    if 'write_time_inc' in tasks:
      
        # definitions
        #folder = 'D:/Data/HGS/SNW/EnKF/TWC/forcing/'
        inc_file = 'output.inc'
        #date_range = ('2017-05-01', '2017-12-31', '1D')
        
        os.chdir(folder)
        # write file
        print('')
        writeTimeTable(inc_file, date_range=date_range, lfeedback=True)
        if not os.path.exists(inc_file): 
            raise IOError(inc_file)
        
        print('\n===\n')

    if 'write_raster_inc' in tasks:
      
        # definitions
        #folder = 'D:/Data/HGS/SNW/EnKF/TWC/forcing/'
        
        inc_files = {'precip.inc':local_forcing_folder+'/pcp_{:s}.asc', 
                     'pet.inc':local_forcing_folder+'/pet_{:s}.asc'}
        #date_range = ('2017-05-01', '2017-12-31', '1D')
        
        os.chdir(folder)
        # loop over file types
        for inc_file,filepattern in inc_files.items():
            # write files
            print('')
            writeTimeRaster(inc_file, date_range=date_range, filepattern=filepattern, lfeedback=True)
            if not os.path.exists(inc_file): 
                raise IOError(inc_file)
            
        print('\n===\n')

    if 'raster_average' in tasks:
      
        # ascii data
        #folder = 'D:/Data/HGS/SNW/EnKF/TWC/forcing/'
        inc_files = {'precip.inc':'precip_values.inc', 'pet.inc':'pet_values.inc'}
        # shape data
        shape_name = 'WholePRW' # Payne River Watershed
        shape_folder = 'C:/Users/aerler/Data/shapes/Basins/Payne River Watershed/'
        
        
        # define projection of raster
        ## parameters for South Nation grids
        projection = "+proj=utm +zone=18 +north +ellps=NAD83 +datum=NAD83 +units=m +no_defs"
        
        # move into work dir
        os.chdir(folder)
    
        # load Shape file
        shape = Shape(name=shape_name, shapefile=shape_name, folder=shape_folder, load=True, ldebug=True)
        
        # loop over include files
        for inc_file,new_inc_file in inc_files.items():
            
            print("\nProcessing '{}':\n".format(inc_file))
            # read inc file
            times, filelist = readTimeRaster(filepath=inc_file, folder=folder, lvalidate=True)
            print(times)
            print(filelist)
            
            # load rasters and average over shape
            print('')
            timeseries = averageRasters(filelist=filelist, projection=projection, shape=shape)
            assert len(timeseries) == len(filelist)
            
            # write new inc file
            writeTimeValue(filepath=new_inc_file, times=times, values=timeseries)
            print(new_inc_file); print('')
            if not os.path.exists(new_inc_file):
                raise IOError(new_inc_file)        
        
        print('\n===\n')
