'''
Created on Dec 30, 2017

A utility program to average rasters over a shape and generate a time value table 
for HGS input, based on a time raster table.

@author: Andre R. Erler, GPL v3
'''

# imports
import os
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


# write time value table
def writeTimeValue(filepath=None, times=None, values=None):
    ''' write a time value file for HGS input '''
    with open(filepath, 'w') as f:
      
        # loop over times/values
        for time,value in zip(times,values):        
            # assemble line
            line = '{:>20s}   {:e}\n'.format(time, value)
            f.write(line)
            

if __name__ == '__main__':
    
    # ascii data
    folder = 'D:/Data/HGS/SNW/EnKF/TWC/forcing/'
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
        
        