'''
Created on Feb. 23, 2019

Utility functions to extract data from xarray Dataset or DataArray classes.

@author: Andre R. Erler, GPL v3
'''

import numpy as np
import xarray as xr
import netCDF4 as nc


# names of valid geographic/projected coordinates
default_x_coords = dict(geo=('lon','long','longitude',), proj=('x','easting') )
default_y_coords = dict(geo=('lat','latitude',),         proj=('y','northing'))


## functions to interface with rasterio

def getGeoCoords(xvar, x_coords=None, y_coords=None, lraise=True):
    '''  helper function to extract geographic/projected coordinates from xarray'''
    
    if x_coords is None: x_coords = default_x_coords
    if y_coords is None: y_coords = default_y_coords
    
    xlon,ylat = None,None # return None, if nothing is found
    
    if isinstance(xvar,(xr.DataArray,xr.Dataset)):
        # test geographic grid and projected grids separately
        for coord_type in x_coords.keys():
            for name,coord in xvar.coords.items():
                if name.lower() in x_coords[coord_type]: 
                    xlon = coord; break
            for name,coord in xvar.coords.items():
                if name.lower() in y_coords[coord_type]: 
                    ylat = coord; break
            if xlon is not None and ylat is not None: break
            else: xlon,ylat = None,None
    elif isinstance(xvar,nc.Variable) and lraise:
        raise TypeError("Cannot infer coordinates from netCDF4 Variable - only Dataset!")
    elif isinstance(xvar,nc.Dataset):
        # test geographic grid and projected grids separately
        for coord_type in x_coords.keys():
            for name in xvar.dimensions:
                if name.lower() in x_coords[coord_type]: 
                    if name in xvar.variables:
                        xlon = xvar.variables[name]; break
            for name in xvar.dimensions:
                if name.lower() in y_coords[coord_type]: 
                    if name in xvar.variables:
                        ylat = xvar.variables[name]; break
            if xlon is not None and ylat is not None: break
            else: xlon,ylat = None,None      
    elif lraise: # optionally check input
        raise TypeError("Can only infer coordinates from xarray - not from {}".format(xvar.__class__))
    else:
        pass # return None,None
        
    # optionally raise error if no coordinates are found, otherwise just return None
    if lraise and (xlon is None or ylat is None):
        raise ValueError("No valid pair of geographic coodinates found:\n {}".format(xvar.dims))
    
    # return a valid pair of geographic or projected coordinate axis
    return xlon,ylat

  
def isGeoVar(xvar, x_coords=None, y_coords=None, lraise=True):
    ''' helper function to identify variables that have geospatial coordinates (geographic or 
        projected), based on xarray or netCDF4 dimension names '''
    
    if x_coords is None: x_coords = default_x_coords
    if y_coords is None: y_coords = default_y_coords

    if isinstance(xvar,(xr.DataArray,xr.Dataset)):
        dims = xvar.coords.keys()
    elif isinstance(xvar,(nc.Dataset,nc.Variable)):
        dims = xvar.dimensions
    elif lraise:
        raise TypeError("Can only infer coordinate system from xarray or netCDF4- not from {}".format(xvar.__class__))
    else:
        return None # evaluates as False, but allows checking
        
    # test geographic grid and projected grids separately
    for coord_type in x_coords.keys():
        xlon,ylat = False,False
        for name in dims:
            if name.lower() in x_coords[coord_type]: 
                xlon = True; break
        for name in dims:
            if name.lower() in y_coords[coord_type]: 
                ylat = True; break
        if xlon and ylat: break
    
    # if it has a valid pair of geographic or projected coordinate axis
    return ( xlon and ylat )

  
def isGeoCRS(xvar, lat_coords=None, lon_coords=None, lraise=True):
    ''' helper function to determine if we have a simple geographic lat/lon CRS (based on xarray dimension names) '''
    lat,lon = False,False
    
    if lon_coords is None: lon_coords = default_x_coords['geo']
    if lat_coords is None: lat_coords = default_y_coords['geo']
    
    if isinstance(xvar,(xr.DataArray,xr.Dataset)):
        dims = xvar.coords.keys()
    elif isinstance(xvar,(nc.Dataset,nc.Variable)):
        dims = xvar.dimensions
    elif lraise:
        raise TypeError("Can only infer coordinate system from xarray or netCDF4- not from {}".format(xvar.__class__))
    else:
        return None # evaluates as False, but allows checking
      
    # check dimension names
    for name in dims:
        if name.lower() in lon_coords: 
            lon = True; break
    for name in dims:
        if name.lower() in lat_coords: 
            lat = True; break
    
    # it is a geographic coordinate system if both, lat & lon are present
    return ( lat and lon )


def getTransform(xvar=None, x=None, y=None, lcheck=True):
    ''' generate an affine transformation from xarray coordinate axes '''
    from rasterio.transform import Affine # to generate Affine transform
    
    if isinstance(xvar,(xr.DataArray,xr.Dataset,nc.Dataset)):
        x,y = getGeoCoords(xvar, lraise=True)
    elif xvar is None and isinstance(x,(xr.DataArray,nc.Variable)) and isinstance(y,(xr.DataArray,nc.Variable)):
        pass # x and y axes are supplied directly
    elif xvar:
        raise TypeError('Can only infer GeoTransform from xarray Dataset or DataArray or netCDF4 Dataset\n - not from {}.'.format(xvar))
    
    # check X-axis
    if isinstance(x,xr.DataArray): x = x.data
    elif isinstance(x,nc.Variable): x = x[:]
    if not isinstance(x,np.ndarray): 
        raise TypeError(x)
    diff_x = np.diff(x); dx = diff_x.min()
    if lcheck and not np.isclose(dx, diff_x.max(), rtol=1.e-2): 
        raise ValueError("X-axis is not regular: {} - {}".format(dx, diff_x.max()))
    
    # check Y-axis
    if isinstance(y,xr.DataArray): y = y.data
    elif isinstance(y,nc.Variable): y = y[:]
    if not isinstance(y,np.ndarray): 
        raise TypeError(y)
    diff_y = np.diff(y); dy = diff_y.min()
    if lcheck and not np.isclose(dy, diff_y.max(), rtol=1.e-2): 
        raise ValueError("Y-axis is not regular. {} - {}".format(dy, diff_y.max()))
    
    # generate transform
    return Affine.from_gdal(x[0]-dx/2.,dx,0.,y[0]-dy/2.,0.,dy), (len(x),len(y))


def getProj(xvar, lraise=True):
    ''' infer projection from a xarray Dataset or DataArray; this function assumes that either a proj4 string or
        an EPSG designation is stored in the attributes of the dataset/variable. '''
    from geospatial.rasterio_tools import genProj # used to generate CRS object
    
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
    
    # return a GDAL/rasterio CRS instance
    return proj


def addGeoReference(xds, proj4_string=None, x_coords=None, y_coords=None):
    ''' helper function to add GDAL/rasterio-style georeferencing to an xarray dataset '''
    xds.attrs['proj4'] = proj4_string
    xlon,ylat = getGeoCoords(xds, x_coords=x_coords, y_coords=y_coords)
    xds.attrs['xlon'] = xlon.name
    xds.attrs['ylat'] = ylat.name
    for xvar in list(xds.data_vars.values()): 
        if isGeoVar(xvar):
            xvar.attrs['proj4'] = proj4_string
            xvar.attrs['xlon'] = xlon.name
            xvar.attrs['ylat'] = ylat.name
            xvar.attrs['dim_order'] = int( xvar.dims[-2:] == (ylat.name, xlon.name) )
            # N.B.: the NetCDF-4 backend does not like Python bools
    return xds


def rechunkTo2Dslices(xvar):
    ''' convenience function to rechunk an xarray so that the horizontal dimensions are contiguous (not chunked)
        N.B.: rechunking in a way that does not simply combine existing chunks seems to cause all chunks/data
              to be loaded into memory (we want to avoid that); also, chunks are defined by their size, not by 
              their number, i.e. the definition for one large 2D chunk is (len(y),len(x)) and *not* (1,1) '''
    if not isinstance(xvar,(xr.DataArray,xr.Dataset)):
        raise TypeError(xvar)
    # find horizontal/map dimentions
    xlon = xvar.coords[xvar.attrs['xlon']]; ylat = xvar.coords[xvar.attrs['ylat']]
    return xvar.chunk(chunks={xlon.name:len(xlon),ylat.name:len(ylat)}) # rechunk x/lon and y/lat
         

if __name__ == '__main__':
    pass