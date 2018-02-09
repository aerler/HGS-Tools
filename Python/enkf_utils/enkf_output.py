'''
Created on Jan 4, 2018

A collection of functions to read EnKF output files.

@author: Andre R. Erler, GPL v3
'''

# imports
import os
import numpy as np
import pandas as pd
from geodata.misc import ArgumentError
# internal imports


## load functions

def loadHydro(filename='discharge_out.mpiio', folder=None, nreal=None, ntime=None, dtype=np.float64):
    ''' function to load hydroggraphs/discharge from EnKF output '''
    if not nreal and not ntime: 
      raise ArgumentError("Please specify number of realizations 'nreal' or number of time steps 'ntime'.")
    filepath = os.path.join(folder,filename)
    if isinstance(dtype, basestring): dtype = getattr(np,dtype)
    # load data
    data = np.fromfile(filepath, dtype=dtype)    
    # reshape (need to know number or realizations or time steps)
    n = data.size
    if nreal:
        nt = int(n/nreal)
        nr = nreal
        if ntime  and nt != ntime:
            raise ValueError("Given number of time steps is not consistent with file size or data type ({} != {}).".format(ntime,nt))
    elif ntime:
        nt = ntime
        nr = int(n/ntime)
    if nt*nr != n: 
        raise ValueError("Given number of realizations ({}) and time steps do not divide number of data points ({}).".format(nr,nt,n))
    data = data.reshape((nt,nr))
    # return timeseries
    return data


def loadObs(filename='obs_h.dat', folder=None, varlist=None, lpandas=False):
    ''' function to load observation/innovation data from EnKF output '''
    filepath = os.path.join(folder,filename)
    if lpandas:
        vardata = pd.read_csv(filepath)
        # mark missing values
        vardata[vardata > 10000.] = np.NaN
        # N.B.: because each time step has information for all observations and 
        #       all realizations, each variable actually has three dimensions
    else:
        # read headers
        with open(filepath, 'r') as f:
            varnames = f.readline().split()
            data = np.loadtxt(f,)
        if varlist is None:
            varlist = [varname for varname in varnames if varname not in ('Timestep', 'NoObs', 'NoReal')]
        else:
            for varname in varlist:
                if varname not in varnames: 
                    raise ValueError(varname)
        # infer dimensions
        ntime = int(data[:,varnames.index('Timestep')].max())
        nobs  = int(data[:,varnames.index('NoObs')].max())
        nreal = int(data[:,varnames.index('NoReal')].max())
        nvars = len(varnames)
        # remove incomplete timestep
        if (ntime*nobs*nreal,nvars) != data.shape:
            ntime -= 1
            data = data[:(ntime*nobs*nreal),:]
        assert (ntime*nobs*nreal,nvars) == data.shape, (ntime,nobs,nreal,nvars)
        # reshape array
        data = data.reshape((ntime,nobs,nreal,nvars))
        # detect missing values
        data[data > 10000.] = np.NaN
        # pick variables
        vardata = dict()
        for varname in varlist:
            vardata[varname] = data[:,:,:,varnames.index(varname)].squeeze()
    # return dataframe or dictionary of variable arrays
    # dimentsions are: time, observation, realization
    return vardata


if __name__ == '__main__':
  
    ## settings
#     folder = 'D:/Data/HGS/SNW/EnKF/TWC/enkf_test_open_dec/' # experiment folder
    folder = 'D:/Data/HGS/SNW/EnKF/TWC/enkf_test_closed_dec/' # experiment folder
    out_folder = os.path.join(folder,'out/') # default output folder
    if not os.path.exists(out_folder): raise IOError(out_folder)
    
#     mode = 'test_load_hydro'
    mode = 'test_load_obs'
    
    # do some tests
    if mode == 'test_load_hydro':
        
        data = loadHydro(folder=out_folder, nreal=90, ntime=31, dtype='float64')
        
#         data = data[0,:].reshape((1,20))
        print(data.mean(),data.std(axis=0).mean()) # mean and ensemble spread 
        print(data.shape)
  
    elif mode == 'test_load_obs':
      
        lpandas = False
        vardata = loadObs(folder=out_folder, lpandas=lpandas, )
        
        if lpandas: 
            print(vardata)
        else:
            for varname,data in vardata.items():
                print(varname, data.shape) 

    elif mode == 'test_load_binary':
      
        raise NotImplementedError
  
    elif mode == 'test_load_stats':
      
        raise NotImplementedError

