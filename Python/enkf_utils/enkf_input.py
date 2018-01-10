'''
Created on Jan 1, 2018

A collection of functions to generate EnKF input files.

@author: Andre R. Erler, GPL v3
'''

# imports 
import os, yaml
import numpy as np
import pandas as pd
from glob import glob
from collections import OrderedDict
from scipy import stats as ss
# internal/local imports
from hgs_output import binary 
from geodata.misc import ArgumentError, isNumber
# Graham's package to read HGS binary data: https://github.com/Aquanty/hgs_output
# This package requires Cython-compiled code; on Windows the easiest way to get 
# this to work is via a Wheel installation file, which can be obtained here:
# \\AQUANTY-NAS\share\resources\scripts\Python\hgs_output (inside Aquanty LAN)


## helper functions

def readKister(filepath=None, period=None, resample='1D', missing=None, bias=None, 
               header=3, name='value', lpad=True, lvalues=True, outliers=None):
    ''' read a Kister csv file and slice and resample timeseries '''
    df = pd.read_csv(filepath, header=header, index_col=0, parse_dates=True, names=('time',name))
    # slice
    if period: 
        begin,end = pd.to_datetime(period[0]),pd.to_datetime(period[1])
        df = df[begin:end]
    if resample:
        df = df.resample(resample).mean()
    if period and resample and lpad:
        # extend time axis/index, if necessary, and pad with missing values
        df = df.reindex(pd.date_range(begin,end, freq=resample))
    if outliers is not None:
        # remove values that are more than 'outliers' x standard deviation away from the mean
        df[( ( df[name] - df[name].mean() ) / df[name].std() ).abs() > outliers] = np.NaN
    if bias is not None:
        df += bias 
    if missing:
        df[np.isnan(df)] = missing 
    if lvalues: data = df.values.squeeze()
    else: data = df 
    # return data as pandas dataframe or as numpy array
    return data
    

## functions to write EnKF input files

def writeEnKFini(enkf_folder=None, prefix=None, input_folders=None, glob_pattern='????', lfeedback=True):
    ''' loop over PM and OLF files using some conventions to write IC files for EnKF '''
    if isinstance(input_folders,basestring): input_folders = [input_folders]
    if not os.path.exists(enkf_folder): raise IOError(enkf_folder)
    prefixo = prefix + 'o'
    # loop over OM and OLF
    pm_file = os.path.join(enkf_folder,'inihead.dat'); pm_data = []
    olf_file = os.path.join(enkf_folder,'iniheadolf.dat'); olf_data = []
    ## load data
    # loop over folders and timesteps
    npm = None; nolf = None
    for folder in input_folders:
        glob_path = os.path.join(folder,prefix+'o.head_pm.'+glob_pattern)
        filelist = glob(glob_path)
        if not filelist: 
            raise IOError(glob_path)
        # loop over file list and load data
        for ic_file in filelist:
            idx = int(ic_file[-4:]) # get index number
            reader = binary.IO(prefixo,os.path.dirname(ic_file),idx)
            # extract data and validate PM data
            coords_pm = reader.read_coordinates_pm()
            tmp = coords_pm.shape[0]
            if npm is None: npm = tmp
            elif npm != tmp: 
                raise ValueError("Total number of nodes does not match in input files: {} != {}".format(npm,tmp))
            head_pm = reader.read_var("head_pm", npm)
            pm_data.append(head_pm.values)
            # extract data and validate PM data
            tmp = reader.read_coordinates_olf(coords_pm).shape[0]
            if nolf is None: nolf = tmp
            elif nolf != tmp: 
                raise ValueError("Total number of nodes does not match in input files: {} != {}".format(nolf,tmp))
            head_olf = reader.read_var("head_olf", nolf)
            olf_data.append(head_olf.values)
            # read number of elements for printing later
            if lfeedback:
                nepm = len(reader.read_elements('pm'))
                neolf = len(reader.read_elements('olf'))
    # print number of elements
    if lfeedback:
        print("Number of PM elements: {}".format(nepm))
        print("Number of OLF elements: {}".format(neolf))
        print('')
    # assemble data into arrays and transpose 
    # N.B.: in the EnKF IC file the rows are nodes and the columns are realisations
    pm_data = np.stack(pm_data).squeeze().transpose()
    assert pm_data.shape[0] == npm, pm_data.shape
    if lfeedback: print("Number of PM nodes: {}".format(npm))
    olf_data = np.stack(olf_data).squeeze().transpose()
    assert olf_data.shape[0] == nolf, olf_data.shape
    if lfeedback: print("Number of OLF nodes: {}".format(nolf))
    assert olf_data.shape[1] == pm_data.shape[1]
    nreal = pm_data.shape[1] 
    if lfeedback: print("Number of realizations: {}".format(nreal))
    ## write output files
    fmt = [' %18.0f ']+[' %.18f ']*nreal # node number (first) should be read as integer!
    if lfeedback: print('')
    pm_table = np.concatenate([np.arange(1,npm+1).reshape((npm,1)),pm_data], axis=1)
    # N.B.: we have to add a columns with the node numbers
    np.savetxt(pm_file, pm_table, fmt=fmt)
    if lfeedback: print("Wrote PM IC data to file:\n '{}'.".format(pm_file))
    olf_table = np.concatenate([np.arange(1,nolf+1).reshape((nolf,1)),olf_data], axis=1)
    # N.B.: we have to add a columns with the node numbers
    # N.B.: log10-transform is only done for hydraulic conductivities
    np.savetxt(olf_file, olf_table, fmt=fmt)
    if lfeedback: print("Wrote OLF IC data to file:\n '{}'.".format(olf_file))
    # return file names
    return pm_file, olf_file, nreal
  
  
def writeEnKFbdy(enkf_folder=None, bdy_files=None, filename='flux_bc.dat', mode='deterministic', 
                 nreal=None, scalefactors=None, intermittency=None, default_factor=0.1, lfeedback=True):
    ''' read flux boundary conditions from HGS/Grok .inc files and write an EnKF broundary condition file '''
    if isinstance(bdy_files, dict): bdy_files = OrderedDict(bdy_files) 
    else: raise TypeError(bdy_files)
    if scalefactors is None: scalefactors = dict()
    elif isinstance(scalefactors,(np.integer,int,np.inexact,float)): 
        default_factor = scalefactors
        scalefactors = dict()
    elif not isinstance(scalefactors,dict): raise TypeError(scalefactors)
    if intermittency is None: intermittency = dict()
    elif not isinstance(intermittency,dict): raise TypeError(intermittency)
    if not os.path.exists(enkf_folder): raise IOError(enkf_folder)
    filepath = os.path.join(enkf_folder,filename) # assemble complete path or trunk
    nbdy = len(bdy_files)
    # read boundary flux data
    bdy_data = None
    for i,bdy_file in enumerate(bdy_files.values()):
        data = np.loadtxt(bdy_file,)
        assert data.shape[1] == 2, data.shape 
        if bdy_data is None: bdy_data = np.zeros((data.shape[0],nbdy))
        bdy_data[:,i] = data[:,1]
    ntime = bdy_data.shape[0]
    # write boundary file(s)
    header = [str(nbdy),] + bdy_files.keys() # assemble header 
    if lfeedback:
        print("Number of flux boundary conditions: {}".format(nbdy))
        for head in header[1:]: print(head)
        print("Number of time steps: {}".format(ntime))
    header = [line+'\n' for line in header] # add line breaks
    fmt = '  '.join(['{:18e}']*nbdy)+'\n' # line format
    # there are two modes: deterministic and stochastic
    if mode.lower() == 'deterministic':
        if lfeedback: print("\nWriting 'deterministic' boundary conditions to single file.")
        # all ensemble members get the same input
        # open file and write header
        with open(filepath, 'w') as f: 
            f.writelines(header)
            # loop over actual values
            for i in range(ntime):
                f.write(fmt.format(*bdy_data[i,:]))
        if lfeedback: print("\nWrote flux boundary condition data to file:\n '{}'".format(filepath))
        filelist = filepath
    elif mode.lower() == 'stochastic':
        # every ensemble member gets a different input
        if lfeedback: print("\nWriting 'stochastic' boundary conditions, one file per timestep:")
        if nreal is None: raise ValueError(nreal)
        # variable-dependent scale factors
        bdy_factors = []; bdy_nooccurence = []
        for bdy_file in bdy_files.keys():
            if bdy_file in scalefactors: bdy_factors.append(scalefactors[bdy_file])
            elif default_factor: bdy_factors.append(default_factor)
            else: raise ValueError(bdy_file)
            if bdy_file in intermittency: bdy_nooccurence.append(intermittency[bdy_file])
            else: bdy_nooccurence.append(1.) # always occur
        bdy_factors = np.asarray(bdy_factors).reshape((1,nbdy)).repeat(nreal, axis=0)
        bf_1 = 1-bdy_factors; bf_2 = 2*bdy_factors # shortcuts used below
        # compute actual occurence
        actual_occurence = ( bdy_data > 0 ).sum(axis=0, dtype=bdy_data.dtype) / ntime
        if lfeedback:
            print('Actual Occurence: {}'.format(actual_occurence))
        actual_occurence = actual_occurence.reshape((1,nbdy)).repeat(nreal, axis=0)
        # probability of no occurence, given actual occurence
        bdy_nooccurence = np.asarray(bdy_nooccurence).reshape((1,nbdy)).repeat(nreal, axis=0)
        if lfeedback:
            print('No Occurence: {}'.format(bdy_nooccurence.mean(axis=0)))
        # probability of occurence, given no actual occurence
        bdy_occurence = actual_occurence * bdy_nooccurence / ( 1. - actual_occurence )
        if lfeedback:
            print('New Occurence: {}'.format(bdy_occurence.mean(axis=0)))
        # parameters for random values        
        mean = bdy_data.mean(axis=0)
        std  = bdy_data.std(axis=0)
        # loop over timesteps
        filetrunk = filepath; filelist = []
        for i in range(ntime):
            filepath = filetrunk + '.{:05d}'.format(i+1)
            # prepare data
            scalefactor = np.random.ranf((nreal,nbdy))*bf_2 + bf_1 # uniform random distribution
            rnd_data = bdy_data[i,:].reshape((1,nbdy)).repeat(nreal, axis=0) * scalefactor
            #random_occurence = fake_data[i,:].reshape((1,nbdy)).repeat(nreal, axis=0) * scalefactor
            fake_data = [ss.expon.rvs(loc=mean[i], scale=std[i], size=nreal) for i in range(nbdy)]
            random_occurence = np.stack(fake_data, axis=1)
            assert random_occurence.shape == (nreal,nbdy), random_occurence.shape
            # make random occurences
            lsetZero = np.logical_and( rnd_data > 0, np.random.ranf((nreal,nbdy)) < bdy_nooccurence )
            lcreateNew = np.logical_and( rnd_data == 0, np.random.ranf((nreal,nbdy)) < bdy_occurence )
            rnd_data[lsetZero] = 0
            rnd_data = np.where(lcreateNew, random_occurence, rnd_data)
            # open file and write header
            with open(filepath, 'w') as f: 
                f.writelines(header)
                # loop over actual values
                for j in range(nreal):
                    f.write(fmt.format(*rnd_data[j,:]))
            if lfeedback: print(" '{}'".format(filepath))     
    else:
        raise ValueError(mode)
    # return filepath (or list of files)
    return filelist


def writeEnKFobs(enkf_folder=None, obs_wells=None, filename='obs_head.dat', lfeedback=True, 
                 yaml_file='obs_meta.yaml', lYAML=True):
    ''' write an EnKF observation file with node number, observation error and time-series '''
    if not isinstance(obs_wells, (list,tuple)): raise TypeError(obs_wells)
    if not os.path.exists(enkf_folder): raise IOError(enkf_folder)
    filepath = os.path.join(enkf_folder,filename) # assemble complete path or trunk
    # prepare header
    header = ''; nobs = len(obs_wells); ntime = 0
    print("Number of boundary conditions: {}".format(nobs))
    for i,obs in enumerate(obs_wells):
        if 'error' in obs: error = obs['error']
        else: raise ValueError(obs)
        header += '{:5d}   {:8d}   {:18f}\n'.format(i+1, obs['node'], error**2)
        ntime = max(ntime,len(obs['data']))
    if lfeedback: print(header)
    # assemble time series
    data = np.stack([obs['data'] for obs in obs_wells], axis=1)    
    assert data.shape == (ntime,nobs), data.shape
    print("Number of time steps: {}".format(ntime))
    # write to file
    with open(filepath, 'w') as f:
        f.write(header)
        np.savetxt(f, data, delimiter=' ', fmt=' %.18f ')
    if lfeedback: print("\nWrote observation well data to file:\n '{}'".format(filepath))
    # write obs meta data
    if lYAML:  
        # remove actual binary data from dictionary (only save meta data)
        obs_meta = []
        for obs in obs_wells:
            meta = obs.copy()
            meta['data'] = filepath
            obs_meta.append(meta)
        # write YAML file
        yaml_path = os.path.join(enkf_folder,yaml_file)
        with open(yaml_path, 'w') as yf:
            yaml.dump(obs_meta, yf)
        if lfeedback: print("Also wrote well meta data to YAML file:\n '{}'".format(yaml_path))
    # return filepath
    return filepath


if __name__ == '__main__':
    
    ## settings
    # available range
    folder = 'D:/Data/HGS/SNW/EnKF/TWC/enkf_may/' # folder where files a written
    date_range = ('2017-05-01', '2017-12-31', '1D') # date range for files
    glob_pattern = '00[0123]?'; nreal = 80 # first 40 x 2
    # just december
#     folder = 'D:/Data/HGS/SNW/EnKF/TWC/enkf_december/' # folder where files a written
#     date_range = ('2017-12-01', '2017-12-31', '1D') # date range for files
# #     glob_pattern = '021?' # output timesteps to use for initial conditions; 0215 is Dec. 1st
#     glob_pattern = '02??' # output timesteps to use for initial conditions; 0215 is Dec. 1st
    
    # work folder setup
    if not os.path.exists(folder): os.mkdir(folder)
    input_folder = 'input_data/'
    enkf_folder = os.path.join(folder,input_folder)
    if not os.path.exists(enkf_folder): os.mkdir(enkf_folder)

    
    # execution taskes
    tasks = []
#     tasks += ['test_read_kister']
    tasks += ['write_ic_file']
    tasks += ['write_bdy_file']
    tasks += ['write_obs_file']

    if 'test_read_kister' in tasks:
      
        # test file
        csv_file = 'D:/Data/HGS/SNW/EnKF/Kister/W268-1.csv'
        time_sampling = ('2017-05-01', '2017-12-31', '1D')
        
        # load data
        data = readKister(filepath=csv_file, period=time_sampling[:2], resample=time_sampling[2])
        # test
        datelist = pd.date_range(pd.to_datetime(time_sampling[0]), pd.to_datetime(time_sampling[1]), 
                                 freq=time_sampling[2])
        assert len(datelist) == len(data), (data.shape,len(datelist))
        print((data.shape,len(datelist)))
        
        print('\n===\n')
    
    if 'write_ic_file' in tasks:
      
        # definitions
        prefix = 'prw'
        input_folders =['D:/Data/HGS/SNW/EnKF/TWC/open_value/', 
                        'D:/Data/HGS/SNW/EnKF/TWC/open_raster/']
        #enkf_folder = 'D:/Data/HGS/SNW/EnKF/TWC/enkf_test/input_deterministic/'
        
        # create input files
        pm_file, olf_file, nreal = writeEnKFini(enkf_folder=enkf_folder, prefix=prefix, 
                                                input_folders=input_folders, glob_pattern=glob_pattern)
        if not os.path.exists(pm_file): raise IOError(pm_file)
        if not os.path.exists(olf_file): raise IOError(olf_file)
        
        # create dummy file for initial K values (not read, but still needed)
        open(os.path.join(enkf_folder,'k_dummy.dat'), 'w').close()
                
        print('\n===\n')
    
    if 'write_bdy_file' in tasks:
      
        # definitions
        bdy_files = {'precip.inc': os.path.join(folder,'precip_values.inc'),
                     'pet.inc'   : os.path.join(folder,'pet_values.inc'),}
        scalefactors = {'precip.inc':0.5, 'pet.inc':0.5,}
        intemittency = {'precip.inc':0.3, 'pet.inc':0.,}
        #enkf_folder = 'D:/Data/HGS/SNW/EnKF/TWC/enkf_test/input_deterministic/'        
        
        # create boundary files
        for mode in ('deterministic','stochastic'):
            filelist = writeEnKFbdy(enkf_folder=enkf_folder, bdy_files=bdy_files, mode=mode, 
                                    nreal=nreal, scalefactors=scalefactors, intermittency=intemittency)
            if isinstance(filelist,(list,tuple)):
                for bdy_file in filelist:
                    if not os.path.exists(bdy_file): raise IOError(bdy_file)
            else:
                if not os.path.exists(filelist): raise IOError(filelist)
            
        print('\n===\n')
    
    if 'write_obs_file' in tasks:
      
        # definitions
        #enkf_folder = 'D:/Data/HGS/SNW/EnKF/TWC/enkf_test/input_deterministic/'
        #date_range = ('2017-05-01', '2017-12-31', '1D')   
        datelist = pd.date_range(pd.to_datetime(date_range[0]), pd.to_datetime(date_range[1]), 
                                 freq=date_range[2]) 
        ntime = len(datelist) 
        stderr = 0.01 # observation error
        missing = 99999 # larger than 10,000 indicates missing value
        # actual observation wells
        obs_wells = [
                     # W268-1, 48.52-61.32m, sheet 2-3, possibly 1 (1-2 according to Omar)
                     dict(name='W268-1', z=-35.0, sheet=1, node= 2617, bias=0.24, 
                          csv='D:/Data/HGS/SNW/EnKF/Kister/W268-1.csv'),
                     dict(name='W268-1', z=57.08, sheet=2, node= 5501, bias=0.24, 
                          csv='D:/Data/HGS/SNW/EnKF/Kister/W268-1.csv'),
                     dict(name='W268-1', z=58.08, sheet=3, node= 8385, bias=0.24,
                          csv='D:/Data/HGS/SNW/EnKF/Kister/W268-1.csv'),
                     # W350-2, 104.13-107.13m, sheet 3, possibly 4 (3-4 according to Omar)
                     dict(name='W350-2', z=106.81, sheet=3, node= 7685, bias=-1.65,
                          csv='D:/Data/HGS/SNW/EnKF/Kister/W350-2.csv'),
                     dict(name='W350-2', z=109.93, sheet=4, node=10569, bias=-1.65, 
                          csv='D:/Data/HGS/SNW/EnKF/Kister/W350-2.csv'),
                     # W350-3, 87.33-96.73m, sheet 2 (2-3 according to Omar)
                     dict(name='W350-3', z=91.67, sheet=2, node= 4801, error=0.05, # very unreliable well 
                           csv='D:/Data/HGS/SNW/EnKF/Kister/W350-3.csv', bias=-1.65,)
                     ]
               
        
        # produce open and closed loop observation files
        mode = {'fake_head.dat':False, 'obs_head.dat':True,}
        for filename,lreal in mode.items():
            
            print('')
            for obs_well in obs_wells:
                # add defaults        
                if 'error' not in obs_well: obs_well['error'] = stderr
                if 'bias' not in obs_well: obs_well['bias'] = None       
                #print(obs_well) # feedback without data  
                if lreal:
                    # load actual observation data
                    print("Reading well observations:\n '{}'".format(obs_well['csv']))
                    obs_well['data'] = readKister(filepath=obs_well['csv'], bias=obs_well['bias'],
                                                  period=date_range[:2], resample=date_range[2], 
                                                  missing=missing, lpad=True, lvalues=True,
                                                  outliers=3)
                else:  
                    # create fake/missing data (for equivalent open loop testing)
                    obs_well['data'] = np.ones((ntime,))*missing
            
            # create boundary files
            obs_file = writeEnKFobs(enkf_folder=enkf_folder, obs_wells=obs_wells, filename=filename, 
                                    lYAML=True, )
            if not os.path.exists(obs_file): raise IOError(obs_file)
        
        print('\n===\n')
    
        