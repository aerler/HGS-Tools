'''
Created on Dec 8, 2017

A module that contains a bunch of helper functions for the HGS wrapper class; 
mainly functions that parse configuration and output files.

@author: Andre R. Erler, GPL v3
'''

import os, subprocess, glob

# filename patterns
grok_file = '{PROBLEM:s}.grok' # the Grok configuration file; default includes problem name
# binary files
restart_file = '{PROBLEM:s}o.hen' # restart file name
out_files = '{PROBLEM:s}o.{FILETYPE}.{IDX}' # general pattern for output files
# N.B.: eventually IDX should have the format '{IDX:04d}'
head_files = dict(pm = 'head_pm', olf = 'head_olf', chan = 'head_Chan')
# timeseries files
newton_file = '{PROBLEM:s}o.newton_info.dat' # newton info time series file
water_file = '{PROBLEM:s}o.water_balance.dat' # water balance time series file
hydro_files = '{PROBLEM:s}o.hydrograph.{TAG}.dat' # hydrograph time series file
well_files = '{PROBLEM:s}o.observation_well_flow.{TAG}.dat' # observation well time series file 


## general utility 

# WindowsError is not defined on Linux - need a dummy
try: 
    lWin = True
    WindowsError
except NameError:
    lWin = False
    WindowsError = None
        
## patch symlink on Windows
# adapted from Stack Overflow (last answer: "Edit 2"):
# https://stackoverflow.com/questions/6260149/os-symlink-support-in-windows
if os.name == "nt":
  def symlink_ms(source, link_name):
    import ctypes
    csl = ctypes.windll.kernel32.CreateSymbolicLinkW
    csl.argtypes = (ctypes.c_wchar_p, ctypes.c_wchar_p, ctypes.c_uint32)
    csl.restype = ctypes.c_ubyte
    flags = 1 if os.path.isdir(source) else 0
    if csl(link_name, source.replace('/', '\\'), flags) == 0:
        raise ctypes.WinError()
    if not os.path.exists(link_name):
        raise WindowsError("Creation of Symbolic Link '{}' failed - do you have sufficient Privileges?".format(link_name))
  
  # replace os symlink with this function
  os.symlink = symlink_ms
else: symlink_ms = None # avoid import errors

  
# some named exceptions
class GrokError(Exception):
  ''' Exception indicating an Error in Grok '''
  pass
class HGSError(Exception):
  ''' Exception indicating an Error in HGS '''
  pass


# a platform-independent solution to clear a folder...
def clearFolder(folder, lWin=lWin, lmkdir=True):
    ''' create a folder; if it already exists, remove it and all files and subfolder in it, and recreate it '''
    if os.path.exists(folder):
        if lWin: cmd = 'rmdir /q /s {}'.format(folder.replace('/', '\\'))
        else: cmd = ['rm','-r',folder]
        # N.B.: don't use shutil.rmtree(folder)! It appears that on Windows the use of soft links causes rmtree to 
        #       follow the links and also delete all source files!
        p = subprocess.Popen(cmd, shell=lWin, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # on windows we also need the shell, because rmdir is not an executable in itself
        stdout, stderr = p.communicate(); ec = p.poll(); del stdout
        if ec > 0: 
            raise IOError("The command '{}' to remove the existing directory failed; exit code {}\n{}".format(cmd,ec,stderr))
    # crease folder
    if lmkdir: os.mkdir(folder)


## file collectors

# function to find numeric values of occurences of a numbered file pattern
def numberedPattern(name_pattern, nidx, folder=None):
    ''' function to find numeric values of occurences of a numbered file pattern;
        works only with indices at the end '''
    if folder: 
        pwd = os.getcwd()
        os.chdir(folder)
    pattern = name_pattern + '[0-9]'*nidx 
    filenames = glob.glob(pattern)
    numbers = [ int(filename[-nidx:]) for filename in filenames ]
    if folder: 
        os.chdir(pwd)
    return numbers


# function to collect all time-dependent binary output files
ignore_files = ['{PROBLEM:s}o.ElemK_pm.0001'] # files to ignore when backing up
def binaryFiles(prefix, folder=None, nidx=4, ldict=True, ignore_list=None):
    ''' function to collect all time-dependent binary output files '''
    if ignore_list is None: ignore_list = ignore_files
    if folder: 
        pwd = os.getcwd()
        os.chdir(folder)
    # find binary file pattern
    pattern = out_files.format(PROBLEM=prefix, FILETYPE='*', IDX='[0-9]'*nidx )
    filenames = glob.glob(pattern)
    if folder: 
        os.chdir(pwd)
    # reorganize files
    ignore_list = [ignore.format(PROBLEM=prefix) for ignore in ignore_list]
    if ldict:
        # split head files and store in dictionary
        binary_files = dict(); head_dict= dict()
        for tag,pattern in head_files.items():
            head_tag = 'head_{}'.format(tag)
            binary_files[head_tag] = []
            head_dict[head_tag] = out_files.format(PROBLEM=prefix,FILETYPE=pattern,IDX='')
        binary_files['others'] = [] # collect the rest here 
        # typically binary_files will contain the following lists: head_pm, head_olf, head_chan, others 
        for filename in filenames:
            lmatch = False
            for head_tag,head_file in head_dict.items():
                if filename.startswith(head_file): 
                  binary_files[head_tag].append(filename); lmatch = True; break
            if not lmatch and not filename in ignore_list: 
                binary_files['others'].append(filename)
    else:
        # only remove unwanted files
        binary_files = [filename for filename in filenames if not filename in ignore_list]
    return binary_files
#TOTO: remove hard-coding of file patterns and add wrapper as class method to Grok

# function to collect all time-series output files
def timeseriesFiles(prefix, folder=None, ldict=True, llogs=True, lcheck=True):
    ''' function to collect all time-series output files '''
    if folder: 
        pwd = os.getcwd()
        os.chdir(folder)
    ts_dict = dict()
    # log files
    if llogs:
        log_files = ['log.grok','log.hgs_run']
        if lcheck: 
            for log_file in log_files:
                if not os.path.exists(log_file): 
                    raise IOError(log_file)
        ts_dict['logs'] = log_files        
    # special files
    special_files = [ name.format(PROBLEM=prefix) for name in [newton_file, water_file,] ]
    if lcheck: 
        for special_file in special_files:
            if not os.path.exists(special_file): 
                raise IOError(special_file)
    ts_dict['special'] = special_files
    # hydrographs
    pattern = hydro_files.format(PROBLEM=prefix, TAG='*')
    ts_dict['hydrographs']= glob.glob(pattern)
    # observation wells
    pattern = well_files.format(PROBLEM=prefix, TAG='*') 
    ts_dict['wells']= glob.glob(pattern)
    if folder: 
        os.chdir(pwd)
    # reorganize
    if not ldict:
        flat_list = []
        for files in ts_dict.values(): flat_list.extend(files)
        return flat_list
    else:
        return ts_dict
#TOTO: remove hard-coding of file patterns and add wrapper as class method to Grok


## file parsere

# helper function to recursively parse a Grok file with includes
def parseGrokFile(filepath, line_list):
    ''' recursively parse files which include other files, starting with 'filepath', 
        and append their lines to 'line_list' '''
    # change working directory locally to resolve relative path'
    pwd = os.getcwd() # need to go back later!
    os.chdir(os.path.dirname(filepath))
    with open(os.path.basename(filepath), 'r') as f:
        lines = f.readlines() # load lines all at once (for performance)
    # loop over lines, clean them and find includes
    for line in lines:
        line = line.strip() # remove which space 
        # skip comments and empty lines
        if not line.startswith('!'):
            line = line.replace('\\','/') # enforce Unix folder convention
            # scan for include statements (except precip.inc and pet.inc) 
            if line.startswith('include') and not ( line.endswith('precip.inc') or line.endswith('pet.inc') ):
                # N.B.: apparently HGS/Grok is case-sensitive... 
                # figure out new file path
                incpath = line[7:].strip()
                if not os.path.lexists(incpath):
                    raise IOError("Unable to open include file '{}'.".format(incpath))
                # initiate recursion
                parseGrokFile(incpath, line_list)
            else:
                # append line to line_list
                line_list.append(line)
                # N.B.: only convert non-path statements to lower
    # now we are done - we don't return anything, since line_list was modified in place
    os.chdir(pwd) # this is important: return to previous working directory!!!


if __name__ == '__main__':

    test = 'filelists'
    
    if test == 'filelists':
        
        from hgs.HGS import root_folder, prefix_file
        exp_folder = 'GRW/grw2/NRCan/clim_30/hgs_run_v3/'
        folder = os.path.join(root_folder,exp_folder)
        if not os.path.exists(folder): raise IOError(folder)
        # read prefix from simulation
        with open(os.path.join(folder,prefix_file), 'r') as pfx:
            prefix = ''.join(pfx.readlines()).strip()
        
        # count numbered patter
        print('\n  Numbered Pattern')
        files = numberedPattern(prefix+'o.head_olf.', folder=folder, nidx=4)
        print(files)
        
        # list binary files
        print('\n  Binary Files')
        files = binaryFiles(prefix=prefix, folder=folder, nidx=4, ldict=True)
        for key,val in files.items():
            print('{}: {}'.format(key,len(val)))
            print(val)
        
        # list time-series files
        print('\n  Time-series Files')
        files = timeseriesFiles(prefix=prefix, folder=folder, llogs=True, lcheck=True, ldict=True)
        for key,val in files.items():
            print('{}: {}'.format(key,len(val)))
            print(val)
        