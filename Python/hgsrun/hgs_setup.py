'''
Created on Jul 31, 2016

A Python class (Grok) that facilitates the configuration of an HGS simulation via grok and a child class
that also handles folder setup (based on template) and execution of HGS.

@author: Andre R. Erler, GPL v3
'''

# external imports
import numpy as np
import os, shutil
import subprocess # launching external programs
import bisect
# internal imports
from hgsrun.input_list import generateInputFilelist, resolveInterval
from hgsrun.misc import lWin, symlink_ms, GrokError, HGSError, timeseriesFiles,\
  binaryFiles
from hgsrun.misc import parseGrokFile, clearFolder, numberedPattern
from geodata.misc import ArgumentError
from utils.misc import tail


## patch symlink on Windows
# adapted from Stack Overflow (last answer: "Edit 2"):
# https://stackoverflow.com/questions/6260149/os-symlink-support-in-windows
if os.name == "nt":
  os.symlink = symlink_ms # replace os symlink with this function


## a class to handle Grok for HGS simulations
class Grok(object):
  '''
    A class that loads a grok configuration file into memory, provides functions for editing,
    saving the file, and running Grok.
  '''
  batchpfx = 'batch.pfx' # a file that defines the HGS problem name for Grok
  grok_dbg = 'grok.dbg' # Grok debug output
  grok_file = '{PROBLEM:s}.grok' # the Grok configuration file; default includes problem name
  restart_file = '{PROBLEM:s}o.hen' # restart file name
  out_files = '{PROBLEM:s}o.{FILETYPE}.{{IDX:04d}}' # porous media head output (for restart)
  pm_tag = 'head_pm'
  olf_tag = 'head_olf'
  chan_tag = 'head_Chan'
  newton_file = '{PROBLEM:s}o.newton_info.dat' # newton info time series file
  water_file = '{PROBLEM:s}o.water_balance.dat' # water balance time series file
  hydro_file = '{PROBLEM:s}o.hydrograph.{{TAG}}.dat' # hydrograph time series file
  well_file = '{PROBLEM:s}o.observation_well_flow.{{TAG}}.dat' # observation well time series file 
  lchannel = None # whether or not we are using 1D channels  
  rundir = None # folder where the experiment is set up and executed
  project = None # a project designator used for file names
  problem = None # the HGS problem name (defaults to project)
  input_mode = None # type of simulation (forcing data): stead-state, periodic, transient
  input_interval = None # update interval for climate forcing
  input_vars = 'PET' # input variable configuration
  input_prefix = None # prefix for input files
  input_folder = '../climate_forcing' # default folder for input data
  output_interval = None # how many evenly spaced restart files to write (supports multiple nested intervals)
  starttime = 0 # model time at initialization
  runtime = None # run time for the simulations (in seconds)
  length = None # run time in multiples of the interval length
  _lines = None # list of lines in file
  _sourcefile = None # file that the configuration was read from
  _targetfile = None # file that the configuration is written to
  
  def __init__(self, rundir=None, project=None, problem=None, starttime=0, runtime=None, length=None, 
               output_interval='default', input_mode=None, input_interval=None, input_vars='PET', 
               input_prefix=None, input_folder='../climate_forcing', pet_folder=None, lcheckdir=True,
               grok_bin='grok.exe'):
    ''' initialize a Grok configuration object with some settings '''
    if lcheckdir and not os.path.isdir(rundir): raise IOError(rundir)
    # determine end time in seconds (begin == 0) or number of intervals (length)
    length, runtime = resolveInterval(length=length, end_time=runtime, interval=input_interval)
    # assign class variables
    self.grok_bin = grok_bin # Grok executable: first try local folder, then $HGSDIR/bin/
    self.rundir = rundir
    self.project = project
    self.problem = project if problem is None else problem
    self.starttime = starttime
    self.runtime = runtime
    self.length = length
    self.grok_file = self.grok_file.format(PROBLEM=self.problem) # default name
    self.restart_file = self.restart_file.format(PROBLEM=self.problem) # default name
    # binary files for restart
    self.pm_files = self.out_files.format(PROBLEM=self.problem, FILETYPE=self.pm_tag) # default name
    self.olf_files = self.out_files.format(PROBLEM=self.problem, FILETYPE=self.olf_tag) # default name
    self.chan_files = self.out_files.format(PROBLEM=self.problem, FILETYPE=self.chan_tag) # default name
    # time-series files (for later concatenatenation)
    self.newton_file = self.newton_file.format(PROBLEM=self.problem)
    self.water_file = self.water_file.format(PROBLEM=self.problem)
    self.hydro_file = self.hydro_file.format(PROBLEM=self.problem)
    self.well_file = self.well_file.format(PROBLEM=self.problem)
    # input configuration
    self.setInputMode(input_mode=input_mode, input_interval=input_interval,
                      input_vars=input_vars, input_prefix=input_prefix,
                      input_folder=input_folder, pet_folder=pet_folder)
    # output configuration
    self.resolveOutput(output_interval=output_interval)
  
  def readConfig(self, filename=None, folder=None):
    ''' Read a grok configuration file into memory (or a template to start from). '''    
    filename = filename or self.grok_file; self.grok_file = filename
    filename = '{:s}/{:s}'.format(folder or self.rundir, filename) # prepend folder
    if not os.path.isfile(filename): 
        raise IOError("Grok configuration file not found: '{}'".format(filename))
    self._sourcefile = filename # use  different file as template
    # read source file (and recurse into includes)
    self._lines = [] # initialize empty list (will be extended in-place)
    parseGrokFile(filename, self._lines)
    # check if we are dealing with channels
    self.lchannel = 'channel' in self._lines    
    # return exit code
    return 0 if isinstance(self._lines,list) else 0
      
  def writeConfig(self, filename=None, folder=None):
    ''' Write the grok configuration to a file in run dir. '''    
    filename = filename or self.grok_file; self.grok_file = filename # or default name    
    filename = '{:s}/{:s}'.format(folder or self.rundir,filename) # prepend run dir
    self._targetfile = filename # use  different file as template
    # move existing file to backup
    if os.path.isfile(filename): shutil.move(filename, '{:s}.backup'.format(filename))
    # write configuration to file
    with open(filename, 'w') as tgt: # with-environment should take care of closing the file
        tgt.write('\n'.join(self._lines)+'\n')
        # N.B.: this is necessary, because our list does not have newlines and Python does not add them...
    # return exit code
    return 0 if os.path.isfile(filename) else 1
      
  def setParam(self, param, value, formatter=None, after=None, start=0):
    ''' edit a single parameter, based on the assumption that the parameter value follows in the 
        line below the one where the parameter name appears (case in-sensitive); format is a
        regular format string used to write the value
        N.B.: the method is not able to detect nested lists/vector-valued parameters; it is only possible to 
              manipulate inner-most lists, because the insertion is terminated at the first 'end' '''
    start = self._lines.index(after, start) if after else start # search offset for primary paramerter
    formatter = str if formatter is None else formatter.format # apply appropriate formatting
    if isinstance(value,(tuple,list)):
      # vector-valued parameter
      values = [formatter(val) for val in value] # apply formatter to list items
      start = self._lines.index(param, start) # find begin of vector statement
      end = self._lines.index('end', start) # 'end' marks the end of a vector statement
      # insert list of vector entries into file, line by line
      self._lines = self._lines[:start+1] + values + self._lines[end:]
      # N.B.: the Grok class is not able to detect nested lists/vector-valued parameters; it is only possible to 
      #       manipulate inner-most lists, because the insertion is terminated at the first 'end'
    else:
      # single-valued parameter
      self._lines[self._lines.index(param, start)+1] = formatter(value) # replace value in list of lines

  def getParam(self, param, dtype=None, llist=None, after=None, start=0, lindex=False, lerror=True):
    ''' read a single parameter, based on the assumption that the parameter value follows in the 
        line below the one where the parameter name appears (case sensitive); dtype is a 
        numpy data type to which the value string is cast; a list can be inferred if values of the 
        proper type follow in consecutive lines and are terminated by and 'end' '''
    if after is not None: start = self._lines.index(after, start) # search offset for primary paramerter
    if isinstance(dtype,basestring): dtype = getattr(np,dtype) # convert to given numpy dtype
    elif dtype is None: dtype = str # read as string data type as default
    # lists of strings cannot be detected properly
    if llist is None and ( np.issubdtype(dtype,np.string_) or np.issubdtype(dtype,np.unicode_) ): llist = False
    # find entry
    try: 
      i = self._lines.index(param, start)+1
    except ValueError:
      if lerror: raise
      else: return None,None if lindex else None
    # different handling for scalars and lists
    if llist is False:
      value = dtype(self._lines[i]) # just a single value
    else:
      lterm = False; values = [] # indicate if list is complete
      while i < len(self._lines) and not lterm:
        value = self._lines[i] # lines should be stripped (but still case-sensitive!)
        if value == 'end': 
          lterm = True # proper termination of list
        elif value == '':
          if not llist and len(values) == 1: # if it is a list, just skip empty lines... 
              lterm = True # proper termination for a scalar
        else:
          try: 
            value = dtype(value) # convert string and append
            values.append(value) # append to list
          except ValueError: 
            if llist: raise ValueError("Illegal list termination '{}' for parameter '{}'".format(value, param))
            else: lterm = True # terminate with invalid value
        i += 1 # increment to next line
      i -= 1 # decrement to pass on position
      # check results
      if value == 'end': value = values # legitimately a list
      elif llist is None and len(values) == 1: value = values[0] # auto-detect scalar
      else: raise ValueError(value)
      # N.B.: this is very fragile, because it assumes there are no empty lines in the list
    # return value and optionally index position
    return (value,i) if lindex else value

  def replaceValue(self, old, new, formatter=None, after=None, start=0):
    ''' repalce a value with a new value (does not work for lists) '''
    if formatter: 
      new = formatter.format(new) # apply appropriate formatting
      old = formatter.format(old) # apply appropriate formatting
    else:
      new = str(new); old = str(old)
    if after is not None: start = self._lines.index(after, start)# search offset for primary paramerter
    self._lines[self._lines.index(old, start)] = new # replace value in list of lines
    
  def setParams(self, **params):
    ''' edit a bunch of parameters, which are defined as key/values pairs using self.editParam;
        note that there is no validation and 'format' is only supported for single edits '''
    for param,value in params.iteritems():
      self.setParam(param, value)
      
  def remParam(self, param, llist=False, after=None, start=0, lerror=True, lall=False):
    ''' find a parameter definition and comment it out, including the value of the parameter; the
        value is assumed to follow directly after the definition; lists cannot be auto-detected as 
        in consecutive values terminated by an 'end'; lall indicates that all occurences should be
        commented out; lall returns a list of occurences, while normally only the last index of 
        the occurence is returned. '''
    # lall indicates wether all occurences or just the first one should be commented out
    if after is not None: start = self._lines.index(after, start) # search offset for primary paramerter
    # find entry
    try: 
      i = self._lines.index(param, start)
    except ValueError:
      if lerror: raise
      else: return None
    # comment out command
    self._lines[i] = '! '+self._lines[i]
    i += 1 # increment index
    # comment out values; different handling for scalars and lists
    if llist is False:
      self._lines[i] = '! '+self._lines[i] # just one value
    else:
      lterm = False # indicate if list is complete
      while i < len(self._lines) and not lterm:
        self._lines[i] = '! '+self._lines[i] # comment line by line
        value = self._lines[i] # lines should be stripped (but still case-sensitive!)
        if value == 'end': lterm = True # proper termination of list
        i += 1 # increment to next line
    # find more occurences recursively
    if lall: 
        ii = self.remParam(param, llist=llist, after=None, start=i, lerror=False, lall=True)
        if ii is None: 
            i = [i] # terminate recursion - return current occurence as last in the list
        else:
            i = [i] + ii # prepend this current occurence to list from recursion
    # return index position 
    return i
  
  def remParams(self, *params, **kwargs):
    ''' comment out a list of parameters (and their values) '''
    for param in params:
        self.remParam(param, **kwargs)
      
  def changeICs(self, ic_pattern=None):
    ''' change the initial condition files in Grok using either .hen/restart or head/output files'''
    os.chdir(self.rundir)
    # check which type of restart file we are dealing with
    if not isinstance(ic_pattern,basestring): 
        raise TypeError(ic_pattern)
    lhenf = ic_pattern.endswith('.hen')
    lheadf = '{FILETYPE}' in ic_pattern
    # expand pattern and make sure files are present
    ic_pattern = ic_pattern.format(FILETYPE='{FILETYPE}',RUNDIR=self.rundir)
    if lheadf and lhenf: 
        raise ArgumentError('Filetype expansion is not supported for .hen files:',ic_pattern)
    elif lhenf:
        ic_file_hen = ic_pattern # for consistency
        if not os.path.exists(ic_file_hen): IOError(ic_file_hen)
    elif lheadf:
        ic_file_pm = ic_pattern.format(FILETYPE=self.pm_tag)
        if not os.path.exists(ic_file_pm): IOError(ic_file_pm)
        ic_file_olf = ic_pattern.format(FILETYPE=self.olf_tag)
        if not os.path.exists(ic_file_olf): IOError(ic_file_olf)
        if self.lchannel:
            ic_file_chan = ic_pattern.format(FILETYPE=self.chan_tag)
            if not os.path.exists(ic_file_chan): IOError(ic_file_chan) 
    else: raise ArgumentError('Restart filetype not recognized:',ic_pattern)
    # see if we are using .hen restart files
    hen_file,i = self.getParam('restart file for heads', dtype=str, llist=False, lindex=True, lerror=False)
    if hen_file:
        # this Grok file/simulation uses .hen files for initialization
        if lhenf:
            self._lines[i] = ic_file_hen
        else:
            # remove all occurrences of 'restart file for heads'
            idx = self.remParam('restart file for heads', llist=False, start=0, lerror=True, lall=True)
            if len(idx) < 1: raise GrokError(idx)
            # append required 'initial head from output file' sections at the end
            domains_files = [('porous media',ic_file_pm),('surface',ic_file_olf),]
            if self.lchannel: domains_files += [('channel',ic_file_chan)]
            for domain,ic_file in domains_files:
                lines = ["! restart file '{}' for domain {} added by restart function".format(ic_file,domain),'',
                         'use domain type', domain, '','clear chosen nodes', 'choose nodes all','',
                         'initial head from output file', ic_file, '',
                         'clear chosen nodes', '',]
                self._lines += lines
    else:
        # this Grok file/simulation uses regular head output files for initialization
        if lheadf:
            lpm = False; lolf = False; lchan = not self.lchannel; i = 0; ini_file = ''
            while ini_file is not None:
                if self.pm_tag in ini_file: 
                    self._lines[i] = ic_file_pm; lpm = True
                elif self.olf_tag in ini_file: 
                    self._lines[i] = ic_file_olf; lolf = True
                elif self.lchannel and self.chan_tag in ini_file: 
                    self._lines[i] = ic_file_chan; lchan = True
                ini_file, i = self.getParam('initial head from output file', llist=False, start=i, lindex=True, lerror=False)
            if not ( lpm and lolf and lchan ): 
                raise GrokError("Not all IC filetypes have been found/changed (PM,OLF,Chan): ",lpm,lolf,lchan)
        else:
            # remove all occurrences of 'initial head from output file'
            idx = self.remParam('initial head from output file', llist=False, start=0, lerror=True, lall=True)
            if len(idx) < 2: raise GrokError(idx)
            # append 'restart file for heads' section at the end
            lines = ["! restart file '{}' added by restart function".format(ic_file_hen),'',
                     'clear chosen nodes', 'choose nodes all', '',
                     'restart file for heads', ic_file_hen, '',
                     'clear chosen nodes', '',]
            self._lines += lines
      
  def resolveOutput(self, output_interval):
    ''' method to check and determine default output intervals '''
    # detect output interval based on length and input_interval
    if isinstance(output_interval, basestring):
      if output_interval.lower() == 'default':
        # monthly for the last year and yearly otherwise
        if self.input_interval == 'monthly':
            self.output_interval = (int(self.length/12.),12)
        elif self.input_interval == 'daily':
            self.output_interval = (int(self.length/365.),12)
        else: raise NotImplementedError(self.input_interval)
        # sanity check
        if self.output_interval[0] < 1: 
            raise NotImplementedError("No default 'output_interval' for simulation times less than one year: {}".format(self.output_interval))
      elif output_interval.lower() == 'yearly':
        if self.input_interval == 'monthly':
            self.output_interval = (int(self.length/12.),)
        elif self.input_interval == 'daily':
            self.output_interval = (int(self.length/365.),)
        else: raise NotImplementedError(self.input_interval)
      elif output_interval.lower() == 'monthly':
        if self.input_interval == 'monthly':
            self.output_interval = (self.length,)
        elif self.input_interval == 'daily':
            self.output_interval = (int(self.length/(365./12.)),)
        else: raise NotImplementedError(self.input_interval)
      else: raise NotImplementedError(output_interval)
    elif isinstance(output_interval,(int,np.integer)):
        self.output_interval = (output_interval,)
    elif isinstance(output_interval,(tuple,list)):
        self.output_interval = output_interval
    else:
        raise TypeError(output_interval)
      
  def setRuntime(self, runtime=None, starttime=0, output_interval=None):
    ''' set the run time of the simulations, as well as the initial time and output times (model time in seconds) '''
    self.runtime = runtime if runtime is not None else self.runtime
    self.starttime = starttime if starttime is not None else self.starttime
    if output_interval is not None: self.resolveOutput(output_interval)
    if self._lines:
      # set start and run times
      self.setParam('initial time', self.starttime, formatter='{:e}', )
      # N.B.: there is no actual "run time" parameter
      # figure out output/restart times
      outtimes = []; outinit = self.starttime
      for nout in self.output_interval:
          if not isinstance(nout,(int,np.integer)): raise TypeError(nout)
          if nout < 1: raise ValueError("An 'output_interval' of {} < 1 is illegal!".format(nout))
          timedelta = self.runtime - outinit
          if nout == 1: tmp = [outinit + timedelta] 
          else: tmp = [ outinit + ( timedelta * float(r) / float(nout) ) for r in xrange(1,nout)]
          outtimes.extend(tmp) # append new/refined list
          outinit = tmp[-1] # use last value as starting point for next refinement interval
      outtimes.append(self.runtime) # append final time for output
      # change grok file
      self.setParam('output times', outtimes, formatter='{:e}', )
    else: 
        raise GrokError("No Grok file loaded.")
    # check
    return 0 if self.getParam('initial time', float, llist=False) == self.starttime else 1
  
  
  def setInputMode(self, input_mode=None, input_interval=None, input_vars='PET', input_prefix=None, 
                   input_folder='../climate_forcing', pet_folder=None):
    ''' set the type of the simulation: mean/steady-state, climatology/periodic, time-series/transient '''
    # resolve type of input data 
    input_mode = input_mode.lower()
    if input_mode in ('mean','steady','steady-state') or input_mode[-5:] == '-mean': input_mode = 'steady-state'
    elif input_mode in ('clim','peri','climatology','periodic'): input_mode = 'periodic'
    elif input_mode in ('time-series','timeseries','trans','transient'): input_mode = 'transient'
    elif input_mode in ('quasi-trans','quasi-transient'): input_mode = 'quasi-transient'
    self.input_mode = input_mode
    # set input interval
    input_interval = input_interval.lower()
    if input_interval[:5].lower() == 'month': input_interval = 'monthly'
    elif input_interval[:3].lower() == 'day': input_interval = 'daily'
    else: raise NotImplementedError(input_interval)
    self.input_interval = input_interval
    # set other variables
    if input_vars.upper() not in ('PET','WRFPET','NET'): raise NotImplementedError(input_vars)
    if input_mode == 'quasi-transient' and input_vars.upper() not in ('PET','WRFPET'): 
        raise ArgumentError("Quasi-transient input requires a PET variable.") 
    # N.B.: the input config NET actually requires changes to the grok file, which is currently not done
    self.input_vars = input_vars
    self.input_prefix = input_prefix
    self.input_folder = input_folder
    self.pet_folder = input_folder if pet_folder is None else pet_folder
  
  def generateInputLists(self, input_vars=None, input_prefix=None, input_folder=None, pet_folder=None,
                         lvalidate=True, axis='iTime', lcenter=True, l365=True, lFortran=True):
    ''' generate and validate lists of input files and write to appropriate files '''
    input_vars = self.input_vars if input_vars is None else input_vars
    input_prefix = self.input_prefix if input_prefix is None else input_prefix
    input_folder = self.input_folder if input_folder is None else input_folder
    pet_folder = self.pet_folder if pet_folder is None else pet_folder
    # generate default input vars
    if isinstance(input_vars,basestring):
      if input_vars.upper() == 'PET': # liquid water + snowmelt & PET as input
        input_vars = dict(precip=('rainfall','rain','liqwatflx'),  
                          pet=('pet','potential evapotranspiration','pet'))
      elif input_vars.upper() == 'WRFPET': # liquid water + snowmelt & PET from WRF as input
        input_vars = dict(precip=('rainfall','rain','liqwatflx'),  
                          pet=('pet','potential evapotranspiration','pet_wrf'))
      elif input_vars.upper() == 'NET': # liquid water + snowmelt - ET as input
        input_vars = dict(precip=('rainfall','rain','waterflx'),)
        if self.getParam('name', after='potential evapotranspiration', llist=False).upper() == 'PET':
            raise NotImplementedError("Cannot use 'NET' input type if 'PET' input is defined in grok file.")
      else: raise ArgumentError("Invalid or missing input_vars: {}".format(input_vars))
    elif not isinstance(input_vars, dict): raise TypeError(input_vars)
    # iterate over variables and generate corresponding input lists
    ec = 0 # cumulative exit code
    for varname,val in input_vars.iteritems():
      grokname,vartype,wrfvar = val
      filename = '{}.inc'.format(varname)
      if self.getParam('name', after=vartype, llist=False).lower() != grokname:
            raise GrokError("No entry for boundary condition type '{}'/'{}' found in grok file!".format(vartype,grokname))
      self.setParam('time raster table', 'include {}'.format(filename), after=vartype)      
      # special handlign for quasi-transient forcing based on variable
      if self.input_mode == 'quasi-transient':
          input_mode = 'periodic' if varname == 'pet' else 'transient' 
          actual_input_folder = pet_folder if varname == 'pet' else input_folder
      else: 
          input_mode = self.input_mode; actual_input_folder = input_folder
      # select interval and output format
      if self.input_interval == 'monthly':
        if input_mode == 'steady-state': input_pattern = 'iTime_{IDX:d}' # IDX will be substituted
        elif input_mode == 'periodic': input_pattern = 'iTime_{IDX:02d}' # IDX will be substituted
        elif input_mode == 'transient': input_pattern = 'iTime_{IDX:03d}' # IDX will be substituted
        else: raise GrokError(self.input_mode)
        # N.B.: this is a very ugly hack - I don'e have a better idea at the moment, since
        #       we don't know the original length of the time series
      else:
        length = self.length + 1 if lFortran else self.length
        input_pattern = '{0:s}_{{IDX:0{1:d}d}}'.format(axis, int(np.ceil(np.log10(length)))) # number of digits
      input_pattern = '{0}_{1}.asc'.format(wrfvar,input_pattern)
      if input_prefix is not None: input_pattern = '{0}_{1}'.format(input_prefix,input_pattern)
      # write file list
      lec = generateInputFilelist(filename=filename, folder=self.rundir, input_folder=actual_input_folder, 
                                  input_pattern=input_pattern, lcenter=lcenter, lvalidate=lvalidate, 
                                  units='seconds', l365=l365, lFortran=lFortran, interval=self.input_interval, 
                                  end_time=self.runtime, mode=input_mode)
      ec += 0 if lec else 1          
    # return exit code
    return ec
  
  def runGrok(self, executable=None, logfile='log.grok', batchpfx=None, lerror=True, lcompress=True, ldryrun=False):
    ''' run the Grok executable in the run directory '''
    pwd = os.getcwd() # save present workign directory to return later
    os.chdir(self.rundir) # go into run Grok/HGS folder
    self.grok_bin = executable if executable is not None else self.grok_bin
    self.batchpfx = batchpfx if batchpfx is not None else self.batchpfx
    if not os.path.lexists(self.grok_bin): 
      raise IOError("Grok executable '{}' not found.\n".format(self.grok_bin))
    # create batch.pfx file with problem name for batch processing
    with open(self.batchpfx, 'w+') as bp: bp.write(self.problem) # just one line...
    # run executable while logging output
    with open(logfile, 'w+') as lf: # output and error log
      if ldryrun:
        lf.write('Dry-run --- no execution')
        lec = True # pretend everything works
      else:
        # run Grok
        subprocess.call([os.path.abspath(self.grok_bin)], stdout=lf, stderr=lf)
        # parse log file for errors
        lec = ( tail(lf, n=3)[0].strip() == '---- Normal exit ----' )
        # i.e. -3, third line from the end (different from HGS)
    # compress grok debug output
    if lcompress and lec:
      with open(logfile, 'a') as lf: # output and error log
        try:
          if os.path.exists(self.grok_dbg):
            # compress using gzip (single file; no shell expansion necessary)
            subprocess.call(['gzip',self.grok_dbg], stdout=lf, stderr=lf)
            if os.path.exists(self.grok_dbg+'.gz'): 
              lf.write('\nCompressed Grok debug output ({}.gz).\n'.format(self.grok_dbg))
            else: raise IOError # just trigger exception (see below)
          else:
            lf.write('\nNo Grok debug output found ({}).\n'.format(self.grok_dbg))
        except:
          lf.write('\nGrok debug output compression failed ({}).\n'.format(self.grok_dbg))
    os.chdir(pwd) # return to previous working directory
    if lerror and not lec: 
      raise GrokError("Grok failed; inspect log-file: {}\n  ('{}')\n".format(logfile,self.rundir))
    return 0 if lec else 1
            
      
class HGS(Grok):
  '''
    A child class of Grok that can also set up the entire run folder and launch an HGS instance;
    otherwise the same as Grok.
  '''
  template_folder = None # temlate for rundir setup
  linked_folders = None # list of folders that should be linked rather than copied
  rundirOK  = None # indicate if rundir setup was successfule
  configOK  = None # indicate if Grok configuration was successful
  GrokOK    = None # indicate if Grok ran successfully
  pidxOK    = None # indicate if parallel index configuration was successful
  HGSOK     = None # indicate if HGS ran successfully
  pidx_file = 'parallelindx.dat' # file with parallel execution settings 
  lindicators = True # use indicator files (default: True)
  lrestart  = False # whether or not this is a restart run (to complete an interrupted run)
  ic_files  = None # pattern for initial condition files (path can be expanded)
  
  def __init__(self, rundir=None, project=None, problem=None, runtime=None, length=None, output_interval='default',
               input_mode=None, input_interval=None, input_vars='PET', input_prefix=None, pet_folder=None,
               input_folder='../climate_forcing', template_folder=None, linked_folders=None, NP=1, lindicator=True,
               grok_bin='grok.exe', hgs_bin='phgs.exe', lrestart=False):
    ''' initialize HGS instance with a few more parameters: number of processors... '''
    # call parent constructor (Grok)
    super(HGS,self).__init__(rundir=rundir, project=project, problem=problem, runtime=runtime, 
                             output_interval=output_interval, input_vars=input_vars, input_prefix=input_prefix,
                             input_mode=input_mode, input_interval=input_interval, pet_folder=pet_folder,
                             input_folder=input_folder, length=length, lcheckdir=False, grok_bin=grok_bin,)
    self.hgs_bin = hgs_bin # HGS executable: first try local folder, then $HGSDIR
    self.template_folder = template_folder # where to get the templates
    # prepare linked folders
    if linked_folders is None: linked_folders = ('etprop', 'gb', 'icbc', 'prop', 'soil', # original 
                                                 'grid', 'init_con', 'K_maps', 'landcover', 'mprops', 'node_lists', # extended
                                                 'pks_table', 'retent_tables', 'inc') # even more...
    linked_folders = tuple(lf[:-1] if lf[-1] == '/' else lf for lf in linked_folders) # trim slash
    self.linked_folders = linked_folders
    # N.B.: these folders just contain static data and do not need to be replicated
    self.NP = NP # number of processors
    self.lindicators = lindicator # use indicator files
    self.lrestart = lrestart # complete an interrupted run
    
  def setupRundir(self, template_folder=None, bin_folder='{HGSDIR:s}', loverwrite=None, lschedule=True):
    ''' copy entire run folder from a template folder and link executables '''
    template_folder = self.template_folder if template_folder is None else template_folder
    if template_folder is None: raise ValueError("Need to specify a template path.")
    if not os.path.exists(template_folder): raise IOError(template_folder)
    if bin_folder: bin_folder = bin_folder.format(HGSDIR=os.getenv('HGSDIR'))
    # check for restart
    if self.lrestart and self.lindicators:
      if os.path.exists('{}/SCHEDULED'.format(self.rundir)) or os.path.exists('{}/ERROR'.format(self.rundir)):
        self.lrestart = False # this mean it probably never ran!
    # clear existing directory
    if loverwrite is None: 
      loverwrite = not self.lrestart # i.e. don't clean for a restart
    if loverwrite and os.path.isdir(self.rundir):
      clearFolder(self.rundir, lWin=lWin, lmkdir=False)
      # N.B.: rmtree is dangerous, because if follows symbolic links and deletes contents of target directories!
    # copy folder tree
    if not os.path.isdir(self.rundir): shutil.copytree(template_folder, self.rundir, symlinks=True,
                                                       ignore=shutil.ignore_patterns('*.grok*',*self.linked_folders))
    # place link to template
    os.symlink(template_folder, os.path.join(self.rundir,'template'))
    # symlinks to static folder that are not copied (linked_folders)
    if self.linked_folders:
      for lf in self.linked_folders:
        link_source = os.path.join(template_folder,lf)
        if os.path.exists(link_source):
          #print('linking {}: {}'.format(lf,os.path.join(self.rundir,lf)))
          os.symlink(link_source, os.path.join(self.rundir,lf))
    # copy or put links to executables in place
    for exe in (self.hgs_bin, self.grok_bin):
      local_exe = os.path.join(self.rundir,exe)
      if os.path.lexists(local_exe): os.remove(local_exe)
      bin_path = os.path.join(template_folder,exe)
      if os.path.exists(bin_path):
        os.symlink(bin_path, local_exe)
      elif bin_folder and os.path.exists(bin_folder):
        bin_path = os.path.join(bin_folder,exe)    
        if os.path.exists(bin_path):
          os.symlink(bin_path, local_exe)
        else:
          raise IOError("Executable file '{}' not found in bin folder.\n ('{}') ".format(exe,bin_folder))
      else:
        raise IOError("Executable file '{}' not found in template folder (no bin folder found).\n ('{}') ".format(exe,template_folder))
      # check executables
      if os.path.islink(local_exe):
        if not os.path.exists(local_exe): 
          raise IOError("Link to executable '{}' in run folder is broken.\n ('{}') ".format(exe,self.rundir))
      elif not os.path.isfile(local_exe): 
        raise IOError("Executable file '{}' not found in run folder.\n ('{}') ".format(exe,self.rundir)) 
    # set rundir status
    self.rundirOK = os.path.isdir(self.rundir)
    if self.lindicators and lschedule: 
      if self.rundirOK: open('{}/SCHEDULED'.format(self.rundir),'a').close()
      else: open('{}/ERROR'.format(self.rundir),'a').close()
    return 0 if self.rundirOK else 1
    
  def rewriteRestart(self, backup_folder='restart_', lerror=True, nidx=4):
    ''' rewrite grok file for a restart based on existing output files '''
    os.chdir(self.rundir) # move into directory and work with relative path
    # extract end time from time series to find restart time   
    with open(self.newton_file, 'r') as nf:
        lines = nf.readlines()
    last_time = float(lines[-1].split()[0])
    out_times = self.getParam('output times', dtype='float', llist=True)
    restart_index = bisect.bisect(out_times, last_time)
    times_done = out_times[:restart_index]
    if len(times_done)==0:        
        self.lrestart = False
        return # no restart necessary - just start from beginning
    times_todo = out_times[restart_index:]
    if len(times_todo) == 0:
        raise HGSError("Simulations appears to be complete --- no restart necesary/possible.")          
    restart_time = times_done[-1]
    # set new time parameters in Grok file
    self.setParam('output times', times_todo, formatter='{:e}', )
    self.setParam('initial time', restart_time, formatter='{:e}', )
    # find last head files
    filetypes = [self.pm_files,self.olf_files]
    if self.lchannel: filetypes += [self.chan_files]
    indices = [] # last indices out head output files
    for filetype in filetypes: 
      # nidx: number of digits at the end
      name_pattern = filetype.format(IDX=0)[:nidx] # cut off last four digits
      file_inidces = numberedPattern(name_pattern, nidx=nidx, folder=None)
      indices.append(max(file_inidces))
    if min(indices) < max(indices): # should all be the same
        raise IOError("Available head output files (PM,OLF,Chan) are not numbered consistently -- cannot restart!")    
    # now we know that we are actually restarting, so set RESTART indicator
    open('RESTARTED','a').close()
    # assemble name for restart file pattern
    tmp = self.out_files.format(PROBLEM=self.problem, FILETYPE='{FILETYPE}')
    # N.B.: the double-format is necessary, because IDX is enclosed in double-braces (see Grok.__init__)      
    restart_pattern = tmp.format(IDX=indices[0], FILETYPE='{FILETYPE}')
    # determine new restart backup folder to store time-dependent output
    # N.B.: to prevent data loss, a new folder with a 4-digit running number is created
    idx = max(numberedPattern(backup_folder, nidx=nidx, folder=None)) # absolute path
    backup_folder = backup_folder + '{:04d}'.format(idx)
    os.mkdir(backup_folder)
    # backup grok files
    shutil.copy2(self._targetfile, backup_folder)
    # move time-series and binary files
    ts_list = timeseriesFiles(prefix=self.problem, folder=None, ldict=False, llogs=True, lcheck=True)
    binary_list = binaryFiles(prefix=self.problem, folder=None, nidx=nidx, ldict=False)
    for backup_file in ts_list + binary_list:
        shutil.move(backup_file,backup_folder) # relative path
    # add backup folder to restart file pattern (will be checked when restart files are updated)
    restart_pattern = os.path.join(backup_folder,restart_pattern)
    # return name of restart file
    return restart_pattern
    
  def setupConfig(self, template_folder=None, linput=True, lpidx=True, runtime_override=None):
    ''' load config file from template and write configuration to rundir '''
    if template_folder is None:
      template_folder = self.rundir if self.template_folder is None else self.template_folder
    ec = 0 # cumulative exit code
    # load config file from template folder
    ec += self.readConfig(folder=template_folder)
    # apply time setting
    ec += self.setRuntime()
    # write input lists to run folder
    if linput: ec += self.generateInputLists(lvalidate=True)
    # rewrite Grok file for restarts
    if self.lrestart:
      self.ic_files = self.rewriteRestart()      
    # change initial condition files in Grok
    if self.ic_files:
      self.changeICs(ic_pattern=self.ic_files)
    # after input lists are written, apply runtime override (just for testing)
    if runtime_override is not None: 
        ec += self.setRuntime(runtime=runtime_override, output_interval=1)
    # write config file to run folder
    ec += self.writeConfig()
    # write parallelindex with default settings
    if lpidx: ec += self.writeParallelIndex() # can also run just before HGS
    # set config status
    self.configOK = True if ec == 0 else False
    return ec
    
  def runGrok(self, executable=None, logfile='log.grok', lerror=False, lconfig=True, linput=True, ldryrun=False, lcompress=True):
    ''' run the Grok executable in the run directory and set flag indicating success '''
    if lconfig: self.writeConfig()
    ec = super(HGS,self).runGrok(executable=executable, logfile=logfile, lerror=lerror, ldryrun=ldryrun, lcompress=lcompress)
    self.GrokOK = True if ec == 0 else False # set Grok flag
    return ec
  
  def writeParallelIndex(self, NP=None, dom_parts=None, solver=None, input_coloring=False, 
                         run_time=-1., restart=1, parallelindex=None):
    ''' write the parallelindex.dat input file for HGS execution (executed by runHGS) '''
    pwd = os.getcwd() # save present workign directory to return later
    os.chdir(self.rundir) # go into run Grok/HGS folder
    # fix up arguments
    self.pidx_file = parallelindex if parallelindex is not None else self.pidx_file
    if NP is None: NP = self.NP
    if dom_parts is None: dom_parts = NP
    if solver is None: solver = 1 if NP == 1 else 2
    input_coloring = 'T' if input_coloring else 'F'
    # the other two are just numbers (float and int)
    # insert arguments
    contents  = '__Number_of_CPU\n  {:d}\n'.format(NP)
    contents += '__Num_Domain_Partitiong\n  {:d}\n'.format(dom_parts)
    contents += '__Solver_Type\n  {:d}\n'.format(solver)
    contents += '__Coloring_Input\n  {:s}\n'.format(input_coloring)
    contents += '__Wrting_Output_Time\n  {:g}\n'.format(run_time)
    contents += '__Simulation_Restart\n  {:d}\n'.format(restart)
    # write file
    with open(self.pidx_file, 'w+') as pi: pi.writelines(contents)
    if os.path.isfile(self.pidx_file): self.pidxOK = True # make sure file was written
    else: raise IOError(self.pidx_file)
    os.chdir(pwd) # return to previous working directory
    return 0 if self.pidxOK else 1
    
  def runHGS(self, executable=None, logfile='log.hgs_run', lerror=True, lcompress=True,
             skip_config=False, skip_grok=False, skip_pidx=False, ldryrun=False):
    ''' check if all inputs are in place and run the HGS executable in the run directory '''
    pwd = os.getcwd() # save present workign directory to return later    
    os.chdir(self.rundir) # go into run Grok/HGS folder
    self.hgs_bin = executable if executable is not None else self.hgs_bin
    if not os.path.isfile(self.hgs_bin): 
      raise IOError("HGS executable '{}' not found.".format(self.hgs_bin))
    ## check prerequisites and run, if necessary
    # Grok configuration
    if not skip_config and not self.configOK:
      ec = self.setupConfig() # will run with defaults, assuming template is already defined
      if lerror and ec != 0: raise GrokError('Grok configuration did not complete properly.')
    # Grok run
    if not skip_grok and not self.GrokOK: 
      ec = self.runGrok(lconfig=not skip_config, lerror=lerror, ldryrun=ldryrun, lcompress=lcompress) # run grok (will raise exception if failed)
      if lerror and ec != 0: raise GrokError('Grok did not run or complete properly.')
    # parallelindex configuration
    if not skip_pidx and not self.pidxOK:
      self.writeParallelIndex() # write parallel index with default settings
      if not os.path.isfile(self.pidx_file): 
        raise HGSError('Parallel index file was not written properly.')  
    # set indicator file to 'in progress'
    if self.lindicators: 
      shutil.move('{}/SCHEDULED'.format(self.rundir),'{}/IN_PROGRESS'.format(self.rundir))  
    ## run executable while logging output
    with open(logfile, 'w+') as lf: # output and error log
      if ldryrun:
        lf.write('\nDry-run --- no execution\n')
        lec = True # pretend everything works
      else:
        # run HGS as subprocess
        subprocess.call([os.path.abspath(self.hgs_bin)], stdout=lf, stderr=lf)
        # parse log file for errors
        lec = ( tail(lf, n=2)[0].strip() == '---- NORMAL EXIT ----' )
        # i.e. -2, second line from the end (and different capitalization from Grok!)
    # compress binary 3D output fields
    if lcompress and lec:
      with open(logfile, 'a') as lf: # output and error log
        try:  
          # compress, using tar; appending to HGS log
          tar_file = 'binary_fields.tgz'; bin_regex = '*.[0-9][0-9][0-9][0-9]'
          ec = subprocess.call('tar czf {} {}'.format(tar_file,bin_regex), shell=True, stdout=lf, stderr=lf)
          if ec == 0 and os.path.isfile(tar_file):
            lf.write('\nBinary 3D output has been compressed: \'{:s}\'\n'.format(tar_file))
            # if tarring was, remove the binary fields 
            ec = subprocess.call('rm {}'.format(bin_regex), shell=True, stdout=lf, stderr=lf)
            if ec == 0: lf.write('All binary 3D output files (\'{:s}\') have been removed.\n'.format(bin_regex))
            else: lf.write('Cleanup of binary 3D output files (\'{:s}\') failed.\n'.format(bin_regex))
          else:
            lf.write('\nBinary output compression failed; tar exit code: {}\n'.format(ec))
        except:
          lf.write('\nBinary output compression failed for unknown reasons; is \'tar\' availalbe?.\n'.format(ec))
          raise
    os.chdir(pwd) # return to previous working directory
    self.HGSOK = lec # set Grok flag
    # set indicator file to indicate result
    if self.lindicators:
      if lec: shutil.move('{}/IN_PROGRESS'.format(self.rundir),'{}/COMPLETED'.format(self.rundir))  
      else: shutil.move('{}/IN_PROGRESS'.format(self.rundir),'{}/FAILED'.format(self.rundir))
    # after indicators are set, we can raise an error  
    if lerror and not lec: 
      raise HGSError("HGS failed; inspect log-file: {}\n  ('{}')".format(logfile,self.rundir))
    # return a regular (POSIX) exit code
    return 0 if lec else 1
  
  def concatOutput(self):
    ''' a function to concatenate HGS timeseries files after a restart; the following file types are 
        concatenated: hydrograph's, observation_well_flow's, water_balance, and newton_info '''
    # hydrograph's and observation_wll_flow's require globbing expressions to find all relevant files
    # hydrograph's, water_balance, and newton_info essentially follow the same pattern,
    # but observation_well_flow's require somewhat different method
    # in all cases search from end backwards and drop everything after the last/initial time step
    # finally, we also need to renumber 
    raise NotImplementedError