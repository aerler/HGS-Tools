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
# internal imports
from input_list import generateInputFilelist, resolveInterval
from geodata.misc import ArgumentError
from utils.misc import tail

# some named exceptions
class GrokError(Exception):
  ''' Exception indicating an Error in Grok '''
  pass
class HGSError(Exception):
  ''' Exception indicating an Error in HGS '''
  pass

## a class to handle Grok for HGS simulations
class Grok(object):
  '''
    A class that loads a grok configuration file into memory, provides functions for editing,
    saving the file, and running Grok.
  '''
  grok_bin = './grok_premium.x' # default Grok executable
  batchpfx = 'batch.pfx' # a file that defines the HGS problem name for Grok
  rundir = None # folder where the experiment is set up and executed
  project = None # a project designator used for file names
  problem = None # the HGS problem name (defaults to project)
  input_mode = None # type of simulation (forcing data): stead-state, periodic, transient
  input_interval = None # update interval for climate forcing
  input_vars = 'PET' # input variable configuration
  input_prefix = None # prefix for input files
  input_folder = '../climate_forcing' # default folder for input data
  runtime = None # run time for the simulations (in seconds)
  length = None # run time in multiples of the interval length
  _lines = None # list of lines in file
  _sourcefile = None # file that the configuration was read from
  _targetfile = None # file that the configuration is written to
  
  def __init__(self, rundir=None, project=None, problem=None, runtime=None, length=None, 
               input_mode=None, input_interval=None, input_vars='PET', input_prefix=None, 
               input_folder='../climate_forcing', lcheckdir=True):
    ''' initialize a Grok configuration object with some settings '''
    if lcheckdir and not os.path.isdir(rundir): raise IOError(rundir)
    # determine end time in seconds (begin == 0) or number of intervals (length)
    length, runtime = resolveInterval(length=length, end_time=runtime, interval=input_interval)
    # assign class variables
    self.rundir = rundir
    self.project = project
    self.problem = project if problem is None else problem
    self.runtime = runtime
    self.length = length
    # input configuration
    if input_mode or input_interval:
      self.setInputMode(input_mode=input_mode, input_interval=input_interval,
                        input_vars=input_vars, input_prefix=input_prefix,
                        input_folder=input_folder)
  
  def readConfig(self, filename=None, folder=None):
    ''' Read a grok configuration file into memory (or a template to start from). '''    
    filename = filename or '{:s}.grok'.format(self.project) # or default name
    filename = '{:s}/{:s}'.format(folder or self.rundir, filename) # prepend folder
    if not os.path.isfile(filename): raise IOError(filename)
    self._sourcefile = filename # use  different file as template
    # read source file
    with open(filename, 'r') as src: # with-environment should take care of closing the file
      self._lines = src.readlines() # read all lines into list
    # strip white spaces and convert to lower case
    self._lines = [line.strip() for line in self._lines]
    # N.B.: converting to lower case creates problems with file/folder path
    # apply time setting
    if self.runtime is not None: self.setRuntime(time=self.runtime)
    # return exit code
    return 0 if isinstance(self._lines,list) else 0
      
  def writeConfig(self, filename=None):
    ''' Write the grok configuration to a file in run dir. '''    
    filename = filename or '{:s}.grok'.format(self.project) # or default name    
    filename = '{:s}/{:s}'.format(self.rundir,filename) # prepend run dir
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
        regular format string used to write the value '''
    if formatter: value = formatter.format(value) # apply appropriate formatting
    else: value = str(value)
    start = self._lines.index(after, start) if after else start # search offset for primary paramerter
    self._lines[self._lines.index(param.lower())+1] = value # replace value in list of lines

  def getParam(self, param, dtype=None, after=None, start=0):
    ''' read a single parameter, based on the assumption that the parameter value follows in the 
        line below the one where the parameter name appears (case in-sensitive); dtype is a 
        numpy data type to which the value string is cast '''
    start = self._lines.index(after, start) if after else start # search offset for primary paramerter
    value = self._lines[self._lines.index(param.lower(), start)+1]
    if isinstance(dtype,basestring): value = getattr(np,dtype)(value) # convert to given numpy dtype
    elif dtype: value = dtype(value) # try conversion this way...
    return value

  def replaceParam(self, old, new, formatter=None, after=None, start=0):
    ''' repalce a parameter value with a new value '''
    if formatter: 
      new = formatter.format(new) # apply appropriate formatting
      old = formatter.format(old) # apply appropriate formatting
    else:
      new = str(new); old = str(old)
    start = self._lines.index(after, start) if after else start # search offset for primary paramerter
    self._lines[self._lines.index(old.lower(), start)] = new # replace value in list of lines
    
  def editParams(self, **params):
    ''' edit a bunch of parameters, which are defined as key/values pairs using self.editParam;
        note that there is no validation and 'format' is only supported for single edits '''
    for param,value in params.iteritems():
      self.editParam(param, value, format=None)
      
  def setRuntime(self, time):
    ''' set the run time of the simulations (model time in seconds) '''
    self.runtime = time
    if self._lines:
      self.setParam('output times', time, formatter='{:.3e}', )
  
  def setInputMode(self, input_mode=None, input_interval=None, input_vars='PET', input_prefix=None, 
                   input_folder='../climate_forcing'):
    ''' set the type of the simulation: mean/steady-state, climatology/periodic, time-series/transient '''
    # resolve type of input data 
    input_mode = input_mode.lower()
    if input_mode in ('mean','steady','steady-state') or input_mode[-5:] == '-mean': input_mode = 'steady-state'
    elif input_mode[:4] in ('clim','peri','climatology','periodic'): input_mode = 'periodic'
    elif input_mode in ('time-series','timeseries','trans','transient'): input_mode = 'transient'
    self.input_mode = input_mode
    # set input interval
    input_interval = input_interval.lower()
    if input_interval[:5].lower() == 'month': input_interval = 'monthly'
    elif input_interval[:3].lower() == 'day': input_interval = 'daily'
    else: raise NotImplementedError(input_interval)
    self.input_interval = input_interval
    # set other variables
    self.input_vars = input_vars
    self.input_prefix = input_prefix
    self.input_folder = input_folder
  
  def generateInputLists(self, input_vars=None, input_prefix=None, input_folder=None,
                         lvalidate=True, axis='iTime', lcenter=True, l365=True, lFortran=True):
    ''' generate and validate lists of input files and write to appropriate files '''
    input_vars = self.input_vars if input_vars is None else input_vars
    input_prefix = self.input_prefix if input_prefix is None else input_prefix
    input_folder = self.input_folder if input_folder is None else input_folder
    # generate default input vars
    if isinstance(input_vars,basestring):
      if input_vars.upper() == 'PET': # liquid water + snowmelt & PET as input
        input_vars = dict(precip=('rain','liqwatflx'),  
                          pet=('potential evapotranspiration','pet'))
      elif input_vars.upper() == 'NET': # liquid water + snowmelt - ET as input
        input_vars = dict(precip=('rain','waterflx'),)
      else: raise ArgumentError("Invalid or missing input_vars: {}".format(input_vars))
    elif not isinstance(input_vars, dict): raise TypeError(input_vars)
    # iterate over variables and generate corresponding input lists
    ec = 0 # cumulative exit code
    for varname,val in input_vars.iteritems():
      vartype,wrfvar = val
      filename = '{0}.inc'.format(varname)
      self.setParam('time raster table', 'include {}'.format(filename), after=vartype)
      length = self.length + 1 if lFortran else self.length
      input_pattern = '{0:s}_{{IDX:0{1:d}d}}'.format(axis, int(np.ceil(np.log10(length)))) # number of digits
      input_pattern = '{0}_{1}.asc'.format(wrfvar,input_pattern)
      if input_prefix is not None: input_pattern = '{0}_{1}'.format(input_prefix,input_pattern)
      # write file list
      generateInputFilelist(filename=filename, folder=self.rundir, input_folder=self.input_folder, 
                            input_pattern=input_pattern, lcenter=lcenter, 
                            lvalidate=lvalidate, units='seconds', l365=l365, lFortran=lFortran, 
                            interval=self.input_interval, end_time=self.runtime, mode=self.input_mode)
      ec += 0 if os.path.isfile(filename) else 1    
    # return exit code
    return ec
  
  def runGrok(self, executable=None, logfile='log.grok', batchpfx=None, lerror=True):
    ''' run the Grok executable in the run directory '''
    pwd = os.getcwd() # save present workign directory to return later
    os.chdir(self.rundir) # go into run Grok/HGS folder
    self.grok_bin = executable if executable is not None else self.grok_bin
    self.batchpfx = batchpfx if batchpfx is not None else self.batchpfx
    if not os.path.isfile(self.grok_bin): 
      raise IOError("Grok executable '{}' not found.".format(self.grok_bin))
    # create batch.pfx file with problem name for batch processing
    with open(self.batchpfx, 'w+') as bp: bp.write(self.problem) # just one line...
    # run executable while logging output
    with open(logfile, 'w+') as lf: # output and error log
      # run Grok
      subprocess.call([self.grok_bin], stdout=lf, stderr=lf)
      # parse log file for errors
      lec = ( tail(lf, n=3)[0].strip() == '---- Normal exit ----' )
      # i.e. -3, third line from the end (different from HGS)
    os.chdir(pwd) # return to previous working directory
    if lerror and not lec: 
      raise GrokError("Grok failed; inspect log-file: {}\n  ('{}')".format(logfile,self.rundir))
    return 0 if lec else 1
            
      
class HGS(Grok):
  '''
    A child class of Grok that can also set up the entire run folder and launch an HGS instance;
    otherwise the same as Grok.
  '''
  template_folder = None # temlate for rundir setup
  hgs_bin   = './hgs_premium.x' # default HGS executable
  rundirOK  = None # indicate if rundir setup was successfule
  configOK  = None # indicate if Grok configuration was successful
  grokOK    = None # indicate if Grok ran successfully
  pidxOK    = None # indicate if parallel index configuration was successful
  pidx_file = 'parallelindx.dat' # file with parallel execution settings 
  
  def __init__(self, rundir=None, project=None, problem=None, runtime=None, length=None, 
               input_mode=None, input_interval=None, input_vars='PET', input_prefix=None, 
               input_folder='../climate_forcing', template_folder=None, NP=1):
    ''' initialize HGS instance with a few more parameters: number of processors... '''
    # call parent constructor (Grok)
    super(HGS,self).__init__(rundir=rundir, project=project, problem=problem, runtime=runtime, 
                             input_mode=input_mode, input_interval=input_interval,
                             input_vars=input_vars, input_prefix=input_prefix, 
                             input_folder=input_folder, length=length, lcheckdir=False)
    self.template_folder = template_folder # where to get the templates
    self.NP = NP # number of processors
    
  def setupRundir(self, template=None, bin_folder=None, loverwrite=True):
    ''' copy entire run folder from a template and link executables '''
    template = self.template_folder if template is None else template
    if template is None: raise ValueError("Need to specify a template path.")
    if not os.path.isdir(template): raise IOError(template)
    # clear existing directory
    if loverwrite and os.path.isdir(self.rundir): shutil.rmtree(self.rundir)
    # copy folder tree
    if not os.path.isdir(self.rundir): shutil.copytree(template, self.rundir, symlinks=True)
    # put links to executables in place
    for exe in (self.hgs_bin, self.grok_bin):
      local_exe = '{}/{}'.format(self.rundir,exe)
      if bin_folder is not None:
        os.symlink('{}/{}'.format(bin_folder,exe), local_exe)
      # check executables
      if os.path.islink(local_exe):
        if not os.path.exists(local_exe): 
          raise IOError("Link to executable '{}' in run folder is broken.\n ('{}') ".format(exe,self.rundir))
      elif not os.path.isfile(local_exe): 
        raise IOError("Executable file '{}' not found in run folder.\n ('{}') ".format(exe,self.rundir)) 
    # set rundir status
    self.rundirOK = True if os.path.isdir(self.rundir) else False
    return 0 if self.rundirOK else 1
    
  def setupConfig(self, template=None, linput=True, lpidx=True):
    ''' load config file from template and write configuration to rundir '''
    template = self.rundir if template is None else template
    ec = 0 # cumulative exit code
    # load config file from template folder
    ec += self.readConfig(folder=template)
    # N.B.: runtime is already set by readConfig
    # write input lists to run folder
    if linput: ec += self.generateInputLists(lvalidate=True)
    # write config file to run folder
    ec += self.writeConfig()
    # write parallelindex with default settings
    if lpidx: ec += self.writeParallelIndex() # can also run just before HGS
    # set config status
    self.configOK = True if ec == 0 else False
    return ec
    
  def runGrok(self, executable=None, logfile='log.grok', lerror=False, lconfig=True, linput=True):
    ''' run the Grok executable in the run directory and set flag indicating success '''
    if lconfig: self.writeConfig()
    ec = super(HGS,self).runGrok(executable=executable, logfile=logfile, lerror=lerror)
    self.grokOK = True if ec == 0 else False # set Grok flag
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
    
  def runHGS(self, executable=None, logfile='log.hgs', lerror=True,
             skip_config=False, skip_grok=False, skip_pidx=False):
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
    if not skip_grok and not self.grokOK: 
      ec = self.runGrok(lerror=lerror) # run grok (will raise exception if failed)
      if lerror and ec != 0: raise GrokError('Grok did not run or complete properly.')
    # parallelindex configuration
    if not skip_pidx and not self.pidxOK:
      self.writeParallelIndex() # write parallel index with default settings
      if not os.path.isfile(self.pidx_file): 
        raise HGSError('Parallel index file was not written properly.')    
    ## run executable while logging output
    with open(logfile, 'w+') as lf: # output and error log
      # run HGS as subprocess
      subprocess.call([self.hgs_bin], stdout=lf, stderr=lf)
      # parse log file for errors
      ec = ( tail(lf, n=2)[0].strip() == '---- Normal exit ----' )
      # i.e. -2, second line from the end (different from Grok)
    os.chdir(pwd) # return to previous working directory
    if lerror and not ec: 
      raise HGSError("HGS failed; inspect log-file: {}\n  ('{}')".format(logfile,self.rundir))
    return 0 if ec else 1
  
  