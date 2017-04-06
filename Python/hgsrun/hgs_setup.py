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
from hgsrun.input_list import generateInputFilelist, resolveInterval
from geodata.misc import ArgumentError
from utils.misc import tail

## patch symlink on Windows
# adapted from Stak Overflow (last answer: "Edit 2"):
# https://stackoverflow.com/questions/6260149/os-symlink-support-in-windows
if os.name == "nt":
  def symlink_ms(source, link_name):
    import ctypes
    csl = ctypes.windll.kernel32.CreateSymbolicLinkW
    csl.argtypes = (ctypes.c_wchar_p, ctypes.c_wchar_p, ctypes.c_uint32)
    csl.restype = ctypes.c_ubyte
    flags = 1 if os.path.isdir(source) else 0
    try:
      if csl(link_name, source.replace('/', '\\'), flags) == 0:
          raise ctypes.WinError()
    except:
      pass
  os.symlink = symlink_ms

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
  grok_dbg = 'grok.dbg' # Grok debug output
  grok_file = '{PROBLEM:s}.grok' # the Grok configuration file; default includes problem name
  pm_files = '{PROBLEM:s}o.head_pm.{{IDX:04d}}' # porous media head output (for restart)
  olf_files = '{PROBLEM:s}o.head_olf.{{IDX:04d}}' # over-land flow head output (for restart)
  restart_file = '{PROBLEM:s}o.hen' # restart file name
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
  restarts = None # how many evenly spaced restart files to write
  _lines = None # list of lines in file
  _sourcefile = None # file that the configuration was read from
  _targetfile = None # file that the configuration is written to
  
  def __init__(self, rundir=None, project=None, problem=None, runtime=None, length=None, 
               restarts=10, input_mode=None, input_interval=None, input_vars='PET', 
               input_prefix=None, input_folder='../climate_forcing', lcheckdir=True):
    ''' initialize a Grok configuration object with some settings '''
    if lcheckdir and not os.path.isdir(rundir): raise IOError(rundir)
    # determine end time in seconds (begin == 0) or number of intervals (length)
    length, runtime = resolveInterval(length=length, end_time=runtime, interval=input_interval)
    # assign class variables
    self.rundir = rundir
    self.project = project
    self.problem = project if problem is None else problem
    self.grok_file = self.grok_file.format(PROBLEM=self.problem) # default name
    self.pm_files = self.pm_files.format(PROBLEM=self.problem) # default name
    self.olf_files = self.olf_files.format(PROBLEM=self.problem) # default name
    self.restart_file = self.restart_file.format(PROBLEM=self.problem) # default name
    self.runtime = runtime
    self.length = length
    self.restarts = restarts
    # input configuration
    self.setInputMode(input_mode=input_mode, input_interval=input_interval,
                      input_vars=input_vars, input_prefix=input_prefix,
                      input_folder=input_folder)
  
  def readConfig(self, filename=None, folder=None):
    ''' Read a grok configuration file into memory (or a template to start from). '''    
    filename = filename or self.grok_file; self.grok_file = filename
    filename = '{:s}/{:s}'.format(folder or self.rundir, filename) # prepend folder
    if not os.path.isfile(filename): raise IOError(filename)
    self._sourcefile = filename # use  different file as template
    # read source file
    with open(filename, 'r') as src: # with-environment should take care of closing the file
      self._lines = src.readlines() # read all lines into list
    # strip white spaces and convert to lower case
    self._lines = [line.strip() for line in self._lines]
    # N.B.: converting to lower case creates problems with file/folder path
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
      start = self._lines.index(param.lower(), start) # find begin of vector statement
      end = self._lines.index('end', start) # 'end' marks the end of a vector statement
      # insert list of vector entries into file, line by line
      self._lines = self._lines[:start+1] + values + self._lines[end:]
      # N.B.: the Grok class is not able to detect nested lists/vector-valued parameters; it is only possible to 
      #       manipulate inner-most lists, because the insertion is terminated at the first 'end'
    else:
      # single-valued parameter
      self._lines[self._lines.index(param.lower())+1] = formatter(value) # replace value in list of lines

  def getParam(self, param, dtype=None, llist=None, after=None, start=0):
    ''' read a single parameter, based on the assumption that the parameter value follows in the 
        line below the one where the parameter name appears (case in-sensitive); dtype is a 
        numpy data type to which the value string is cast; a list can be inferred if values of the 
        proper type follow in consecutive lines and are terminated by and 'end' '''
    if after is not None: start = self._lines.index(after, start) # search offset for primary paramerter
    if isinstance(dtype,basestring): dtype = getattr(np,dtype) # convert to given numpy dtype
    elif dtype is None: dtype = str # read as string data type as default
    i = self._lines.index(param.lower(), start)+1
    # different handling for scalars and lists
    if llist is False:
      value = dtype(self._lines[i]) # just a single value
    else:
      lterm = False; values = [] # indicate if list is complete
      while i < len(self._lines) and not lterm:
        value = self._lines[i].strip().lower()
        if value == 'end': 
          lterm = True # proper termination of list
        elif value == '' and llist is None: 
          lterm = True # proper termination for a scalar
        else:
          try: 
            value = dtype(value) # convert string and append
            values.append(value) # append to list
          except ValueError: 
            if llist: raise ValueError("Illegal list termination: '{}'".format(value))
            else: lterm = True # terminate with invalid value
        i += 1 # increment to next line
      # check results
      if value == 'end': value = values # legitimately a list
      elif llist is None and len(values) == 1: value = values[0] # auto-detect scalar
      else: raise ValueError()
    return value

  def replaceParam(self, old, new, formatter=None, after=None, start=0):
    ''' repalce a parameter value with a new value '''
    if formatter: 
      new = formatter.format(new) # apply appropriate formatting
      old = formatter.format(old) # apply appropriate formatting
    else:
      new = str(new); old = str(old)
    if after is not None: start = self._lines.index(after, start)# search offset for primary paramerter
    self._lines[self._lines.index(old.lower(), start)] = new # replace value in list of lines
    
  def editParams(self, **params):
    ''' edit a bunch of parameters, which are defined as key/values pairs using self.editParam;
        note that there is no validation and 'format' is only supported for single edits '''
    for param,value in params.iteritems():
      self.editParam(param, value, format=None)
      
  def setRuntime(self, time, restarts=None):
    ''' set the run time of the simulations (model time in seconds) '''
    if restarts is not None: self.restarts = restarts
    else: restarts = self.restarts
    self.runtime = time
    if self._lines:
      # figure out restart times
      times = [time * float(r) / float(restarts) for r in xrange(1,restarts+1)]
      # change grok file
      self.setParam('output times', times, formatter='{:e}', )
  
  def rewriteRestart(self, backup_folder=None, lnoOut=True):
    ''' rewrite grok file for a restart based on existing output files '''
    if not backup_folder: backup_folder = '{}/backup/'.format(self.rundir)
    if not os.path.isdir(backup_folder): os.mkdir(backup_folder)
    backup_file = '{}/{}'.format(backup_folder,os.path.basename(self._targetfile))
    shutil.copy2(self._targetfile, backup_file)
    # read output times and check presence of files
    out_times = self.getParam('output times', dtype='float', llist=True)
    last_pm = last_olf = last_time = None
    for i,out_time in enumerate(out_times):
        ii = i+1
        pm_file = '{}/{}'.format(self.rundir,self.pm_files.format(IDX=ii))
        olf_file = '{}/{}'.format(self.rundir,self.olf_files.format(IDX=ii))
        if os.path.isfile(pm_file) and os.path.isfile(olf_file):
            last_time = out_time # candidate for last completed time step
            # move output to backup folder
            last_pm = '{}/{}'.format(backup_folder,self.pm_files.format(IDX=ii))
            last_olf = '{}/{}'.format(backup_folder,self.olf_files.format(IDX=ii))
            shutil.move(pm_file, last_pm); shutil.move(olf_file, last_olf)                
        elif os.path.isfile(pm_file):
            raise IOError("Matching OLF file for output time step {} not found!\n('{}')".format(ii,olf_file))
        elif os.path.isfile(olf_file):
            raise IOError("Matching PM file for output time step {} not found!\n('{}')".format(ii,pm_file))
        else: break # terminate loop
    # modify grok file accordingly
    if last_time is None:
        # same restart point as last time... no need to change anything
        if lnoOut: 
            raise IOError("No output files found in run folder!\n('{}')".format(self.rundir))
    else:
      if i == len(out_times)-1: 
          new_times = [] # no more output if this was the last output step
      else:
          assert i > 0, i 
          new_times = out_times[i:] # use the remaining output times
      self.setParam('output times', new_times, formatter='{:e}', )
      self.setParam('initial time', last_time, formatter='{:e}', )
      # merge pm and olf file
      restart_file = '{}/{}'.format(self.rundir, self.restart_file) 
      if os.path.exists(restart_file):     
          shutil.move(restart_file,backup_folder)
      # open new restart file for writing
      with open(restart_file,'wb') as rf:
          # open porous media file and dump into restart file
          with open(last_pm, 'rb') as pm: shutil.copyfileobj(pm, rf)
          # open over land flow file and dump into restart file
          with open(last_olf, 'rb') as olf: shutil.copyfileobj(olf, rf)
    # return name of restart file
    return restart_file
    
  
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
    if input_vars.upper() not in ('PET','WRFPET','NET'): raise NotImplementedError(input_vars)
    # N.B.: the input config NET actually requires changes to the grok file, which is currently not done
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
      filename = '{0}.inc'.format(varname)
      if self.getParam('name', after=vartype, llist=False).lower() != grokname:
            raise GrokError("No entry for boundary condition type '{}'/'{}' found in grok file!".format(vartype,grokname))
      self.setParam('time raster table', 'include {}'.format(filename), after=vartype)
      if self.input_interval == 'monthly':
        if self.input_mode == 'steady-state': input_pattern = 'iTime_{IDX:d}' # IDX will be substituted
        elif self.input_mode == 'periodic': input_pattern = 'iTime_{IDX:02d}' # IDX will be substituted
        elif self.input_mode == 'transient': input_pattern = 'iTime_{IDX:03d}' # IDX will be substituted
        else: raise GrokError(self.input_mode)
        # N.B.: this is a very ugly hack - I don'e have a better idea at the moment, since
        #       we don't know the original length of the time series
      else:
        length = self.length + 1 if lFortran else self.length
        input_pattern = '{0:s}_{{IDX:0{1:d}d}}'.format(axis, int(np.ceil(np.log10(length)))) # number of digits
      input_pattern = '{0}_{1}.asc'.format(wrfvar,input_pattern)
      if input_prefix is not None: input_pattern = '{0}_{1}'.format(input_prefix,input_pattern)
      # write file list
      lec = generateInputFilelist(filename=filename, folder=self.rundir, input_folder=self.input_folder, 
                                  input_pattern=input_pattern, lcenter=lcenter, lvalidate=lvalidate, 
                                  units='seconds', l365=l365, lFortran=lFortran, interval=self.input_interval, 
                                  end_time=self.runtime, mode=self.input_mode)
      ec += 0 if lec else 1          
    # return exit code
    return ec
  
  def runGrok(self, executable=None, logfile='log.grok', batchpfx=None, lerror=True, lcompress=True, ldryrun=False):
    ''' run the Grok executable in the run directory '''
    pwd = os.getcwd() # save present workign directory to return later
    os.chdir(self.rundir) # go into run Grok/HGS folder
    self.grok_bin = executable if executable is not None else self.grok_bin
    self.batchpfx = batchpfx if batchpfx is not None else self.batchpfx
    if not os.path.isfile(self.grok_bin): 
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
        subprocess.call([self.grok_bin], stdout=lf, stderr=lf)
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
  hgs_bin   = './hgs_premium.x' # default HGS executable
  rundirOK  = None # indicate if rundir setup was successfule
  configOK  = None # indicate if Grok configuration was successful
  GrokOK    = None # indicate if Grok ran successfully
  pidxOK    = None # indicate if parallel index configuration was successful
  HGSOK     = None # indicate if HGS ran successfully
  pidx_file = 'parallelindx.dat' # file with parallel execution settings 
  lindicators = True # use indicator files (default: True)
  
  def __init__(self, rundir=None, project=None, problem=None, runtime=None, length=None, restarts=10,
               input_mode=None, input_interval=None, input_vars='PET', input_prefix=None, 
               input_folder='../climate_forcing', template_folder=None, linked_folders=None, NP=1, lindicator=True):
    ''' initialize HGS instance with a few more parameters: number of processors... '''
    # call parent constructor (Grok)
    super(HGS,self).__init__(rundir=rundir, project=project, problem=problem, runtime=runtime, 
                             restarts=restarts, input_vars=input_vars, input_prefix=input_prefix,
                             input_mode=input_mode, input_interval=input_interval,
                             input_folder=input_folder, length=length, lcheckdir=False)
    self.template_folder = template_folder # where to get the templates
    # prepare linked folders
    if linked_folders is None: linked_folders = ('etprop','gb','icbc','prop','soil', # original 
                                                 'grid','init_con','K_maps','landcover','mprops') # extended
    linked_folders = tuple(lf[:-1] if lf[-1] == '/' else lf for lf in linked_folders) # trim slash
    self.linked_folders = linked_folders
    # N.B.: these folders just contain static data and do not need to be replicated
    self.NP = NP # number of processors
    self.lindicators = lindicator # use indicator files
    
  def setupRundir(self, template_folder=None, bin_folder=None, loverwrite=True):
    ''' copy entire run folder from a template folder and link executables '''
    template_folder = self.template_folder if template_folder is None else template_folder
    if template_folder is None: raise ValueError("Need to specify a template path.")
    if not os.path.isdir(template_folder): raise IOError(template_folder)
    # clear existing directory
    if loverwrite and os.path.isdir(self.rundir):
#       ec = subprocess.call(['rm','-r',self.rundir], stdout=subprocess.STDOUT, stderr=subprocess.STDOUT)
      p = subprocess.Popen(['rm','-r',self.rundir], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      stdout, stderr = p.communicate(); ec = p.poll(); del stdout
      if ec > 0: 
        raise IOError("The command rm -r {} to remove the existing directory failed; exit code {}\n{}".format(self.rundir,ec,stderr))
      #shutil.rmtree(self.rundir) 
      # N.B.: rmtree is dangerous, because if follows symbolic links and deletes contents of target directories!
    # copy folder tree
    if not os.path.isdir(self.rundir): shutil.copytree(template_folder, self.rundir, symlinks=True,
                                                       ignore=shutil.ignore_patterns(*self.linked_folders))
    # place link to template
    os.symlink(template_folder, '{}/template'.format(self.rundir))
    # symlinks to static folder that are not copied (linked_folders)
    if self.linked_folders:
      for lf in self.linked_folders:
        link_source = '{}/{}/'.format(template_folder,lf)
        if os.path.exists(link_source):
          #print('linking {}: {}'.format(lf,'{}/{}'.format(self.rundir,lf)))
          os.symlink(link_source, '{}/{}'.format(self.rundir,lf))
    # copy or put links to executables in place
    for exe in (self.hgs_bin, self.grok_bin):
      local_exe = '{}/{}'.format(self.rundir,exe)
      if os.path.lexists(local_exe): os.remove(local_exe)
      if bin_folder is not None:
        os.symlink('{}/{}'.format(bin_folder,exe), local_exe)
      else: 
        os.symlink('{}/{}'.format(template_folder,exe), local_exe)
      # check executables
      if os.path.islink(local_exe):
        if not os.path.exists(local_exe): 
          raise IOError("Link to executable '{}' in run folder is broken.\n ('{}') ".format(exe,self.rundir))
      elif not os.path.isfile(local_exe): 
        raise IOError("Executable file '{}' not found in run folder.\n ('{}') ".format(exe,self.rundir)) 
    # set rundir status
    self.rundirOK = os.path.isdir(self.rundir)
    if self.lindicators: 
      if self.rundirOK: open('{}/SCHEDULED'.format(self.rundir),'a').close()
      else: open('{}/ERROR'.format(self.rundir),'a').close()
    return 0 if self.rundirOK else 1
    
  def setupConfig(self, template_folder=None, linput=True, lpidx=True):
    ''' load config file from template and write configuration to rundir '''
    if template_folder is None:
      template_folder = self.rundir if self.template_folder is None else self.template_folder
    ec = 0 # cumulative exit code
    # load config file from template folder
    ec += self.readConfig(folder=template_folder)
    # apply time setting
    if self.runtime is not None: self.setRuntime(time=self.runtime)
    # write input lists to run folder
    if linput: ec += self.generateInputLists(lvalidate=True)
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
      ec = self.runGrok(lerror=lerror, ldryrun=ldryrun, lcompress=lcompress) # run grok (will raise exception if failed)
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
        subprocess.call([self.hgs_bin], stdout=lf, stderr=lf)
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
  
  