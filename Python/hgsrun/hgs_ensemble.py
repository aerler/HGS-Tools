'''
Created on Aug 13, 2016

A class that manages an ensemble of HGS simulations from creating the run folder over running Grok
to actually running HGS.

@author: Andre R. Erler, GPL v3
'''

# external imports
import os, inspect, multiprocessing
# internal imports
from utils.misc import expandArgumentList
from geodata.misc import ArgumentError
from hgsrun.hgs_setup import HGS, HGSError, GrokError


# named exception
class EnsembleError(Exception):
  ''' Exception indicating an Error with the HGS Ensemble '''
  pass

# a function that executes a class/instance method for use in apply_async
def apply_method(member, attr, **kwargs): 
  ''' execute the method 'attr' of instance 'member' with keyword arguments 'kwargs';
      return a tuple containing member instance (possibly changed) and method results '''
  return member,getattr(member, attr)(**kwargs) # returns a TUPLE!!!

# callback function to print reports of completed simulations
def reportBack(result):
  ''' function that prints the results of a simulations from a multiprocessing batch;
      N.B.: the callback function is passed a result from apply_method, which is a tuple '''
  member, ec = result
  if ec == 0: # simulation completed successfully
    print("The simulation in folder '{:s}' completed successfully!".format(member.rundir))
  else: # simulation failed
    print("FAILURE: The simulation in folder '{:s}' terminated with exit code {:d}!".format(member.rundir,ec)) 

## define ensemble wrapper class
class EnsembleWrapper(object):
  ''' 
    A class that applies an attribute or method call on the ensemble class to all of its members
    and returns a list of the results. The ensemble class is assumed to have an attribute 
    'members', which is a list of the ensemble member.
    The EnsembleWrapper class is instantiated with the ensemble class and the called attibure or
    method in the __getattr__ method and is returned instead of the class attribute.
  '''

  def __init__(self, klass, attr):
    ''' the object has to be initialized with the ensemlbe class 'klass' and the called attribute
        'attr' in the __getattr__ method '''
    self.klass = klass # the object that the attribute is called on
    self.attr = attr # the attribute name that is called
    
  def __call__(self, lparallel=False, NP=None, inner_list=None, outer_list=None, callback=None, **kwargs):
    ''' this method is called instead of a class or instance method; it applies the arguments 
        'kwargs' to each ensemble member; it also supports argument expansion with inner and 
        outer product (prior to application to ensemble) and parallelization using multiprocessing '''
    # expand kwargs to ensemble list
    kwargs_list = expandArgumentList(inner_list=inner_list, outer_list=outer_list, **kwargs)
    if len(kwargs_list) == 1: kwargs_list = kwargs_list * len(self.klass.members)
    elif len(kwargs_list) != len(self.klass.members): 
      raise ArgumentError('Length of expanded argument list does not match ensemble size! {} ~= {}'.format(
                          len(kwargs_list),len(self.klass.members)))
    # loop over ensemble members and execute function
    if lparallel:
      # parallelize method execution using multiprocessing
      pool = multiprocessing.Pool(processes=NP) # initialize worker pool
      if callback is not None and not callable(callback): raise TypeError(callback)
      # N.B.: the callback function is passed a result from the apply_method function, 
      #       which returns a tuple of the form (member, exit_code)
      # define work loads (function and its arguments) and start tasks      
      results = [pool.apply_async(apply_method, (member,self.attr), kwargs, callback=callback) 
                                      for member,kwargs in zip(self.klass.members,kwargs_list)]          
      # N.B.: Beware Pickling!!!
      pool.close(); pool.join() # wait to finish
      # retrieve and assemble results 
      results = [result.get() for result in results]
      # divide members and results (apply_method returns both, in case members were modified)
      self.klass.members = [result[0] for result in results]
      results = [result[1] for result in results]
    else:
      # get instance methods
      methods = [getattr(member,self.attr) for member in self.klass.members]
      # just apply sequentially
      results = [method(**kwargs) for method,kwargs in zip(methods,kwargs_list)]
    if len(results) != len(self.klass.members): 
      raise ArgumentError('Length of results list does not match ensemble size! {} ~= {}'.format(
                          len(results),len(self.klass.members)))
    return tuple(results)
  
#   def __get__(self, ):
#     ''' get attribute values of all ensemble members and return as list '''
#     print('\nGETTER\n')
#     return [getattr(member, self.attr) for member in self.klass.members]

#   def setter(self, value):
#     ''' set attribute of all ensemble members to the same value (argument expansion not supported) '''
#     for member in self.klass.members: setattr(member, self.attr, value)
# N.B.: setting attributes is done via the __setattr__ method, which does not involve this wrapper...

#   def __iter__(self):
#     ''' return an iterator over the attribute values of all ensemble members '''
#     return iter([getattr(member, self.attr) for member in self.klass.members])
  

## HGS Ensemble manager class
class EnsHGS(object):
  '''
    A class that manages an ensemble of HGS simulations from creating the run folder over running Grok
    to actually running HGS.
  '''
  members = None # list of ensemble members
  rundirs = None # list of HGS run dirs
  hgsargs = None # list of kwargs used to instantiate ensemble members
  lreport = True # print short simulation summaries and other info (using callback)
  lindicator = True # use indicator files
  loverwrite = False # overwrite existing folders
  lrunfailed = False # rerun failed experiments
  
  def __init__(self, inner_list=None, outer_list=None, **kwargs):
    ''' initialize an ensemble of HGS simulations based on HGS arguments and project descriptors;
        all keyword arguments are automatically expanded based on inner/outer product rules, defined
        using the inner_list/outer_list arguments; the expanded argument lists are used to initialize
        the individual ensemble members; note that a string substitution is applied to all folder 
        variables (incl. 'rundir') prior to constructing the HGS instance, i.e. rundir.format(**kwargs) '''
    self.lreport    = kwargs.get('lreport',self.lreport)
    self.loverwrite = kwargs.get('loverwrite',self.loverwrite)
    self.lindicator = kwargs.get('lindicator',self.lindicator)
    self.lrunfailed = kwargs.get('lrunfailed',self.lrunfailed)
    # expand argument list (plain, nothing special)
    kwargs_list = expandArgumentList(inner_list=inner_list, outer_list=outer_list, **kwargs)
    # loop over ensemble members
    self.members = []; self.rundirs = []; self.hgsargs = [] # ensemble lists
    for kwargs in kwargs_list:
      # isolate folder variables and perform variable substitution
      for folder_type in ('template_folder','input_folder','pet_folder','rundir'):
        if folder_type in kwargs:
          folder = kwargs[folder_type]
          if not isinstance(folder,basestring): raise TypeError(folder)
          # perform keyword substitution with all available arguments
          kwargs[folder_type] = folder.format(**kwargs)
      # check rundir
      rundir = kwargs['rundir']
      if rundir in self.rundirs:
        raise ArgumentError("Multiple occurence of run directory:\n '{}'".format(rundir))
      # figure out skipping      
      if os.path.exists(rundir):
        if self.loverwrite:
          if self.lreport: print("Overwriting existing experiment folder '{:s}'.".format(rundir))
          lskip = False
        elif self.lindicator and os.path.exists('{}/SCHEDULED'.format(rundir)):
          if self.lreport: print("Skipping experiment folder '{:s}' (scheduled).".format(rundir))
          lskip = True
        elif self.lindicator and os.path.exists('{}/IN_PROGRESS'.format(rundir)):
          if self.lreport: print("Skipping experiment folder '{:s}' (in progress).".format(rundir))
          lskip = True
        elif self.lindicator and os.path.exists('{}/COMPLETED'.format(rundir)):
          if self.lreport: print("Skipping experiment folder '{:s}' (completed).".format(rundir))
          lskip = True
        elif self.lindicator and os.path.exists('{}/FAILED'.format(rundir)):
          # this should be the last option, so as to prevent overwriting data
          if self.lrunfailed:            
            if self.lreport: print("Overwriting failed experiment folder '{:s}'.".format(rundir))
            lskip = False # rundir will be deleted
          else: 
            if self.lreport: print("Skipping experiment folder '{:s}' (failed).".format(rundir))
            lskip = True
        else: # no/unknown indicator file
          if self.lreport: print("Overwriting existing experiment folder '{:s}'.".format(rundir))
          lskip = False # rundir will be deleted
      else:
        if self.lreport: print("Creating new experiment folder '{:s}'.".format(rundir))
        lskip = False
      if not lskip:
        self.rundirs.append(rundir)
        # isolate HGS constructor arguments
        hgsargs = inspect.getargspec(HGS.__init__).args # returns args, varargs, kwargs, defaults
        hgsargs = {arg:kwargs[arg] for arg in hgsargs if arg in kwargs} 
        self.hgsargs.append(hgsargs)
        # initialize HGS instance      
        hgs = HGS(**hgsargs)
        self.members.append(hgs)
    # final check
    if len(self.members) == 0: 
      raise EnsembleError("No experiments to run (empty list).")
    
  @property
  def size(self):
    assert len(self.members) == len(self.rundirs) == len(self.hgsargs) 
    return len(self.members) 
  
  def __len__(self): 
    assert len(self.members) == len(self.rundirs) == len(self.hgsargs)
    return len(self.members)
  
  def __iter__(self):
    return self.members.__iter__()
  
  def __setattr__(self, attr, value):
    ''' redirect setting of attributes to ensemble members if the ensemble class does not have it '''
    if attr in self.__dict__ or attr in EnsHGS.__dict__: 
      super(EnsHGS, self).__setattr__(attr, value)
    else:
      for member in self.members: setattr(member, attr, value)
  
  def __getattr__(self, attr):
    ''' execute function call on ensemble members, using the same arguments; list expansion with 
        inner_list/outer_list is also supported'''
    # N.B.: this method returns an on-the-fly EnsembleWrapper instance which expands the argument       
    #       list and applies it over the list of methods from all ensemble members
    # N.B.: this method is only called as a fallback, if no class/instance attribute exists,
    #       i.e. Variable methods and attributes will always have precedent 
    # determine attribute type
    attrs = [callable(getattr(member, attr)) for member in self.members]
    if not any(attrs):
      # treat as regular attributes and return list of attibutes of all members
      return [getattr(member, attr) for member in self.members]
    elif all(attrs):
      # instantiate new wrapper with current arguments and return wrapper instance
      return EnsembleWrapper(self,attr)
    else: raise EnsembleError("Inconsistent attribute type '{}'".format(attr))
        
  def setupExperiments(self, inner_list=None, outer_list=None, lgrok=False, lparallel=True, NP=None, 
                       runtime_override=None, **allargs):
    ''' set up run dirs as execute Grok for each member; check setup and report results '''
    ec = 0 # cumulative exit code (sum of all members)
    # create run folders and copy data
    kwargs = {arg:allargs[arg] for arg in inspect.getargspec(HGS.setupRundir).args if arg in allargs}
    ecs = self.setupRundir(inner_list=inner_list, outer_list=outer_list, lparallel=lparallel, NP=NP, **kwargs)
    if any(ecs) or not all(self.rundirOK): 
      raise GrokError("Run folder setup failed in {0} cases:\n{1}".format(sum(ecs),self.rundirs[ecs]))
    ec += sum(ecs)
    # write configuration
    kwargs = {arg:allargs[arg] for arg in inspect.getargspec(HGS.setupConfig).args if arg in allargs}
    ecs = self.setupConfig(inner_list=inner_list, outer_list=outer_list, lparallel=lparallel, NP=NP, 
                           runtime_override=runtime_override, **kwargs)
    if any(ecs) or not all(self.configOK): 
      raise GrokError("Grok configuration failed in {0} cases:\n{1}".format(sum(ecs),
                                                                            [rd for rd,e in zip(self.rundirs,ecs) if e > 0]))
    ec += sum(ecs)
    # run Grok
    if lgrok: 
      kwargs = {arg:allargs[arg] for arg in inspect.getargspec(HGS.runGrok).args if arg in allargs}
      ecs = self.runGrok(inner_list=inner_list, outer_list=outer_list, lparallel=lparallel, NP=NP, **kwargs)
      if any(ecs) or not all(self.GrokOK): 
        raise GrokError("Grok execution failed in {0} cases:\n{1}".format(sum(ecs),self.rundirs[ecs]))
      ec += sum(ecs)
    # return sum of all exit codes
    return ec
    
  def runSimulations(self, inner_list=None, outer_list=None, lsetup=True, lgrok=False, 
                     lparallel=True, NP=None, runtime_override=None, callback=reportBack, **allargs):
    ''' execute HGS for each ensemble member and report results; set up run dirs as execute Grok,
        if necessary; note that Grok will be executed in runHGS, just prior to HGS '''
    if not self.lreport: callback = None # suppress output
    ec = 0 # cumulative exit code (sum of all members)
    # check and run setup and configuration
    if lsetup:
      arglist = set(inspect.getargspec(HGS.setupRundir).args)
      arglist.union(inspect.getargspec(HGS.setupConfig).args)
      if lgrok: arglist.union(inspect.getargspec(HGS.runGrok).args)
      kwargs = {arg:allargs[arg] for arg in arglist if arg in allargs}
      ec = self.setupExperiments(inner_list=inner_list, outer_list=outer_list, lgrok=lgrok, lparallel=lparallel, NP=NP, 
                                 runtime_override=runtime_override, **kwargs)
      if ec > 0 or not all(self.configOK): 
        rundirs = [rundir for rundir,OK in zip(self.rundirs,self.configOK) if not OK]
        raise GrokError("Experiment setup failed in {0} cases:\n{1}".format(ec,rundirs))
    # run HGS
    kwargs = {arg:allargs[arg] for arg in inspect.getargspec(HGS.runHGS).args if arg in allargs}
    ecs = self.runHGS(inner_list=inner_list, outer_list=outer_list, lparallel=lparallel, NP=NP, callback=callback, 
                      skip_config=True, **kwargs) # setup already ran (or was skipped intentionally)
    if any(ecs) or not all(self.HGSOK): 
      raise HGSError("Grok configuration failed in {0} cases:\n{1}".format(sum(ecs),self.rundirs[ecs]))
    ec += sum(ecs)
    # return sum of all exit codes
    return ec
    