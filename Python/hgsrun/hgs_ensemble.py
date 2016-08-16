'''
Created on Aug 13, 2016

A class that manages an ensemble of HGS simulations from creating the run folder over running Grok
to actually running HGS.

@author: Andre R. Erler, GPL v3
'''

# external imports
import inspect, multiprocessing
# internal imports
from utils.misc import expandArgumentList
from geodata.misc import ArgumentError
from hgs_setup import HGS, HGSError, GrokError


# named exception
class EnsembleError(Exception):
  ''' Exception indicating an Error with the HGS Ensemble '''
  pass

# a function that executes a class/instance method for use in apply_async
def apply_method(member, attr, **kwargs): 
  return member,getattr(member, attr)(**kwargs)

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
    
  def __call__(self, lparallel=False, NP=None, inner_list=None, outer_list=None, **kwargs):
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
      # define work loads: functions and their arguments
      results = [pool.apply_async(apply_method, (member,self.attr), kwargs) for member,kwargs in zip(self.klass.members,kwargs_list)]
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
  
  def getter(self, value):
    ''' get attribute values of all ensemble members and return as list '''
    return [getattr(member, self.attr) for member in self.klass.members]

#   def setter(self, value):
#     ''' set attribute of all ensemble members to the same value (argument expansion not supported) '''
#     for member in self.klass.members: setattr(member, self.attr, value)
# N.B.: setting attributes is done via the __setattr__ method, which does not involve this wrapper...

  def __iter__(self):
    ''' return an iterator over the attribute values of all ensemble members '''
    return iter([getattr(member, self.attr) for member in self.klass.members])
  

## HGS Ensemble manager class
class EnsHGS(object):
  '''
    A class that manages an ensemble of HGS simulations from creating the run folder over running Grok
    to actually running HGS.
  '''
  members = None # list of ensemble members
  rundirs = None # list of HGS run dirs
  hgsargs  = None # list of kwargs used to instantiate ensemble members
  
  def __init__(self, inner_list=None, outer_list=None, **kwargs):
    ''' initialize an ensemble of HGS simulations based on HGS arguments and project descriptors;
        all keyword arguments are automatically expanded based on inner/outer product rules, defined
        using the inner_list/outer_list arguments; the expanded argument lists are used to initialize
        the individual ensemble members; note that a string substitution is applied to 'rundir'
        prior to constructing the HGS instance, i.e. rundir.format(**kwargs) '''
    # expand argument list (plain, nothing special)
    kwargs_list = expandArgumentList(inner_list=inner_list, outer_list=outer_list, **kwargs)
    # loop over ensemble members
    self.members = []; self.rundirs = []; self.hgsargs = [] # ensemble lists
    for kwargs in kwargs_list:
      # expand run folder 
      rundir = kwargs.pop('rundir')
      if not isinstance(rundir,basestring): raise TypeError(rundir)
      rundir = rundir.format(**kwargs) # make all arguments available
      if rundir in self.rundirs:
        raise ArgumentError("Multiple occurence of run directory:\n '{}'".format(rundir))
      self.rundirs.append(rundir)
      # isolate HGS constructor arguments
      hgsargs = inspect.getargspec(HGS.__init__).args
      hgsargs = {arg:kwargs[arg] for arg in hgsargs if arg in kwargs} # returns args, varargs, kwargs, defaults
      self.hgsargs.append(hgsargs)
      # initialize HGS instance      
      hgs = HGS(rundir=rundir, **hgsargs)
      self.members.append(hgs)
    
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
    # instantiate new wrapper with current arguments and return wrapper instance
    return EnsembleWrapper(self,attr)      
        
  def setupExperiments(self, inner_list=None, outer_list=None, lgrok=False, lparallel=True, NP=None, **allargs):
    ''' set up run dirs as execute Grok for each member; check setup and report results '''
    ec = 0 # cumulative exit code (sum of all members)
    # create run folders and copy data
    kwargs = {arg:allargs[arg] for arg in inspect.getargspec(HGS.setupRundir).args if arg in allargs}
    ecs = self.setupRundir(inner_list=None, outer_list=None, lparallel=lparallel, NP=NP, **kwargs)
    if any(ecs) or not all(self.rundirOK): 
      raise GrokError("Run folder setup failed in {0} cases:\n{1}".format(sum(ecs),self.rundirs[ecs]))
    ec += sum(ecs)
    # write configuration
    kwargs = {arg:allargs[arg] for arg in inspect.getargspec(HGS.setupConfig).args if arg in allargs}
    ecs = self.setupConfig(inner_list=None, outer_list=None, lparallel=lparallel, NP=NP, **kwargs)
    if any(ecs) or not all(self.configOK): 
      raise GrokError("Grok configuration failed in {0} cases:\n{1}".format(sum(ecs),self.rundirs[ecs]))
    ec += sum(ecs)
    # run Grok
    if lgrok: 
      kwargs = {arg:allargs[arg] for arg in inspect.getargspec(HGS.runGrok).args if arg in allargs}
      ecs = self.runGrok(inner_list=None, outer_list=None, lparallel=lparallel, NP=NP, **kwargs)
      if any(ecs) or not all(self.GrokOK): 
        raise GrokError("Grok execution failed in {0} cases:\n{1}".format(sum(ecs),self.rundirs[ecs]))
      ec += sum(ecs)
    # return sum of all exit codes
    return ec
    
  def runSimulations(self, inner_list=None, outer_list=None, lsetup=True, lgrok=False, 
                     lparallel=True, NP=None, **allargs):
    ''' execute HGS for each ensemble member and report results; set up run dirs as execute Grok,
        if necessary; note that Grok will be executed in runHGS, just prior to HGS '''
    ec = 0 # cumulative exit code (sum of all members)
    # check and run setup and configuration
    if lsetup:
      ec = self.setupExperiments(inner_list=inner_list, outer_list=outer_list, lgrok=lgrok, lparallel=lparallel, NP=NP)
      if ec > 0 or not all(self.configOK): 
        rundirs = [rundir for rundir,OK in zip(self.rundirs,self.configOK) if not OK]
        raise GrokError("Experiment setup failed in {0} cases:\n{1}".format(ec,rundirs))
    # run HGS
    kwargs = {arg:allargs[arg] for arg in inspect.getargspec(HGS.runHGS).args if arg in allargs}
    ecs = self.runHGS(inner_list=None, outer_list=None, lparallel=lparallel, NP=NP, **kwargs)
    if any(ecs) or not all(self.HGSOK): 
      raise HGSError("Grok configuration failed in {0} cases:\n{1}".format(sum(ecs),self.rundirs[ecs]))
    ec += sum(ecs)
    # return sum of all exit codes
    return ec
    