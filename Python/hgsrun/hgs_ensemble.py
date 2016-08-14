'''
Created on Aug 13, 2016

A class that manages an ensemble of HGS simulations from creating the run folder over running Grok
to actually running HGS.

@author: Andre R. Erler, GPL v3
'''

# external imports
import inspect
# internal imports
from utils.misc import expandArgumentList
from hgs_setup import HGS
from geodata.misc import ArgumentError


# named exception
class EnsembleError(Exception):
  ''' Exception indicating an Error with the HGS Ensemble '''
  pass


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
  
  def __getattr__(self, attr):
    ''' execute function call on ensemble members, using the same arguments; list expansion with 
        inner_list/outer_list is also supported'''
    # N.B.: this method returns an on-the-fly wrapper function which expands the argument list and
    #       applies it over the list of methods from all ensemble members
    # N.B.: this method is only called as a fallback, if no class/instance attribute exists,
    #       i.e. Variable methods and attributes will always have precedent 
    # get list of member methods
    methods = [getattr(member,attr) for member in self.members]
    if all([callable(method) for method in methods]):
      ## define wrapper function
      def wrapper(inner_list=None, outer_list=None, **kwargs):
        # expand kwargs to ensemble list
        kwargs_list = expandArgumentList(inner_list=inner_list, outer_list=outer_list, **kwargs)
        if len(kwargs_list) == 1: kwargs_list = kwargs_list * len(self.members)
        elif len(kwargs_list) != len(self.members): 
          raise ArgumentError('Length of expanded argument list does not match ensemble size! {} ~= {}'.format(
                              len(kwargs_list),len(self.members)))
        # loop over ensemble members and execute function
        results = [method(**kwargs) for method,kwargs in zip(methods,kwargs_list)]
        if len(results) != len(self.members): 
          raise ArgumentError('Length of results list does not match ensemble size! {} ~= {}'.format(
                              len(results),len(self.members)))
        return results
    else:
      return methods # in this case, these are just class/instance variables
        
  def setupExperiments(self):
    ''' set up run dirs as execute Grok for each member; check setup and report results '''
    raise NotImplementedError
  
  def runExperiments(self, NP):
    ''' run multiple experiments in parallel using multi-processing and report results '''
    raise NotImplementedError