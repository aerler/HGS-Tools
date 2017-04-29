'''
Created on Jul 31, 2016

A Python module to generate the input file list that contains the atmospheric foring for HGS.

@author: Andre R. Erler, GPL v3
'''

import numpy as np
import os
from geodata.misc import days_per_month, days_per_month_365, seconds_per_month, seconds_per_month_365, abbr_of_month, ArgumentError
list_format = '{T:18.3f} {F:s}'


# iterator for monthly intervals
class MonthlyIter(object):
  ''' iterator to return the the elapsed time at the middle or end of each month '''
  _len = None
  _i = 0
  _sum = 0
  _time_per_month = None
  
  def __init__(self, length, start=0, l365=True, units='seconds', lctr=True):
    ''' initialize with number of month (length) and start month (start) '''
    # interprete first month
    if isinstance(start,basestring):
      try: start = abbr_of_month.index(start[:3].lower())
      except ValueError: raise ValueError("Invalid name of month: {}".format(start))
    # figure out units and data convention (l365 means no leap days)
    if units[:6].lower() == 'second':
      self._time_per_month = seconds_per_month_365 if l365 else seconds_per_month
    elif units[:3].lower() == 'day':
      self._time_per_month = days_per_month_365 if l365 else days_per_month
    elif units[:5].lower() == 'month':
      self._time_per_month = np.ones(12)
    else:
      raise NotImplementedError("Unknown units: '{:s}'".format(units))
    self._i = start # start counting here
    self._len = length + start
    self._sum = 0 # cumulative sum (initialize with 0)
    self._lctr = lctr # return mid-point or boundary values
    
  def __iter__(self):
    ''' make iterator iterable'''
    return self
  
  def next(self):
    ''' return cumulative elapsed time for this month '''
    if self._i > self._len:
      raise StopIteration
    else:
      cumsum = self._sum # return the previous value
      self._sum += self._time_per_month[self._i%12] # cyclical
      self._i += 1 # increase for next month
      if self._lctr: cumsum = (cumsum+self._sum)/2. # mid-point or boundary values
      # N.B.: the values in self._sum are the end points/boundaries, not the centers/mid-points
      return cumsum

# iterator for daily intervals
class DailyIter(object):
  ''' iterator to return the the elapsed time at the middle or end of each day '''
  _len = None
  _i = 0
  _time_per_day = None
  
  def __init__(self, length, start=0, units='seconds', lctr=True, l365=True):
    ''' initialize with number of days (length) and start day (start) '''
    if not l365: raise NotImplementedError('Need to think about how to handle leap years')
    # figure out units and data convention (l365 means no leap days)
    if units[:6].lower() == 'second': 
      self._time_per_day = 86400. if l365 else NotImplemented
    elif units[:3].lower() == 'day':
      self._time_per_day = 1 if l365 else NotImplemented
    else:
      raise NotImplementedError("Unknown units: '{:s}'".format(units))
    self._i = start # start counting here
    self._len = length + start
    self._sum = 0 # cumulative sum (initialize with 0)
    self._lctr = lctr # return mid-point or boundary values
    
  def __iter__(self):
    ''' make iterator iterable'''
    return self
  
  def next(self):
    ''' return cumulative elapsed time for this month '''
    if self._i > self._len:
      raise StopIteration
    else:
      i = self._i + 0.5 if self._lctr else self._i # mid-point or boundary values
      self._i += 1 # increase count of days
      # N.B.: the computer with self._i alone would be the end points/boundaries, not the centers/mid-points
      return self._time_per_day * i

# function to resolve length, end_time, and interval
def resolveInterval(length=None, end_time=None, interval=None):
  ''' determine length or end_time based on interval and the other variable '''
  if length and end_time: raise ArgumentError
  if isinstance(interval,basestring):
    # interval is a string
    if interval[:5].lower() == 'month': 
      if length: # determing end_time from length and interval
        end_time = 86400. * 365. * length/12.
      else: # determine length from end_time and interval 
        length = end_time * 12. / 86400. / 365.  
    elif interval.lower() in ('day','daily'): 
      if length: end_time = 86400. * length # determing end_time from length and interval        
      else: length = end_time / 86400. # determine length from end_time and interval 
    else:
      raise NotImplementedError(interval)
  else: 
    # interval is a number (in 'units')
    if length: end_time = interval * length
    else: length = end_time//interval
  # return consistent length and end_time
  return length, end_time


# function to determine period based on 
def resolvePeriod(period=None, interval=None, units=None, l365=False):
  ''' determine the period length based on input variables '''
  if period.lower() in ('annual','yearly'):
    # period length in right units
    if units[:6].lower() == 'second': period = 365.*86400. if l365 else 365.2425*86400.
    elif units[:3].lower() == 'day': period = 365. if l365 else 365.2425
    elif units[:5].lower() == 'month': period = 12
    else: raise NotImplementedError("Unknown units: '{:s}.".format(units))
    # period length as multiples of interval
    if isinstance(interval,basestring):
      if interval[:5].lower() == 'month': idxprd = 12 # one year for monthly input
      elif interval.lower() in ('day','daily'): idxprd = 365. if l365 else 365.2425
      else: raise NotImplementedError("Unknown interval: '{:s}.".format(interval))
    else: idxprd = period // interval # period should be a number by now
  # return period in correct units and as multiples of interval
  return period, idxprd


# function to write input file list
def generateInputFilelist(filename=None, folder=None, input_folder=None, input_pattern=None, lcenter=True,
                          listformat=list_format, lvalidate=True, units='seconds', l365=True, lFortran=True, 
                          interval='monthly', length=0, end_time=0, mode='climatology', period='yearly'):
  ''' a function to generate a list of climate data input files for HGS '''
  if folder: filename = '{}/{}'.format(folder,filename) 
  # determine end time in seconds (begin == 0) or number of intervals (length)
  length, end_time = resolveInterval(length=length, end_time=end_time, interval=interval)
  # construct time and file name lists
  if mode[-5:] == '-mean' or mode in ('mean','steady','steady-state'):
    time_iter = iter([0,end_time])
    lperiodic = False
  else:
    # determine mode
    if mode in ('clim','peri','climatology','periodic'): lperiodic = True      
    elif mode in ('time-series','timeseries','transient','trans'): lperiodic = False
    else: raise NotImplementedError(mode)
    # determine length of period (always one year, but different units)
    if lperiodic:
      period, idxprd = resolvePeriod(period=period, interval=interval, units=units)      
    # initialize time iterator
    if interval[:5].lower() == 'month':
      time_iter = MonthlyIter(length=length-1 if lcenter else length, start=0, l365=l365, 
                              units=units, lctr=lcenter)
    elif interval.lower() in ('day','daily'): 
      time_iter = DailyIter(length=length-1 if lcenter else length, start=0, l365=l365, 
                              units=units, lctr=lcenter)
    else:
      raise NotImplementedError(interval)
  # write time/filepath list based on iterators
  listformat = listformat+'\n' # need *two* end of line character
  if os.path.exists(filename): os.remove(filename) # remove old file
  with open(filename, 'w') as openfile: # open file list
    list_time = 0 # in case there is no iteration (needed below)
    for idx,time in enumerate(time_iter):
      #print idx,time
      list_time = time; list_idx = idx # cumulative time/index in list
      # construct filename
      if lperiodic:
        time %= period; idx %= idxprd
      if lFortran: idx += 1 # Fortran index starts at 1, not 0    
      input_file = input_pattern.format(TIME=time,IDX=idx)
      if input_folder is None: filepath = input_file
      else: filepath = '{:s}/{:s}'.format(input_folder,input_file) # assemble current file name
      # check if file actually exists
      abspath = filepath if os.path.isabs(filepath) else '{:s}/{:s}'.format(folder,filepath) 
      if lvalidate and not os.path.exists(abspath): 
        raise IOError("The input file '{:s}' does not exist.\n(run folder: '{:s}')".format(filepath,folder))
      # write list entry(s)
      if list_idx == 0 and list_time != 0: # add a first line starting at t=0
        openfile.write(listformat.format(T=0,F=filepath))
      openfile.write(listformat.format(T=list_time,F=filepath))
    if list_time < end_time: # add another line until the end of the simulation
      openfile.write(listformat.format(T=end_time,F=filepath))
  # check results and return file status
  return os.path.isfile(filename)
    
    
# abuse for testing
if __name__ == '__main__':
    
    # test cases
#     test_case = 'simple_mean'
#     test_case = 'climatology'
    test_case = 'time-series'
    
    ## file settings
    # work directory settings ("global" variable)
    # the environment variable RAMDISK contains the path to the RAM disk
    RAM = bool(os.getenv('RAMDISK', '')) # whether or not to use a RAM disk
    # either RAM disk or data directory
    workdir = os.getenv('RAMDISK', '') if RAM else '{:s}/test/'.format(os.getenv('DATA_ROOT', '')) 
    if not os.path.isdir(workdir): raise IOError(workdir)
    # path to test file
    testfolder = '{:s}/input_list/'.format(workdir)
    if not os.path.exists(testfolder):os.mkdir(testfolder)
    testfile = 'test.inc'; testpattern = 'test_file'; inputfolder = '../test_folder'
    
    # for production
    grid = 'brd1'; project = 'ASB'; length = 360
    varname = 'liqwatlfx'; testfile = 'precip.inc'; 
#     varname = 'pet'; testfile = 'pet.inc'      
    inputfolder = '../climate_forcing/na12_maritime/'; testpattern = '{:s}_{:s}_iTime'.format(grid,varname)
    testfolder = 'E:/Data/HGS/{PRJ}/{GRD}/NRCan/'.format(PRJ=project,GRD=grid)
      
    ## write test file
    if test_case == 'simple_mean':
      # test simple mean
      generateInputFilelist(filename=testfile, folder=testfolder,
                            input_folder=inputfolder, input_pattern=testpattern+'.asc', 
                            length=100, interval='daily', mode='mean', lvalidate=False)
    
    elif test_case == 'climatology':
      # test simple mean
      testfolder += '/clim_30/'
      generateInputFilelist(filename=testfile, folder=testfolder,
                            input_folder=inputfolder, input_pattern=testpattern+'_{IDX:02d}.asc', 
                            length=length, mode='climatology', lvalidate=False)
    
    elif test_case == 'time-series':
      # test simple mean
      testfolder += '/timeseries/'
      generateInputFilelist(filename=testfile, folder=testfolder,
                            input_folder=inputfolder, input_pattern=testpattern+'_{IDX:02d}.asc', 
                            length=length, mode='time-series', lvalidate=False)
#                             length=10, interval='daily', mode='time-series', lvalidate=False)


    ## read and print test file
    testfilepath = testfolder+testfile
    openfile = open(testfilepath, 'r')
    for line in openfile: print(line)
    openfile.close()
    print('\nFilepath: \'{:s}\''.format(testfilepath))
    
      