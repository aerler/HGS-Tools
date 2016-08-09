'''
Created on Jul 31, 2016

A Python module to generate the input file list that contains the atmospheric foring for HGS.

@author: Andre R. Erler, GPL v3
'''

import numpy as np
import os
from geodata.misc import days_per_month, days_per_month_365, seconds_per_month, seconds_per_month_365, abbr_of_month
list_format = '{T:18.3f} {F:s}'


# iterator for monthly intervals
class MonthlyIter(object):
  ''' iterator to return the cumulative at the end of each month '''
  _len = None
  _i = 0
  _sum = 0
  _time_per_month = None
  
  def __init__(self, length, start=0, l365=True, units='seconds'):
    ''' initialize with number of month (length) and start month (start) '''
    # interprete first month
    if isinstance(start,basestring):
      try: start = abbr_of_month.index(start[:3].lower())
      except ValueError: raise ValueError, "Invalid name of month: {}".format(start)
    # figure out units and data convention (l365 means no leap days)
    if units[:6].lower() == 'second':
      self._time_per_month = seconds_per_month_365 if l365 else seconds_per_month
    elif units[:3].lower() == 'day':
      self._time_per_month = days_per_month_365 if l365 else days_per_month
    elif units[:5].lower() == 'month':
      self._time_per_month = np.ones(12)
    else:
      raise NotImplementedError, "Unknown units: '{:s}'".format(units)
    self._i = start # start counting here
    self._len = length + start
    self._sum = 0 # cumulative sum (initialize with 0)
    
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
      return cumsum
        

# function to write input file list
def generateInputFilelist(filename=None, folder=None, input_folder=None, input_pattern=None, 
                          listformat=list_format, lvalidate=True, units='seconds', l365=True,
                          lFortran=True, interval='monthly', length=0, mode='climatology'):
  ''' a function to generate a list of climate data input files for HGS '''
  # construct time and file name lists
  if mode[-5:] == '-mean' or mode in ('mean','steady-state'):
    # determin begin and end times in seconds (begin == 0)
    if interval[:5].lower() == 'month': end_time = 86400. * 365. * length/12. 
    else: end_time = interval * length
    time_iter = iter([0,end_time])
  else:
    # determine mode
    if mode in ('climatology','periodic'): lperiodic = True      
    elif mode in ('time-series','transient'): lperiodic = False
    else: raise NotImplementedError, mode
    # determine length of period (always one year, but different units)
    if lperiodic:
        if units[:6].lower() == 'second': period = 365.*86400. if l365 else 365.2425*86400.
        elif units[:3].lower() == 'day': period = 365. if l365 else 365.2425
        elif units[:5].lower() == 'month': period = 12
        else: raise NotImplementedError, "Unknown units: '{:s}.".format(units)
    # initialize time iterator
    if interval[:5].lower() == 'month':
      idxprd = 12 # one year for monthly input
      time_iter = MonthlyIter(length=length, start=0, l365=l365, units=units)
    elif interval[:3].lower() == 'day': 
      raise NotImplementedError, interval
  # write time/filepath list based on iterators
  listformat = listformat+'\n' # need *two* end of line character
  os.chdir(folder) # move into run folder
  if os.path.exists(filename): os.remove(filename) # remove old file
  openfile = open(filename, 'w') # open file list
  for idx,time in enumerate(time_iter):
    list_time = time # cumulative time in list
    # construct filename
    if lperiodic:
      time %= period; idx %= idxprd
    if lFortran: idx += 1 # Fortran index starts at 1, not 0    
    input_file = input_pattern.format(TIME=time,IDX=idx)
    filepath = '{:s}/{:s}'.format(input_folder,input_file) # assemble current file name
    # check if file actually exists
    if lvalidate and not os.path.exists(filepath): 
      raise IOError, "The input file '{:s}' does not exist.\n(run folder: '{:s}')".format(filepath,folder)
    # write list entry
    openfile.write(listformat.format(T=list_time,F=filepath))
  openfile.close()
    
    
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
    if not os.path.isdir(workdir): raise IOError, workdir
    # path to test file
    testfolder = '{:s}/input_list/'.format(workdir)
    if not os.path.exists(testfolder):os.mkdir(testfolder)
    testfile = 'test.inc'
    testfilepath = testfolder+testfile
    
    ## write test file
    if test_case == 'simple_mean':
      # test simple mean
      generateInputFilelist(filename=testfile, folder=testfolder,
                            input_folder='../test_folder', input_pattern='test_file.asc', 
                            length=180, mode='mean', lvalidate=False)
    
    elif test_case == 'climatology':
      # test simple mean
      generateInputFilelist(filename=testfile, folder=testfolder,
                            input_folder='../test_folder', input_pattern='test_file_{IDX:02d}.asc', 
                            length=24, mode='climatology', lvalidate=False)
    
    elif test_case == 'time-series':
      # test simple mean
      generateInputFilelist(filename=testfile, folder=testfolder,
                            input_folder='../test_folder', input_pattern='test_file_{IDX:02d}.asc', 
                            length=24, mode='time-series', lvalidate=False)

    ## read and print test file
    openfile = open(testfilepath, 'r')
    for line in openfile: print(line)
    openfile.close()
    
      