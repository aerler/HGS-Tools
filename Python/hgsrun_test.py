'''
Created on Aug 11, 2016

Unittests for hgsrun components.

@author: Andre R. Erler, GPL v3
'''

import unittest
import numpy as np
import os, sys, gc, shutil


# import modules to be tested
from setup import Grok

# work directory settings ("global" variable)
data_root = os.getenv('DATA_ROOT', '')
# test folder either RAM disk or data directory
RAM = bool(os.getenv('RAMDISK', '')) # whether or not to use a RAM disk
# N.B.: the environment variable RAMDISK contains the path to the RAM disk
workdir = os.getenv('RAMDISK', '') if RAM else '{:s}/test/'.format(data_root)
if not os.path.isdir(workdir): raise IOError(workdir)
# other settings
NP = 2
ldebug = False


## tests for multiprocess module
class GrokTest(unittest.TestCase):  
  # some HGS test data
  hgs_template = data_root+'/HGS/Templates/GRW-test/' 
  hgs_testcase = 'grw_omafra' # name of test project (for file names)
   
  def setUp(self):
    ''' create two test variables '''
    if not os.path.isdir(self.hgs_template): 
      raise IOError("HGS Template for testing not found:\n '{}'".format(self.hgs_template))
    # test folder
    self.rundir = '{}/grok_test/'.format(workdir,)
    if os.path.isdir(self.rundir): shutil.rmtree(self.rundir)
    os.mkdir(self.rundir)
    # grok test files
    self.grok_input  = '{}/{}.grok'.format(self.hgs_template,self.hgs_testcase)
    self.grok_output = '{}/{}.grok'.format(self.rundir,self.hgs_testcase)
    # create Grok instance
    self.grok = Grok(rundir=self.rundir, project=self.hgs_testcase)
    # load a config file from template
    if not os.path.isfile(self.grok_input):
      raise IOError("Grok configuration file for testing not found:\n '{}'".format(self.grok_input))
    self.grok.read(folder=self.hgs_template)
    assert isinstance(self.grok._lines, list), self.grok._lines
      
  def tearDown(self):
    ''' clean up '''
    self.grok.write()
    del self.grok
    gc.collect()
    #shutil.rmtree(self.rundir)
 
  def testClass(self):
    ''' test instantiation of class '''    
    # instantiation done in self.setUp()
    assert self.grok.rundir, self.grok.rundir
    assert os.path.isdir(self.grok.rundir), self.grok.rundir
    
  def testSetTime(self):
    ''' test setting the run time variable in Grok config file '''
    grok = self.grok
    time = 24*60*60 # in seconds, i.e. one day
    # set run time
    grok.setRuntime(time)
    # test
    assert grok.runtime == time
    # read config file into in-memory file and verify time
    output = ''.join(grok._lines) # don't need newlines 
    #print(output)
    assert '{:.3e}'.format(time) in output, '{:.3e}'.format(time)
    
    
    
  def testWrite(self):
    ''' load config file from template and write to rundir (on disk) '''
    grok = self.grok
    grok._lines = None # delete already loaded file contents
    # read from template
    assert os.path.isfile(self.grok_input), self.grok_input 
    grok.read(folder=self.hgs_template)
    assert isinstance(grok._lines, list), grok._lines
    # write to rundir
    grok.write() 
    assert os.path.isfile(self.grok_output), self.grok_output
    
    
if __name__ == "__main__":

    
    specific_tests = []
#     specific_tests += ['SetTime']


    # list of tests to be performed
    tests = [] 
    # list of variable tests
    tests += ['Grok']
    

    # construct dictionary of test classes defined above
    test_classes = dict()
    local_values = locals().copy()
    for key,val in local_values.iteritems():
      if key[-4:] == 'Test':
        test_classes[key[:-4]] = val


    # run tests
    report = []
    for test in tests: # test+'.test'+specific_test
      if len(specific_tests) > 0: 
        test_names = ['hgsrun_test.'+test+'Test.test'+s_t for s_t in specific_tests]
        s = unittest.TestLoader().loadTestsFromNames(test_names)
      else: s = unittest.TestLoader().loadTestsFromTestCase(test_classes[test])
      report.append(unittest.TextTestRunner(verbosity=2).run(s))
      
    # print summary
    runs = 0; errs = 0; fails = 0
    for name,test in zip(tests,report):
      #print test, dir(test)
      runs += test.testsRun
      e = len(test.errors)
      errs += e
      f = len(test.failures)
      fails += f
      if e+ f != 0: print("\nErrors in '{:s}' Tests: {:s}".format(name,str(test)))
    if errs + fails == 0:
      print("\n   ***   All {:d} Test(s) successfull!!!   ***   \n".format(runs))
    else:
      print("\n   ###     Test Summary:      ###   \n" + 
            "   ###     Ran {:2d} Test(s)     ###   \n".format(runs) + 
            "   ###      {:2d} Failure(s)     ###   \n".format(fails)+ 
            "   ###      {:2d} Error(s)       ###   \n".format(errs))
    