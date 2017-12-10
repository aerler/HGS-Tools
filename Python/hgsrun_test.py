'''
Created on Aug 11, 2016

Unittests for hgsrun components.

@author: Andre R. Erler, GPL v3
'''

import unittest
import numpy as np
import os, sys, gc, shutil
import subprocess
# from subprocess import STDOUT

# WindowsError is not defined on Linux - need a dummy
try: 
    lWin = True
    WindowsError
except NameError:
    lWin = False
    WindowsError = None

# import modules to be tested
from hgsrun.hgs_setup import lWin, clearFolder
from hgsrun.hgs_setup import Grok, GrokError, HGS, HGSError
from hgsrun.hgs_ensemble import EnsHGS, EnsembleError

# work directory settings ("global" variable)
data_root = os.getenv('HGS_ROOT', '')
# test folder either RAM disk or data directory
RAM = bool(os.getenv('RAMDISK', '')) # whether or not to use a RAM disk
# RAM = WindowsError is None # RAMDISK may not exist on Windows...
loverwrite = RAM # copying folders on disk takes long...

# N.B.: the environment variable RAMDISK contains the path to the RAM disk
workdir = os.getenv('RAMDISK', '') if RAM else '{:s}/test/'.format(data_root)
if not os.path.isdir(workdir): raise IOError(workdir)
# other settings
NP = 2 # for parallel tests
ldebug = False # additional debug output
lbin = True # execute binaries to test runner methods
# lbin = False # don't execute binaries (takes very long)

# executables
grok_bin = 'grok.exe' # Grok executable
hgs_bin  = 'phgs.exe' # HGS executable
hgsdir   = os.getenv('HGSDIR',) # HGS license file    

# # use old GRW model for testing
# hgs_testcase = 'grw_omafra' # name of test project (for file names)
# hgs_template = data_root+'/Test/Templates/GRW-test/' 
# test_prefix  = 'grw2' # pefix for climate input
# clim_data    = data_root+'/Test/grw2/test-run/clim/climate_forcing/'
# ts_data      = data_root+'/Test/grw2/test-run/timeseries/climate_forcing/'

# use small Payne River model for testing
hgs_testcase = 'test' # name of test project (for file names)
hgs_template = data_root+'/Test/Templates/PRC-test/' 
test_prefix  = 'snw1' # pefix for climate input
clim_data    = data_root+'/Test/snw1/test-run/clim/climate_forcing/'
ts_data      = data_root+'/Test/snw1/test-run/timeseries/climate_forcing/'

# # use new hires GRW model for testing
# hgs_testcase = 'grw_omafra' # name of test project (for file names)
# hgs_template = data_root+'/GRW/Templates/GRW-V3-hires/' 
# test_prefix  = 'grw2' # pefix for climate input
# clim_data    = data_root+'/Test/grw2/test-run/clim/climate_forcing/'
# ts_data      = data_root+'/Test/grw2/test-run/timeseries/climate_forcing/'

## tests for Grok class
class GrokTest(unittest.TestCase):  
  # some Grok test data
  hgs_template = hgs_template  
  clim_data    = clim_data 
  ts_data      = ts_data
  hgs_testcase = hgs_testcase
  test_prefix  = test_prefix  
  grok_bin     = grok_bin
  lvalidate    = False # validate input data
  rundir       = '{}/grok_test/'.format(workdir,) # test folder
  # some grok settings
  runtime = 5*365*24*60*60 # two years in seconds
  input_interval = 'monthly'
  input_mode = 'periodic'
   
  def setUp(self):
    ''' initialize a Grok instance '''
    if not os.path.isdir(self.hgs_template): 
      raise IOError("HGS Template for testing not found:\n '{}'".format(self.hgs_template))
    clearFolder(self.rundir) # make a clean folder
    # grok test files
    self.grok_input  = '{}/{}.grok'.format(self.hgs_template,self.hgs_testcase)
    self.grok_output = '{}/{}.grok'.format(self.rundir,self.hgs_testcase)
    # create Grok instance
    self.grok = Grok(rundir=self.rundir, project=self.hgs_testcase, runtime=self.runtime,
                     input_mode=self.input_mode, input_interval=self.input_interval)
    # load a config file from template
    if not os.path.isfile(self.grok_input):
      raise IOError("Grok configuration file for testing not found:\n '{}'".format(self.grok_input))
    self.grok.readConfig(folder=self.hgs_template)
    if self.grok.runtime is not None: # need to set runtime manually here 
      self.grok.setRuntime(runtime=self.grok.runtime)
    assert isinstance(self.grok._lines, list), self.grok._lines
      
  def tearDown(self):
    ''' clean up '''
    self.grok.writeConfig()
    del self.grok
    gc.collect()
    # leave directory, so we can delete folder
    os.chdir(os.path.expanduser('~'))

  def testClass(self):
    ''' test instantiation of class '''    
    # instantiation done in self.setUp()
    assert self.grok.rundir, self.grok.rundir
    assert os.path.isdir(self.grok.rundir), self.grok.rundir
    
  def testInputLists(self):
    ''' test writing of input list files with climate forcings '''
    grok = self.grok  
    # write lists for fictional scenario
    grok.generateInputLists(input_vars='WRFPET', input_prefix=self.test_prefix, pet_folder=self.clim_data,
                            input_folder=self.clim_data, lvalidate=self.lvalidate,)
    # convert config file list into string and verify
    output = ''.join(grok._lines) # don't need newlines 
    assert 'precip.inc' in output
    assert 'pet.inc' in output    
    
  def testRunGrok(self):
    ''' test the Grok runner command (will fail, because other inputs are missing) '''
    grok = self.grok  
    # write to rundir
    grok.writeConfig() 
    exe = os.path.join(self.hgs_template,self.grok_bin)
    logfile = '{}/log.grok'.format(grok.rundir)
    # run Grok
    if not os.path.isfile(exe): raise IOError(exe)
    try: 
      ec = grok.runGrok(executable=exe, logfile=logfile, ldryrun=not lbin)
      if not lbin: ec = 1 # because the test should fail...
    except GrokError: ec = 1
    # check output
    batchpfx = '{}/batch.pfx'.format(self.rundir)
    assert os.path.isfile(batchpfx), batchpfx
    assert os.path.isfile(logfile), logfile
    assert ec > 0, ec
    
  def testSetTime(self):
    ''' test setting the run time variable in Grok config file '''
    grok = self.grok
    time = 24*60*60 # in seconds, i.e. one day
    # set run time
    ec = grok.setRuntime(time)
    # test
    assert grok.runtime == time and ec == 0
    # test output times
    outtimes = grok.getParam('output times', dtype='float', llist=None)
    assert all(np.diff(outtimes) > 0), outtimes
    lenot = np.sum(grok.output_interval)-len(grok.output_interval)+1
    assert len(outtimes) == lenot, grok.output_interval  
    # convert config file list into string and verify
    output = ''.join(grok._lines) # don't need newlines 
    #print(output)
    assert '{:e}'.format(time) in output, '{:e}'.format(time)
    
  def testWrite(self):
    ''' load config file from template and write to rundir (on disk) '''
    grok = self.grok
    grok._lines = None # delete already loaded file contents
    # read from template
    assert os.path.isfile(self.grok_input), self.grok_input 
    grok.readConfig(folder=self.hgs_template)
    assert isinstance(grok._lines, list), grok._lines
    # write to rundir
    grok.writeConfig() 
    assert os.path.isfile(self.grok_output), self.grok_output
    

## tests for HGS class
class HGSTest(GrokTest):  
  # some HGS test data
  hgs_bin      = hgs_bin
  lvalidate    = True # validate input data
  rundir       = '{}/hgs_test/'.format(workdir,) # test folder
  hgsdir       = os.getenv('HGSDIR',) # HGS license file
  # some grok settings
  runtime = 5*365*24*60*60 # five years in seconds
  input_interval = 'monthly'
#   input_mode = 'periodic'
  input_mode = 'quasi-transient'
   
  def setUp(self):
    ''' initialize an HGS intance '''
    if not os.path.isdir(self.hgs_template): 
      raise IOError("HGS Template for testing not found:\n '{}'".format(self.hgs_template))
    if lbin and not os.path.exists(hgsdir+'/hgs.lic'):
      raise HGSError("License file for HGS not found:\n '{}/hgs.lic'".format(hgsdir))
    clearFolder(self.rundir) # provide clean folder for testing 
    # grok test files
    self.grok_input  = '{}/{}.grok'.format(self.hgs_template,self.hgs_testcase)
    self.grok_output = '{}/{}.grok'.format(self.rundir,self.hgs_testcase)
    # data sources
    if self.input_mode == 'periodic':
        input_folder = self.clim_data
        pet_folder   = self.clim_data
    elif self.input_mode == 'quasi-transient': 
        input_folder = self.ts_data
        pet_folder   = self.clim_data
    elif self.input_mode == 'transient': 
        input_folder = self.ts_data
        pet_folder   = self.ts_data
    else:
        raise ValueError(self.input_mode)
    # HGS settings
    self.NP = NP
    # create Grok instance
    self.hgs = HGS(rundir=self.rundir, project=self.hgs_testcase, runtime=self.runtime,
                   input_mode=self.input_mode, input_interval=self.input_interval, 
                   input_prefix=self.test_prefix, input_folder=input_folder, pet_folder=pet_folder,
                   template_folder=self.hgs_template, NP=self.NP)
    self.grok = self.hgs
    # load a config file from template
    if not os.path.isfile(self.grok_input):
      raise IOError("Grok configuration file for testing not found:\n '{}'".format(self.grok_input))
    self.grok.readConfig(folder=self.hgs_template)
    assert isinstance(self.hgs._lines, list), self.hgs._lines

  def tearDown(self):
    ''' clean up '''
    self.grok.writeConfig()
    del self.grok, self.hgs
    gc.collect()
    # leave directory, so we can delete folder
    os.chdir(os.path.expanduser('~'))

  def testInputLists(self):
    ''' test writing of input list files with climate forcings '''
    grok = self.grok  
    # write lists for fictional scenario
    grok.generateInputLists(lvalidate=self.lvalidate,)
    # convert config file list into string and verify
    output = ''.join(grok._lines) # don't need newlines 
    assert 'precip.inc' in output
    assert 'pet.inc' in output    

  def testParallelIndex(self):
    ''' test writing of parallelindex file '''
    hgs = self.hgs
    pidx_file = self.rundir+hgs.pidx_file
    # write parallelindex file to rundir
    hgs.writeParallelIndex(NP=1, parallelindex=pidx_file) 
    assert os.path.isfile(pidx_file), pidx_file
  
  def testRestart(self):
    ''' load config file from rundir and modify restart time etc. '''
    hgs = self.hgs
    os.chdir(hgs.rundir)
    # write to rundir
    hgs.writeConfig() 
    assert os.path.isfile(self.grok_output), self.grok_output
    old_times = hgs.getParam('output times', dtype='float', llist=True)
    hgs._lines = None # delete already loaded file contents
    # create fake output files for restart
    for i in range(5):
        pm_file = os.path.join(hgs.rundir,hgs.pm_files.format(IDX=i+1))
        open(pm_file,'w').close()
        olf_file = os.path.join(hgs.rundir,hgs.olf_files.format(IDX=i+1))
        open(olf_file,'w').close()
        if hgs.lchannel:
            chan_file = os.path.join(hgs.rundir,hgs.chan_files.format(IDX=i+1))
            open(chan_file,'w').close()
    # create fake newton_info and water_balance files with time-stamps
    with open(hgs.newton_file, 'w') as nf:
        nf.write('Title\nVARIABLES\nzone\n')
        nf.writelines(['{:e}\n'.format(t) for t in old_times[:5]])
    #with open(hgs.newton_file, 'r') as nf:
    #    for line in nf.readlines(): print(line)
    with open(hgs.water_file, 'w') as wf:
        wf.write('Title\nVARIABLES\nzone\n')
        wf.writelines(['{:e}\n'.format(t) for t in old_times[:5]])
    # create fake log files, too
    open(hgs.grok_log, 'a').close(); open(hgs.hgs_log, 'a').close() 
    # read from template
    hgs.readConfig(folder=self.rundir)
    assert isinstance(hgs._lines, list), hgs._lines
    new_times = hgs.getParam('output times', dtype='float', llist=None)
    assert all([old == new for old,new in zip(old_times,new_times)]), old_times
    assert all(np.diff(new_times) > 0), np.diff(new_times)
    # apply modifications for restart
    restart_file = hgs.rewriteRestart()
    test_file = os.path.join(self.rundir,restart_file.format(FILETYPE=hgs.olf_tag))
    assert os.path.isfile(test_file), test_file
    # write new restart files into Grok file
    hgs.changeICs(ic_pattern=restart_file)
    # write modified file to rundir
    hgs.writeConfig()
    assert os.path.isfile(self.grok_output), self.grok_output
    # verify grok file
    hgs._lines = None # delete already loaded file contents
    hgs.readConfig(folder=self.rundir)
    assert isinstance(hgs._lines, list), hgs._lines
    new_times = hgs.getParam('output times', dtype='float', llist=None)
    assert len(new_times) == max(len(old_times)-5,0), (len(new_times),len(old_times))
    assert all([old == new for old,new in zip(old_times[5:],new_times)]), old_times
    assert all(np.diff(new_times) > 0), np.diff(new_times)
    
  def testRunGrok(self):
    ''' test the Grok runner command and check if flag is set correctly '''
    hgs = self.hgs  
    exe = os.path.join(self.hgs_template,self.grok_bin) # run folder not set up
    logfile = '{}/log.grok'.format(hgs.rundir)
    assert hgs.GrokOK is None, hgs.GrokOK
    # climate data
    if not os.path.isdir(self.clim_data): raise IOError(self.clim_data)
    if not os.path.isdir(self.ts_data): raise IOError(self.ts_data)
    # run Grok
    if not os.path.isfile(exe): raise IOError(exe)
    ec = hgs.runGrok(executable=exe, logfile=logfile, lerror=False, linput=lbin, ldryrun=not lbin)
    assert ec == ( 1 if lbin else 0 ), ec # since we did not set up the run folder, an actual run will fail
    # check output
    batchpfx = self.rundir+hgs.batchpfx
    assert os.path.isfile(batchpfx), batchpfx
    assert os.path.isfile(logfile), logfile
    # check flag
    assert hgs.GrokOK is not lbin, hgs.GrokOK

  def testRunHGS(self):
    ''' test the HGS runner command (will fail without proper folder setup) '''
    hgs = self.hgs  
    exe = os.path.join(self.hgs_template,self.hgs_bin) # run folder not set up
    logfile = '{}/log.hgs_run'.format(hgs.rundir)
    assert hgs.GrokOK is None, hgs.GrokOK
    # print environment variable for license file
    print('\nHGSDIR: {}'.format(self.hgsdir))
    # attempt to run HGS
    if not os.path.isfile(exe): raise IOError(exe)
    # run setup
    ec = hgs.setupRundir(template_folder=self.hgs_template, loverwrite=loverwrite, bin_folder=None)
    assert ec == 0, ec
    assert hgs.rundirOK, hgs.rundir
    indicator = '{}/SCHEDULED'.format(self.rundir)
    assert os.path.exists(indicator), indicator
    # create fake input file to cause crash
    open(os.path.join(hgs.rundir,hgs.problem+'o.gen'),'w').close()
    # create some required files without running full setup
    ec = hgs.runHGS(executable=exe, logfile=logfile, skip_grok=True, ldryrun=not lbin, lerror=False)
    assert ec == ( 1 if lbin else 0 ), ec # since we did not run Grok etc., an actual run will fail
    # check output
    pidx_file = self.rundir+hgs.pidx_file
    assert os.path.isfile(pidx_file), pidx_file
    assert os.path.isfile(logfile), logfile
    if hgs.HGSOK: indicator = '{}/COMPLETED'.format(self.rundir)
    else: indicator = '{}/FAILED'.format(self.rundir)
    assert os.path.isfile(indicator), indicator
    # check flag
    assert hgs.GrokOK is None, hgs.GrokOK # we skipped Grok

  def testSetup(self):
    ''' test copying of a run folder from a template '''
    hgs = self.hgs
    if not os.path.isdir(self.hgs_template): raise IOError(self.hgs_template)
    # run setup
    hgs.setupRundir(template_folder=self.hgs_template, loverwrite=loverwrite, bin_folder=None)
    # check that all items are there
    assert os.path.isdir(self.rundir), self.rundir
    for exe in (self.hgs_bin, self.grok_bin):
      local_exe = os.path.join(self.rundir,exe)
      assert os.path.exists(local_exe), local_exe
    indicator = '{}/SCHEDULED'.format(self.rundir)
    assert os.path.exists(indicator), indicator
#     print('\nRundir: {}'.format(self.rundir))
    assert hgs.rundirOK, hgs.rundir
   
    
## tests for EnsHGS class
class EnsHGSTest(unittest.TestCase):  
  # some HGS test data
  hgs_template = hgs_template 
  hgs_testcase = hgs_testcase
  clim_data    = clim_data
  ts_data      = ts_data
  test_prefix  = test_prefix
  rundir       = '{}/enshgs_test/'.format(workdir,) # test folder
  grok_bin     = grok_bin
  hgs_bin      = hgs_bin
  hgsdir       = hgsdir  
  # some grok settings
  runtime = 5*365*24*60*60 # two years in seconds
  input_interval = 'monthly'
  input_mode = 'periodic'
#   input_mode = 'quasi-transient'
   
  def setUp(self):
    ''' initialize an HGS ensemble '''
    if not os.path.isdir(self.hgs_template): 
      raise IOError("HGS Template for testing not found:\n '{}'".format(self.hgs_template))
    if lbin and not os.path.exists(hgsdir+'/hgs.lic'):
      raise HGSError("License file for HGS not found:\n '{}/hgs.lic'".format(hgsdir))    
    clearFolder(self.rundir) # provide clean folder for testing
    # grok test files
    self.grok_input  = '{}/{}.grok'.format(self.hgs_template,self.hgs_testcase)
    self.grok_output = '{}/{}.grok'.format(self.rundir,self.hgs_testcase)
    # data sources
    if self.input_mode == 'periodic':
        input_folder = self.clim_data
        pet_folder   = self.clim_data
    elif self.input_mode == 'quasi-transient': 
        input_folder = self.ts_data
        pet_folder   = self.clim_data
    elif self.input_mode == 'transient': 
        input_folder = self.ts_data
        pet_folder   = self.ts_data
    else:
        raise ValueError(self.input_mode)
    # HGS settings
    self.NP = NP
    # create 2-member ensemble
    self.enshgs = EnsHGS(rundir=self.rundir + "/{A}/", project=self.hgs_testcase, runtime=self.runtime,
                         output_interval=(2,12), input_mode=self.input_mode, input_interval=self.input_interval, 
                         input_prefix=self.test_prefix, input_folder=input_folder, pet_folder=pet_folder,
                         NP=self.NP, A=['A1','A2'], outer_list=['A'], template_folder=self.hgs_template,
                         loverwrite=True)
    # load a config file from template
    if not os.path.isfile(self.grok_input):
      raise IOError("Grok configuration file for testing not found:\n '{}'".format(self.grok_input))

  def tearDown(self):
    ''' clean up '''
    del self.enshgs
    gc.collect()
    
  def testInitEns(self):
    ''' initialize the an HGS ensemble using list expansion '''
    # define simple rundir pattern
    rundir = self.rundir + "/{A}/{B}/{C}/"
    rundir_args = dict(A=['A1','A2'], B=['B'], C=['C1','C2'], outer_list=['A','B',('C','input_mode')])
    # initialize ensemble with general and rundir arguments
    enshgs = EnsHGS(rundir=rundir, project=self.hgs_testcase, runtime=self.runtime,
                    input_mode=['steady-state','periodic'], input_interval=self.input_interval, 
                    NP=self.NP, **rundir_args)
    assert len(enshgs) == 4, len(enshgs)
    # test expansion of folder arguments
    for As in rundir_args['A']:
      for Bs in rundir_args['B']:
        for Cs in rundir_args['C']:        
          tmpdir = rundir.format(A=As, B=Bs, C=Cs)
          assert tmpdir in enshgs.rundirs, tmpdir
          # test concurrent expansion of input_mode with folder argument C
          i = enshgs.rundirs.index(tmpdir)
          if Cs == 'C1': assert enshgs.members[i].input_mode == 'steady-state'
          if Cs == 'C2': assert enshgs.members[i].input_mode == 'periodic'
    
  def testSetTime(self):
    ''' mainly just test the method application mechanism '''
    enshgs = self.enshgs
    assert all(rt == self.runtime for rt in enshgs.runtime), enshgs.runtime
    # load Grok file
    enshgs.readConfig(folder=self.hgs_template)
    time = 1 # new runtime
    # setter method
    enshgs.setRuntime(runtime=time)
    assert all(hgs.runtime == time for hgs in enshgs), enshgs.members # EnsHGS supports iteration over members
    assert all(rt == time for rt in enshgs.runtime), enshgs.runtime # EnsHGS supports iteration over members
    #print(enshgs.runtime)
    # using new ensemble wrapper
    time = 10
    enshgs.runtime = time
    assert all(hgs.runtime == time for hgs in enshgs), enshgs.members # EnsHGS supports iteration over members
    assert all(rt == time for rt in enshgs.runtime), enshgs.runtime # EnsHGS supports iteration over members
    
  def testRunEns(self):
    ''' test running the ensemble; the is the primary application test '''
    enshgs = self.enshgs
    assert all(not g for g in enshgs.HGSOK), enshgs.HGSOK
    # check license
    print('\nHGSDIR: {}'.format(self.hgsdir))
    # setup run folders and run Grok
    enshgs.runSimulations(lsetup=True, lgrok=True, loverwrite=loverwrite, skip_grok=True, lparallel=True, NP=NP, 
                          runtime_override=120, ldryrun=not lbin) # set runtime to 2 minutes
    assert not lbin or all(g for g in enshgs.HGSOK), enshgs.HGSOK
    for rundir in enshgs.rundirs:
      assert os.path.isdir(rundir), rundir
      if lbin:
        hgslog = '{0}/log.hgs_run'.format(rundir)
        assert os.path.isfile(hgslog), hgslog

  def testSetupExp(self):
    ''' test experiment setup '''
    enshgs = self.enshgs
    assert all(not g for g in enshgs.GrokOK), enshgs.GrokOK
    # setup run folders and run Grok
    enshgs.setupExperiments(lgrok=lbin, loverwrite=loverwrite, lparallel=True)
    assert not lbin or all(g for g in enshgs.GrokOK), enshgs.GrokOK
    for rundir in enshgs.rundirs:
      assert os.path.isdir(rundir), rundir
      if lbin:
        grokfile = '{0}/{1}.grok'.format(rundir,self.hgs_testcase)
        assert os.path.isfile(grokfile), grokfile

  def testSetupRundir(self):
    ''' test rundir setup '''
    enshgs = self.enshgs
    assert all(not g for g in enshgs.GrokOK), enshgs.GrokOK
    # setup run folders
    enshgs.setupExperiments(lgrok=False, loverwrite=loverwrite, lparallel=True)
    for rundir in enshgs.rundirs:
      assert os.path.isdir(rundir), rundir
      grok_bin = '{0}/{1}'.format(rundir,self.grok_bin)
      assert os.path.exists(grok_bin), grok_bin


if __name__ == "__main__":

    
    specific_tests = []
#     specific_tests += ['Class']
#     specific_tests += ['InitEns']
#     specific_tests += ['InputLists']
#     specific_tests += ['ParallelIndex']
#     specific_tests += ['Restart']
#     specific_tests += ['RunEns']
#     specific_tests += ['RunGrok']
#     specific_tests += ['RunHGS']
#     specific_tests += ['SetTime']
#     specific_tests += ['Setup']
#     specific_tests += ['SetupExp']
#     specific_tests += ['SetupRundir']
#     specific_tests += ['Write']


    # list of tests to be performed
    tests = [] 
    # list of variable tests
    tests += ['Grok']
    tests += ['HGS']    
#     tests += ['EnsHGS']

    # construct dictionary of test classes defined above
    test_classes = dict()
    local_values = locals().copy()
    for key,val in local_values.items():
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
    