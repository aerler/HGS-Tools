#!/usr/bin/python
# encoding: utf-8
'''
This is the command line interface of HGSrun, a program to manage HGS simulations.
HGSrun is capable of running single simulations or entire ensembles of simulations.

@author:     Andre R. Erler

@copyright:  2016 Aquanty Inc.. All rights reserved.

@license:    GPL v3

@contact:    aerler@aquanty.com
@deffield    updated: 06/03/2019
'''

# external imports
import sys, os, yaml
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
# hgsrun imports
from hgsrun.hgs_ensemble import EnsHGS

# meta data
__all__ = []
__version__ = 0.10
__date__ = '2016-08-16'
__updated__ = '2019-03-06'

DEBUG = 0
TESTRUN = 0
PROFILE = 0

class CLIError(Exception):
    '''Generic exception to raise and log different fatal errors.'''
    def __init__(self, msg):
        super(CLIError).__init__(type(self))
        self.msg = "E: %s" % msg
    def __str__(self):
        return self.msg
    def __unicode__(self):
        return self.msg

def main(argv=None): # IGNORE:C0111
    ''' parse ommand line options '''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    program_name = os.path.basename(sys.argv[0])
    program_version = "v{:}".format(__version__)
    program_build_date = str(__updated__)
    program_version_message = '%(prog)s {:s} ({:s})'.format(program_version, program_build_date)
    program_shortdesc = __import__('__main__').__doc__
    program_license = '''  This is the command line interface of HGSrun, a program to manage HGS simulations.
  HGSrun is capable of running single simulations or entire ensembles of simulations.
  The program is controlled through a YAML configuration file which has two major sections:
    HGS_parameters: parameters to initialize the HGS Ensemble (passed to EnsHGS.__init__)
    batch_config: run-time parameters for parallel batch execution (passed to EnsHGS.runSimulations)

  Created by Andre R. Erler on 16/08/2016
  Copyright 2016 Aquanty Inc. All rights reserved.

  Licensed under the GNU General Public License 3.0
  https://www.gnu.org/licenses/gpl.html

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

  Type python {:s} --help for usage information.
'''.format(program_name)

#     try:
    # Setup argument parser
    parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('config', metavar='config.yaml', type=str, help="the ensemble configuration file (in YAML format)")
    parser.add_argument('-V', '-v', '--version', action='version', version=program_version_message)
    parser.add_argument("--debug", dest="debug", action="store_true", help="print debug output [default: %(default)s]")
    parser.add_argument("-q", "--quiet", dest="quiet", action="store_true", help="do not print status reports [default: %(default)s]")
    #parser.add_argument("-v", "--verbose", dest="verbose", action="count", help="set verbosity level [default: %(default)s]")
    parser.add_argument("--ignore-indicator", dest="noindicator", action='store_true', help="run all simulations, ignoring indicator files [default: %(default)s]")
    parser.add_argument("--overwrite", dest="overwrite", action='store_true', help="overwrite all run folders, ignoring indicator files [default: %(default)s]")
    parser.add_argument("--rerun-failed", dest="runfailed", action='store_true', help="rerun failed experiments, ignore completed [default: %(default)s]")
    parser.add_argument('--data-root', nargs='?', const=None, default=None, type=str, 
                        help="Override the root folder of the data archive [default: $DATA_ROOT or config.yaml]")
    parser.add_argument('--hgs-root', nargs='?', const=None, default=None, type=str, 
                        help="Override the HGS data root folder [default: $HGS_ROOT or config.yaml]")
    parser.add_argument('--hgsdir', nargs='?', const=None, default=None, type=str, 
                        help="Override the root folder for the HGS license and fall-back executables [default: $HGSDIR or config.yaml]")
    parser.add_argument("--skip-setup", dest="nosetup", action='store_true', help="skip run folder setup; start simulations immediately [default: %(default)s]")
    parser.add_argument("--only-setup", dest="nosim", action='store_true', help="only set up run folder; don't start simulation [default: %(default)s]")
    parser.add_argument("--grok-first", dest="grok", action='store_true', 
                        help="run Grok for all folders during setup [default: %(default)s]")
    parser.add_argument("--skip-grok", dest="skipgrok", action='store_true', help="do not run Grok at all [default: %(default)s]")
    parser.add_argument("--restart", dest="restart", action='store_true', help="complete an ensemble, restarting simulations 'in progress' [default: %(default)s]")
    parser.add_argument("--dry-run", dest="dryrun", action='store_true', 
                        help="do not actually run simulations [default: %(default)s]")
    parser.add_argument("-np", "-n", "--processes", dest="NP", default=None, type=int, 
                        help="number of concurrent simulations to run [default: number of available CPUs]")
    parser.add_argument("--serial", dest="serial", action='store_true', 
                        help="run batch execution in serial mode [default: %(default)s]")
    parser.add_argument("--runtime", dest="runtime", default=None, type=int, 
                        help="override runtime setting [default: set in Grok configuration file]")
    
    # Process arguments
    args = parser.parse_args()
    
    ldebug       = args.debug
    lquiet       = args.quiet
    lnoindicator = args.noindicator
    loverwrite   = args.overwrite
    lrunfailed   = args.runfailed
    yamlfile     = args.config
    data_root    = args.data_root
    hgs_root     = args.hgs_root
    hgsdir       = args.hgsdir
    lnosetup     = args.nosetup
    lnosim       = args.nosim
    lgrok        = args.grok
    lskipgrok    = args.skipgrok
    lrestart     = args.restart
    ldryrun      = args.dryrun
    NP           = args.NP
    lserial      = args.serial
    runtime      = args.runtime
    

    if ldebug: print("Printing Debug Output")

    ## load arguments from config file
    if not os.path.exists(yamlfile):
        raise CLIError("Configuration file '{:s}' not found!".format(yamlfile))          
    
    # read YAML configuration file 
    with open(yamlfile) as f: 
      config = yaml.load(f, Loader=yaml.Loader)
    # parameters to initialize the HGS Ensemble (passed to EnsHGS.__init__)
    hgs_config = config['HGS_parameters']
    # run-time parameters for parallel batch execution (passed to EnsHGS.runSimulations)
    batch_config = config['batch_config']
    
    # add some variables based on cli override, YAML file, and environment variables
    for envval,envvar in [(data_root,'DATA_ROOT'),(hgs_root,'HGS_ROOT'),(hgsdir,'HGSDIR')]:
        tmpvar = os.getenv(envvar, None)
        if envval is not None: # override
            hgs_config[envvar] = envval 
            os.environ[envvar] = envval
        elif envvar in hgs_config: # use that value
            os.environ[envvar] =  hgs_config[envvar]
        elif tmpvar: # use environment variable 
            hgs_config[envvar] = tmpvar
    
    # override some settings with command-line arguments
    if lnoindicator: hgs_config['lindicator'] = False
    if loverwrite: hgs_config['loverwrite'] = True
    if lrunfailed: hgs_config['lrunfailed'] = True
    if lrestart: hgs_config['lrestart'] = True
    
    # instantiate ensemble
    if not lquiet:
      print('\n\n   ---              Setting up Run Folders               ---\n') # two newlines
    if ldebug:
      print("\nYAML Configuration that will be passed to EnsHGS class:")
      print(hgs_config)
      print('')
    enshgs = EnsHGS(**hgs_config)
    if lquiet: enshgs.lreport = False

    ## here we are actually running the program
    if 'lsetup' in batch_config:
        print("Run folder setup is handled via command line arguments; removing batch option 'lsetup'.")
        del  batch_config['lsetup']
    if 'lgrok' in batch_config: del batch_config['lgrok']
    if lskipgrok: batch_config['skip_grok'] = True
    if ldryrun: batch_config['ldryrun'] = True
    if NP is not None: batch_config['NP'] = NP
    if lserial or NP == 1: batch_config['lparallel'] = False
    elif NP > 1: batch_config['lparallel'] = True
    if runtime is not None: batch_config['runtime_override'] = runtime
    
    # run setup
    if lnosetup:
      
        # mark experiments as scheduled, before we begin batch execution
        if not lnoindicator:
          for m in enshgs:
              open('{}/SCHEDULED'.format(m.rundir),'a').close()
              
    else:
        
        # run setup for experiments
        ec = enshgs.setupExperiments(lschedule=not (lnosim or lnoindicator), lgrok=(lgrok and not lskipgrok), **batch_config) 
        # don't use indicators, if we don't actuall run a batch execution
      
        # report results
        if ec == 0 and all(g for g in enshgs.rundirOK):
          if lnosim:
              print("\n   ***   All HGS Run Folder Setups Completed Successfully!!!   ***\n")
        elif not any(g for g in enshgs.rundirOK):
          print("\n   ###           All HGS Setups Failed!!!           ###\n")
        else:
          print("\n   ===           Some HGS Setups Failed!!!          ===\n")
          print((  "         {} out of {} completed ".format(sum(g for g in enshgs.rundirOK), len(enshgs))))
      
    
    # run simulations
    if not lnosim:
        if not lquiet:
          print('\n\n   ---             Beginning Batch Execution             ---\n') # two newlines
        if ldebug:
          print("\nYAML Configuration for batch execution:")
          print(batch_config)
          print('')
        
        # begin actual batch execution
        ec = enshgs.runSimulations(lsetup=False, **batch_config) # setup handled above
    
        # check results
        if not lquiet:
          print('\n') # two newlines
          if ec == 0 and all(g for g in enshgs.HGSOK):
            print("\n   ***   All HGS Simulations Completed Successfully!!!   ***\n")
          elif not any(g for g in enshgs.HGSOK):
            print("\n   ###           All HGS Simulations Failed!!!           ###\n")
          else:
            print("\n   ===           Some HGS Simulations Failed!!!          ===\n")
            print((  "         {} out of {} completed ".format(sum(g for g in enshgs.HGSOK), len(enshgs))))
          print('') # one newline
    
    # return exit code of HGS ensemble
    return ec
      
#     except KeyboardInterrupt:
#         ### handle keyboard interrupt ###
#         return 0
#     except Exception as e:
#         if DEBUG or TESTRUN:
#             raise(e)
#         indent = len(program_name) * " "
#         sys.stderr.write(program_name + ": " + repr(e) + "\n")
#         sys.stderr.write(indent + "  for help use --help")
#         return 2


if __name__ == "__main__":
  
    if DEBUG:
        #sys.argv.append("-h")
        sys.argv.append("-d")
    if TESTRUN:
        import doctest
        doctest.testmod()
    if PROFILE:
        import cProfile
        import pstats
        profile_filename = 'hgsrun.run_hgs_ensemble_profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    sys.exit(main())
