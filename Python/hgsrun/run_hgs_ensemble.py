#!/usr/local/bin/python
# encoding: utf-8
'''
hgsrun.run_hgs_ensemble -- shortdesc

hgsrun.run_hgs_ensemble is a description

It defines classes_and_methods

@author:     Andre R. Erler

@copyright:  2016 Aquanty Inc.. All rights reserved.

@license:    GPL v3

@contact:    aerler@aquanty.com
@deffield    updated: 16/08/2016
'''

# external imports
import sys, os, yaml
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
# hgsrun imports
from hgs_ensemble import EnsHGS

# meta data
__all__ = []
__version__ = 0.1
__date__ = '2016-08-16'
__updated__ = '2016-08-16'

DEBUG = 1
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
    '''Command line options.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    #program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    program_license = '''  This is the command line interface of HGSrun, a program to manage HGS simulations.
  HGSrun is capable of running single simulations or entire ensembles of simulations.

  Created by Andre R. Erler on 16/08/2016
  Copyright 2016 Aquanty Inc. All rights reserved.

  Licensed under the GNU General Public License 3.0
  http://www.apache.org/licenses/LICENSE-2.0

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
'''

    try:
        # Setup argument parser
        parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
        parser.add_argument("-v", "--verbose", dest="verbose", action="count", help="set verbosity level [default: %(default)s]")
        parser.add_argument('-V', '--version', action='version', version=program_version_message)
        parser.add_argument('config', help="the ensemble configuration file (in YAML format)", type=str )

        # Process arguments
        args = parser.parse_args()

        yamlfile = args.config
        verbose = args.verbose

        if verbose > 0:
            print("Verbose mode on")

        ## load arguments from config file
        if not os.path.exists(yamlfile):
            raise CLIError("Configuration file '{:s}' not found!".format(yamlfile))          
        # read file 
        with open(yamlfile) as f: 
          config = yaml.load(f, Loader=yaml.Loader)
        # instantiate ensemble
        if verbose > 0:
          print("\nYAML Configuration that will be passed to EnsHGS class:")
          print(config['HGS_parameters'])
        enshgs = EnsHGS(**config['HGS_parameters'])

        ## here we are actually running the program
        if verbose > 0:
          print("\nYAML Configuration for batch execution:")
          print(config['batch_config'])
        ec = enshgs.runSimulations(**config['batch_config'])
        # check results
        if ec == 0 and all(g for g in enshgs.HGSOK):
          print("\n   ***   All HGS Simulations Completed Successfully!!!   ***\n")
        elif not any(g for g in enshgs.HGSOK):
          print("\n   ###   All HGS Simulations Failed!!!   ###\n")
        else:
          print("\n   ===   Some HGS Simulations Failed!!!   ===\n")
          print(  "         {} out of {} completed ".format(sum(g for g in enshgs.HGSOK), len(enshgs)))
        
        # return exit code of HGS ensemble
        return ec
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
    except Exception as e:
        if DEBUG or TESTRUN:
            raise(e)
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        return 2


if __name__ == "__main__":
  
    if DEBUG:
        #sys.argv.append("-h")
        sys.argv.append("-v")
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