#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue May 23 14:28:22 2017, updated June 23 2017

A script designed for comand line execution, which utilizes the MPPROPlines class to modify the mprops file.
There are two modes:
  1) Add commands to the mprops file to cause Grok to output tabulated values of the van Genuchten (pks) curve 
  2) Modify the mprops file to use the tabulated values for the van Genuchten curve and prepare the table files

@author: Fan Yang and Andre R. Erler

@copyright:  2017 Aquanty Inc.. All rights reserved.

@license:    GPL v3

@contact:    aerler@aquanty.com

@deffield    updated: 24/06/2017
"""

# external imports
import os, sys
from argparse import ArgumentParser, RawDescriptionHelpFormatter
# class module
from mprop import MPROPlines

# meta data
__all__ = []
__version__ = 0.9
__date__ = '2017-05-23'
__updated__ = '2017-06-24'

program_name = os.path.basename(sys.argv[0])
program_version = "v{:}".format(__version__)
program_build_date = str(__updated__)
program_version_message = '%(prog)s {:s} ({:s})'.format(program_version, program_build_date)
program_shortdesc = __import__('__main__').__doc__
program_license = '''  This is the command line interface of the mprop module, a program designed to edit the HGS 
  mprops file for the purpose of generating and using van Genuchten function tables.

  Created by Fan Yang on 23/05/2017, updated by Fan Yang and Andre R. Erler.
  Copyright 2016 Aquanty Inc. All rights reserved.

  Licensed under the GNU General Public License 3.0
  https://www.gnu.org/licenses/gpl.html

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

  Type python {:s} --help for usage information.
'''.format(program_name)


## parse arguments

# Setup argument parser
parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
# mprops file (only positional argument, since it is required)
parser.add_argument('mprops', metavar='mprops', type=str, help="the mprops file to be used/modified")
# execution mode
parser.add_argument("-g", "--generate", dest="gentab", action="store_true", 
                    help="prepare the mprops file to generate van Genuchten tables")
parser.add_argument("-i", "--include", dest="inctab", action="store_true", 
                    help="include table files in mprops file (for HGS execution)")
parser.add_argument("--no-tables", dest="notabs", action="store_true", 
                    help="do not assemble actual tables files when including tables in mprops file")
# problem parameters
parser.add_argument("--prefix","--problem", dest="prefix", nargs=1, default=None, type=str, 
                    help="the problem name/prefix [default: inferred from batch.pfx]")
parser.add_argument("--root-folder","--grok-folder", dest="folder", nargs=1, default=None, type=str, 
                    help="the main execution folder for Grok and HGS [default: current working directory]")
# misc settings
parser.add_argument("--debug", dest="debug", action="store_true", help="print debug output [default: %(default)s]")
parser.add_argument('-V', '-v', '--version', action='version', version=program_version_message)


## assign parameters
args = parser.parse_args()

# command options
lgentab = args.gentab
linctab = args.inctab
lnotabs = args.notabs

# other options
ldebug = args.debug

# grok folder first
grok_dirc = args.folder
if grok_dirc is None: grok_dirc = os.getcwd()
if not os.path.exists(grok_dirc): raise IOError(grok_dirc)
# grok_dirc = 'C:/Users/fyang/Desktop/pcs/data/GRW-V3'
# grok_dirc = r'D:\Data\HGS\Templates\GRW-V3-test'

# now mprops file
mprop_file = args.mprops
if not os.path.exists(mprop_file): raise IOError(mprop_file)
mprop_path = os.path.dirname(mprop_file)
# mprop_path = r'C:\Users\fyang\Desktop\pcs\data\GRW-V3\mprops'
# mprop_path = grok_dirc+'/mprops/'
mprop_name = os.path.basename(mprop_file)
# mprop_name = 'GRB.mprops'

# now prefix/problem name
grok_name = args.prefix
if grok_name is None:
    # infer from batch.pfx file
    batchpfx = os.path.join(grok_dirc,'batch.pfx')
    if os.path.exists(batchpfx):
        with open(batchpfx) as bf:
            grok_name = bf.readline().strip()
        print("\nUsing problem name/prefix from file '{}': {}".format(batchpfx,grok_name))
    else:
        raise IOError("No problem name/prefix specified and no 'batch.pfx' file found in root folder.")
# grok_name = 'grw_omafra'


## create MPROPlines class instance
mp = MPROPlines(grok_dirc = grok_dirc, mprop_dirc = mprop_path, mprop_name = mprop_name, grok_name = grok_name)


## begin execution based on command option

if lgentab:
    print("\nModifying mprops file '{}' to output van Genuchten tables. Run Grok to generate table files.".format(mprop_file)) 
    mp.walk_mprop(lgentab=True)
elif linctab: 
    print("\nModifying mprops file '{}' to use tabulated values for van Genuchten function.".format(mprop_file)) 
    mp.walk_mprop(lgentab=False)
    if not lnotabs:
        print("\nAssembling van Genuchten table files in folder '{}'.".format(mp.pks_path))
        mp.combine_pkstable(ldebug=ldebug)
else: 
    print("\nTo modify the mprops file, specify either the --generate or --include command option - exiting.")

print('')


### This section is for updating the grok to generate pks tables ### 
#get the index where str need to be added 
# mp.gen_tables()

#add in the extra str and create a new file
# mp.write_mprop()

### This section reads in the tables and link them to properties in mprops
#list all the p_c table generated from grok
# a=mp.list_pstable(file_pattern = 'p_s_table', file_path = grok_dirc)
# #print(a)

# mp.combine_pkstable()

#get the name of all the zones in mprops
#zone names is used to refer to the pks table
# mp.get_zone_name()

# mp.comment_out_unsat_func()

# mp.walk_mprop(lgentab=True)


