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
from pks_curve.mprop import MPROPlines

def read_mprops(grok_dirc, grok_name, mprop_dirc, mprop_name, 
                lgentab = False, linctab = False, lnotabs = False, ldebug = False,
                sr = 0.13, kr_min = 1.072e-12, tsf = 1e-4, p_min = -1.0e3 ):

## create MPROPlines class instance
  mp = MPROPlines(grok_dirc = grok_dirc, mprop_dirc = mprop_dirc, mprop_name = mprop_name, grok_name = grok_name,
                  sr = sr,
                  kr_min = kr_min,
                  tsf = tsf,
                  p_min = p_min,
                  ldebug = ldebug)

  ## begin execution based on command option

  if lgentab:
      print(("\nModifying mprops file '{}' to output van Genuchten tables. Run Grok to generate table files.".format(mprop_name))) 
      mp.walk_mprop(lgentab=True, ldebug=ldebug)
  elif linctab: 
      print(("\nModifying mprops file '{}' to use tabulated values for van Genuchten function.".format(mprop_name))) 
      mp.walk_mprop(lgentab=False, ldebug=ldebug)
      if not lnotabs:
          print(("\nAssembling van Genuchten table files in folder '{}'.".format(mp.pks_path)))
          mp.combine_pkstable(ldebug=ldebug)
  else: 
      print("\nTo modify the mprops file, specify either the --generate or --include command option - exiting.")


if __name__ == "__main__":
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
  parser.add_argument("--prefix","--problem", dest="prefix", default=None, type=str, 
                      help="the problem name/prefix [default: inferred from batch.pfx]")
  parser.add_argument("--root-folder","--grok-folder", dest="folder", default=None, type=str, 
                      help="the main execution folder for Grok and HGS [default: current working directory]")
  # parameter for plotting van genuchten curve
  parser.add_argument("-s","--sr", dest="sr", type=float, default=0.13,
                      help="Residual Saturation, default as 0.13")
  parser.add_argument("-k","--kr_min", dest="kr_min", type=float, default=1.072e-12,
                      help="minimum relative permeability, default as 1.072e-12")
  parser.add_argument("-t","--tsf", dest="tsf", type=float, default=1e-4,
                      help="table smooth factor, default as 1e-4")
  parser.add_argument("-p","--p_min", dest="p_min", type=float, default=-1.0e3,
                      help="table smooth factor, default as -1.0e3")
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

  # parameters for plotting unsat curve (van genuchten)
  tsf    = args.tsf
  sr     = args.sr
  kr_min = args.kr_min
  p_min  = args.p_min

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
          print(("\nUsing problem name/prefix from file '{}': {}".format(batchpfx,grok_name)))
      else:
          raise IOError("No problem name/prefix specified and no 'batch.pfx' file found in root folder.")

  read_mprops(grok_dirc = grok_dirc, mprop_dirc = mprop_path, mprop_name = mprop_name, grok_name = grok_name,
              lgentab = lgentab, linctab = linctab, ldebug = ldebug,
              sr = sr, kr_min = kr_min, tsf = tsf, p_min = p_min)


