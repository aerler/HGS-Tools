# -*- coding: utf-8 -*-
"""
Created on Tue May 23 14:28:22 2017, updated June 23 2017

A script designed for comand line execution, which utilizes the MPPROPlines class to modify the mprops file.
There are two modes:
  1) Add commands to the mprops file to cause Grok to output tabulated values of the van Genuchten (pks) curve 
  2) Modify the mprops file to use the tabulated values for the van Genuchten curve and prepare the table files

@author: Fan Yang and Andre R. Erler
"""

from mprop import MPROPlines

#TODO: add command line arguments


grok_name = 'grw_omafra'

# grok_dirc = 'C:/Users/fyang/Desktop/pcs/data/GRW-V3'
grok_dirc = r'D:\Data\HGS\Templates\GRW-V3-test'

# mprop_path = r'C:\Users\fyang\Desktop\pcs\data\GRW-V3\mprops'
mprop_path = grok_dirc+'/mprops/'

mprop_name = 'GRB.mprops'


mp = MPROPlines(grok_dirc = grok_dirc, mprop_dirc = mprop_path, mprop_name = mprop_name, grok_name = grok_name)

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

mp.walk_mprop(lgentab=True)


