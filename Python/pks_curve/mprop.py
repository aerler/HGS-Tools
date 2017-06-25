# -*- coding: utf-8 -*-
"""
Created on Tue May 23 14:28:22 2017, updated June 23 2017

A module defining the MPPROPlines class. This class provides functions to edit the mprops file for the purpose 
of outputting tabulated values of the van Genuchten function or to remove the van Genuchten function definitions
in order to read them from input tables. It also assembles the tabulated values written by Grok in the form 
that is required to read them as file inputs.

@author: Fan Yang and Andre R. Erler
"""


import os, glob, warnings, shutil


# TODO: clean up, test table file concatenation

class MPROPlines():
    """ This class provides functions to edit the mprops file for the purpose 
        of outputting tabulated values of the van Genuchten function or to remove the van Genuchten function definitions
        in order to read them from input tables. It also assembles the tabulated values written by Grok in the form 
        that is required to read them as file inputs. """

    def __init__ (self,grok_dirc, grok_name, mprop_dirc, mprop_name):
        self.grok_dirc = grok_dirc
        self.mprop_path = os.path.join(mprop_dirc,mprop_name)
        self.grok_name = grok_name
        self.pks_folder = 'pks_table'
        self.pks_path = os.path.join(self.grok_dirc, self.pks_folder)
        
        
    def get_index_addstr(self, add_str, ref_str):
        """ get the line index of ref_str """
        # assume "add_str" need to be added before the "ref_str"
        # this function get the index of the "ref_str", and index can be used by other function to add in "add_str"
        # if "add_str" already exist "self.rd" lines before the "ref_str", the index of that "ref_str" will be ignored since
        # no str need to be added

        # self.end_func_index :    index of the ref_str
        # self.end_func_index_2 :  index of the add_str which already exist in the mprop file

    
        self.end_func_index = []
        
        self.end_func_index_2 = []
        
        add_str_num=-999
     
        with open(self.mprop_path) as fhand:
        
            for i, line in enumerate(fhand):
                
                # add_str record the position of the str that is going to be added
                # if the str is too close to the reference str, no str need to be added
                # therefore, not index is recorded
                if add_str in line:
                    add_str_num=i
                    self.end_func_index_2.append(i)
                    #
                
                # find the line that contains the reference string
                # filter out the sections that already has the target string
                if ref_str in line and add_str_num < i-self.rd: 
                    self.end_func_index.append(i)
                        
                        
    def gen_tables(self):
        ''' wrapper for get_index_addstr with default commands to generate table '''
        ref_str = 'end'
        add_str = 'generate tables from unsaturated functions'
        self.get_index_addstr(add_str, ref_str)
        
    
#     def write_mprop(self):
#         """ based on the index in end_func_index, add string to mprop"""
#         
#         with open(self.mprop_path) as fhand, open('gen_table.mprops','w') as fcopy:
#                 for i, line in enumerate(fhand):
#                     if i in self.end_func_index:
#                         fcopy.write('residual saturation \n')
#                         fcopy.write('0.13 \n\n')
#                         fcopy.write('minimum relative permeability \n')
#                         fcopy.write('1.072e-12 \n\n')
#                         fcopy.write('table smoothness factor \n')
#                         fcopy.write('1e-4 \n\n')
#                         fcopy.write('table minimum pressure \n')
#                         fcopy.write('-1.0e3 \n\n')
#                         fcopy.write('generate tables from unsaturated functions \n')
#                     
#                     if i in self.end_func_index_2:
#                         fcopy.write('residual saturation \n')
#                         fcopy.write('0.13 \n\n')
#                         fcopy.write('minimum relative permeability \n')
#                         fcopy.write('1.072e-12 \n\n')
#                         fcopy.write('table smoothness factor \n')
#                         fcopy.write('1e-4 \n\n')
#                         fcopy.write('table minimum pressure \n')
#                         fcopy.write('-1.0e3 \n\n')
#                         
#                     fcopy.write(line)                        
       
    
    def get_gen_table(self,):
        ''' get the code section that controls table generation ''' 
        code = ( ' ! Parameters to generate van Genuchten table output'
                 ' residual saturation \n'
                 ' 0.13 \n\n'
                 ' minimum relative permeability \n'
                 ' 1.072e-12 \n\n'
                 ' table smoothness factor \n'
                 ' 1e-4 \n\n'
                 ' table minimum pressure \n'
                 ' -1.0e3 \n\n'
                 ' generate tables from unsaturated functions \n\n' )
        return code
#         code = '''residual saturation \n 0.13 \n\n minimum relative permeability \n 1.072e-12 \n\n 
#                   table smoothness factor \n 1e-4 \n\n table minimum pressure \n -1.0e3 \n\n 
#                   generate tables from unsaturated functions \n'''

    
    def combine_pkstable(self, ldebug=False):
        """ take the p_s and Kr_s table and combine them together"""
        
        combine_pkstable_psfiles = glob.glob( '{}/{}o.p_s_table.*.dat'.format(self.grok_dirc, self.grok_name))
        if ldebug:
            print('')
            print('{}/{}.p_s_table.*.dat'.format(self.grok_dirc, self.grok_name))
            print(combine_pkstable_psfiles)
        
        if len(combine_pkstable_psfiles) == 0:
            print("\nERROR: No p_s and Kr_s tables found in folder '{}'. Did you run Grok?".format(self.grok_dirc))
            return
        
        # file name of the combined table
        pksfiles = [w.replace('p_s_table', 'p_k_s_table') for w in combine_pkstable_psfiles]
        # path for the pks tables
        pksfiles = list(map(os.path.basename,pksfiles))
        pksfiles = list(map(lambda x:os.path.join(self.grok_dirc, self.pks_folder,x),pksfiles))
        if ldebug: 
          print(''); print(pksfiles)
        
        # the combined table is saved under directory 'pks_table'
        if not os.path.exists( self.pks_path):
            os.makedirs(self.pks_path)
        else:
            warnings.warn('folder ' + self.pks_path + ' already exist')
            
        for i, f in enumerate(combine_pkstable_psfiles):

            #corresponding Kr_S table
            ksfiles = f.replace('p_s_table','s_k_table')
            
            with open(f) as psf, open(pksfiles[i],'w') as fcopy, open(ksfiles) as ksf:
                
                #above 14 are the headers for TECPLOT
                for i, line in enumerate(psf,1):
                    if i >=14:
                        #remove # from the table
                        if line.startswith('#'):
                            line = line.replace('#','')
                            
                        fcopy.write(line)
                        
                for i, line in enumerate(ksf,1):
                    if i >=14:
                        #remove # from the table
                        if line.startswith('#'):
                            line = line.replace('#','')
                            
                        fcopy.write(line)
    
#     def comment_out_unsat_func(self):
#         
# 
#         start_func = 'unsaturated van genuchten functions'
#         end_func = 'end ! functions'
# 
#         #couf_pks_table = self.list_pstable(file_pattern = 'p_k_s_table', file_path=os.path.join(self.grok_dirc,'pks_table'))
#         
#         j = 0 #
#         
#         with open(self.mprop_path) as fhand, open('comment_out_unsat_func.mprops','w') as fcopy:
#             for i,line in enumerate(fhand):
#                 
#                 if start_func in line:
#                     while True:
#                         fcopy.write('!' + line)
#                         
#                         if end_func in line:
#                             fcopy.write('\n include ./pks_table/' + self.grok_name + '.p_k_s_table.' + self.zn[j] + '.dat' + '\n')
#                             j=j+1
#                             break
#                         
#                         line = fhand.readline()
#                 else: 
#                     fcopy.write(line)
#                         
#                     
#             
# 
#     def get_zone_name(self):
#             """ get the line index at where HGS should be added"""
#         
#             self.zn = []
#             
#             ref_str = 'end ! material'
#             
#             a= '!'
#          
#             with open(self.mprop_path) as fhand:
#             
#                 for line in fhand:
#                     
#                     # add_str record the position of the str that is going to be added
#                     # if the str is too close to the reference str, no thr need to be added
#                     # therefore, not index is recorded
#                     if ref_str in line:
#                         line = fhand.readline()
#                         while True:
#                             if a in line or line in [' \n', ' \r\n', '\n']:
# 
#                                 line = fhand.readline()
#                                 
#                             else:
#                                 self.zn.append(line.strip())
#                                 break
#                 
#             self.zn = ['1'] + self.zn


    def insert_pks(self, material_name):
        ''' construct string for command to read values from table (incl. table file path) '''        
        pks_path =  os.path.join('.', self.pks_folder, '{}o.p_k_s_table.{}.dat'.format(self.grok_name, material_name)) 
        return ' include ' + pks_path


    def walk_mprop(self, lgentab=False):
            """ walk through mprops file, identify van Genuchten function definitions, and edit file """
            
            material_flag = False
            unsat_fun_flag = False
            unsat_fun_cmd = 'unsaturated van genuchten functions'
    
         
            # make backup copy 
            # N.B.: this should be persistent, i.e. should be the actual original for all subsequent executions
            backup_copy = self.mprop_path + '.preunsat_original'
            if not os.path.exists(backup_copy): 
                shutil.copy(self.mprop_path, backup_copy)
            # remove original and write new file (creat a one-time backup)
            shutil.move(self.mprop_path, self.mprop_path + '.backup')            
         
            # read file from persistent original and write to location of original mprops file
            with open(backup_copy) as fhand, open(self.mprop_path,'w') as fcopy:
                
                # loop over lines in persistent original
                line = fhand.readline() 
                
                # N.B.: readline preserves the trailing newline; an empty string indicates EOF
                while line:
                    
                    if not line.strip().startswith('!') and line.strip():

                        material_name = line.strip()  
                        material_flag = True # enter material block mode
                        lhasunsat = False # flag indicating if the material has a unsat function definition

                        # Material Mode
                        while material_flag:

                            # unsat/van Genuchten function mode
                            if unsat_fun_cmd in line.lower():

                                unsat_fun_flag = True # enter van Genuchten function block mode
                                lhasunsat = True # indicate that this material has a custom van Genuchten function definition 

                                while unsat_fun_flag:

                                    # warn from duplicates
                                    if 'generate tables from unsaturated functions' in line.lower() and lgentab:
                                        warnings.warn("Potentially duplicate table generation command detected!")
                                    
                                    # turn off unsat func mode
                                    if 'end' in line.lower():
                                        unsat_fun_flag = False # exit van Genuchten function block mode
                                    
                                    # decide what to do with van Genuchten definition
                                    if lgentab:
                                        # if in table-generation mode, include table directives just before 'end'
                                        if 'end' in line.lower(): 
                                            fcopy.write(self.get_gen_table())
                                        fcopy.write(line) # and print original line, too                                        
                                    else:
                                        # if in include-table mode, just comment out
                                        fcopy.write('! '+line)
                                    
                                    line = fhand.readline()
                                        
                            # turn off material mode
                            elif 'end' in line.lower():
                              
                                # end of the material
                                material_flag = False # exit material block mode
                                
                                # if in include-table mode, insert include command                                
                                if lgentab:
                                  
                                    # if the material does not have a unsaturated function definition, add
                                    # default values, so that a default table is output
                                    if not lhasunsat:
                                         
                                        # add unsat function block header
                                        tmp_cmd = ' {} ! added to output default tables \n\n'.format(unsat_fun_cmd)
                                        fcopy.write(tmp_cmd)
                                        fcopy.write(self.get_gen_table())
                                        fcopy.write(' end ! function \n\n')
                                  
                                else:
                                  
                                    # add the command to load tabulated van Genuchten values from a file
                                    fcopy.write(self.insert_pks(material_name)+'\n\n')
                                
                                # write out end of the material line
                                fcopy.write(line)
                                fcopy.write('\n')

                            else:
                              
                                # write out material commands
                                fcopy.write(line)
                                line = fhand.readline()

                    else:
                        fcopy.write(line)
                    
                    # advance loop
                    line = fhand.readline()