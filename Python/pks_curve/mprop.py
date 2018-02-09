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

    def __init__ (self,grok_dirc, grok_name, mprop_dirc, mprop_name, sr, kr_min, tsf, p_min, ldebug = False):
        self.grok_dirc  = grok_dirc
        self.mprop_path = os.path.join(mprop_dirc,mprop_name)
        self.grok_name  = grok_name
        self.pks_folder = 'pks_table'
        self.pks_path   = os.path.join(self.grok_dirc, self.pks_folder)
        self.sr         = sr
        self.kr_min     = kr_min
        self.tsf        = tsf
        self.p_min      = p_min

        if ldebug:

            print (" van Genuchten parameters \n ")
            print({"sr":sr, "kr_min":kr_min, "tsf":tsf, "p_min":p_min})

    def get_gen_table(self, sr, kr_min, tsf, p_min):
        ''' get the code section that controls table generation ''' 

        code = ( ' ! Parameters to generate van Genuchten table output \n\n'
                 ' residual saturation \n' 
                 '{} \n\n'.format(self.sr) + 
                 ' minimum relative permeability \n' 
                 '{} \n\n'.format(self.kr_min) + 
                 ' table smoothness factor \n'
                 '{} \n\n'.format(self.tsf) +
                 ' table minimum pressure \n'
                 '{} \n\n'.format(self.p_min) +
                 ' generate tables from unsaturated functions \n\n' )

        return code
#         code = '''residual saturation \n 0.13 \n\n minimum relative permeability \n 1.072e-12 \n\n 
#                   table smoothness factor \n 1e-4 \n\n table minimum pressure \n -1.0e3 \n\n 
#                   generate tables from unsaturated functions \n'''

    
    def combine_pkstable(self, ldebug=False):
        """ take the p_s and Kr_s table and combine them together"""
        
        psfiles = glob.glob( '{}/{}o.p_s_table.*.dat'.format(self.grok_dirc, self.grok_name))
        if ldebug:
            print('')
            print('{}/{}o.p_s_table.*.dat'.format(self.grok_dirc, self.grok_name))
            print(psfiles)
        
        if len(psfiles) == 0:
            print("\nERROR: No p_s and Kr_s tables found in folder '{}'. Did you run Grok?".format(self.grok_dirc))
            return
        
        # file name of the combined table
        pksfiles = [w.replace('p_s_table', 'p_k_s_table') for w in psfiles]
        # path for the pks tables
        pksfiles = list(map(os.path.basename,pksfiles))
        pksfiles = list(map(lambda x:os.path.join(self.grok_dirc, self.pks_folder,x),pksfiles))

        #get list of files for the KrS table
        ksfile = [w.replace('p_s_table', 's_k_table') for w in psfiles]

        if ldebug: 
          print(''); print(pksfiles)
        
        # the combined table is saved under directory 'pks_table'
        if not os.path.exists( self.pks_path):
            os.makedirs(self.pks_path)
        else:
            warnings.warn('folder ' + self.pks_path + ' already exist')
            
        # loop over pressure-saturation, relative permeability, and combined tables (f, g, h)
        for f, g, h in zip(psfiles, ksfile, pksfiles):
            
            with open(f) as psf, open(h,'w') as pksf, open(g) as ksf:
                
                #above 14 are the headers for TECPLOT
                for i, line in enumerate(psf,1):
                    if i >=14:
                        #remove # from the table
                        if line.startswith('#'):
                            line = line.replace('#','')
                            if line.strip().lower() != 'unsaturated tables' and i==14:
                                raise ValueError( 'unrecognized table header: {} \n {} '.format(line, f))
                            if line.strip().lower() != 'pressure-saturation' and i==15:
                                raise ValueError( 'unrecognized table header: {} \n {} '.format(line, f))
                            
                        pksf.write(line)
                        
                for i, line in enumerate(ksf,1):
                    if i >=14:
                        #remove # from the table
                        if line.startswith('#'):
                            line = line.replace('#','')
                            if line.strip().lower() != 'saturation-relative k' and i==14:
                                raise ValueError( 'unrecognized table header: {} \n {} '.format(line, g))
                            
                        pksf.write(line)


    def insert_pks(self, material_name):
        ''' construct string for command to read values from table (incl. table file path) '''      
        material_name = material_name.lower() # this is necessary, because Grok makes all output filenames lowercase
        # N.B.: I'm not sure if Grok also makes the problem prefix lowercase or not...
        pks_path =  os.path.join('.', self.pks_folder, '{}o.p_k_s_table.{}.dat'.format(self.grok_name, material_name)) 
        return ' include ' + pks_path


    def walk_mprop(self, lgentab=False, ldebug=False):
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
                                            fcopy.write(self.get_gen_table(sr = self.sr, 
                                                                           kr_min = self.kr_min, 
                                                                           tsf = self.tsf, 
                                                                           p_min = self.p_min))
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
                                        fcopy.write(unsat_fun_cmd+' \n\n')
                                        fcopy.write(self.get_gen_table(sr = self.sr, 
                                                                        kr_min = self.kr_min, 
                                                                        tsf = self.tsf, 
                                                                        p_min = self.p_min))
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