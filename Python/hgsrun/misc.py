'''
Created on Dec 8, 2017

A module that contains a bunch of helper functions for the HGS wrapper class; 
mainly functions that parse configuration and output files.

@author: Andre R. Erler, GPL v3
'''

import os, subprocess

## general utility 

# WindowsError is not defined on Linux - need a dummy
try: 
    lWin = True
    WindowsError
except NameError:
    lWin = False
    WindowsError = None
        
## patch symlink on Windows
# adapted from Stack Overflow (last answer: "Edit 2"):
# https://stackoverflow.com/questions/6260149/os-symlink-support-in-windows
if os.name == "nt":
  def symlink_ms(source, link_name):
    import ctypes
    csl = ctypes.windll.kernel32.CreateSymbolicLinkW
    csl.argtypes = (ctypes.c_wchar_p, ctypes.c_wchar_p, ctypes.c_uint32)
    csl.restype = ctypes.c_ubyte
    flags = 1 if os.path.isdir(source) else 0
    if csl(link_name, source.replace('/', '\\'), flags) == 0:
        raise ctypes.WinError()
    if not os.path.exists(link_name):
        raise WindowsError("Creation of Symbolic Link '{}' failed - do you have sufficient Privileges?".format(link_name))
  
  # replace os symlink with this function
  os.symlink = symlink_ms
  
# some named exceptions
class GrokError(Exception):
  ''' Exception indicating an Error in Grok '''
  pass
class HGSError(Exception):
  ''' Exception indicating an Error in HGS '''
  pass

# a platform-independent solution to clear a folder...
def clearFolder(folder, lWin=lWin, lmkdir=True):
    ''' create a folder; if it already exists, remove it and all files and subfolder in it, and recreate it '''
    if os.path.exists(folder):
        if lWin: cmd = 'rmdir /q /s {}'.format(folder.replace('/', '\\'))
        else: cmd = ['rm','-r',folder]
        # N.B.: don't use shutil.rmtree(folder)! It appears that on Windows the use of soft links causes rmtree to 
        #       follow the links and also delete all source files!
        p = subprocess.Popen(cmd, shell=lWin, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # on windows we also need the shell, because rmdir is not an executable in itself
        stdout, stderr = p.communicate(); ec = p.poll(); del stdout
        if ec > 0: 
            raise IOError("The command '{}' to remove the existing directory failed; exit code {}\n{}".format(cmd,ec,stderr))
    # crease folder
    if lmkdir: os.mkdir(folder)


## file parsere

# helper function to recursively parse a Grok file with includes
def parseGrokFile(filepath, line_list):
    ''' recursively parse files which include other files, starting with 'filepath', 
        and append their lines to 'line_list' '''
    # change working directory locally to resolve relative path'
    pwd = os.getcwd() # need to go back later!
    os.chdir(os.path.dirname(filepath))
    with open(os.path.basename(filepath), 'r') as f:
        lines = f.readlines() # load lines all at once (for performance)
    # loop over lines, clean them and find includes
    for line in lines:
        line = line.strip() # remove which space 
        # skip comments and empty lines
        if not line.startswith('!'):
            line = line.replace('\\','/') # enforce Unix folder convention
            # scan for include statements (except precip.inc and pet.inc) 
            if line.startswith('include') and not ( line.endswith('precip.inc') or line.endswith('pet.inc') ):
                # N.B.: apparently HGS/Grok is case-sensitive... 
                # figure out new file path
                incpath = line[7:].strip()
                if not os.path.lexists(incpath):
                    raise IOError("Unable to open include file '{}'.".format(incpath))
                # initiate recursion
                parseGrokFile(incpath, line_list)
            else:
                # append line to line_list
                line_list.append(line)
                # N.B.: only convert non-path statements to lower
    # now we are done - we don't return anything, since line_list was modified in place
    os.chdir(pwd) # this is important: return to previous working directory!!!


if __name__ == '__main__':
    pass