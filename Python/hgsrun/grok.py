'''
Created on Jul 31, 2016

A Python class that provides access to elements of the grok configuration file.

@author: Andre R. Erler, GPL v3
'''

class Grok(object):
  '''
    A class that provides access to elements of the grok configuration file.
  '''


  def __init__(self, params):
    '''
      Read the grok configuration file or create a new file if not present.
    '''
        