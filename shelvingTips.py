# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 10:34:14 2019

@author: Azanca
"""

# This code can be used for saving and loading workspace variables at different points in the program
# Has to be run in the IPython console because it doesn't load local variables properly when run from .py file (for whatever reason).

import shelve

#use code below to save all local variables at a given point in the program to a file

shelf = shelve.open('localsave2.out', 'n') 

for k in list(locals()):
    try:
        shelf[k] = locals()[k]
    except:
        print('ERROR shelving: {0}'.format(k))

shelf.close()

#use code below to load variables from file as locals

shelf = shelve.open('localsave.out')

for k in shelf:
    print(k)
    globals()[k] = shelf[k]
    
shelf.close()