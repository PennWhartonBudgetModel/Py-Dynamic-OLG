# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 13:33:25 2019

@author: Azanca
"""
import numpy as np

#collection of helper functions

def makeIterable(x):
    
    if isinstance(x, (int,float)):
        return np.array([x])
    elif isinstance(x, list):
        return np.array(x)
    elif isinstance(x, np.ndarray):
        return np.copy(x)
    else:
        raise Exception("Input could not be cast to iterable array")
        
#checks if the input is an int/float/array whose Matlab equivalent would have had just 1 column
        
def checkOneColumn(x):
    
    try:
        y = makeIterable(x)
    
        if len(y.shape) == 1:
            if y.shape[0] == 1:
                return True
        elif len(y.shape) > 1:
            if y.shape[1] == 1:
                return True
        else:
            return False
    except:
        return False
    
#replicates Matlab round (rounds .5 up instead of to nearest even like in Python)

import math
def normal_round(n):
    if n - math.floor(n) < 0.5:
        return math.floor(n)
    return math.ceil(n)