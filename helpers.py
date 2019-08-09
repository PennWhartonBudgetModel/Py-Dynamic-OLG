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
        return x
    else:
        raise Exception("Input could not be cast to iterable array")