# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 09:58:42 2019

@author: Azanca
"""

#Run Solver with test parameters

from scenarioModule import Scenario
from modelTesterModule import ModelTester
from modelSolverModule import ModelSolver

t = ModelTester.test_params
s = Scenario(t)

ModelSolver.solve(s)