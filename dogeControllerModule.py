##
# DOGEController is the interface to run the PWBM OLG model
#  from scripted requests -- eventually from Luigi type calls
#

from modelSolverModule import ModelSolver
from scenarioModule import Scenario
from outputWriterModule import OutputWriter

import os
import pandas as pd
import numpy as np

class DOGEController:
     
    # Solve a list of scenarios
    #    Inputs:
    #      workListFileName : full path to CSV file with Scenarios
    #      index (optional) : line number in the worklist to solve
    #                               ignoring the other items
    def solve(workListFileName, index = None):
        
        # Get name of worklist from name of file
        workListName = os.path.split(workListFileName)[1]
        workListName = os.path.splitext(workListName)[0]
        scenarios = DOGEController.getScenarios(workListFileName)
           
        # If just running one, collapse list to one
        isSingleItem = False
        numScenarios = len(scenarios)
        if index != None:
            isSingleItem = True
            numScenarios = 1
        
        for i in range(numScenarios):
            
            idx = i
            if isSingleItem:
                idx = index
                
            print('DOGEController.solve: Solving worklist item %d \n' % idx)
            
            # Tag this solution to avoid collisions with other model runs
            callerTag = '%s%u' % (workListName, idx)
            taggedDir = ModelSolver.solve(scenarios[idx], callerTag)
        
    # Export scenarios in worklist. 
    # These should already be cached or it will warn.
    # Inputs:
    # workListFileName : full path to CSV file with Scenarios
    def export(workListFileName):
        
        # Get name of worklist from name of file
        workListName = os.path.split(workListFileName)[1]
        workListName = os.path.splitext(workListName)[0]
        scenarios = DOGEController.getScenarios(workListFileName)
        
        numWritten = OutputWriter.writeScenarios(scenarios)
        print( 'DOGEController.export: Wrote %u scenarios in worklist "%s".\n' % (numWritten, workListName) )

        # TBD: Add report on terminations    
    
    ## Make Scenarios list from worklist. 
    #  Inputs:
    #  workListFileName : full path to CSV file with Scenarios
    def getScenarios(workListFileName):
        
        # Get name of worklist from name of file
        workListName = os.path.split(workListFileName)[1]
        workListName = os.path.splitext(workListName)[0]

        # Force read CSV -- readtable gets confused and skips first
        # rows sometimes
        workList        = (pd.read_csv(workListFileName, header = None)).to_dict(orient='list')
        numScenarios    = len(workList)
        
        print('DOGEController.getScenarios: Worklist <%s> size = %d \n' % (workListName, numScenarios))
        
        scenarios = np.array([])
        for i in range(numScenarios):
            scenarios.append(Scenario(workList[i]))
        return scenarios