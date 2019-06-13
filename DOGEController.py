##
# DOGEController is the interface to run the PWBM OLG model
#  from scripted requests -- eventually from Luigi type calls
#

import os as path
import pandas as pd
import Scenario
import ModelSolver
import OutputWriter

class DOGEController:
	 
	# Solve a list of scenarios
	#    Inputs:
	#      workListFileName : full path to CSV file with Scenarios
	#      index (optional) : line number in the worklist to solve
	#                               ignoring the other items
	def solve(workListFileName, index):
		
		# Get name of worklist from name of file
		workListName = path.splitext(workListFileName)[0]
		scenarios = DOGEController.getScenarios(workListFileName)
		   
		# If just running one, collapse list to one
		isSingleItem = 'index' in globals() and len(index) != 0
		if isSingleItem:
			scenarios   = scenarios(index)
		numScenarios    = len(scenarios)
		
		for i in range(numScenarios):
			scenario = scenarios[i]
			
			if isSingleItem:
				idx = index
			else:
				idx = i
				
			print('DOGEController.solve: Solving worklist item %d \n' % idx)
			
			# Tag this solution to avoid collisions with other model runs
			callerTag = '%s%u' % (workListName, idx)
			taggedDir = ModelSolver.solve(scenario, callerTag)
        
	# Export scenarios in worklist. 
	# These should already be cached or it will warn.
	# Inputs:
	# workListFileName : full path to CSV file with Scenarios
	def export(workListFileName):
		
		# Get name of worklist from name of file
		workListName = path.splitext(workListFileName)[0]
		
		scenarios = DOGEController.getScenarios(workListFileName)
		
		numWritten = OutputWriter.writeScenarios(scenarios)
		print( 'DOGEController.export: Wrote %u scenarios in worklist "%s".\n' % (numWritten, workListName) )

		# TBD: Add report on terminations	
	
	## Make Scenarios list from worklist. 
	#  Inputs:
	#  workListFileName : full path to CSV file with Scenarios
	def getScenarios(workListFileName):
		
		# Get name of worklist from name of file
		workListName = path.splitext(workListFileName)[0]

		# Force read CSV -- readtable gets confused and skips first
		# rows sometimes
		workList        = (pd.read_csv(workListFileName)).to_dict(orient='dict')
		numScenarios    = len(workList)
		
		print('DOGEController.getScenarios: Worklist <%s> size = %d \n' % (workListName, numScenarios))
		
		scenarios = Scenario.Scenario(workList)
		return scenarios