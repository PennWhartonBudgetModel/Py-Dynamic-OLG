# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 13:33:04 2019

@author: Azanca
"""

#Test ParamGenerator.labinc_discretization
from modelTesterModule import ModelTester
from scenarioModule import Scenario
from pathFinderModule import PathFinder
from paramGeneratorModule import ParamGenerator
from inputReaderModule import InputReader
from socialSecurityModule import SocialSecurity
from initialGuessModule import InitialGuess

import pandas as pd
import numpy as np
import os
import pickle

t = ModelTester.test_params
scenario = Scenario(t)

scenario = scenario.currentPolicy().steady()

s = ParamGenerator.labinc_discretization(scenario)


'''
#Test ParamGenerator.social_security
from modelTesterModule import ModelTester
from scenarioModule import Scenario
from pathFinderModule import PathFinder
from paramGeneratorModule import ParamGenerator
from inputReaderModule import InputReader
from socialSecurityModule import SocialSecurity
from initialGuessModule import InitialGuess

import pandas as pd
import numpy as np
import os
import pickle

t = ModelTester.test_params
scenario = Scenario(t)

scenario = scenario.currentPolicy().steady()

params = ParamGenerator.social_security(scenario)
'''

'''
#Test some Scenario methods
from modelTesterModule import ModelTester
from scenarioModule import Scenario
from pathFinderModule import PathFinder
from paramGeneratorModule import ParamGenerator
from inputReaderModule import InputReader
from socialSecurityModule import SocialSecurity
from initialGuessModule import InitialGuess

import pandas as pd
import numpy as np
import os
import pickle

t = ModelTester.test_params
scenario = Scenario(t)



scenario = scenario.currentPolicy().steady()

print(vars(scenario))
'''

'''
#Test InitialGuess
from modelTesterModule import ModelTester
from scenarioModule import Scenario
from pathFinderModule import PathFinder
from paramGeneratorModule import ParamGenerator
from inputReaderModule import InputReader
from socialSecurityModule import SocialSecurity
from initialGuessModule import InitialGuess

import pandas as pd
import numpy as np
import os
import pickle

t = ModelTester.test_params
scenario = Scenario(t)

filename = os.path.join(PathFinder.getSourceDir(), 'InitialGuess.pkl')
         
with open(filename, 'rb') as f:
    dict = pickle.load(f)
    
scenario2 = dict[0,0]

print(scenario.isEquivalent(scenario2))
print(scenario.isEquivalentIgnoreVersion(scenario2))
'''

'''
#Test Paramgenerator.social_security()
from modelTesterModule import ModelTester
from scenarioModule import Scenario
from pathFinderModule import PathFinder
from paramGeneratorModule import ParamGenerator
from inputReaderModule import InputReader
from socialSecurityModule import SocialSecurity
from initialGuessDataModule import InitialGuessData

import pandas as pd
import numpy as np

t = ModelTester.test_params
scenario = Scenario(t)

timing                  = ParamGenerator.timing(scenario)
T_model                 = timing['T_model']
T_life                  = timing['T_life']
first_year              = timing['yearFirst']
last_year               = timing['yearLast']
startyears              = timing['startyears']
realage_entry           = timing['realage_entry']

params = ParamGenerator.social_security(scenario)
bv = ParamGenerator.grids(scenario)['bv']
Market = InitialGuessData.market
priceIndices = Market['priceindices']

#theSocialSecurity = SocialSecurity(params, bv, Market['priceindices'], startyears, realage_entry, T_model)

t_1                 = priceIndices['T_modelStart']
t_end               = priceIndices['T_modelEnd']
wage_inflations     = priceIndices['wage_inflations'][t_1-1:t_end]
incomeFloor    = (params['ssincmins'] * wage_inflations)
incomeCeiling  = (params['ssincmaxs'] * wage_inflations)
'''

'''
#Test SocialSecurity.getBenefitsForCohort(...)
from modelTesterModule import ModelTester
from scenarioModule import Scenario
from pathFinderModule import PathFinder
from paramGeneratorModule import ParamGenerator
from inputReaderModule import InputReader
from socialSecurityModule import SocialSecurity
from initialGuessDataModule import InitialGuessData

import pandas as pd
import numpy as np

t = ModelTester.test_params
scenario = Scenario(t)
scenario.IsSteadyState = True

timing                  = ParamGenerator.timing(scenario)
T_model                 = timing['T_model']
T_life                  = timing['T_life']
first_year              = timing['yearFirst']
last_year               = timing['yearLast']
startyears              = timing['startyears']
realage_entry           = timing['realage_entry']

params = ParamGenerator.social_security(scenario)
bv = ParamGenerator.grids(scenario)['bv']
priceindices = InitialGuessData.market['priceindices']
                
socialSecurity = SocialSecurity(params, bv, priceindices, startyears, realage_entry, T_model)

i = -1

if i == -1:
    # Note: retire_year = 1 so that ssbenefits is calculated 
    # for T_model = 1 and indexed for first startyear
    startYear   = socialSecurity.startyears[0]
    retireYear  = max(socialSecurity.T_model - 1, 1)
    brackets    = socialSecurity.params['benefitParamsByCohort'][-1]['brackets']
    rates       = socialSecurity.params['benefitParamsByCohort'][-1]['rates']
else:
    startYear   = socialSecurity.startyears[i]
    retireYear  = socialSecurity.params['retire_years'][i]
    brackets    = socialSecurity.params['benefitParamsByCohort'][i]['brackets']
    rates       = socialSecurity.params['benefitParamsByCohort'][i]['rates']
    
# Fetch index used to grow brackets
wage_index = socialSecurity.priceIndices['wage_inflations']  # full long index (has all cohorts)

# Year cohort turns 62 in current model. 
# Index benefits by wage index of that year.
year62 = startYear + (62 - socialSecurity.realage_entry)

# Cohort based benefits by year
#   REM: Every household is at a grid point, so can
#        calculate benefits directly here (outside of
#        solve_cohort).
#   NOTE: We set ssbenefits to -Inf otherwise -- it should not be
#   used in solve_cohort (this will blow it up just in case)
benefits  = np.ones((socialSecurity.T_model, (socialSecurity.bv).shape[0])) * (- float('Inf'))

#Build cohort-specific bracket cutoffs for each year indexed by
# wage_index of when it turns 62.
w_year62    = year62 + socialSecurity.priceIndices['T_modelStart'] # year in index vector when cohort turns 62
bracket_idx = wage_index[int(w_year62 - 1)] * np.ones((socialSecurity.T_model,1))
adjbrackets = brackets * np.tile(bracket_idx, (1, brackets.shape[1]))

# A possible step to implement here would be to deflate the adjusted
# brackets by another index, e.g. CPI old people,  but currently
# it's only deflated by CPI, which is already done in ParamGenerator

# Build cumulative benefits matrix to aid benefit calculation
adjtotben   = np.cumsum(np.diff(adjbrackets, axis = 1) * rates[:, 0:-1], axis = 1) 
adjtotben   = np.hstack((np.zeros((adjbrackets.shape[0], 1)), adjtotben))  # rem: first column is zero

peek = socialSecurity.bv

# Calculate benefits for each year of retirement until end of time.
for t in range((max(retireYear,1) -  1) , socialSecurity.T_model):
    for ib in range((socialSecurity.bv).shape[0]):
        thebracket = np.where(adjbrackets[t,:] <= socialSecurity.bv[ib])[0][-1]
        benefits[t,ib] = adjtotben[t,thebracket] + rates[t,thebracket]*(socialSecurity.bv[ib] - adjbrackets[t,thebracket])
'''

'''
#Test ParamGenerator.budget(scenario)
from modelTesterModule import ModelTester
from scenarioModule import Scenario
from pathFinderModule import PathFinder
from paramGeneratorModule import ParamGenerator
from inputReaderModule import InputReader
import pandas as pd
import numpy as np

t = ModelTester.test_params
scenario = Scenario(t)
scenario.IsSteadyState = True

timing                  = ParamGenerator.timing( scenario )
T_model                 = timing['T_model']
T_life                  = timing['T_life']
first_year              = timing['yearFirst']
last_year               = timing['yearLast']

pathFinder = PathFinder(scenario)
                
projections_file   = pathFinder.getProjectionsInputPath('Projections')
projections_series = InputReader.read_series(projections_file, 'Year', first_year, first_year)
        
taxcalculator_file      = pathFinder.getTaxCalculatorInputPath( 'Aggregates' )
taxcalculator_series    = InputReader.read_series(taxcalculator_file, 'Year', first_year, first_year )
        
oasicalculator_file     = pathFinder.getOASIcalculatorInputPath( 'aggregates' )
oasicalculator_series   = InputReader.read_series(oasicalculator_file, 'Year', first_year, first_year)

deepHistory_file   = pathFinder.getProjectionsInputPath('DeepHistory')
history_series     = InputReader.read_series(deepHistory_file, 'Year', first_year - T_life + 1, last_year)
'''

'''
#Test InputReader.read_cohort_brackets_rates_indices
from modelTesterModule import ModelTester
from scenarioModule import Scenario
from pathFinderModule import PathFinder
from paramGeneratorModule import ParamGenerator
from inputReaderModule import InputReader
import pandas as pd

t = ModelTester.test_params
scenario = Scenario(t)
scenario.IsSteadyState = True

timing                  = ParamGenerator.timing(scenario)
T_model                 = timing['T_model']
T_life                  = timing['T_life']
first_year              = timing['yearFirst']
last_year               = timing['yearLast']
nstartyears             = len(timing['startyears'])
realage_entry           = timing['realage_entry']

first_birthyear = first_year - (T_life + realage_entry)
last_birthyear  = first_birthyear + nstartyears - 1

filename = r'Z:OASIcalculator\Interfaces\2019-05-29-00-00-efraim-transition\oasicalculator1\WriteSeriesInterface_0_OASIcalculator_oasicalculator_4b46e74754\PIAParameters.csv'

T = pd.read_csv(filename)

birthYear = first_birthyear
T_cohort = T[T.BirthYear == birthYear]

#(brackets, rates, indices) = InputReader.parse_brackets_rates_indices(T_cohort, firstYear, Tyears)

tableData = T_cohort
firstYear = first_year
Tyears = T_model

yearStart = tableData[tableData.Year == firstYear].index[0]

brackets = tableData[tableData.columns[['Bracket' in s for s in tableData.columns.values]]]
brackets = brackets.loc[yearStart:,]
'''

'''
# Test InputReader.read_brackets_rates_indices
from modelTesterModule import ModelTester
from scenarioModule import Scenario
from pathFinderModule import PathFinder
from paramGeneratorModule import ParamGenerator
from inputReaderModule import InputReader
import pandas as pd

t = ModelTester.test_params
scenario = Scenario(t)

timing                  = ParamGenerator.timing(scenario)
T_model                 = timing['T_model']
T_life                  = timing['T_life']
first_year              = timing['yearFirst']
last_year               = timing['yearLast']
nstartyears             = len(timing['startyears'])
realage_entry           = timing['realage_entry']

pathFinder = PathFinder(scenario)
file       = pathFinder.getOASIcalculatorInputPath( 'BenefitParameters' )
(brackets, _, _)   = InputReader.read_brackets_rates_indices(file, first_year, T_model)
'''

'''
#Test InputReader.readDemographics(...)
from modelTesterModule import ModelTester
from scenarioModule import Scenario
from pathFinderModule import PathFinder
from paramGeneratorModule import ParamGenerator
from inputReaderModule import InputReader

import numpy as np
import pandas as pd

t = ModelTester.test_params
scenario = Scenario (t)
scenario.IsSteadyState = 1


pathFinder = PathFinder(scenario)

        
# Initializing the timing variables and the grids 
timing    = ParamGenerator.timing( scenario )

file      = pathFinder.getMicrosimInputPath( 'SurvivalProbability' )
series    = InputReader.read_series( file, 'Age', 1 + timing['realage_entry'] )
survival  = series['SurvivalProbability']
s = {}
s['surv']             = np.tile(np.reshape(survival, (1, timing['T_life'])), [timing['T_model'], 1])
'''

'''
#Test ParamGenerator.demographicsNew(scenario) and InputReader.demographicsNew(...)
from modelTesterModule import ModelTester
from scenarioModule import Scenario
from pathFinderModule import PathFinder
from paramGeneratorModule import ParamGenerator

import numpy as np
import pandas as pd

t = ModelTester.test_params
scenario = Scenario (t)
scenario.IsSteadyState = 1
pathFinder = PathFinder(scenario)
timing = ParamGenerator.timing(scenario)

# Initializing the timing variables and the grids 
timing     = ParamGenerator.timing( scenario )
modelYears = timing['T_model'] + 1
grids      = ParamGenerator.grids(  scenario )
firstYear  = timing['TransitionFirstYear']
if scenario.isSteady():
    firstYear = timing['TransitionFirstYear'] - 1
    
# Set initial distribution of individuals by age-productivity-group
# Note: at some point we could make DISTz time-varying.
newDISTz = grids['DISTz']

file       = pathFinder.getDemographicsInputPath( 'rates' )

microSimStatus = {}
microSimStatus['naturalized']      = 0
microSimStatus['legal']            = 1
microSimStatus['illegal']          = 2
microSimStatus['native']           = 3 
        
statusMap = [['naturalized', 'citizen'],   # Naturalized citizens in Microsim are citizens in the OLG
             ['legal',       'legal'],     # Lawful noncitizen is "legal" in OLG
             ['illegal',     'illegal'],   # Unathorized noncitizen is "illegal" in OLG
             ['native',      'citizen']]   # Native born citizen is "citizen" in the OLG

#demographics = InputReader.read_demographics(file, firstYear, modelYears, timing['realage_entry'], timing['T_life'], grids['g'], microSimStatus, statusMap, scenario)
        
filename = file
T_model = modelYears
realage_entry = timing['realage_entry']
T_life = timing['T_life']
g = grids['g']

T = pd.read_csv(filename)

#print(T)
                
groups    = {}
groupList = list(g.keys())
microsimGroups = {}
microsimGroupList = list(microSimStatus.keys())

varlist   = ['EmigratedThisYear', 'ImmigratedThisYear', 'BornThisYear', 'DiedThisYear', 'Population']

for group in range(len(microsimGroupList)):
    # Extract sub-table for each group
    
    T_group = T.loc[T.LegalStatus == group]
    #print(T_group)
    microsimGroups[microsimGroupList[group]] = {
            'EmigratedThisYear': np.zeros((T_model, T_life)),
            'ImmigratedThisYear': np.zeros((T_model, T_life)),
            'BornThisYear': np.zeros((T_model, T_life)),
            'DiedThisYear': np.zeros((T_model, T_life)),
            'Population': np.zeros((T_model, T_life))
            }
    
    for age in range(realage_entry + 1, realage_entry + T_life + 1):

        # Find age
        if sum(T_group.Age.values == age) == 0:
            raise Exception('Cannot find age %u for group %s in input table.' % (age, group))
            
        # Extract sub-sub-table for each group & age
        T_group_age = T_group.loc[T_group.Age == age]
        
        #print(T_group_age) 
        

        # Find first year & check for missing years
        if not firstYear in T_group_age.Year.values:
            raise Exception('Cannot find first index (%u) for age %u of group %s in input table.' % (firstYear, age, microsimGroupList[group]))
        years = (T_group_age.Year >= firstYear) 
        gap = np.diff(T_group_age.loc[years, 'Year'].values)
        if any(gap>1):
            missing_years = '%d ' % (T_group_age.Year[gap > 1] + 1)
            raise Exception('Missing year(s) %sfor age %u of group %s in input table.' % (missing_years, age, microsimGroupList[group]))
        
        # Truncating the data if microsim is longer than OLG.
        # Throwing an error if OLG is longer than microsim.
        
        for var in varlist:
                    
            temp_array = np.array(T_group_age.loc[years,var].values)
                
            # Pad series if not long enough, truncate if too long
            numYears = len(temp_array)
            if T_model - numYears <= 0:
                temp_array = temp_array[0:T_model]
            else:
                raise Exception('read_demographics:PADDING Microsim does not extend far enough into the future.')
                
            microsimGroups[microsimGroupList[group]][var][:,age-realage_entry-1] = temp_array
'''

'''
#Test ParamGenerator.demographicsNew(scenario) and ParamGenerator.truncateDISTz(...)
from modelTesterModule import ModelTester
from scenarioModule import Scenario
from pathFinderModule import PathFinder
from paramGeneratorModule import ParamGenerator

import numpy as np
import pandas as pd

t = ModelTester.test_params
scenario = Scenario (t)
pathFinder = PathFinder(scenario)
timing = ParamGenerator.timing(scenario)

# Initializing the timing variables and the grids 
timing     = ParamGenerator.timing( scenario )
modelYears = timing['T_model'] + 1
grids      = ParamGenerator.grids(  scenario )
firstYear  = timing['TransitionFirstYear']
if scenario.isSteady():
    firstYear = timing['TransitionFirstYear'] - 1
    
# Set initial distribution of individuals by age-productivity-group
# Note: at some point we could make DISTz time-varying.
newDISTz = grids['DISTz']

#newDISTz[:,:,grids['g']['illegal']] = ParamGenerator.truncateDISTz(scenario.prem_illegal, newDISTz[:,:,grids['g']['citizen']], timing['T_life'], np.squeeze(grids['zs'][0,:,:]) )

premium = scenario.prem_illegal
baseDISTz = newDISTz[:,:,grids['g']['citizen']]
T_life = timing['T_life']
zs = np.squeeze(grids['zs'][0,:,:])

groupDISTz     = baseDISTz           # Initializing the output
dampingFactor  = 0.05                # If the algorithm is not converging to a solution, try smaller values of this to make sure it's not "jumping" over the solution.
shareTruncated = 1                  # One minus this number is what is removed from the top of the distribution; this is the initial guess
        
for age in range(T_life):
            
    # Checking to make sure that this is a non-zero productivity (non-retirement)
    if (sum(np.squeeze(zs[age,:]) * np.squeeze(baseDISTz[:,age])) > .0001):
                
        for iter in range(1000):
                    
            productivityArray = np.sort(np.squeeze(zs[age,:]))
            sortIndex = np.argsort(np.squeeze(zs[age,:]))
                    
            # Initialize the truncated distribution for one age
            truncAgeDISTz = groupDISTz[sortIndex,age]

            # This loop goes through each grid point in that productivity
            # group and reallocates the distribution.
                   
            totProb       = 0
            upperBound    = max(productivityArray.shape)
            prodThreshold = sum(truncAgeDISTz[0:upperBound]) * shareTruncated
            for jz in range(upperBound):
                # Walk through the distribution from the bottom, making sure
                # that we're not at the threshold yet
                if totProb < prodThreshold:
                    totProb = totProb + truncAgeDISTz[jz]
                    # If this puts the probability over the threshold, stop
                    # and only take away enough from the top grid point so
                    # that it hits, but does not exceed the threshold
                    if totProb > prodThreshold:
                        truncAgeDISTz[jz] = truncAgeDISTz[jz] - (totProb - prodThreshold)
                    else:
                        truncAgeDISTz[jz] = 0
                    
            # Reallocate the truncated part proportionally to the remaining
            # part of the distribution.
            truncAgeDISTz = truncAgeDISTz / shareTruncated
            
            # This is the inverse sort: it's taking the density for
            # the sorted productivity and unwinding back to the
            # original way that the productivity shocks were
            # ordered.
            truncAgeDISTz[sortIndex] = truncAgeDISTz
            
            # Update the share to be truncated
            groupRatio     = (sum(np.squeeze(zs[age,:]) * np.squeeze(truncAgeDISTz)) / 
                              sum(np.squeeze(zs[age,:]) * np.squeeze(baseDISTz[:,age])))
            shareTruncated = min(max(shareTruncated + dampingFactor * (premium - groupRatio),0),1)
                
        groupDISTz[:,age] = truncAgeDISTz
        shareTruncated  = 1
        
    # Check to make sure that the ratio is correct and throw a warning,
    # skipping those with NaN (the denominator is 0 b/c there is no
    # productivity for retired people)
    groupRatioCohortError = max(np.abs(premium - np.sum(zs * np.transpose(groupDISTz),axis=1) / np.sum(zs * np.transpose(baseDISTz),axis=1)))
    if np.isnan(groupRatioCohortError):
        groupRatioCohortError = 0
    if np.abs(groupRatioCohortError) > 0.0001:
        print( 'WARNING! Did not hit the target income ratio for a citizenship status group.\n')
'''

'''
#Test InputReader.read_transitions
from scenarioModule import Scenario
from modelTesterModule import ModelTester
from pathFinderModule import PathFinder
from paramGeneratorModule import ParamGenerator
from inputReaderModule import InputReader
import pandas as pd
import numpy as np

t = ModelTester.test_params
s = Scenario(t)

pathFinder = PathFinder(s)

timing = ParamGenerator.timing(s)
T_life    = timing['T_life']
firstYear = timing['yearFirst']
lastYear  = timing['yearLast']

file        = pathFinder.getLaborShockInputPath( 'transitions' )
#transitions = InputReader.read_transitions(file, yearFirst,(yearLast - yearFirst + 1), T_life)

T = pd.read_csv(file)
firstAge    = min(T.Age.values)
lastAge     = max(T.Age.values)
numStates   = np.sum([True if 'Transition_1_' in x else False for x in T.columns.values])

for year in range(firstYear,lastYear):
    
    T_year_rows = (T.Year == 2018)
    
    for initState in range(1,numStates+1):
        
        T_columns = [True if ('Transition_%u_' % initState) in x else False for x in T.columns.values]
        trans = T.values[T_year_rows,:][:, T_columns]
        s = np.sum(trans, axis = 1) - 1
        if any(s > 1e-15) or any(s < -1e-15) :
            print('year' + str(year) + ' and state ' + str(initState))
'''

'''
#Test InputReader.read_initDIST
from scenarioModule import Scenario
from modelTesterModule import ModelTester
from pathFinderModule import PathFinder
from paramGeneratorModule import ParamGenerator
from inputReaderModule import InputReader
import pandas as pd
import numpy as np

t = ModelTester.test_params
s = Scenario(t)

pathFinder = PathFinder(s)

timing = ParamGenerator.timing(s)
T_life    = timing['T_life']
firstYear = timing['yearFirst']
lastYear  = timing['yearLast']

file        = pathFinder.getLaborShockInputPath( 'initial_DIST' )
#DISTz       = InputReader.read_initDIST(file, yearFirst, T_life)

T = pd.read_csv(file)
T_year_rows = (T.Year.values == 2018)
        
# Select sub-table with the appropriate year
T = T.iloc[T_year_rows,:]
        
# Find ranges of ages and states
firstAge    = min(T.Age.values)
lastAge     = max(T.Age.values)
numStates   = sum(['Mass_' in s for s in T.columns.values])
numGroups   = sum(T.Age.values == firstAge)

groupState = 1

rows_group = (T.Group.values == groupState + 1)
dist_group = T.iloc[rows_group,:]
dist_group = dist_group.drop(dist_group.columns.values[0:3], axis = 1)
dist_group = dist_group.values

#print(dist_group)

s = np.sum(dist_group, axis = 1)
dif = np.tile(s[:,np.newaxis], [1, numStates])

print(s.shape)
print(dist_group.shape)
print(dif.shape)
'''


'''
#Test Firm:
from firmTesterModule import FirmTester
FirmTester.doGELoop()
'''