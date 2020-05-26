
from scenarioModule import Scenario
from pathFinderModule import PathFinder
from outputWriterModule import OutputWriter
from modelSolverModule import ModelSolver

import numpy as np
import numbers
import pickle
import os
import scipy.io as sio

class ModelTester:
    
    '''
    test_params = {
        'IsSteadyState'         : 0                     ,
           'beta'               : 0.9945000000000           ,
           'gamma'              : 0.7210000000000           ,
           'sigma'              : 1.6000000000000           ,
           'bequest_phi_1'      : 0                         ,
           'modelunit_dollar'   : 3.600000000e-05           ,
        'IsLowReturn'           : 0                     ,
        'TransitionFirstYear'   : 2019                      ,
        'OpennessPath'          : 'closed'                  ,
        'AllowBusinessDebt'     : 1                      ,
        'UseStaticDebt'         : 0                     ,
        'CapitalAdjustmentCost' : 0.0                       ,
        'LeverageSensitivity'   : 3.5                       ,
        'TransitionLastYear'    : (2019+10)                 ,
        'ClosureYear'           : (2019+5)                  ,
        'PolicyShockYear'       : (2019)                    ,
       # POLICY PARAMS 
        'Description'           : '0_0_0_0_0_67_0_2031_0_0_0_0_0_0_3_0_1'           ,
       # VERSIONS
        'VersionOpenness'       : '2019-05-09-16-49-efraim-76b302a'                 ,
        'VersionTaxCalculator'  : '2019-05-15-12-00-efraim-LatestBaseline'          ,
        'VersionOASIcalculator' : '2019-03-18-13-20-ses-2ab2306'                    ,
        'VersionMicrosim'       : '2019-01-08-09-11-njanetos-transition'            ,
        'VersionEntryExitRates' : '2019-01-31-07-19-njanetos-0c2b55d'               ,
        'VersionProjections'    : '2019-03-15-00-00-jricco-LatestProjectionsBaseline' ,   
        'VersionOutOfModel'     : '2019-05-03-13-49-jhuntley-7524391'               }
    '''
    
    #Use 0 and 1 for boolean values instead of False and true. Important for consistency when generating comparisontags later on
    test_params = {
       'IsSteadyState'         : 0,
       'LaborElasticity'       : 0.5,                       
       'LaborShock'            : 'Kent9',                   
       'IsLowReturn'           : 0  ,                   
       'TransitionFirstYear'   : 2019    ,                  
       'OpennessPath'          : 'closed'   ,               
       'AllowBusinessDebt'     : 1       ,               
       'UseStaticDebt'         : 0       ,              
       'CapitalAdjustmentCost' : 0           ,            
       'LeverageSensitivity'   : 3.5          ,             
       'TransitionLastYear'    : (2019+10)      ,           
       'ClosureYear'           : (2019+5)         ,         
       'PolicyShockYear'       : (2019)           ,         
       'UseNewDemographics'    : 1                ,      
       'Demographics'          : 'baseline'           ,    
     # POLICY PARAMS %
       'OASIPolicy'            : '0_0_0_0_0_67_0_2031_0_0_0_0_0_0_3_0_1' ,          
     # VERSIONS
       'VersionOpenness'       : '2019-05-09-16-49-efraim-76b302a'       ,          
       'VersionTaxCalculator'  : '2019-05-15-12-00-efraim-LatestBaseline'   ,      
       'VersionOASIcalculator' : '2019-05-29-00-00-efraim-transition'     ,         
       'VersionMicrosim'       : '2019-05-20-04-17-danielav_6daf726'       ,        
       'VersionEntryExitRates' : '2019-01-31-07-19-njanetos-0c2b55d'       ,        
       'VersionProjections'    : '2019-03-15-00-00-jricco-LatestProjectionsBaseline' ,   
       'VersionOutOfModel'     : '2019-05-03-13-49-jhuntley-7524391'        ,       
       'VersionLaborShock'     : '2019-05-03-13-49-jhuntley-7524391'              }

    policyShockShift = 3    # Num years after which to do unanticipated shock


    # Levels of matching
    DEVIATION_NONE          = 0    # Perfect match
    DEVIATION_TINY          = 1    # Less than 1e-4% deviation or < 1e-12 numerical
    DEVIATION_SMALL         = 1    # Less than 0.01% deviation but more than TINY
    DEVIATION_FATAL         = 2    # More than 0.01% deviation, or missing field, or extra field

 
    # Test steady state solution and elasticities
    @staticmethod
    def steady():
        scenario = Scenario(ModelTester.test_params).currentPolicy().steady()
        testName = 'steady'
        ModelTester.testOutput(scenario, testName, True)    
    
    # Test open economy baseline solution, dynamic aggregates, and static aggregates
    @staticmethod
    def open_base():
        scenario = Scenario(ModelTester.test_params).currentPolicy().open()
        testName = 'open_base'
        ModelTester.testOutput(scenario, testName, True)
    
    # Test open economy counterfactual solution, dynamic aggregates, and static aggregates
    @staticmethod
    def open_counter():
        scenario = Scenario(ModelTester.test_params).open()
        testName = 'open_counter'
        ModelTester.testOutput(scenario, testName, True)   
    
    # Test closed economy baseline solution, dynamic aggregates, and static aggregates
    @staticmethod
    def closed_base():
        scenario = Scenario(ModelTester.test_params).currentPolicy().closed()
        testName = 'closed_base'
        ModelTester.testOutput(scenario, testName, True)
    
    # Test closed economy counterfactual solution, dynamic aggregates, and static aggregates
    @staticmethod
    def closed_counter():
        scenario = Scenario(ModelTester.test_params).closed()
        testName = 'closed_counter'
        ModelTester.testOutput( scenario, testName, True )
    
    # Test all basic cases, for when expecting no change
    def all_cases():
        
        testNames   = ['steady', 'open_base', 'open_counter', 'closed_base', 'closed_counter']
        for o in testNames:
            if o == 'steady':
                scenario = Scenario(ModelTester.test_params).currentPolicy().steady()
            elif o == 'open_base':
                scenario = Scenario(ModelTester.test_params).currentPolicy().open()
            elif o == 'open_counter':
                scenario = Scenario(ModelTester.test_params).open()
            elif o == 'closed_base':
                scenario = Scenario(ModelTester.test_params).currentPolicy().closed()
            elif o == 'closed_counter':
                scenario = Scenario(ModelTester.test_params).closed()
            else:
                scenario = []
            
            if (ModelTester.testOutput( scenario, o, 0 ) != ModelTester.DEVIATION_NONE):
                return
    
    
    # Test unanticipated shock. Non-shock should be within tolerances
    @staticmethod
    def unanticipated_shock():
        
        # Make the baseline scenario and "non-shock" version
        t                   = ModelTester.test_params
        
        # baseline scenario is not shocked
        s_baseline          = Scenario(t).currentPolicy().baseline()
        
        # Make "non-shock" shock baseline
        t                   = s_baseline.getParams()
        t.PolicyShockYear   = t.TransitionFirstYear + ModelTester.policyShockShift
        s_next              = Scenario(t)

        # Get baseline Market, Dynamic
        ModelSolver.removeCached(s_baseline)                 # Clear cached Scenario
        
        tagged_dir      = ModelSolver.solve(s_baseline)
        baseline_dir    = PathFinder.getCacheDir(s_baseline)
        with open(os.path.join(baseline_dir, 'market.pkl'), 'rb') as handle:
            baseMarket      = pickle.load(handle)
        with open(os.path.join(baseline_dir, 'dynamics.pkl'), 'rb') as handle:
            baseDynamic     = pickle.load(handle)   
        
        # Get shocked Market, Dynamic
        ModelSolver.removeCached(s_next)                     # Clear cached scenario
        
        tagged_dir      = ModelSolver.solve(s_next)
        x_dir           = PathFinder.getCacheDir(s_next)
        with open(os.path.join(x_dir, 'market.pkl'), 'rb') as handle:
            xMarket         = pickle.load(handle)
        with open(os.path.join(x_dir, 'dynamics.pkl'), 'rb') as handle:
            xDynamic        = pickle.load(handle)
        
        # Compare baseline and shocked path
        print( '\n' )
        
        def do_check (baseD, xD, dName):
            passed = 1
            for p in baseD.keys():
                valuename = p
                if (not isinstance(baseD[valuename], numbers.Number) or ('_next' in valuename)):
                    continue

                # Check for within percent tolerance, also check 
                #    within numerical deviation (this is in case div by
                #    zero or close to zero)
                # TBD: Standardize deviations and tolerances
                percentDeviation    = abs((xD[valuename] - baseD[valuename]) / baseD[valuename])
                absoluteDeviation   = abs(baseD[valuename] - xD[valuename])
                if not np.all(np.array(percentDeviation) < 1e-4):
                    if not np.all(np.array(absoluteDeviation) < 1e-13):
                        m1 = print( 'Max percentdev = %f' % max(percentDeviation) )
                        m2 = print( 'Max abs dev = %0.14f' % max(absoluteDeviation) )
                        print( '%s.%s outside tolerance;\t\t %s; %s \n' % (dName, valuename, m1, m2))
                        passed = 0
                
            return passed
        
        passed = do_check( baseMarket , xMarket , 'Market'  )
        passed = do_check( baseDynamic, xDynamic, 'Dynamic' )
        if passed:
            print( 'All values within convergence tolerances.\n' )
        
        return passed
    
    
    ## Run Test suite on Jenkins server
    @staticmethod
    def jenkinsTests():
        
        try:
            isHPCC      = PathFinder.isHPCCRun()

            # Run just the matching cases for now
            testNames   = ['steady', 'open_base', 'open_counter', 'closed_base', 'closed_counter']
            for o in testNames:
                if o == 'steady':
                    scenario = Scenario(ModelTester.test_params).currentPolicy().steady()
                elif o == 'open_base':
                    scenario = Scenario(ModelTester.test_params).currentPolicy().open()
                elif o == 'open_counter':
                    scenario = Scenario(ModelTester.test_params).open()
                elif o == 'closed_base':
                    scenario = Scenario(ModelTester.test_params).currentPolicy().closed()
                elif o == 'closed_counter':
                    scenario = Scenario(ModelTester.test_params).closed()
                else:
                    scenario = []

                typeDeviation = ModelTester.testOutput( scenario, o, 0 )

                if typeDeviation != ModelTester.DEVIATION_NONE:
                    if typeDeviation == ModelTester.DEVIATION_TINY and isHPCC:
                        continue
                    else:
                        exit(1)

            # Test writing the 'series' interface with the last scenario
            # Requires that 'baseline' scenario exists
            PathFinder.setToTestingMode()
            print( 'TESTING OutputWriter.writeScenarios\n' )
            ModelSolver.solve( scenario.baseline() )
            OutputWriter.writeScenarios( [scenario] )
            PathFinder.setToDevelopmentMode()

            print( 'ALL TESTS PASSED.\n' )
            exit(0)
        except:
            exit(1)
    
    # Test solver output against target values
    @staticmethod
    def testOutput(scenario, testName, isInteractive):

        # Set to testing environment
        PathFinder.setToTestingMode()
        
        # Clear the old results and solve
        ModelSolver.removeCached(scenario)
        taggedDir = ModelSolver.solve(scenario)
        cacheDir  = PathFinder.getCacheDir(scenario)
        
        # Set to development environment 
        #   TBD: Set back to original environment?
        PathFinder.setToDevelopmentMode()

        # testSet depends on type of scenario
        if( scenario.isSteady() ):
            setNames = ['market', 'dynamics']
        elif( scenario.isCurrentPolicy() ):
            setNames = ['market', 'dynamics' ]
        else:
            setNames = ['market', 'dynamics', 'statics']
        
        # Load target values
        targetfile = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ModelTester.pkl')
        with open(targetfile, 'rb') as handle:
            s = pickle.load(handle)
        target = s.target

        # Initialize match flag
        typeDeviation = ModelTester.DEVIATION_NONE

        # Define function to flag issues
        # NOTE: Relies on severity of deviation to be increasing
        def flag(str, deviation):
            print('\t%-15s%-20s%s\n' % (setname, valuename, str))
            global typeDeviation
            if deviation > typeDeviation:
                typeDeviation = deviation        

        print('\n[Test results]\n')
        for i in range(len(setNames)):

            # Extract output and target values by set
            setname = setNames[i]
            output = {}
            with open(os.path.join(cacheDir, ('%s.pkl' % setname)), 'rb') as handle:
                output[testName][setname] = pickle.load(handle)
            outputset = output[testName][setname]
            targetset = target[testName][setname]

            # Iterate over target values
            targetvaluenames = targetset.keys()

            for j in range(len(targetvaluenames)):

                valuename = targetvaluenames[j]

                if not valuename in outputset.keys():

                    # Flag missing value
                    flag('Not found', ModelTester.DEVIATION_FATAL)
                    continue

                if isinstance(outputset[valuename], dict):

                    # Skip checking of structs -- it is currently just
                    # priceindex which does not need to be checked
                    print('\tSkipping %s because it is a struct.\n' % valuename)
                    continue

                if np.any(np.isnan(outputset[valuename][:])):

                    # Flag NaN value
                    flag('NaN value', ModelTester.DEVIATION_FATAL)
                    continue

                if np.any(outputset[valuename].shape != targetset[valuename].shape):

                    # Flag for size mismatch
                    flag('Size mismatch', ModelTester.DEVIATION_FATAL)
                    continue

                # Classify deviation
                deviation = ModelTester.calculateDeviation(outputset[valuename][:], targetset[valuename][:])
                if deviation > 0:
                    if (deviation < 1e-6): 
                        msg = 'TINY : %06.16f%% deviation' % deviation*100
                        flag(msg, ModelTester.DEVIATION_TINY)
                    elif deviation < 1e-4:
                        msg = 'SMALL: %06.16f%% deviation' % deviation*100
                        flag( msg, ModelTester.DEVIATION_SMALL )
                    else:
                        msg = 'LARGE: %06.4f%% deviation' % deviation*100
                        flag( msg, ModelTester.DEVIATION_FATAL )

            # Identify new values, if any
            outputvaluenames = outputset.keys()

            for j in range(len(outputvaluenames)):

                valuename = outputvaluenames[j]

                if not valuename in targetset.keys():
                    flag('New', ModelTester.DEVIATION_FATAL)

        # Check for match
        if typeDeviation == ModelTester.DEVIATION_NONE:
            print('\tTarget matched.\n\n')
        else:

            if not isInteractive: 
                print( '\tTarget not matched.\n\n' )
                return
            
            # Query user for target update
            ans = input('\n\tUpdate test target with new values? Y/[N]: ')
            if ans == 'Y':
                target[testName] = output[testName]
                with open(targetfile) as f:
                    pickle.dump(target, f)
                print('\tTarget updated.\n\n')
            else:
                print('\tTarget retained.\n\n')

        return typeDeviation

        
    ##
    #  Calculation of deviation (as scalar) in a vector
    #      Inputs: output and target vectors
    def calculateDeviation(test, target):
        # First check for equality
        if not np.any(test - target):
            deviation = 0
            return deviation
        
        smallCutoff         = 1e-10
        smalls              = (abs(target) < smallCutoff )
        notSmalls           = not smalls
        
        # Function for calculating moving average
        def movmean(x, N):
            cumsum = np.cumsum(np.insert(x, 0, 0)) 
            return (cumsum[N:] - cumsum[:-N]) / float(N)

        # Check deviations for not small elements as percentage
        if np.any(notSmalls):
            # Find deviations using moving average AND
            #  element-by-element. Take smallest deviation.
            MAdeviation      = abs(movmean(test[notSmalls], 5) / movmean(target[notSmalls], 5) - 1)
            elementDeviation = abs(test[notSmalls]             / target[notSmalls]             - 1)
            percentDeviation = max( min(MAdeviation, elementDeviation) )
        else:
            percentDeviation = 0
        
        # Check small elements and apply custom deviation
        if( any(smalls) ):
            
            smallDeviation  = 0
            testSmalls      = abs(test[smalls]) < smallCutoff
            
            # For similar sized elements which are small,
            # set any deviations to be a set small size
            if( any(testSmalls) ):
                if( test[testSmalls] != target[testSmalls]):
                    smallDeviation = smallCutoff
        
            # For any test elements which are not small,
            # apply a percent deviation off of smallCutoff
            if( any(testSmalls) ):
                smallDeviation = max( abs((test[testSmalls] / smallCutoff) - 1) )
            
        else: # no smalls
            smallDeviation = 0
        
        deviation = max(percentDeviation, smallDeviation)
        
        return deviation

    
