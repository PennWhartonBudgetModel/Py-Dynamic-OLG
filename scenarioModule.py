##
# Scenario definition for dynamic model execution.
#

import numpy as np
import warnings
import os
import scipy.io as sio
import csv
import shutil

class Scenario:
    
    # Identifier tags
    basedeftag = None         # level 1 in dir structure
    counterdeftag = None      # level 2 in dir structure
    transitiontag = None      # level 3 in dir structure
    comparisontag = None      # REM: This is for isEquivalent()
    nonversioncomparisontag = None  # This is for isEquivalentExceptVersion()
    
    # REQUIRED parameters
    
    # Economy assumption: steady-state or transition
    
    IsSteadyState = None
    
    # Preference parameters
    beta = None
    gamma = None
    sigma = None
    bequest_phi_1 = None
    
    # Conversion to match US dollar amounts
    modelunit_dollar = None
    
    # Core economy parameters
    IsLowReturn = None               # Add risk-premium reduction to MPK return
    AllowBusinessDebt = None         # Whether business debt is on in Firm class
    LeverageSensitivity = None        # Parameter in Firm leverage cost function
    CapitalAdjustmentCost = None      # Capital Adjustment Cost (eta)
    
    # International 
    OpennessPath = None               # string ID of capital and debt takeups by foreigners
    
    # Timing
    TransitionFirstYear = None        # Actual year to start transition path
    TransitionLastYear = None         # Year transition path ends
    ClosureYear = None                # Year to fix D/Y onward
    PolicyShockYear = None            # Year unanticipated policy shock occurs
    
    # Simulation Options
    UseStaticDebt = None              # Whether using static projections of debt interest rates
    UseNewDemographics = None         # Use the original module for demographics
    
    # OPTIONAL policy parameters
        
    # Immigration parameters
    prem_legal = None
    prem_illegal = None
    amnesty = None
    # TBD: manage immigrantToCitizenRatio;   
    
    # Tax parameters
    TaxCode = None
    
    # Government expenditures
    OutlaysPolicy = None
    
    # Microsim parameters
    Microsim = None
    
    # Demographics parameters
    Demographics = None
        
    # Social Security parameters
    OASIPolicy = None
        
    # Out of model experiment
    OutOfModel = None
        
    # Labor shock process
    LaborShock = None
    
    # OPTIONAL calibration target params
    LaborElasticity = None
    Description = None
    
    # Default VERSIONS for input interfaces
    VersionOpenness = None
    VersionTaxCalculator = None
    VersionOASIcalculator = None
    VersionMicrosim = None
    VersionEntryExitRates = None
    VersionProjections = None
    VersionOutOfModel = None
    VersionDemographics = None
    VersionLaborShock = None
    
    #Property
    ConstructedDescription = None
    
    # Define list of required initial state parameters
    initial_params = ['beta','gamma','sigma','bequest_phi_1',
                            'modelunit_dollar','IsLowReturn',
                            'AllowBusinessDebt','CapitalAdjustmentCost',
                            'LeverageSensitivity','TransitionFirstYear', 'UseNewDemographics']
    
    # Define list of required transition path parameters
    transition_params = ['IsSteadyState','OpennessPath',
                               'TransitionLastYear','ClosureYear',
                               'PolicyShockYear','UseStaticDebt'] 
    
    # Specify default values for optional parameters
    policy_params = {'prem_legal': 1, 'prem_illegal': 0.62043,
                     'amnesty': 0,'TaxCode': 'Baseline',
                  'OutlaysPolicy': 'CurrentPolicy','Microsim': 'baseline',
                  'Demographics'              : 'baseline'            , 
                  'OASIPolicy'                : 'baseline'            , 
                  'OutOfModel'                : 'baseline'            , 
                  'LaborShock'                : 'Kent9' }
    
    # Specify default input versions
    version_params = {'VersionOpenness' : '2019-05-09-16-49-efraim-76b302a', 
            'VersionTaxCalculator'      : '2019-05-15-12-00-efraim-LatestBaseline'  ,  
            'VersionOASIcalculator'     : '2019-05-29-00-00-efraim-transition'      ,  
            'VersionMicrosim'           : '2019-05-20-04-17-danielav_6daf726'       ,  
            'VersionEntryExitRates'     : '2019-01-31-07-19-njanetos-0c2b55d'       ,  
            'VersionProjections'        : '2019-03-15-00-00-jricco-LatestProjectionsBaseline', 
            'VersionOutOfModel'         : '2019-05-03-13-49-jhuntley-7524391',         
            'VersionDemographics'       : '2019-05-10-15-00-jhuntley-6f1d22',          
            'VersionLaborShock'         : '2019-05-03-13-49-jhuntley-7524391'}

    # Parameters on which to match in Map files
    mapfile_params = [
            'OpennessPath',  
            'TaxCode',       
            'OutlaysPolicy', 
            'OutOfModel',    
            'Demographics',   
            'LaborShock',    
            'OASIPolicy']
                          
    # Constructor
    def __init__(self, params):
        
        assert isinstance(params, dict), 'Scenario constructor expects dictionary of parameter values.'
        
        # Set certain required missing params to defaults
        if not 'IsSteadyState' in params.keys():
            params['IsSteadyState'] = 0
        if not 'LaborShock' in params.keys():
            params['LaborShock'] = Scenario.policy_params['LaborShock']
            
        # If required parameters are missing, 
        # find them from calibrator (which looks at other params)
        check = [1 if not v in params.keys() else 0 for v in Scenario.initial_params]

        if any(check):
            print('[INFO] Scenario core parameters missing. Fetching from Calibrator.\n' )
            from paramGeneratorModule import ParamGenerator
            x = ParamGenerator.invert(params)
            params['beta']             = x['beta']
            params['gamma']            = x['gamma']
            params['sigma']            = x['sigma']
            params['modelunit_dollar'] = x['modelunit_dollar']
            params['bequest_phi_1']    = x['bequest_phi_1']
            
            
        # Assign fields to Scenario, warn if any extras
        for k in params.keys():
            if hasattr(self, k):
                setattr(self, k, params[k]) 
            else:
                warnings.warn('Field <%s> does not match any Scenario fields.' % k)
                
        # Check for required initial parameters
        # TBD: leave this here or take it out from needCalibrate check?
        for i in Scenario.initial_params: 
            assert i in params.keys() and params[i] != None, 'Scenario constructor requires nonempty <%s> parameter.' % i
            
        # Check for required transition path parameters
        for i in Scenario.transition_params:
            assert i in params.keys() and params[i] != None, 'Scenario constructor requires nonempty <%s> parameter.' % i
            
        # Set optional policy parameters defaults where unspecified
        for i in Scenario.policy_params:
            if not i in params.keys():
                setattr(self, i, Scenario.policy_params[i])
                
        # Construct Scenario Description if it is not present
        if self.Description == None:
            if self.isCurrentPolicy():
                policy = 'CurrentPolicy'
            else:
                policy = 'Counterfactual'
            self.Description = self.OpennessPath + '-' + policy
            if self.IsSteadyState:
                self.Description = 'Steady-state'
            self.ConstructedDescription = 1
        else:
            self.ConstructedDescription = 0
            # TBD: Check validity of Description, since sometimes used
            # as filename.
        
        # Set version defaults where unspecified
        for i in Scenario.version_params:
            if not i in params.keys():
                setattr(self, i, Scenario.version_params[i])
                
        # Fix timing inconsistencies, if any
        self.ClosureYear = min(max(self.ClosureYear, self.TransitionFirstYear), self.TransitionLastYear)
        self.PolicyShockYear = min(max(self.PolicyShockYear, self.TransitionFirstYear), self.TransitionLastYear)
        
        # Generate identifier tags for baseline and counterfactual definitions
        #   1. Make string of concatenated params
        #   2. Hash the string down to 120 chars
        # NOTE: comparisontag is built for isEquivalent 
        tags = {}
        tags['initial'] = ''
        for i in Scenario.initial_params:
            tags['initial'] += '_' + str(getattr(self,i))
        
        tags['initialExVersion'] = tags['initial']
        
        for i in Scenario.version_params.keys():
            tags['initial'] += '_' + getattr(self,i)
        
        self.basedeftag = Scenario.compactifyTag(tags['initial'])
        
        tags['policy'] = ''
        for i in Scenario.policy_params.keys():
            tags['policy'] += '_' + str(getattr(self,i))
        
        if self.isCurrentPolicy():
            self.counterdeftag = 'currentpolicy'
        else:
            self.counterdeftag = Scenario.compactifyTag(tags['policy'])
        
        tags['transition'] = ''
        if self.IsSteadyState:
            tags['transition'] = 'steady'
            self.transitiontag = 'steady'
        else:
            for i in Scenario.transition_params:
                tags['transition'] += '_' + str(getattr(self,i))
            self.transitiontag = Scenario.compactifyTag(tags['transition'])
            
        self.comparisontag              = tags['initial'] + tags['policy'] + tags['transition']
        self.nonversioncomparisontag    = tags['initialExVersion'] + tags['policy'] + tags['transition']
        
    ##
    # Identify if scenario is equivalent to another scenario
    # Parameter representations in tags determine precision for equivalency evaluation
    def isEquivalent(self, scenario):
        flag = (self.comparisontag == scenario.comparisontag)
        return flag

    ##
    # Identify if scenario is equal to another scenario
    #   on everything except versions.
    def isEquivalentIgnoreVersion(self, scenario):
        flag = (self.nonversioncomparisontag == scenario.nonversioncomparisontag)
        return flag
    
    ##
    # Identify if scenario represents current policy
    def isCurrentPolicy(self):
        #For clarity, current policy should not have a "non-shock"
        if self.isPolicyShock():
            flag = False
            return flag
        
        # Current policy identified by default values for all optional parameters
        for i in Scenario.policy_params.keys():
            if getattr(self, i) != Scenario.policy_params[i]:
                flag = False
                return flag
        
        flag = True
        return flag
    
    ##
    # Identify if scenario represents steady state
    def isSteady(self):
        flag = self.IsSteadyState
        return flag
    
    ##
    # Identify if scenario represents open economy
    def isOpen(self):
        flag = not self.IsSteadyState and self.OpennessPath == 'open'
        return flag
    
    ##
    # Identify if scenario represents closed economy
    def isClosed(self):
        flag = not self.IsSteadyState and self.OpennessPath == 'closed'
        return flag
    
    ##
    # Identify if scenario represents no-policy shock economy
    def isPolicyShock(self):
        flag = self.PolicyShockYear > self.TransitionFirstYear
        return flag
    
    ##
    # Check if the Scenario has been solved and stored to files
    def isSolved(self):
        from pathFinderModule import PathFinder
        flag = os.path.exists(os.path.join(PathFinder.getCacheDir(self), 'solved'))
        if flag :
            # Check that hashed location is correct scenario
            s = sio.loadmat(os.path.join(PathFinder.getCacheDir(self), 'scenario.mat' ))
            flag = self.isEquivalent(s['scenario'])
            if not flag :
                # TBD: Until cache-ing is revised.
                raise Exception('Scenario:BAD_CACHEING - WARNING! Cached scenario at location is not this scenario. ')
        return flag
    
    # Get params for matching in Map files
    def getMapFileMatchParams(self):
        params = {}
        for o in Scenario.mapfile_params:
            params[o] = getattr(self,o)
        return params
    
    ##
    # Get all scenario parameters
    def getParams(self):
        params = {}
        for f in vars(self):
            # Do not return 'tag' properties. 
            # TBD: Do this more elegantly.
            if not f.endswith('tag'):
                params[f] = getattr(self, f)
        return params
    
    # Human-readable description
    def shortDescription(self):
            
        from paramGeneratorModule import ParamGenerator   
        T_model = ParamGenerator.timing(self)['T_model']
        desc = (('[ %s ]' % self.Description) + 
            ('\n \t%-25s= %u' % ('T_model', T_model)) +
            ('\n \t%-25s= %u' % ('IsLowReturn', self.IsLowReturn)) +
            ('\n \t%-25s= %7.8f' % ('Beta', self.beta)) +
            ('\n \t%-25s= %7.8f' % ('Gamma', self.gamma)) + 
            ('\n \t%-25s= %7.8f' % ('Sigma', self.sigma)) +
            ('\n \t%-25s= %e' % ('Model$', self.modelunit_dollar)))
        
        return desc
    
    # Full list of all params in a single text line
    def longDescription(self):
        params = self.getParams()
        desc = 'Scenario: '
        for o in params.keys():
            vName = o
            vVal = str(getattr(self,vName))
            desc = desc + vName + '=' + vVal + ';'
        return desc
    
    # Generate corresponding current policy scenario
    def currentPolicy(self):
        params = self.getParams()
        for f in Scenario.policy_params.keys():
            params[f] = Scenario.policy_params[f]
            
        # Make non-shock
        params['PolicyShockYear'] = self.TransitionFirstYear
        if self.ConstructedDescription :
            params.pop('Description')
        scenario = Scenario(params)
        return scenario
    
    # Generate corresponding steady state scenario
    def steady(self):
        params = self.getParams()
        params['IsSteadyState'] = 1
        if self.ConstructedDescription :
            params.pop('Description')
        scenario = Scenario(params)
        return scenario
    
    # Generate corresponding open economy scenario
    def open(self):
        params = self.getParams()
        params['IsSteadyState'] = 0
        params['OpennessPath'] = 'open'
        if self.ConstructedDescription :
            params.pop('Description')
        scenario = Scenario(params)
        return scenario
    
    # Generate corresponding closed economy scenario
    def closed(self):
        params = self.getParams()
        params['IsSteadyState'] = 0
        params['OpennessPath'] = 'closed'
        if self.ConstructedDescription :
            params.pop('Description')
        scenario = Scenario(params)
        return scenario
    
    # Generate corresponding 'baseline' scenario
    def baseline(self):
        params = self.getParams()
        params['IsSteadyState'] = 0
        params['OpennessPath'] = 'baseline'
        if self.ConstructedDescription :
            params.pop('Description')
        scenario = Scenario(params).currentPolicy()
        return scenario
    
    # Generate corresponding post-shock policy continuation scenario
    def postShock(self):
        params = self.getParams()
        params['TransitionFirstYear'] = params['PolicyShockYear']
        if self.ConstructedDescription :
            params.pop('Description')
        scenario = Scenario(params)
        return scenario
    
    ##
    #       Writes an already-solved scenario's optimal decision rules
    #       and productivity transitions to file.
    @staticmethod
    def writeTransitionMatrix(scenario):
        
        # load solution objects
        from pathFinderModule import PathFinder
        cacheDir = PathFinder.getCacheDir(scenario)
        OPTs = sio.loadmat(os.path.join(cacheDir, 'decisions.mat'))
        
        # get the base output directory
        baseOutputDir = PathFinder.getTransitionMatrixOutputDir()
        
        # create output folder if it does not exist
        if not os.path.exists(baseOutputDir):
            os.path.mkdir(baseOutputDir)

        # get the tagged subfolder output directory
        outputDir = os.path.join(baseOutputDir, PathFinder.getScenarioPathTag(scenario))
        
        # check for whether scenario output subfolder exists
        # if it does, then this is a duplicate writing out
        if os.path.exists(outputDir):
            return None
        
        # check if map file exists, create it if it does not
        if not os.path.exists(os.path.join(baseOutputDir, 'map.csv')):
            fileHandle = open(os.path.join(baseOutputDir, 'map.csv'), 'w')
            for k in scenario:
                fileHandle.write(k + ',')
            fileHandle.write('\n')
            fileHandle.close()

        # append scenario info to map file by writing out to text file
        # then loading text file back in
        with open('.temp.txt', 'w') as f:
            values = scenario.getParams()
            w = csv.DictWriter(f, values.keys())
            w.writerow(values)
        f = open('.temp.txt','r')
        text = f.read()
        f.close()
        os.path.remove('.temp.txt')
        fileHandle = open(os.path.join(baseOutputDir, 'map.csv'), 'a+')
        print(fileHandle, scenario.basedeftag + ',' + scenario.counterdeftag + ',' + text)
        fileHandle.close()
         
        # create a folder to store output
        os.path.mkdir(outputDir)
        
        # converts policy function into discretized transition matrix
        # if policy doesn't fall neatly into grid, averages between two
        # nearest points proportionally to distance from that point
        def convertToTransitionMatrix(policy, values, dim):
            discrete = np.digitize(policy, values)
            distanceToBinEdge = policy - values(discrete)
            distanceToBinEdgeUpper = policy - values(discrete + 1)
            upperProbability = distanceToBinEdge / (distanceToBinEdge - distanceToBinEdgeUpper)
            transition = np.zeros((len(discrete), dim))
            transition[np.ravel_multi_index((np.array(range(grids['nz']*grids['nk']*grids['nb'])), (discrete+1)), transition.shape)] = upperProbability
            transition[np.ravel_multi_index((np.array(range(grids['nz']*grids['nk']*grids['nb'])), discrete), transition.shape)] = 1 - upperProbability
            return transition
            
        # for a given age, year, discretize assets and lifetime earning
        # average transitions. store output in `transitions` variable.
        transitions = {}
        
        # store grids for easy access
        from paramGeneratorModule import ParamGenerator
        grids = ParamGenerator.grids(scenario)
        
        for age in range(OPTs['SAVINGS'].shape[3]):
            for year in range(OPTs['SAVINGS'].shape[4]):

                # compute transition matrices for full state -> assets,
                # earnings grid
                assetsTransition   = convertToTransitionMatrix( 
                    OPTs['SAVINGS'][:, :, :, age, year],           
                    grids['kv'],                                   
                    grids['nk']                                    )
                    
                earningsTransition = convertToTransitionMatrix( 
                        OPTs['AVG_EARNINGS'][:, :, :, age, year],    
                        grids['bv'],                                  
                        grids['nb']                               )
                
                # compute joint transition of assets and earnings
                assetEarningsTransition = (                
                        np.kron(np.ones((1, grids['nb'])), assetsTransition)       
                        * np.kron(earningsTransition, np.ones((1, grids['nk']))))
                
                # expand joint transition of asset and earnings to full
                # state space size
                assetEarningsTransition = np.kron(np.ones((1, grids['nz'])), assetEarningsTransition)
                
                # get the productivity transition matrix
                productivityTransition = grids['transz']
                productivityTransition = np.squeeze(productivityTransition[age, :, :])
                
                # expand it to the full state space size
                productivityTransition = np.kron(                   
                        productivityTransition,                      
                        np.ones(grids['nb'] * grids['nk'], grids['nb'] * grids['nk']))
                    
                # multiply to get full transition matrix
                transitionMatrix = productivityTransition * assetEarningsTransition
                
                # save transition matrix into struct
                transitions['age' + str(age) + 'year' + str(year)] = transitionMatrix
            
        sio.savemat(os.join.path(outputDir, 'data.mat'), transitions)
        
    # Export scenario to named folder
    def export(self, outputName=None):
        # If no outputName, create one from Scenario             
        if outputName == None:
            outputName = self.Description
            
        if not self.isSolved() :
            from modelSolverModule import ModelSolver
            ModelSolver.solve(self)
            
        from pathFinderModule import PathFinder
        cacheDir    = PathFinder.getCacheDir( self )
        outDir      = PathFinder(self).getNamedOutputPath( outputName )
            
        print( 'Exporting scenario to %s \n' % outDir )
        if os.path.exists(outDir):
            shutil.rmtree(outDir)
        shutil.copyfile( cacheDir, outDir )

    @staticmethod
    def compactifyTag( tag ):
            
        TAG_SIZE    = 80
        # Allow chars ASCII 65-90 and 97-122 only
        # Remap ones that fall between into numbers (0-6)
        MIN_CHAR1   = 65
        MAX_CHAR1   = 90
        MIN_CHAR2   = 97
        MAX_CHAR2   = 122
        CHAR0       = 48 
        newtag      = list(tag[0:min(len(tag), TAG_SIZE)])
        for i in range(TAG_SIZE - 1, len(tag)):
            d = (i + 1) % (TAG_SIZE) 
            c = (ord(newtag[d]) + ord(tag[i])) % (MAX_CHAR2 - MIN_CHAR1) + MIN_CHAR1
            if (c > MAX_CHAR1) and (c < MIN_CHAR2) :
                c = c - (MAX_CHAR1 - CHAR0)

            newtag[d] = chr(c)
            
        newtag = "".join(newtag)
        return newtag
    
    @staticmethod
    def exportList( l ):
        for i in range(len(l)):
            scenario = l[i]
            scenario.export()
    