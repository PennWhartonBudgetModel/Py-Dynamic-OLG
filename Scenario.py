##
# Scenario definition for dynamic model execution.
#
import numpy as np
import pandas as pd
import warnings
from os import path
import scipy.io as sio
import csv

class Scenario:
	
	# Identifier tags
	basedeftag         # level 1 in dir structure
	counterdeftag      # level 2 in dir structure
	transitiontag      # level 3 in dir structure
	comparisontag      # REM: This is for isEquivalent()
	nonversioncomparisontag  # This is for isEquivalentExceptVersion()
	
	# REQUIRED parameters
	
	# Economy assumption: steady-state or transition
	
	IsSteadyState
	
	# Preference parameters
	beta
	gamma
	sigma
	bequest_phi_1
	
	# Conversion to match US dollar amounts
	modelunit_dollar
	
	# Core economy parameters
	IsLowReturn               # Add risk-premium reduction to MPK return
	AllowBusinessDebt         # Whether business debt is on in Firm class
	LeverageSensitivity        # Parameter in Firm leverage cost function
	CapitalAdjustmentCost      # Capital Adjustment Cost (eta)
	
	# International 
	OpennessPath               # string ID of capital and debt takeups by foreigners
    
	# Timing
	TransitionFirstYear        # Actual year to start transition path
	TransitionLastYear         # Year transition path ends
	ClosureYear                # Year to fix D/Y onward
	PolicyShockYear            # Year unanticipated policy shock occurs
	UseStaticDebt              # Whether using static projections of debt interest rates
	
	# OPTIONAL policy parameters
        
	# Immigration parameters
	prem_legal
	amnesty
        
	# Tax parameters
	TaxCode
	
	# Government expenditures
	OutlaysPolicy
	
	# Microsim parameters
	Microsim
	
	# Social Security parameters
	Description
	
	# Out of model experiment
	OutOfModel
	
	# OPTIONAL calibration target params
	LaborElasticity
	
	# Default VERSIONS for input interfaces
	VersionOpenness
	VersionTaxCalculator
	VersionOASIcalculator
	VersionMicrosim
	VersionEntryExitRates
	VersionProjections
	VersionOutOfModel
	
	# Define list of required initial state parameters
	initial_params = ['beta','gamma','sigma','bequest_phi_1',
							'modelunit_dollar','IsLowReturn',
							'AllowBusinessDebt','CapitalAdjustmentCost',
							'LeverageSensitivity','TransitionFirstYear']
	
	# Define list of required transition path parameters
	transition_params = ['IsSteadyState','OpennessPath',
							   'TransitionLastYear','ClosureYear',
							   'PolicyShockYear','UseStaticDebt'] 
	
	# Define list of calibration target params
	calibration_params = ['LaborElasticity','IsLowReturn']
	
	# Specify default values for optional parameters
	policy_params = {'prem_legal': 1,'amnesty': 0,'TaxCode': 'Baseline',
				  'OutlaysPolicy': 'CurrentPolicy','Microsim': 'baseline',
				  'Description': 'baseline','OutOfModel': 'baseline'}
	
	# Specify default input versions
	version_params = {'VersionOpenness': '2019-05-09-16-49-efraim-76b302a', 
				   'VersionTaxCalculator': '2019-05-15-12-00-efraim-LatestBaseline',
				   'VersionOASIcalculator': '2019-03-18-13-20-ses-2ab2306',
				   'VersionMicrosim': '2019-03-21-03-24-jricco-GaleRuns',
				   'VersionEntryExitRates': '2019-01-31-07-19-njanetos-0c2b55d',
				   'VersionProjections': '2019-03-15-00-00-jricco-LatestProjectionsBaseline',
				   'VersionOutOfModel': '2019-05-03-13-49-jhuntley-7524391'}

	# Constructor
	def __init__(self, params):
		
		assert (isinstance(params, dict),
		 'Scenario constructor expects dictionary of parameter values.')
		
		# If calibration parameters are missing, find them from calibrator 
		needCalibrate = true
		for i in range(len(Scenario.calibration_params)):
			o = Scenario.calibration_params[i]
			condition1 = not params.has_key(o)
			condition2 = params[o] == NULL
			if condition1 or condition2:
				needCalibrate = false
			else:
				self[o] = params[o]
		
		if needCalibrate:
			x = ParamGenerator.invert(params)
			params['beta']             = x['beta']
			params['gamma']            = x['gamma']
			params['sigma']            = x['sigma']
			params['modelunit_dollar'] = x['modelunit_dollar']
			params['bequest_phi_1']    = x['bequest_phi_1']
			
		# If IsSteadyState is missing, assume it is false
		if not params.has_key('IsSteadyState'):
			params['IsSteadyState'] = false
			
		# Assign fields to Scenario, warn if any extras
		for k in params:
			if k in self:
				self[k] = params[k] #how can it have fields if we are just constractucting self?
			else:
				warnings.warn('Field <%s> does not match any Scenario fields.' % k)
				
		# Check for required initial parameters
		for i in Scenario.initial_params: 
			assert(params.has_key(i) and params[i] != NULL, 'Scenario constructor requires nonempty <%s> parameter.' % i)
			
		# Check for required transition path parameters
		for i in Scenario.transition_params:
			assert(params.has_key(i) and params[i] != NULL, 'Scenario constructor requires nonempty <%s> parameter.' % i)
            
		# Set optional policy parameters defaults where unspecified
		for i in Scenario.policy_params:
			if not params.has_key(i):
				self[i] = Scenario.policy_params[i]
				
		# Set version defaults where unspecified
		for i in Scenario.version_params:
			if not params.has_key(i):
				self[i] = Scenario.version_params[i]
				
		# Fix timing inconsistencies, if any
		self.ClosureYear = min(max(self.ClosureYear, self.TransitionFirstYear), self.TransitionLastYear)
		self.PolicyShockYear = min(max(self.PolicyShockYear, self.TransitionFirstYear), self.TransitionLastYear)
		
		# Generate identifier tags for baseline and counterfactual definitions
		#   1. Make string of concatenated params
		#   2. Hash the string down to 120 chars
		# NOTE: comparisontag is built for isEquivalent 
		tag = ''
		for i in Scenario.initial_params:
			tag += '_' + str(self.i)
		
		tagExVersions = tag
		for i in Scenario.version_params:
			tag += '_' + self.i
		
		self.basedeftag = Scenario.compactifyTag(tag)
		
		if self.isCurrentPolicy():
			self.counterdeftag = 'currentpolicy'
		else:
			tag = ''
			for i in Scenario.policy_params:
				tag += '_' + str(self.i)
			self.counterdeftag = Scenario.compactifyTag(tag)
			
		if self.IsSteadyState:
			self.transitiontag = 'steady'
		else:
			tag = ''
			for i in Scenario.transition_params:
				tag += '_' + str(i)
			self.transitiontag = Scenario.compactifyTag(tag)
			
		self.comparisontag = self.basedeftag + self.counterdeftag + self.transitiontag
		
		self.nonversioncomparisontag = tagExVersions + self.counterdeftag + self.transitiontag
		
		return self
		# Scenario constructor
		
	##
	# Identify if scenario is equivalent to another scenario
	# Parameter representations in tags determine precision for equivalency evaluation
	def isEquivalent(self, scenario):
		flag = self.comparisontag == scenario.comparisontag
		return flag

	##
	# Identify if scenario is equal to another scenario
	#   on everything except versions.
	def isEquivalentIgnoreVersion(self, scenario):
		flag = self.nonversioncomparisontag == scenario.nonversioncomparisontag
		return flag
	
	##
	# Identify if scenario represents current policy
	def isCurrentPolicy(self):
		#For clarity, current policy should not have a "non-shock"
		if self.isPolicyShock():
			flag = false
			return flag
		
		# Current policy identified by default values for all optional parameters
		for i in Scenario.policy_params:
			if self[i] != Scenario.policy_params[i]:
				flag = false
				return flag
		
		flag = true
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
		flag = path.exists(path.join(PathFinder.getCacheDir(self), 'solved'))
		return flag
	
	##
	# Get all scenario parameters
	def getParams(self):
		params = {}
		for f in vars(this):
			# Do not return 'tag' properties. 
			# TBD: Do this more elegantly.
			if not f.endswith('tag'):
				params[f] = self[f]
		return params
	
	# Human-readable description
	def shortDescription(self):
		if self.IsSteadyState:
			sType = 'Steady-state'
		else:
			sType = self.OpennessPath
		if self.isCurrentPolicy():
			policy = 'current policy'
		else:
			policy = 'counterfactual'
            
		T_model = ParamGenerator.timing(self).T_model
		desc = print('[%s - %s]' % (sType, policy))
		desc = print('%s\n \t%-25s= %u' % (desc, 'T_model', T_model))
		desc = print('%s\n \t%-25s= %u' % (desc, 'IsLowReturn', self.IsLowReturn))
		desc = print('%s\n \t%-25s= %7.8f' % (desc, 'Beta', self.beta))
		desc = print('%s\n \t%-25s= %7.8f' % (desc, 'Gamma', self.gamma)) 
		desc = print('%s\n \t%-25s= %7.8f' % (desc, 'Sigma', self.sigma))
		desc = print('%s\n \t%-25s= %e' % (desc, 'Model$', self.modelunit_dollar))
		
		return desc
	
	# Generate corresponding current policy scenario
	def currentPolicy(self):
		params = self.getParams()
		for f in Scenario.policy_params:
			params[f] = Scenario.policy_params[f]
			
		# Make non-shock
		params.PolicyShockYear = self.TransitionFirstYear
		scenario = Scenario(params)
		return scenario
	
	# Generate corresponding steady state scenario
	def steady(self):
		params = self.getParams()
		params['IsSteadyState'] = true
		scenario = Scenario(params)
		return scenario
	
	# Generate corresponding open economy scenario
	def open(self):
		params = self.getParams()
		params['IsSteadyState'] = false
		params['OpennessPath'] = 'open'
		scenario = Scenario(params)
		return scenario
	
	# Generate corresponding closed economy scenario
	def closed(self):
		params = self.getParams()
		params['IsSteadyState'] = false
		params['OpennessPath'] = 'closed'
		scenario = Scenario(params)
		return scenario
	
	# Generate corresponding 'baseline' scenario
	def baseline(self):
		params = self.getParams()
		params['IsSteadyState'] = false
		params['OpennessPath'] = 'baseline'
		scenario = Scenario(params).currentPolicy()
		return scenario
	
	# Generate corresponding post-shock policy continuation scenario
	def postShock(self):
		params = self.getParams()
		params['TransitionFirstYear'] = params['PolicyShockYear']
		scenario = Scenario(params)
		return scenario
	
	##
	#       Writes an already-solved scenario's optimal decision rules
	#       and productivity transitions to file.
	def writeTransitionMatrix(scenario):
		
		# load solution objects
		cacheDir = PathFinder.getCacheDir(scenario)
		globals().update(
				sio.loadmat(path.join(cacheDir, 'decisions.mat'), appendmat = False, variable_names = ['OPTs']))
		
		# get the base output directory
		baseOutputDir = PathFinder.getTransitionMatrixOutputDir()
		
		# create output folder if it does not exist
		if not path.exists(baseOutputDir):
			path.mkdir(baseOutputDir)

		# get the tagged subfolder output directory
		outputDir = path.join(baseOutputDir, PathFinder.getScenarioPathTag(scenario))
		
		# check for whether scenario output subfolder exists
		# if it does, then this is a duplicate writing out
		if path.exists(outputDir):
			return None
		
		# check if map file exists, create it if it does not
		if not path.exists(path.join(baseOutputDir, 'map.csv')):
			fileHandle = open(path.join(baseOutputDir, 'map.csv'), 'w')
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
		path.remove('.temp.txt')
		fileHandle = open(path.join(baseOutputDir, 'map.csv'), 'a+')
		print(fileHandle, scenario.basedeftag + ',' + scenario.counterdeftag + ',' + text)
		close(fileHandle)
		 
		# create a folder to store output
		path.mkdir(outputDir)
		
		# converts policy function into discretized transition matrix
		# if policy doesn't fall neatly into grid, averages between two
		# nearest points proportionally to distance from that point
		def convertToTransitionMatrix(policy, values, dim):
			#? discrete = discretize(policy, values)
			#? discrete = discrete(:)
			#? distanceToBinEdge = policy(:) - values(discrete)
			#? distanceToBinEdgeUpper = policy(:) - values(discrete + 1)
			#? upperProbability = distanceToBinEdge ./ (distanceToBinEdge - distanceToBinEdgeUpper)
			#? transition = zeros(size(discrete, 1), dim)
			#? transition(sub2ind(size(transition), 1:(grids.nz*grids.nk*grids.nb), (discrete+1)')) = upperProbability
			#? transition(sub2ind(size(transition), 1:(grids.nz*grids.nk*grids.nb), discrete')) = 1 - upperProbability
			return transition
			
		# for a given age, year, discretize assets and lifetime earning
		# average transitions. store output in `transitions` variable.
		transitions = {}
		
		# store grids for easy access
		grids = ParamGenerator.grids(scenario)
		
		for age = 1:size(OPTs.SAVINGS, 4)
                for year = 1:size(OPTs.SAVINGS, 5)

                    % compute transition matrices for full state -> assets,
                    % earnings grid
                    assetsTransition   = convertToTransitionMatrix( ...
                        OPTs.SAVINGS(:, :, :, age, year),           ...
                        grids.kv,                                   ...
                        grids.nk                                    ...
                    );
                    
                    earningsTransition = convertToTransitionMatrix( ...
                        OPTs.AVG_EARNINGS(:, :, :, age, year),      ...
                        grids.bv,                                   ...
                        grids.nb                                    ...
                    );

                    % compute joint transition of assets and earnings
                    assetEarningsTransition =                           ...
                        kron(ones(1, grids.nb), assetsTransition)       ...
                        .* kron(earningsTransition, ones(1, grids.nk));

                    % expand joint transition of asset and earnings to full
                    % state space size
                    assetEarningsTransition = kron(ones(1, grids.nz), assetEarningsTransition);

                    % get the productivity transition matrix
                    productivityTransition = grids.transz;
                    productivityTransition = squeeze(productivityTransition(age, :, :));

                    % expand it to the full state space size
                    productivityTransition = kron(                   ...
                        productivityTransition,                      ...
                        ones(grids.nb*grids.nk, grids.nb*grids.nk)   ...
                    );
                    
                    % multiply to get full transition matrix
                    transitionMatrix = productivityTransition ...
                        .* assetEarningsTransition;
                
                    % save transition matrix into struct
                    transitions = setfield(                                  ...
                        transitions,                                         ...
                        strcat('age', int2str(age), 'year', int2str(year)),  ...
                        sparse(transitionMatrix)                             ...
                    );
                end
            end
            
            save(fullfile(outputDir, strcat('data.mat')), 'transitions');

        end % writeTransitionMatrix
        
    end % instance methods, public
    
    
    
    methods (Static, Access = private)
        
        function [newtag] = compactifyTag( tag )
            
            TAG_SIZE    = 120;
            % Allow chars ASCII 65-90 and 97-122 only
            % Remap ones that fall between into numbers (0-6)
            MIN_CHAR1   = 65;   MAX_CHAR1   = 90;
            MIN_CHAR2   = 97;   MAX_CHAR2   = 122;
            CHAR0       = 48; 
            newtag      = tag( 1:min(length(tag), TAG_SIZE) );
            for i=TAG_SIZE:length(tag)
                d = mod(i, TAG_SIZE) + 1;
                c = mod(newtag(d) + tag(i), MAX_CHAR2 - MIN_CHAR1) + MIN_CHAR1;
                if ( (c > MAX_CHAR1) && (c < MIN_CHAR2) )
                    c = c - (MAX_CHAR1 - CHAR0);
                end
                newtag(d) = c;
            end
            
        end % compactifyTag
        
    end

    
end % Scenario


