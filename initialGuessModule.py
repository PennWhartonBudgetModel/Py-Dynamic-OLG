# Used to generate model run guesses.

from pathFinderModule import PathFinder
from scenarioModule import Scenario

import os
import numpy as np
import pickle
import scipy.io as sio
    
class InitialGuess:

    scenario = None
    Dynamic = None
    Market = None
    
    def __init__(self, scenario):
        
        self.scenario = scenario
        self.fetch()
    
    # Add guess to DB
    def save(self):
            
        filename = os.path.join(PathFinder.getSourceDir(), 'InitialGuess.pkl')
         
        dict = np.array([[]])
        
        if os.path.isfile(filename):
            with open(filename, 'rb') as f:
                dict = pickle.load(f)
            # dict is an array with a scenario, Dynamic, Market values
            
        # Find exact match (if any), otherwise add
        found = False
        if dict.shape[1] == 3:
            for i in range(dict.shape[0]):
                scenario = dict[i,0]
                if self.scenario.isEquivalent(scenario):
                    dict[i,:] = [scenario, self.Dynamic, self.Market]
                    found = True
                    break
            
        if not found:
            np.append(dict, [self.scenario, self.Dynamic, self.Market])

        with open(filename, 'wb') as f:
            pickle.dump(dict, f, protocol=pickle.HIGHEST_PROTOCOL)

        if found :
            source = 'Overwrote existing guess.'
        else:
            source = 'Added as new guess.'
            
        print('[INFO] Saved initial guess to file. %s \n' % source)
        
    # Get the guess from DB or generate if not there
    def fetch(self):
        filename = os.path.join(PathFinder.getSourceDir(), 'InitialGuess.pkl')
        if not os.path.isfile(filename):
            self.generate()
            return
        
        # Load the DB
        with open(filename, 'rb') as f:
            dict = pickle.load(f)
        #  Rem: Dict is n x 3 cell array {scenario, Dynamic, Market}
        
        # First, find exact match (if any)
        if dict.shape[1] == 3 :
            for i in range(dict.shape[0]):
                scenario = dict[i, 0]
                if self.scenario.isEquivalent(scenario):
                    self.Dynamic = dict[i,1]
                    self.Market = dict[i,2]
                    return
        
        # Check for non-version math as a back-up
        if dict.shape[1] == 3 :
            for i in range(dict.shape[0]):
                scenario = dict[i, 0]
                if self.scenario.isEquivalentIgnoreVersion(scenario):
                    self.Dynamic    = dict[i,1]
                    self.Market     = dict[i,2]
                    print( '[INFO] Using initial guess from different scenario version.\n' )   
                    return
        
        # If not found, generate
        self.generate()
        
        
    # Make a guess either from a cached file or made-up numbers
    def generate(self):
            
        Market = {}
        Dynamic = {}
        steadyScenario  = self.scenario.currentPolicy().steady()
        steady_dir      = PathFinder.getCacheDir(steadyScenario)

        source   = 'cached scenario'
        if os.path.isfile(os.path.join(steady_dir, 'market.pkl')):
            with open(os.path.join(steady_dir, 'market.pkl'), 'rb') as handle:
                Market = pickle.load(handle)
            
        if os.path.isfile(os.path.join(steady_dir, 'dynamics.pkl')):
            with open(os.path.join(steady_dir, 'dynamics.pkl'), 'rb') as handle:
                Dynamic  = pickle.load(handle)

        if len(Market) == 0 or len(Dynamic) == 0:
                
            source = 'made-up numbers'
                
            # Load initial guesses (values come from some steady state results)
            Dynamic['outs']        = np.array([3.1980566])
            Dynamic['caps']        = np.array([9.1898354]) 
            Dynamic['labs']        = 0.5235                               
            captoout            = Dynamic['caps'] / Dynamic['outs']
            debttoout           = np.array([0.75])

            Market['beqs']         = np.array([0.153155])                                  
            Market['capsharesAM']  = captoout / (captoout + debttoout)        # capshare = (K/Y / (K/Y + D/Y)), where K/Y = captoout = 3 and D/Y = debttoout.
            Market['capsharesPM']  = Market['capsharesAM']                      # capshare = (K/Y / (K/Y + D/Y)), where K/Y = captoout = 3 and D/Y = debttoout.
            Market['rhos']         = 4.94974                                  
            Market['invtocaps']    = 0.0078 + 0.056                           # I/K = pop growth rate 0.0078 + depreciation

            Market['investmentToCapital0']    = 0.16
            Market['equityDividendRates']     = 0.05
            Market['worldAfterTaxReturn']     = 0.05
            Market['corpLeverageCost']        = 2
            Market['passLeverageCost']        = 2

            Dynamic['debts']      = Dynamic['outs'] * debttoout
            Dynamic['assetsAM']   = Dynamic['caps'] + Dynamic['debts']           # Assume p_K(0)=1
            Dynamic['assetsPM']   = Dynamic['assetsAM']
            Dynamic['labeffs']    = Dynamic['caps'] / Market['rhos'] 
            Dynamic['investment'] = Dynamic['caps'] * Market['invtocaps']

            Dynamic['caps_foreign'] = 0

        setattr(self, 'Market', Market)
        setattr(self, 'Dynamic', Dynamic)
            
        print( '[INFO] Generated new initial guess from %s. \n' % source )