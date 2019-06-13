# Used to generate model run guesses.

import os
import scipy.io as sio
import PathFinder
import Scenario
    
class InitialGuess:

    scenario = None
    Dynamic = None
    Market = None
    
    def __init__(self, scenario):
        
        self.scenario = scenario
        self.fetch()
        
        
    # Refresh guess from cache
    def refresh(self):
        self.generate()
        self.save();    
        
    
    # Add guess to DB
    def save(self):
            
        filename = os.join.path(PathFinder.getSourceDir(), 'InitialGuess.mat')
         
        dict = [[]]
        
        if os.path.isfile(filename):
            s = sio.loadmat(filename) 
            dict = s['dict']
            # map is a cell array with a scenario (as Key)
            # and Dynamic, Market as values
            
        # Find exact match (if any), otherwise add
        found = 0
        for i in range(dict.shape[0]):
            scenario = dict[i,0]
            if self.scenario.isEquivalent(scenario):
                dict[i,:] = [scenario, self.Dynamic, self.Market]
                found = 1
                break
            
        if not found:
            dict[0,:] = [self.scenario, self.Dynamic, self.Market]

        sio.savemat(filename, dict)

        
        # Get the guess from DB or generate if not there
        def fetch(self):
            filename = os.join.path(PathFinder.getSourceDir(), 'InitialGuess.mat')
            if not os.join.isfile(filename):
                self.generate();
                self.save();
                return
            
            # Load the DB
            s       = sio.loadmat(filename)
            dict    = s['dict']
            #  Rem: Dict is n x 3 cell array {scenario, Dynamic, Market}

            # First, find exact match (if any)
            for i in range(dict.shape[0]):
                scenario = dict[i, 0]
                if self.scenario.isEquivalent(scenario):
                    self.Dynamic = dict[i,1]
                    self.Market = dict[i,2]
            
            # Check for non-version math as a back-up
            for i in range(dict.shape[0]):
                scenario = dict[i, 0]
                if self.scenario.isEquivalentIgnoreVersion(scenario):
                    self.Dynamic    = dict[i,1]
                    self.Market     = dict[i,2]
                    print( '[INFO] Using initial guess from different scenario version.\n' )           
            
            # If not found, generate
            self.generate()
            self.save()
        
        # Make a guess either from a cached file or made-up numbers
        def generate(self):
            
            Market = []
            Dynamic = []
            steadyScenario  = self.scenario.currentPolicy().steady()
            steady_dir      = PathFinder.getCacheDir(steadyScenario)

            source   = 'cached scenario'
            if os.path.isfile(os.path.join(steady_dir, 'market.mat')):
                Market = sio.loadmat(os.path.join(steady_dir, 'market.mat'))
            
            if os.path.isfile(os.path.join(steady_dir, 'dynamics.mat')):
                Dynamic  = sio.loadmat(os.path.join(steady_dir, 'dynamics.mat'))

            if len(Market) == 0 or len(Dynamic) == 0:
                
                source = 'made-up numbers'
                
                # Load initial guesses (values come from some steady state results)
                Dynamic['outs']        = 3.1980566
                Dynamic['caps']        = 9.1898354 
                Dynamic['labs']        = 0.5235                               
                captoout            = Dynamic['caps'] / Dynamic['outs']
                debttoout           = 0.75

                Market['beqs']         = 0.153155                                  
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

            self['Market']     = Market
            self['Dynamic']    = Dynamic
            
            print( '[INFO] Generated new initial guess from %s. \n' % source )