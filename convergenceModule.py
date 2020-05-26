# Class to manage Convergence of GE loop

from scenarioModule import Scenario

class Convergence :
    
    Market = None
    Dynamic = None
    prevMarket = None
    prevDynamic = None
    tolerance = None
    clearing = {'rhos': None, 'beqs': None, 'invtocaps': None, 'capshares': None}
    scenarioType = None

    isConverged = None
    iteration = None
    maxIterations = 80
    
    def __init__(self, scenario):
        
        self.iteration      = 1
        self.isConverged    = 0
        
        if scenario.isSteady() :
            self.scenarioType = 'steady'
        else:
            self.scenarioType = 'transition'
        
        # Define marketing clearing tolerance 
        self.tolerance = {}
        self.tolerance['rhos']      = 1e-5
        self.tolerance['beqs']      = 1e-4
        self.tolerance['invtocaps'] = 5e-4    

    # Update market price guesses
    def iterate(self, Market, prevMarket):
       
        newMarket = Market

        # Measure errors of iteration prices
        self.clearing = {}
        try:
            self.clearing['rhos']      = max(abs((prevMarket['rhos']           - Market['rhos'])          / Market['rhos']     ))
        except:
            self.clearing['rhos']      = abs((prevMarket['rhos']           - Market['rhos'])          / Market['rhos']     )
        try:
            self.clearing['beqs']      = max(abs((prevMarket['beqs']           - Market['beqs'])          / Market['beqs']     ))
        except:
            self.clearing['beqs']      = abs((prevMarket['beqs']           - Market['beqs'])          / Market['beqs']     )
        try:
            self.clearing['invtocaps'] = max(abs((prevMarket['invtocaps']      - Market['invtocaps'])     / Market['invtocaps']))
        except:
            self.clearing['invtocaps'] = abs((prevMarket['invtocaps']      - Market['invtocaps'])     / Market['invtocaps'])
        try:
            self.clearing['capshares'] = max(abs((prevMarket['capsharesPM']    - Market['capsharesPM'])   / Market['capsharesPM']))
        except:
            self.clearing['capshares'] = abs((prevMarket['capsharesPM']    - Market['capsharesPM'])   / Market['capsharesPM'])
                    
        # Check convergence
        isConverged = ((self.clearing['rhos']      < self.tolerance['rhos']     ) and 
                      (self.clearing['beqs']      < self.tolerance['beqs']     ) and 
                      (self.clearing['invtocaps'] < self.tolerance['invtocaps']))

        self.isConverged = int(isConverged)
        if( isConverged ) :
            return newMarket

        # Update guesses and continue iterations
        
        self.iteration = self.iteration + 1

        # Set damper to update guesses
        #    0 = not dampened, i.e., completely updated to new value
        # 	 1 = fully dampened, i.e., stays the same
        damper = {}
        if self.scenarioType == 'steady':
            damper['rhos']      = 0.75
            damper['beqs']      = 0.75
            damper['capshares'] = 0.75
        elif self.scenarioType == 'transition':
            damper['rhos']      = 0.5
            damper['beqs']      = 0.5
            damper['capshares'] = 0.5
        
        newMarket['beqs']        = damper['beqs']     *prevMarket['beqs']        + (1 - damper['beqs'])      * Market['beqs']
        newMarket['rhos']        = damper['rhos']     *prevMarket['rhos']        + (1 - damper['rhos'])      * Market['rhos']
        newMarket['capsharesAM'] = damper['capshares']*prevMarket['capsharesAM'] + (1 - damper['capshares']) * Market['capsharesAM']
        newMarket['capsharesPM'] = damper['capshares']*prevMarket['capsharesPM'] + (1 - damper['capshares']) * Market['capsharesPM']
 
        return newMarket
    
    # Check whether the convergence should continue
    def stillConverging(self):
        flag = not self.isConverged and (self.iteration < self.maxIterations)
        return flag
    
    # Display convergence status
    def status(self):
        
        status = ('Errors: K/L = %7.6f beqs = %7.6f I/K = %7.6f (K)/A = %7.6f' %
                    (self.clearing['rhos'], self.clearing['beqs'], self.clearing['invtocaps'], self.clearing['capshares']))
        return status
    