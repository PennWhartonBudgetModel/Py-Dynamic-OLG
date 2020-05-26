## Tools to test dynamic firm

from scenarioModule import Scenario
from modelTesterModule import ModelTester
from paramGeneratorModule import ParamGenerator
from dynamicFirmGEModule import DynamicFirmGE

import numpy as np
import time

class FirmTester:
        
    def doGELoop():
            
        # Build baseline scenario
        scenario = Scenario(ModelTester.test_params).currentPolicy()
            
        # Get production parameters, remap some parameters -- this is
        # temporary
        paramsProduction = ParamGenerator.production( scenario )
        paramsProduction['capitalShare']   = paramsProduction['alpha']
        paramsProduction['laborShare']     = 1 - paramsProduction['alpha']
            
        Initial          = {'capital0': 10}  
        wage             = 1.5*np.ones(100)
        discountRate     = 0.04*np.ones(100)
        T                = len(wage)
            
            
        # Specify the tolerance level
        tolerance = {}
        tolerance['wage']   = 1e-5
        tolerance['discountRate']   = 1e-5
        excessLabor      = np.ones(100)
        excessShare      = np.ones(100)
            
        # Start the iterative method
        iteration        = 0
        t = time.time()
        while (max(excessLabor[:]) > tolerance['wage'] or max(excessShare[1:T]) > tolerance['discountRate']): 
                
            iteration        = iteration + 1
            Market = {}
            Market['wage']      = wage
            Market['discountRate'] = discountRate
           
            # Initialize firm
            theFirm          = DynamicFirmGE( Initial, Market, paramsProduction )
            laborSupply      = theFirm.getLS()
            laborDemand      = theFirm.getLD()
            excessLabor      = abs(laborSupply - laborDemand)
            shareDemand      = theFirm.getShare()
            excessShare      = abs(shareDemand-1)
            value            = theFirm.getValue()
                     
            # Update guesses
            wage             = (1 / (1 + ((laborSupply - laborDemand) / laborDemand) * 0.1)) * Market['wage']
            for t in range(T-1):
                discountRate[t] = (1 / (1 + ((1 - value[t+1]) / abs(value[t+1])) * 0.1)) * Market['discountRate'][t]
            
            discountRate[T-1] = discountRate[T-2] 
            
        print(time.time() - t)

      
 

    
  
