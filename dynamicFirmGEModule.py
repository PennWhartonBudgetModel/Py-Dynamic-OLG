##
# Dynamic firms
#    -- firms own the stock of physical capital and undertake the
#    intertemporal investment decisions
#

import numpy as np

class DynamicFirmGE:
    
    T = None                         # total number of time periods
    capitalShare = None              # capital share in production
    laborShare = None                # labor share in production
    depreciationRate = None          # depreciation rate of capital
    wage = None                      # wage 
    discountRate = None              # firm's discount rate
    capital0 = None                  # capital pre-transition
        
    ##
    #  Constructor
    #     INPUTS:   Aggregate = Dynamic or Static aggregates
    #               Market    = the prices of stuff
    #               paramsTax = ParamGenerator.tax()
    #               paramsProduction = ParamGenerator.production()
    def __init__(self, Initial, Market, paramsProduction ):
   
        self.T                          = len(Market['wage'])            
        self.capitalShare               = paramsProduction['capitalShare']
        self.laborShare                 = paramsProduction['laborShare'] 
        self.depreciationRate           = paramsProduction['depreciation']
        self.wage                       = Market['wage']
        self.discountRate               = Market['discountRate']
        self.capital0                   = Initial['capital0']
        
    ##
    #   Solve for labor supply at wage
    def getLS(self):
        laborSupply = np.exp(0.75 * np.log(self.wage))
        return laborSupply
        
    ##
    #   Solve for capital assuming that firm uses the above labor supply
    def getCapital(self):
        laborSupply    = self.getLS()
        capital        = np.zeros(self.T)
        capital[0]     = self.capital0 # capital in the first transition period   
        for t in range(1, self.T):
            capital[t] = ((self.discountRate[t-1] + self.depreciationRate) / (self.capitalShare * laborSupply[t] ** self.laborShare)) ** (1/(self.capitalShare-1))
        
        return capital
        
    ##
    #  Solve for labor demand such that w=MPL
    def getLD(self):
        capital            = self.getCapital()
        # capital[0]         = self.capital0 # capital in the first transition period 
        laborDemand = (self.wage/((1-self.capitalShare) * capital ** self.capitalShare)) ** (-1/self.capitalShare)
        return laborDemand
        
    ##
    #   Calculate the firm's profits given labor and capital
    def getProfit(self):
        profit        = np.zeros(self.T)
        capital       = self.getCapital() # getCapital
        laborSupply   = self.getLS()
        for t in range(self.T-1):
            profit[t] = (capital[t] ** self.capitalShare * laborSupply[t] ** self.laborShare 
                          - capital[t+1] + (1 - self.depreciationRate) * capital[t]
                          - self.wage[t] * laborSupply[t])
        return profit

    ##
    #   Solve for value function backwards
    def getValue(self):
        value        = np.zeros(self.T)
        profit       = self.getProfit() # getProfit
        value[self.T-1]= ((1 + self.discountRate[self.T-1]) / self.discountRate[self.T-1]) * profit[self.T-2]
        for t in range(self.T-2,-1,-1):
            value[t] = profit[t] + (1/(1 + self.discountRate[t])) * value[t+1]
        return value
        
    ##
    #   Solve for household demand for share
    def getShare(self):
        value              = self.getValue() # getValue  
        shareDemand = 1./value # assume a demand function of 1/p where ex-dividend price p=value
        
        return shareDemand
        
    ##
    #   Flag negative profits and value
    def checkMe(self):
        negativeProfit = np.zeros(self.T)
        negativeValue  = np.zeros(self.T)
        profit         = self.getProfit() # getProfit
        value          = self.getValue()  # getValue
        for t in range(self.T):
            if profit[t] < 0:
                negativeProfit[t] = 1
            if value[t] < 0:
                negativeValue[t]  = 1
        
        return (negativeProfit, negativeValue)
