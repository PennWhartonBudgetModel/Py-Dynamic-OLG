##
# Law of motion on debt distribution and effective rate on federal debt.
#
##

import numpy as np

from helpers import makeIterable

class GovtDebt:
    
    #issueYears indexes correspond to years 1,2,3,5,7,10,30
    issueYears = np.array([0,1,2,4,6,9,29])
    maxDuration = 30
    
    # These properties are pulled from the "Notes" tab in the Interest Rate
    # spreadsheet.  Look on the "Notes" tab on the Interest Rate
    # spreadsheet for instructions on how to get the values from the
    # Treasury's website and compute these numbers.
    
    # The debtDistributionbyMaturity are all marketable issues excluding
    # matured issued and the Federal Financing Bank
    # EffectiveRatesByMaturity are all computed from marketable issues excluding matured
    # issues, Federal Financing bank, floating rate bonds, and TIPS
    
    # This is the distribution as of Dec. 31, YYYY, where YYYY is the year
    # prior to the transition year (Currently Dec. 31, 2018)

    initEffectiveRatesByMaturity = np.array([0.017144911,    0.02157609,    0.023331417,
                                    0.019313077,    0.023684762,    0.02334536,
                                    0.026541761,    0.021956698,    0.028035407,
                                    0.030821047,    0.056841793,    0.0625,    0.05375,
                                    0,    0,    0,    0,    0.045,    0.048908676,
                                    0.044413717,    0.042267492,    0.042847726,
                                    0.039903304,    0.029138061,    0.03343739,
                                    0.032812104,    0.028442139,    0.025357788,
                                    0.028752312,    0.031129337])
    initDebtDistributionByMaturity = np.array([0.265260895,    0.124796689,    0.108126796,
                                      0.078548502,    0.083992989,    0.046798614,
                                      0.051795508,    0.026327696,    0.02646769,
                                      0.027054252,    0.004432692,    0.001091947,
                                      0.001052544,    0.000457525,    0,    0,    0,
                                      0.001691259,    0.002434793,    0.003076968,
                                      0.009654279,    0.012316841,    0.01280166,
                                      0.012878672,    0.012384115,    0.010764614,
                                      0.010812539,    0.010694009,    0.010878851,
                                      0.043407062])
    
    ## Helper function to calculate debt
    #     Inputs: (1) Aggregate (Static or Dynamic), 
    #             (2) Current iteration residuals: Gtilde, Ctilde, Ttilde
    #             (3) Initial debt -- i.e. t=1 debt
    #             (4) Rates on debt
    #             (5) T_model
    #     Outputs: new debts, and 'pure' deficits (excluding Ctilde) 
    # This is the old version of the debt aggregator
    @staticmethod
    def calculateStaticDebt(Aggregate, debt_1, budget, T_model): 
        debts = np.append(debt_1, np.zeros(T_model-1))
        
        #cast
        Aggregate['Ctilde'] = makeIterable(Aggregate['Ctilde'])
        
        deficits = (Aggregate['Gtilde'] + Aggregate['bens'] + (budget['infraSpending'] * Aggregate['pops'])
                    - Aggregate['Ttilde'] - Aggregate['revs'] - Aggregate['lumpSumTaxes'])
    
        #cast
        deficits = makeIterable(deficits)
    
        for t in range(T_model-1):
            debts[t+1] = debts[t]*(1 + budget['debtrates'][t]) + deficits[t] - Aggregate['Ctilde'][t] 
        return (debts, deficits)
    
    # Helper function to calculate debt
    #     Inputs: (1) Aggregate (Static or Dynamic)
    #             (2) Current iteration residuals: Gtilde, Ctilde, Ttilde
    #             (3) Initial debt -- i.e. t=1 debt
    #             (4) Rates on debt
    #             (5) T_model
    #     Outputs: new debts, and 'pure' deficits (excluding Ctilde)
    @staticmethod
    def calculateDynamicDebts(Aggregate, debtPrev, T_model, budget, effectiveRatesByMaturity, debtDistribution, initRate):
        
        # Initializing the debt and the yearly aggregate effective rate
        # Computing the stream of primary deficits (excluding the closure
        # rule)
        debts = np.append(debtPrev, np.zeros(T_model-1))
        effRate = np.zeros(T_model)
        deficits = (Aggregate['Gtilde'] + Aggregate['bens'] + (budget['infraSpending'] * Aggregate['pops'])
                    - Aggregate['Ttilde'] - Aggregate['revs'] - Aggregate['lumpSumTaxes'])
        
        # This is the distribution in year one of the simulation.  This is
        # going to be the small distribution target for every subsequent
        # year
        targetDistribution = GovtDebt.getSmallDistribution(debtDistribution[0,:])
        
        # Initializing the density of the debt and the disaggregated
        # effective rate. The effective rate will be updated.
        # effectiveRatesByMaturity[0] = GovtDebt.initEffectiveRatesByMaturity
        # debtDistribution[0] = GovtDebt.initDebtDistributionByMaturity
        
        # Loop through and compute debt at each period given the interest
        # rates. In each period, we get the new density of debt across all
        # 30 years, and the new effective rates across all 30 years, but
        # these aren't saved.
        for t in range(T_model):
            yieldCurve = np.array([budget['treas1yr'][t], budget['treas2yr'][t], budget['treas3yr'][t], budget['treas5yr'][t], budget['treas7yr'][t], budget['treas10yr'][t], budget['treas30yr'][t]])
            
            # Computes the debts and the effective rates in each period
            # Note that the effective rate (effRate(t)) is the rate of
            # interest paid during time (t).  Given debt tomorrow, debt
            # today, and the deficits, we can back out the coupon rate paid
            # today.  Also, as dynamic.debt ends at time T, and we need
            # debts at tome T+1 in order to compute the effective rate at
            # time T, we have a special case for the last period where debt
            # at time T+1 is saved in tmpDebt.
            
            # Note: don't update effective rate in first period; that is
            # computed externally and is a state variable. 
            
            if t < T_model - 1:
                (debts[t+1], _, effectiveRatesByMaturity[t+1,:], debtDistribution[t+1,:]) = (
                    GovtDebt.GovtDebtDistribution(debts[t], debtDistribution[t,:], deficits[t] - Aggregate['Ctilde'][t], effectiveRatesByMaturity[t,:], yieldCurve, budget['deflator'][t+1]/budget['deflator'][t], targetDistribution))
                if (t > 1):
                    effRate[t] = ((debts[t+1] - deficits[t] + Aggregate['Ctilde'][t])/debts[t] - 1)
                else:
                    effRate[0] = initRate
            else: # for the last year, we change the inflation rate
                 # There is no deflator t+1, so we have to adjust to use
                 # the previous period's inflation (t/t-1).  In the event
                 # that T_model = 1 (one period transition), we hard code
                 # inflation to be 2 percent.
                if ((t-1) < 0): # Setting inflation rate equal to 2 percent if T_model = 1      
                    inflationRate = 1.02 
                else:
                    inflationRate = budget['deflator'][t]/budget['deflator'][t-1]
                
                (tmpDebt, _, _, _) = (
                    GovtDebt.GovtDebtDistribution(debts[t], debtDistribution[t,:], deficits[t] - Aggregate['Ctilde'][t], effectiveRatesByMaturity[t,:], yieldCurve, inflationRate, targetDistribution))
                if (t > 0):
                    effRate[t] = ((tmpDebt - deficits[t] + Aggregate['Ctilde'][t])/debts[t] - 1)
                else:
                    effRate[0] = initRate
                    
        return (debts, deficits, effRate, effectiveRatesByMaturity, debtDistribution)


    ## Helper function to calculate residual needed to fix D/Y 
    #     Inputs: (1) D/Y ratio to target
    #             (2) NIS (as deficit)
    #             (3) outputs (Y)
    #             (4) Rates on debt
    #     Outputs: Residual to be subtracted from debt to make D/Y match
    @staticmethod
    def calculateFixedDebts(closure_debttoout, deficits, outs, debt1, debtrates):
        T = len(deficits)
        # Set debt targets; rem: GDP(T) is not known, so duplicate as kludge
        new_debts = np.append(debt1, np.zeros(T))
        outs = np.append(outs, outs[T])
        residual = np.zeros(T)
        for t in range(T):
            # NOTICE new_debt(t+1) = new_debts(t)*(1+debtrates(t)) + deficit(t) - residual(t)
            tmp_residual = new_debts[t]*(1+debtrates[t]) + deficits[t] - closure_debttoout * outs[t+1]
            # new_debts(t+1)/outs(t+1) > closure_debttoout implies tmp_residual > 0
            if tmp_residual > 0:
                residual[t] = tmp_residual
            new_debts[t+1] = new_debts[t]*(1+debtrates[t]) + deficits[t] - residual[t]
        return residual
    
    @staticmethod
    def GovtDebtDistribution(oldDebt, oldDebtDistribution, deficit, oldEffectiveRates, yieldCurve, inflation, targetDistribution):
        
        # This function takes the existing level of debt, the distribution
        # of that debt over 30 years, the deficit, the effective coupon
        # rate on federal debt, and the yield curve at 1-yr, 2-yrs, 3-yrs,
        # 5-yrs, 7-years, 10-years, and 30-years (7x1 vector)
        
        # The function returns the new level of debt, the effective rate,
        # and the new distribution of the debt
        
        # Start by calculating the levels of debt at each maturity
        oldDebtLevels = oldDebtDistribution * oldDebt
        
        # This is the interest paid this period on the debt inherited from
        # last period.
        interestPaid = np.sum(oldDebtLevels * oldEffectiveRates)
        
        # New issues of debt are the retiring debt plus additional deficits
        # (assuming they're positive, if this is negative more than the
        # size of the maturing debt, we have to think about that.  I doubt
        # it'll ever be a problem). The retiring debt is the debt that's
        # less than one year.  
        newIssuesTotal = oldDebtLevels[0] + deficit + interestPaid
        
        # Advance the age of all old debt and add the new issues
        existingDebtLevels = np.hstack((oldDebtLevels[1:GovtDebt.maxDuration], 0))
        existingDebtRates = np.hstack((oldEffectiveRates[1:GovtDebt.maxDuration], 0))
        
        # Determine the distribution of the new issues of debt.
        # newIssuesDebt is a vector the size of issueYears (usually 7
        # elements: 1, 2, 3, 5, 7, 10, 30). Computnig how much of the debt
        # at each maturity is new and how much is existing debt.  Also
        # loading the yield-curve into a 30-element (for 30 years) array;
        # only the 7 elements at which new debt is issued are non-zero. 

        # This is the distribution of new issues
        if newIssuesTotal < 0 and abs(newIssuesTotal) < sum(existingDebtLevels):
            newIssuesLevels = newIssuesTotal / sum(existingDebtLevels) * existingDebtLevels
        else:
            if (newIssuesTotal < 0):
                newIssuesLevels = -existingDebtLevels
                newIssuesLevels[0] = newIssuesLevels[0] + sum(existingDebtLevels) + newIssuesTotal
            else:
                newIssuesLevels = GovtDebt.allocateNewDebt(existingDebtLevels, newIssuesTotal, targetDistribution)
        
        # The share each which is new and existing debt

        newDebtShare = newIssuesLevels/(existingDebtLevels + newIssuesLevels)
        existingDebtShare = existingDebtLevels/(existingDebtLevels + newIssuesLevels)

        # Checking to make sure that they aren't nan (which happens if
        # there is no debt for a particular year).
        newDebtShare[np.isnan(newDebtShare)] = 0
        existingDebtShare[np.isnan(existingDebtShare)] = 0
        newDebtShare[np.isinf(newDebtShare)] = 0
        existingDebtShare[np.isinf(existingDebtShare)] = 0
        
        # Setting the coupon rates at those 7 issue maturities to the yield
        # curve (note: the coupons are structured to be as close to the
        # expected yield as possible at the time of auction)
        newDebtRates = np.zeros(GovtDebt.maxDuration)
        newDebtRates[GovtDebt.issueYears] = yieldCurve
        
        # Update the effective rate: it's a weighted average of the old
        # debt and the yields on the new debt (the yields tend to roughly
        # equal the coupons on new issues). The effective rate is computed
        # at each maturity year.
        if (newIssuesTotal < 0) and (abs(newIssuesTotal) < np.sum(existingDebtLevels)):
            effectiveRates = existingDebtRates
        else:
            if newIssuesTotal < 0:
                effectiveRates = 0 * existingDebtRates
                effectiveRates[0] = newDebtRates[0]
            else:
                effectiveRates = existingDebtRates * existingDebtShare + newDebtRates * newDebtShare
                
        # Update the debt levels and distribution being carried over the next period.
        debtLevels = existingDebtLevels + newIssuesLevels
        debtDistribution = debtLevels / np.sum(debtLevels)
        
        # Update aggregate debt, deflate the debt back to real terms.
        debt = np.sum(debtLevels)/inflation
        return (debt, interestPaid, effectiveRates, debtDistribution)

    def getSmallDistribution(debtDistribution):
        
        # This function returns the distribution of bonds at the 1, 2, 3,
        # 5, 7, 10, and 30 year intervals. 
        
        smallDistribution= np.zeros(len(GovtDebt.issueYears))
        for x in range(1, len(smallDistribution)+1):
            smallDistribution[x-1] = np.sum(debtDistribution[(GovtDebt.issueYears[max(x-1,1)-1]+1*(x>1)):GovtDebt.issueYears[x-1]+1])
        smallDistribution = smallDistribution / np.sum(smallDistribution)
        return smallDistribution

    def allocateNewDebt(newDebtDist, newIssuesTotal, targetDistribution):
        
        # This function returns the distribution of the new issues at the.  
        # The distribution is
        # calculated to match the distribution of issues for the 1, 2, 3,
        # 5, 7, 10, and 30 year intervals (smallDistribution).
        newIssuesDist = np.zeros(GovtDebt.maxDuration)
        
        #smallOldDist is what the distribution looks like
        #after all of the old issues matured by one year. smallNewDist is
        #our "target" distribution. 
        smallOldDist = GovtDebt.getSmallDistribution(newDebtDist) * np.sum(newDebtDist)
        smallNewDist = targetDistribution * (np.sum(newDebtDist) + newIssuesTotal)
        
        # We're starting with the total value of the new issues
        # (newIssuesAvailable) and then distributing that among the 1, 2,
        # 3, 5, 7, 10, and 30 year bonds in that order to hit our targets.
        # We do each in order until we run out of new issues.
        newIssuesAvailable = newIssuesTotal
        for x in range(max(GovtDebt.issueYears.shape)):
            newIssuesDist[GovtDebt.issueYears[x]] = max(min(newIssuesAvailable, smallNewDist[x] - smallOldDist[x]),0)
            newIssuesAvailable = newIssuesAvailable - newIssuesDist[GovtDebt.issueYears[x]]

        # If any is left over, it's dumped in the 1-year bonds. Not sure
        # that this is a theoretical possibility, but just error
        # correction.
        if newIssuesAvailable > 1e-5:
            print( 'WARNING! newIssues of debt are left over, value: %12.8f put in 1-year bonds.\n' % newIssuesAvailable )
            newIssuesDist += newIssuesAvailable
        return newIssuesDist

    
