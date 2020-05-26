##
# Social Security Support 
#
import numpy as np




class SocialSecurity:
    
    params = None        # Struct of params from ParamGenerator
    bv = None            # Grid of discretized average earnings
    priceIndices = None   # Indices of various prices for this iteration (wages, CPI, etc)
    startyears = None    # Offsets to t=1 of birthyears of cohorts
    realage_entry = None # Actual age of birth in model
    T_model = None       # time of model
    
    incomeFloor = None           # wage indexed floor    for accepting earnings to count for benefits
    incomeCeiling = None         # wage indexed ceiling  for accepting earnings to count for benefits
    
    payrollTaxBrackets = None    # indexed payroll tax brackets
    payrollTaxRates = None       # indexed payroll tax rates
    payrollTaxBurdens = None     # cumulative payroll tax liability
    
    
    # Constructor
    #    Inputs:
    #       params          : struct from ParamGenerator.social_security
    #       bv              : discretized grid of average earnings
    #       priceIndices    : struct of price indicies
    def __init__(self, params, bv, priceIndices, startyears, realage_entry, T_model):
        
        self.params         = params
        self.bv             = bv
        self.priceIndices   = priceIndices
        self.startyears     = startyears
        self.realage_entry  = realage_entry
        self.T_model        = T_model
        
        # Calculate indexed policy variables, use only model periodindices
        t_1                 = priceIndices['T_modelStart']
        t_end               = priceIndices['T_modelEnd']
        wage_inflations     = priceIndices['wage_inflations'][t_1-1:t_end]
        self.incomeFloor    = (params['ssincmins'] * wage_inflations)
        self.incomeCeiling  = (params['ssincmaxs'] * wage_inflations)
        
        # Index the payroll tax structures and put in local vars
        self.indexPayrollTax()

    
    # Get time-varying, indexed Social Security benefits by cohort
    #   Inputs:
    #       i   : cohort number, if -1, then "steady" cohort
    def getBenefitsForCohort(self, i):
        
        if i == -1:
            # Note: retire_year = 1 so that ssbenefits is calculated 
            # for T_model = 1 and indexed for first startyear
            startYear   = self.startyears[0]
            retireYear  = max(self.T_model - 1, 1)
            brackets    = self.params['benefitParamsByCohort'][-1]['brackets']
            rates       = self.params['benefitParamsByCohort'][-1]['rates']
        else:
            startYear   = self.startyears[i]
            retireYear  = self.params['retire_years'][i]
            brackets    = self.params['benefitParamsByCohort'][i]['brackets']
            rates       = self.params['benefitParamsByCohort'][i]['rates']
            
        # Fetch index used to grow brackets
        wage_index = self.priceIndices['wage_inflations']  # full long index (has all cohorts)
        
        # Year cohort turns 62 in current model. 
        # Index benefits by wage index of that year.
        year62 = startYear + (62 - self.realage_entry)
        
        # Cohort based benefits by year
        #   REM: Every household is at a grid point, so can
        #        calculate benefits directly here (outside of
        #        solve_cohort).
        #   NOTE: We set ssbenefits to -Inf otherwise -- it should not be
        #   used in solve_cohort (this will blow it up just in case)
        benefits  = np.ones((self.T_model, (self.bv).shape[0])) * (- float('Inf'))
        
        #Build cohort-specific bracket cutoffs for each year indexed by
        # wage_index of when it turns 62.
        w_year62    = year62 + self.priceIndices['T_modelStart'] # year in index vector when cohort turns 62
        bracket_idx = wage_index[int(w_year62 - 1)] * np.ones((self.T_model,1))
        adjbrackets = brackets * np.tile(bracket_idx, (1, brackets.shape[1]))
        
        # A possible step to implement here would be to deflate the adjusted
        # brackets by another index, e.g. CPI old people,  but currently
        # it's only deflated by CPI, which is already done in ParamGenerator
        
        # Build cumulative benefits matrix to aid benefit calculation
        adjtotben   = np.cumsum(np.diff(adjbrackets, axis = 1) * rates[:, 0:-1], axis = 1) 
        adjtotben   = np.hstack((np.zeros((adjbrackets.shape[0], 1)), adjtotben))  # rem: first column is zero
        
        # Calculate benefits for each year of retirement until end of time.      
        for t in range((max(int(retireYear),1) -  1) , self.T_model):
            for ib in range((self.bv).shape[0]):
                thebracket = np.where(adjbrackets[t,:] <= self.bv[ib])[0][-1]
                benefits[t,ib] = adjtotben[t,thebracket] + rates[t,thebracket]*(self.bv[ib] - adjbrackets[t,thebracket])

        return benefits
  

    ## Calculate SS tax brackets using required indexing
    #    If brackets overlap because of indexing, reorder the brackets 
    def indexPayrollTax( self ):  
            
        # For easier typing, write as local vars
        brackets        = self.params['taxbrackets']
        bracketindices  = self.params['taxindices'][0,:]  # TBD: Make indices time-varying
        rates           = self.params['taxrates']  
            
        if 'cohort_wages' in bracketindices:
            raise Exception('Cannot use index type <cohort_wages>.' )

        indices = np.zeros(brackets.shape)
        for i in range(indices.shape[1]):
            indexname = bracketindices[i]
            # TBD: This should be for all long indices
            if indexname == 'wage_inflations' :
                t_1     = self.priceIndices['T_modelStart']
                t_end   = self.priceIndices['T_modelEnd']
            else :
                t_1     = 1
                t_end   = indices.shape[0]
            
            # Only use model period range of long index
            indices[:,i] = self.priceIndices[indexname][t_1 - 1:t_end]
        
        ssbrackets = brackets * indices

        # Take initial brackets and see where they are after indexing and
        # sorting. Then, apply rates to new topology with rates applying
        # between moved brackets. Resolve any overlaps by taking higher rate.
        new_old_idx = np.argsort(ssbrackets, axis = 1)
        ssbrackets = np.sort(ssbrackets, axis = 1) # rem: sort index from sort() is map from new to old location
        old_new_idx = np.zeros(new_old_idx.shape)
        n_brackets = old_new_idx.shape[1]
        for i in range(n_brackets):          # map sort index from old to new locations 
            old_new_idx[:, new_old_idx[:, i]] = i 

        ssrates = np.full(rates.shape, -float('inf'))
        for year in range(rates.shape[0]): # TBD: Can this be vectorized?
            for r_old in range(n_brackets):
                r_rate      = rates      [year, r_old]
                r_new       = old_new_idx[year, r_old]
                if r_old + 1 == n_brackets : 
                    r_newtop = n_brackets # extends all the way up
                else:
                    r_newtop = old_new_idx[year, r_old+1]
                
                for rr in range(int(r_new),int(r_newtop)):
                    ssrates[year, rr]  = max(r_rate, ssrates[year, rr])
                

        # Calculate tax burdens
        ssburdens = np.cumsum(np.diff(ssbrackets, axis = 1) * ssrates[:, 0:-1], axis = 1) 
        ssburdens = np.hstack((np.zeros((len(ssbrackets), 1)), ssburdens))  # rem: first burden is zero

        self.payrollTaxBrackets = ssbrackets
        self.payrollTaxRates    = ssrates
        self.payrollTaxBurdens  = ssburdens
