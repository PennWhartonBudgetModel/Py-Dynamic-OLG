##
# Dynamic model solver.

import PathFinder
import Scenario
import os as path
import shutil
import numpy as np
import scipy
import math
import time
import warnings

class ModelSolver:
    
    # Remove cached Scenario files
    @staticmethod
    def removeCached(scenario):
        
        cacheDir = PathFinder.getCacheDir(scenario)
        if path.exists(cacheDir):
            shutil.rmtree(cacheDir)
            
        # For policy shock scenario, remove postShock scenario
        # We leave the baseline scenario
        if scenario.isPolicyShock():
            scenario2 = scenario.postShock()
            cacheDir = PathFinder.getCacheDir(scenario)
            if path.exists(cacheDir):
                shutil.rmtree(cacheDir)
                
    # Solve dynamic model
    @staticmethod
    def solve(scenario, callertag, initial_state): ##ok<*FXUP>
        
        if not 'callertag' in globals():
            callertag  = ''
        
        # Return the solved path if scenario already solved 
        if scenario.isSolved():
            print('Scenario already solved. Run ModelSolver.remove() first to re-solve.\n')
            save_dir = PathFinder.getCacheDir(scenario)
            return save_dir
        
        # Re-direct execution to solvePolicyShock() if needed
        if (scenario.isPolicyShock()):
            save_dir = ModelSolver.solvePolicyShock(scenario, callertag)
            return save_dir
        
        # Make a temporary directory for output which will be copied
        # to the scenario's WorkingDir upon completion.
        (save_dir, callingtag) = PathFinder.getTempWorkingDir(scenario, callertag)
        
        ##
        # Aggregate generation function
        
        def generate_aggregates(Market, DIST_steady, LABs_static, savings_static, DIST_static):
            
            # Define dynamic aggregate generation flag
            isdynamic = (len(DIST_static) == 0) or (len(LABs_static) == 0) or (len(savings_static) == 0)
            
            # Set static optimal decision values to empty values for dynamic aggregate generation
            if isdynamic:
                LABs_static = np.empty((nstartyears,nstartyears))
                savings_static = np.empty((nstartyears,nstartyears))
                
            # Initialize optimal decision value arrays
            os = ['V', 'LABOR', 'SAVINGS', 'CONSUMPTION', 'CONSUMPTION_TAX', 'AVG_EARNINGS', 'TAXABLE_INC', 'OASI_BENEFITS', 'ORD_LIABILITY', 'PAYROLL_LIABILITY', 'PREF_LIABILITY']
            for o in os:
                OPTs[o] = np.zeros((nz,nk,nb,T_life,T_model))
            
            # Initialize array of cohort optimal labor values
            LABs    = np.empty((nstartyears,1))
            savings = np.empty((nstartyears,1))
            
            # Initialize population distribution array
            DIST = np.zeros((nz,nk,nb,T_life,ng,T_model))
            if isSteadyEconomy:
                DIST_trans = np.zeros((nz,nk,nb,T_life,ng,T_model))
            else:
                DIST_trans = np.empty((nz,nk,nb,T_life,ng,T_model))
            
            theSocialSecurity = SocialSecurity(socialsecurity, bv, Market.priceindices, startyears, realage_entry, T_model)
            
            # Calculate indexed policy variables, use only model periodindices
            ssincmins_indexed   = theSocialSecurity.incomeFloor
            ssincmaxs_indexed   = theSocialSecurity.incomeCeiling
            
            sstax_brackets_indexed  = theSocialSecurity.payrollTaxBrackets
            sstax_rates_indexed     = theSocialSecurity.payrollTaxRates
            sstax_burdens_indexed   = theSocialSecurity.payrollTaxBurdens
            
            # Package fixed dynamic optimization arguments into anonymous function
            solve_cohort_ = lambda V0, LAB_static, saving_static, T_past, T_shift, T_active, T_works, ssbenefits, cohort_wageindexes: solve_cohort(
                 V0, LAB_static, saving_static, isdynamic,
                 nz, nk, nb, T_past, T_shift, T_active, T_works, T_model,
                 zs, transz, kv, bv, scenario.beta, scenario.gamma, scenario.sigma, surv,
                 bequest_phi_1, bequest_phi_2, bequest_phi_3,
                 ssbenefits, ssincmins_indexed, ssincmaxs_indexed, cohort_wageindexes,
                 sstax_brackets_indexed, sstax_burdens_indexed, sstax_rates_indexed,
                 taxIndividual['sstaxcredit'], taxIndividual['brackets'], taxIndividual['burdens'], taxIndividual['rates'], 
                 taxIndividual['consbrackets'], taxIndividual['consburdens'], taxIndividual['consrates'],
                 budget['lumpSumTaxes'],
                 taxIndividual['sharePreferredCorp'], taxIndividual['shareOrdinaryCorp'], taxIndividual['shareOrdinaryPass'],
                 taxIndividual['prefbrackets'], taxIndividual['prefburdens'], taxIndividual['prefrates'], 
                 Market['capgains'], taxIndividual['rateCapGain'],
                 Market['beqs'],
                 Market['wages'],
                 Market['capsharesAM'], taxBusiness['shareIncomeCorp'],  taxBusiness['shareIncomePass'],
                 Market['equityDividendRates'], Market['bondDividendRates'], # Equity and bond returns
                 Market['corpDividendRates'], Market['passTaxIncRates'],   # Taxable dividends
                 Market['equityPrices'], np.hstack(Market['equityPrice0'], Market['equityPrices'][0:T_model-1])  # Equity prices (pm and am)
                 )
                
            # Initialize series of terminal utility values
            V0s = np.zeros((nz,nk,nb,T_life))
            
            # TBD: Solve steady state / post-transition path cohort
            if isdynamic:
                
                ssbenefits = theSocialSecurity.getBenefitsForCohort(-1)
                
                # Solve dynamic optimization
                # (Note that active time is set to full lifetime)
                OPT = solve_cohort_(V0s[:,:,:,T_life], [], [], T_pasts[-1], T_shifts[-1], T_life, T_works[-1], ssbenefits, Market.priceindices.cohort_wages[:,-1])
                
                # Define series of terminal utility values
                V0s[:,:,:,1:T_life-1] = OPT.V[:,:,:,2:T_life]
                
            if isSteadyEconomy:
                # Store optimal decision values
                for o in os:
                    OPTs[o][:,:,:,:,0] = OPT[o]
                LABs[0]    = OPT['LABOR']
                savings[0] = OPT['SAVINGS']
            else:
                # Solve transition path cohorts
                OPTs_cohort = np.empty((1,nstartyears))
                
                # TBD parallelize for loop
                for i in range(nstartyears): 
                    
                    # Extract terminal utility values
                    V0 = V0s[:,:,:,T_ends[i]] ##ok<PFBNS>
                    
                    # Calculate cohort-based year-varying benefits policy
                    ssbenefits = theSocialSecurity.getBenefitsForCohort(i) 
                    
                    # Solve dynamic optimization
                    OPTs_cohort[i] = solve_cohort_(V0, LABs_static[i], savings_static[i], T_pasts[i], T_shifts[i], T_actives[i], T_works[i], ssbenefits, Market.priceindices.cohort_wages[:,i])
                    
                    LABs[i]    = OPTs_cohort[i].LABOR
                    savings[i] = OPTs_cohort[i].SAVINGS
                    
                # Construct optimal decision value arrays
                for i in range(nstartyears):
                    for t in range(T_actives[i]):
                        age  = t + T_pasts[i]
                        year = t + T_shifts[i]
                        for o in os:
                            OPTs[o][:,:,:,age,year] = OPTs_cohort[i][o][:,:,:,t]
                            
            if isdynamic:
                if isSteadyEconomy:
                    # Determine steady state age distribution without immigration
                    # TBD rewrite eigenvalues formulas in python
                    # (v, _) = eigs([birth_rate(1:T_model)*ones(1,T_life); [diag(surv(1:T_model,1:T_life-1)), zeros(T_life-1,1)]], k=1, which='LR')
                    # DISTage0 = v / sum(v)
                    
                    # Define initial population distribution using steady state distributions over age and productivity without immigration
                    DIST_next = np.zeros((nz,nk,nb,T_life,ng))
                    DIST_next[:,:,:,:,g['citizen']] = np.tile(np.reshape(np.tile(np.reshape(DISTage0, (1,T_life)), [nz,1]) * DISTz[:,:,g['citizen']], (nz,1,1,T_life)), [1,nk,nb,1]) / (nk*nb)

                    # Specify number of distribution generation years as number of years required to achieve steady state
                    nyears = 153
                else:
                    # Define initial population distribution as steady state distribution
                    DIST_next = DIST_steady[:,:,:,:,:,0]

                    # Specify number of distribution generation years as number of model years
                    nyears = T_model

                for year in range(nyears):

                    # Store population distribution for current year
                    DIST_year = DIST_next
                    DIST[:,:,:,:,:,min(year, T_model)] = DIST_year

                    # Extract optimal decision values for current year
                    K = OPTs['SAVINGS'][:,:,:,:,min(year, T_model)-1]
                    B = OPTs['AVG_EARNINGS'][:,:,:,:,min(year, T_model)-1]

                    # Define population growth distribution
                    DIST_grow     = np.zeros((nz,nk,nb,T_life,ng))
                    P             = sum(DIST_year[:])
                    legal_entry   = np.tile(np.reshape(legal_rate_age[min(year, T_model)-1,:], (1,1,1,T_life,1)), [nz,1,1,1,1])
                    illegal_entry = np.tile(np.reshape(illegal_rate_age[min(year, T_model)-1,:], (1,1,1,T_life,1)), [nz,1,1,1,1])
                    f = lambda x: np.reshape(x, [nz,1,1,x.shape[1],1])

                    DIST_grow[:,0,0,0,g['citizen']] = f(DISTz[:,0,g['citizen']]) * P * birth_rate(min(year, T_model))
                    DIST_grow[:,0,0,:,g['legal']] = f(DISTz[:,:,g['legal']]) * P * legal_entry
                    DIST_grow[:,0,0,:,g['illegal']] = f(DISTz[:,:,g['illegal']]) * P * illegal_entry

                    # Generate population distribution for next year
                    DIST_next = generate_distribution(DIST_year, DIST_grow, K, B, nz, nk, nb, T_life, ng,
                                        np.squeeze(transz[min(year, T_model)-1,:,:,:]), kv, bv, surv[min(year, T_model)-1,:])
                    assert all(DIST_next[:]>=0), 'Negative mass of people at DIST_next.'

                    # Increase legal immigrant population for amnesty, maintaining distributions over productivity
                    DISTz_legal = DIST_next[:,:,:,:,g['legal']] / np.tile(np.sum(DIST_next[:,:,:,:,g['legal']], axis = 0), [nz,1,1,1,1])
                    DISTz_legal(np.isnan(DISTz_legal)) = 1/nz

                    DIST_next[:,:,:,:,g['legal']] = DIST_next[:,:,:,:,g['legal']] + np.tile(np.sum(amnesty*DIST_next[:,:,:,:,g['illegal']], axis = 0), [nz,1,1,1,1]) * DISTz_legal

                    # Reduce illegal immigrant population for amnesty
                    DIST_next[:,:,:,:,g['illegal']] = (1-amnesty) * DIST_next[:,:,:,:,g['illegal']]

                    if isSteadyEconomy:
                        DIST_trans[:,:,:,:,:,0] = DIST_next

            else:
                DIST = DIST_static
            
            # Normalize steady state population distribution
            if isSteadyEconomy:
                DIST_trans = DIST_trans / sum(DIST[:])
                DIST = DIST / sum(DIST[:])
            
            # Generate aggregates
            assert all(DIST[:]>=0, 'WARNING! Negative mass of people at DIST.')
            DIST_gs = np.reshape(np.sum(DIST, axis = 4), [nz,nk,nb,T_life,T_model])
            f = lambda F: np.sum(np.sum(np.reshape(DIST_gs * F, ([], T_model)), axis=0), axis=2)
            
            Aggregate['pops']      = f(1)                                                                                        # Population
            Aggregate['bequests']  = f(OPTs['SAVINGS'] * np.tile(np.reshape(1-surv, (1,1,1,T_life,T_model)), [nz,nk,nb,1,1]))    # Bequests
            Aggregate['labs']      = f(OPTs['LABOR'])                                                                            # Labor
            Aggregate['labeffs']   = f(OPTs['LABOR'] * np.tile(np.reshape(np.transpose(zs, [2,1,0]), (nz,1,1,T_life,T_model)), [1,nk,nb,1,1]))         # Effective labor
            Aggregate['lfprs']     = f(OPTs['LABOR'] > 0.01) / f(1)                                                           # Labor force participation rate
            Aggregate['incs']      = f(OPTs['TAXABLE_INC'])                                                                   # Income
            Aggregate['pits']      = f(OPTs['ORD_LIABILITY'] + OPTs['PREF_LIABILITY'])                                        # Personal income tax
            Aggregate['ssts']      = f(OPTs['PAYROLL_LIABILITY'])                                                             # Capital income tax
            Aggregate['bens']      = f(OPTs['OASI_BENEFITS'])                                                                 # Social Security benefits
            Aggregate['cons']      = f(OPTs['CONSUMPTION'])                                                                   # Consumption
            Aggregate['constax']   = f(OPTs['CONSUMPTION_TAX'])                                                               # Consumption tax
            Aggregate['assetsAM']  = f(np.tile(np.reshape(kv, (1,nk,1,1,1)), [nz, 1,nb,T_life,T_model]))                      # Assets before re-pricing
            Aggregate['assetsPM']  = (Aggregate['assetsAM'] * (np.ones((1,T_model)) + Market['capgains'])                      # Assets after re-pricing            
                                    * (Market['capsharesAM']/Market['capsharesPM']))                                           # Note: The definition of assetsPM corresponds to beginning of period assets at new policy prices, that is, accounting for eventual capital gains.
            Aggregate['laborIncomes'] = f(                                                                                  # Total labor income
                OPTs['LABOR']                                           
                    * np.reshape(np.transpose((zs, [2,1,0])), [nz,1,1,T_life,T_model])                
                    * np.reshape(Market.wages, [1,1,1,1,T_model]))
        
            F = lambda x, y: min(x, y)
            Fv = np.vectorize(F)
            Aggregate['laborIncomeSubjectToSocialSecuritys'] = f(                                                            # Total labor income subject to social security payroll tax
                Fv(OPTs['LABOR']
                        * np.reshape(np.transpose(zs, [2,1,0]), (nz,1,1,T_life,T_model))            
                        * np.reshape(Market['wages'], (1,1,1,1,T_model)), 
                    np.reshape(ssincmaxs_indexed, (1,1,1,1,T_model))))

            Aggregate['capitalIncomes'] = Market['equityDividendRates'] * (Aggregate['assetsPM'] * Market['capsharesPM'])
            Aggregate['capitalIncomeSubjectToPITs'] = (1 - taxIndividual['sharePreferredCorp']) * Aggregate['capitalIncomes']
            Aggregate['averageEffectivePITRates'] = Aggregate['pits'] / (Aggregate['capitalIncomeSubjectToPITs'] + Aggregate['laborIncomes'])
            Aggregate['averageEffectivePITRatesProductSocialSecurityBenefits'] = Aggregate['averageEffectivePITRates'] * Aggregate['pits']
            
            return (Aggregate, LABs, savings, DIST, OPTs, DIST_trans)
        
        
        ## 
        # Define special Scenarios and their generators
        #  then load them
        
        steadyBaseScenario = []
        steady_generator = []
        steady_dir = []
        DIST0 = []
        baselineScenario = []
        base_generator = []
        base_dir = []
        openScenario = []
        open_generator = []
        open_dir = []
        
        # Identify baseline run 
        isBaseline          = scenario.isCurrentPolicy()
        
        # Identify steady state run
        isSteadyEconomy     = scenario.isSteady()

        # Steady state 
        steadyBaseScenario = scenario.currentPolicy().steady()
        steady_generator = lambda: ModelSolver.solve(steadyBaseScenario, callingtag)
        steady_dir = PathFinder.getCacheDir(steadyBaseScenario)

        # Identify initial state (Market0, Dynamic0)
        #   it is passed in or taken from steady
        if 'initial_state' in globals() and initial_state != None:
            Market0     = initial_state['Market']
            Dynamic0    = initial_state['Dynamic']
            DIST0       = initial_state['DIST']
            yearSteady  = initial_state['yearSteady']
        else:
            initial_state   = []      # no initial state for scenario generators
            yearSteady      = scenario.TransitionFirstYear - 1
            if not isSteadyEconomy:
                Market0   = hardyload('market.mat'      , steady_generator, steady_dir)
                Dynamic0  = hardyload('dynamics.mat'    , steady_generator, steady_dir)
                s         = hardyload('distribution.mat', steady_generator, steady_dir)
                DIST0     = s.DIST_trans;
            else: # for steady, load initial guess
                guess       = InitialGuess(scenario)
                Market0     = guess['Market']
                Dynamic0    = guess['Dynamic']
                DIST0       = []
        
        # Scenarios for the baseline transition path
        #    and open economy version for closed
        if not isSteadyEconomy:
            baselineScenario = scenario.currentPolicy()
            base_generator = lambda: ModelSolver.solve(baselineScenario, callingtag, initial_state)
            base_dir = PathFinder.getCacheDir(baselineScenario)
            
            openScenario = scenario.open()
            open_generator = lambda: ModelSolver.solve(openScenario, callingtag, initial_state)
            open_dir = PathFinder.getCacheDir(openScenario)
        
        # Load dependent scenarios
        Dynamic_base = []
        Market_base = []
        Dynamic_open = []
        Market_open = []
        
        if not isBaseline:
            # Load baseline market and dynamics conditions
            Market_base     = hardyload('market.mat'  , base_generator, base_dir)
            Dynamic_base    = hardyload('dynamics.mat', base_generator, base_dir)
        if not isSteadyEconomy and not scenario.isOpen():
            # WARNING: This does not work with solvePolicyShock,
            #   ensure these files already exist from solving shocked 'open'
            Market_open     = hardyload('market.mat'  , open_generator, open_dir)
            Dynamic_open    = hardyload('dynamics.mat', open_generator, open_dir)
        
        
        ## PARAMETERS
        
        beta                = scenario.beta 
        gamma               = scenario.gamma
        sigma               = scenario.sigma
        modelunit_dollar    = scenario.modelunit_dollar
        
        # Find Closure year in model years (t=1 is first year)
        closure_year        = scenario.ClosureYear - scenario.TransitionFirstYear
        
        # Immigration policies
        prem_legal          = scenario.prem_legal       
        amnesty             = scenario.amnesty          

        # Define time constants
        s = ParamGenerator.timing(scenario)
        T_life                  = s['T_life']                 # Total life years
        T_model                 = s['T_model']                # Transition path model years
        realage_entry           = s['realage_entry']          # Real age of model age=0
        startyears              = s['startyears']             # Cohort start years as offsets to year 1
        nstartyears             = len(startyears)
        
        T_pasts   = max(-startyears, 0)                            # Life years before first model year
        T_shifts  = max(+startyears, 0)                            # Model years before first life year
        T_actives = min(startyears+T_life, T_model) - T_shifts     # Life years within modeling period
        T_ends    = min(T_model-startyears, T_life)                # Maximum age within modeling period
        T_steady  = scenario.TransitionFirstYear - yearSteady      # Model years after t=0 when the initial state occurred
        
        # Discretized grids, including shock process
        #   ndem is number of permanent types (hi, low)
        #   g    is population subgroup set: 'citizen',
        #           'legal', 'illegal'
        #   NOTE: DISTz and zs are indexed by g
        s = ParamGenerator.grids( scenario )
        g      = s['g']         # groups: citizen, legal, illegal 
        ng     = s['ng']        # num groups
        nz     = s['nz']        # num labor productivity shocks
        transz = s['transz']    # transition matrix of z's
        DISTz  = s['DISTz']     # steady-state distribution of z's
        zs     = s['zs']        # shocks grid (by demographic type and age)
        nk     = s['nk']        # num asset points
        kv     = s['kv']        # assets grid
        nb     = s['nb']        # num avg. earnings points
        bv     = s['bv']        # avg. earnings grid
        
        # Load production parameters
        production      = ParamGenerator.production( scenario )
        depreciation    = production['depreciation']   # Depreciation rate
        
        # Load population growth parameters
        # Load age-dependent parameters
        s = ParamGenerator.demographics( scenario )
        birth_rate       = s['birth_rate']                     # Annual birth rate
        legal_rate_age   = s['legal_rate_age']                 # Annual legal immigration rate
        illegal_rate_age = s['illegal_rate_age']               # Annual illegal immigration rate
        surv             = s['surv']                           # Survival probabilities by age  
        
        # Load Social Security parameters
        socialsecurity  = ParamGenerator.social_security( scenario )
        T_works         = socialsecurity['T_works']      # retirement time per cohort 
        
        
        # Load out-of-model parameters
        s = ParamGenerator.outofmodel( scenario )
        production['TFP'] = s['TFP']
        
        ##  Budget: interest rates, expenditures, and debt
        #       Rem -- budget.debttoout only used in steady economy
        budget = ParamGenerator.budget( scenario )
        budget['infraSpending'] = s['infraSpending']
        budget['lumpSumTaxes'] = s['lumpSumTaxes']
        
        ## Tax parameters
        #    rem: all US dollars have been converted to modelunit_dollars
        taxIndividual           = ParamGenerator.taxIndividual( scenario )
        taxBusiness             = ParamGenerator.taxBusiness  ( scenario )
        
        # International and economy openess
        international           = ParamGenerator.international( scenario )
        
        # Define parameters on residual value of bequest function.
        s = ParamGenerator.bequest_motive( scenario )
        bequest_phi_1 = s['phi1']                 # phi1 reflects parent's concern about leaving bequests to her children (THIS IS THE ONE WE WANT TO CALIBRATE FOR LATER!)
        bequest_phi_2 = s['phi2']                 # phi2 measures the extent to which bequests are a luxury good
        bequest_phi_3 = s['phi3']                 # phi3 is the relative risk aversion coefficient

        
        ##  
        # BEGIN RESULTS GENERATION
        
        # Clear or create save directory
        if os.path.isfile(save_dir):
            shutil.rmtree(save_dir)
        os.mkdir(save_dir)

        
        ## Static aggregate generation
        if not isBaseline and not isSteadyEconomy:
            
            # Load baseline optimal decisions and population distribution
            s      = hardyload('decisions.mat', base_generator, base_dir)
            LABs_static    = s['LABs']
            savings_static = s['savings']
            
            s      = hardyload('distribution.mat', base_generator, base_dir)
            DIST_static = s['DIST']
            
            # Generate static aggregates
            # (Intermediary structure used to filter out extraneous fields)
            (Static, _, _, Static_DIST, Static_OPTs, _) = generate_aggregates(
                    Market_base, [], LABs_static, savings_static, DIST_static)
            
            # Copy additional static aggregates from baseline aggregates
            for series in ['caps', 'caps_domestic', 'caps_foreign', 'capincs', 'labincs', 'outs', 'investment', 'debts_domestic', 'debts_foreign', 'invest_foreign', 'Gtilde', 'Ttilde', 'Ctilde']:
                Static[series] = Dynamic_base[series]
            
            # Make the firm
            theFirm  = Firm( Static, Market_base, taxBusiness, production )
            
            # Calculate static budgetary aggregate variables
            firmDist = theFirm.distributions()
            
            Static['corpTaxs']           = firmDist['corpTaxs']
            Static['corpDividends']      = firmDist['corpDividends']
            Static['corpDebts']          = firmDist['corpDebts']
            Static['corpTaxbase']        = firmDist['corpTaxbase']
            Static['corpCits']           = firmDist['corpCits']
            Static['corpDebtInterest']   = firmDist['corpDebtInterest']
            Static['corpDepreciation']   = firmDist['corpDepreciation']
            Static['corpExpensing']      = firmDist['corpExpensing']
            Static['corpDebtTaxBenefit'] = firmDist['corpDebtTaxBenefit']
            Static['corpDebtCost']       = firmDist['corpDebtCost']           
            Static['passDividends']      = firmDist['passDividends']
            Static['passDebts']          = firmDist['passDebts']
            Static['passTaxableIncome']  = firmDist['passTaxableIncome']
            Static['passDebtInterest']   = firmDist['passDebtInterest']
            Static['passDepreciation']   = firmDist['passDepreciation']
            Static['passExpensing']      = firmDist['passExpensing']
            Static['passDebtTaxBenefit'] = firmDist['passDebtTaxBenefit']
            
            # Notice that firmDist calculates capital income & taxes with
            # equityDividendRates endogenous to the firm class, thus
            # accounting for policy changes
            Static['corpDividendsForeign']   = firmDist['corpDividendsForeign']
            Static['passDividendsForeign']   = firmDist['passDividendsForeign']
            Static['corpForeignWithholding'] = firmDist['corpForeignWithholding']
            Static['passForeignWithholding'] = firmDist['passForeignWithholding']

            Static['revs'] = (Static['pits'] + Static['ssts'] + Static['corpTaxs'] + Static['constax'] +
                          Static['corpForeignWithholding'] + Static['passForeignWithholding'])            
            Static['GNI']  = Static['outs'] - (Static['corpDividendsForeign'] + Static['passDividendsForeign'])
            
            # In general, we have:
            # Static.debts = Static.debts_domestic + Static.debts_foreign;
            # But this is not true in the static economies.
            # Since static decisions are fixed at the baseline, domestic and
            # foreign debts do not change. But total debt changes due to new
            # tax policy, which implies the equality above no longer holds.
            # Notice that debts is a combination of the actual static debts and
            # the residual mismatch from markets not clearing
            #    Rem: Dynamic_base(1) is supposed to be D' debt carried
            #    from steady state (before policy change)
            
            if scenario['UseStaticDebt']:
                (Static['debts'], Static['deficits']) = GovtDebt.calculateStaticDebt( Static, Static['Gtilde'], Static['Ctilde'], Static['Ttilde'], Dynamic_base['debts'][0], budget, budget, T_model)
            else:
                Market0['effectiveRatesByMaturity'][0,:]   = Market0['effectiveRatesByMaturity_next']
                Market0['debtDistributionByMaturity'][0,:] = Market0['debtDistributionByMaturity_next']
                (Static['debts'], Static['deficits'], Market0['bondDividendRates'], Market0['effectiveRatesByMaturity'], Market0['debtDistributionByMaturity']) = GovtDebt.calculateDynamicDebts( Static, Static['Gtilde'], Static['Ctilde'], Static['Ttilde'], Dynamic_base['debts'][0], T_model, budget, Market0['effectiveRatesByMaturity'][0,:], Market0['debtDistributionByMaturity'][0,:], Market0['bondDividendRates'][0]) 
             
            Static['infraSpending'] = budget['infraSpending']
            Static['lumpSumTaxes'] = budget['lumpSumTaxes']
            
            # Total assets
            # Note: tot_assets is a sum of choice variables, those are constant at baseline values
            price               = Market_base['equityPrices']
            Static['tot_assetsAM'] = np.hstack((Market_base['equityPrice0'], price[0:T_model-1])) * Static['caps'] + Static['debts_domestic'] + Static['debts_foreign']
            Static['tot_assetsPM'] = price * Static['caps'] + Static['debts_domestic'] + Static['debts_foreign']
                    
            # Ad field for symmetry w/ Dynamic
            Static['is_converged'] = False
                                
            # Save static aggregates
            sio.savemat(os.join.path(save_dir, 'statics.mat')              , Static)
            sio.savemat(os.join.path(save_dir, 'Static_decisions.mat')     , Static_OPTs)
            sio.savemat(os.join.path(save_dir, 'Static_distribution.mat')  , Static_DIST)        
                
        
        ## Dynamic aggregate generation
        
        # Set initial guesses (Market_init, Dynamic_init)
        # TBD: Revise how we get initial guesses
        guessSource = 'steady'
        if not isSteadyEconomy:
            guessSource = scenario.OpennessPath
            
        if guessSource == 'steady':
            # rem: steady guess is from guess file
            Market_init     = Market0
            Dynamic_init    = Dynamic0
            Market_init['bondDividendRates'] = budget['debtrates'] # Use the long-term average real interest coupon
            Market_init['effectiveRatesByMaturity_next'][0:T_model,0:GovtDebt.maxDuration] = None
            Market_init['debtDistributionByMaturity_next'][0:T_model,0:GovtDebt.maxDuration] =  None
            if not scenario.UseStaticDebt:
                Market_init['bondDividendRates_next'] = (1+sum(GovtDebt.initEffectiveRatesByMaturity * GovtDebt.initDebtDistributionByMaturity)) / budget['steadyStatePriceGrowth'] - 1
                Market_init['effectiveRatesByMaturity_next'][0,:] = GovtDebt.initEffectiveRatesByMaturity
                Market_init['debtDistributionByMaturity_next'][0,:] =  GovtDebt.initDebtDistributionByMaturity
                    
        elif guessSource == 'open':    
            if isBaseline:
                for p in Market0.keys():
                    value0         = Market0[p]
                    Market_init[p] = value0
                    if not isinstance(value0, dict) and value0.shape[1] == 1:
                        Market_init[p] = value0 * np.ones((1,T_model))
                for p in Dynamic0.keys():
                    value0          = Dynamic0[p]
                    Dynamic_init[p] = value0
                    if not isinstance(value0, dict) and value0.shape[1] == 1:
                        Dynamic_init[p] = value0 * np.ones((1,T_model))
            else:
                # Guess is baseline
                Market_init  = Market_base
                Dynamic_init = Dynamic_base
               
        else:                
            if isBaseline:
                # Guess is from open_base
                Market_init  = Market_open
                Dynamic_init = Dynamic_open
            else:
                # Guess from closed_base
                Market_init  = Market_base
                Dynamic_init = Dynamic_base
                
               
        ##
        # Define marketing clearing tolerance and initialize error term
        tolerance['rhos']      = 1e-5
        tolerance['beqs']      = 1e-4
        tolerance['invtocaps'] = 5e-4
        isConverged         = False
        
        # Set damper to update guesses
        #    0 = not dampened, i.e., completely updated to new value
        #      1 = fully dampened, i.e., stays the same
        if isSteadyEconomy:
            damper['rhos']      = 0.75
            damper['beqs']      = 0.75
            damper['capshares'] = 0.75
        else:
            damper['rhos']      = 0.5
            damper['beqs']      = 0.5
            damper['capshares'] = 0.5 
        
        # Initialize iteration count and set maximum number of iterations
        iter    =   0
        itermax = 100
        
        # Display header
        print( 'Started at: %s \n' % datetime )
        print( '%s\n' % scenario.shortDescription )

        # Initialize Market, Dynamic to guesses 
        #   rem: guesses should be right size.
        Market  = Market_init
        Dynamic = Dynamic_init
                
        # After-tax returns to foreign investors (from initial state)
        #   This is used to pin down capital returns when opening
        #   the economy.
        Market['worldAfterTaxReturn']  = Market0['worldAfterTaxReturn'] 
                
        # Leverage cost (nu) parameters from initial state
        Market['corpLeverageCost']     = Market0['corpLeverageCost']
        Market['passLeverageCost']     = Market0['passLeverageCost']

        # Set other vars from initial state
        Market['investmentToCapital0'] = Market0['investmentToCapital0'] # overwritten in 'steady'
        if not isSteadyEconomy:
            # Take some price history from initial state
            Market['averageWagesHistory']  = Market0['averageWagesHistory']   
            Market['userCostCapital0']     = Market0['userCostCapital0']
            Market['priceCapital0']        = Market0['priceCapital0']
        
        # Initialize convergence variables
        beqs            = Market['beqs']
        rhos            = Market['rhos']
        capsharesAM     = Market['capsharesAM']
        capsharesPM     = Market['capsharesPM']
        invtocaps       = Market['invtocaps']
        
        # Iterate over economy until convergence
        while not isConverged and iter < itermax:
            
            # Increment iteration count
            iter = iter + 1
            print(' Iteration %2d  ...  RUNNING' % iter)
            
            # Update guesses
            Market['beqs']        = damper['beqs'] * Market['beqs'] + (1 - damper['beqs']) * beqs
            Market['rhos']        = damper['rhos'] * Market['rhos'] + (1 - damper['rhos']) * rhos
            Market['capsharesAM'] = damper['capshares'] * Market['capsharesAM'] + (1 - damper['capshares']) * capsharesAM
            Market['capsharesPM'] = damper['capshares'] * Market['capsharesPM'] + (1 - damper['capshares']) * capsharesPM

            Market['invtocaps']   = invtocaps
                
            # Calculate g'vt budget residuals
            # TBD: Revise where these come from
            Gtilde      = 0
            Ttilde      = 0
            Ctilde      = 0
            if not isSteadyEconomy:
                if scenario.isOpen():
                    if isBaseline:
                        # Note: These calculations are 1 iteration behind
                        baseline_bens = Dynamic['bens']
                        baseline_revs = Dynamic['revs']
                        baseline_outs = Dynamic['outs']
                    else:
                        baseline_bens = Static['bens']
                        baseline_revs = Static['revs']
                        baseline_outs = Dynamic_base['outs']

                    Gtilde = budget['outlays_by_GDP'] * baseline_outs - baseline_bens
                    Ttilde = budget['tax_revenue_by_GDP'] * baseline_outs - baseline_revs 
                    Ctilde = np.zeros((1,T_model))

                else:

                    # Set g'vt budget residuals from open economy
                    Gtilde      = Dynamic_open['Gtilde']
                    Ttilde      = Dynamic_open['Ttilde']
                    Ctilde      = Dynamic_open['Ctilde']                
                
            ##
            # Firms sector
            theFirm  = Firm( Dynamic, Market, taxBusiness, production )

            # Capital gains are from price of capital (in Firm)
            Market['capgains'] = theFirm.capitalGains()

            # Other prices
            Market['wages']            = theFirm.wageRequired()
            Market['averageWages']     = Market['wages'] * Dynamic['labeffs'] / Dynamic['labs']
            Market['MPKs']             = theFirm.MPK()  # Just for reporting
            
            # If steady-state, calibrate certain things:
            #    * set leverage cost params from new interest rate
            #    * set user cost capital from new I/K (if adjustment costs)
            #    * set the "world return" for use by transition paths
            #    * set average wage 
            #    * set I/K -- this is TEMP: we should be calculating
            #               final steady state I/K directly (so being lazy here)
            if isSteadyEconomy:
                theFirm.resetLeverageCost()
                theFirm.resetPriceCapital()
                Market['priceCapital0']        = theFirm['priceCapital0']
                Market['userCostCapital0']     = theFirm['userCostCapital0']
                Market['corpLeverageCost']     = theFirm['corpLeverageCost']
                Market['passLeverageCost']     = theFirm['passLeverageCost']
                
                Market['averageWagesHistory']  = Market['averageWages'][0]
                Market['worldAfterTaxReturn']  = Market['equityDividendRates'][0] * (1 - taxBusiness['rateForeignBusinessIncome'][0]) + Market['capgains'][0] # Capgains=0 in steady
                Market['investmentToCapital0'] = Dynamic['investment'] / Dynamic['caps']
            
            # Compute price indices
            Market['priceindices']     = ModelSolver.generateIndices(Market,
                  budget, nstartyears, realage_entry, T_steady, T_model, T_life)
            
            # Payouts from firm
            firmDist  = theFirm.distributions()
            
            # 'Price' of assets -- HH own equal shares of both bond & equity funds
            Market['userCostCapital0']     = theFirm['userCostCapital0']
            Market['priceCapital0']        = theFirm['priceCapital0']
            Market['equityPrice0']         = theFirm['priceCapital0']
            Market['equityPrices']         = theFirm['priceCapital']
            Market['equityDividendRates']  = firmDist['equityDividendRates']
            
            Market['corpDividendRates']    = firmDist['corpDividendRates']
            Market['passDividendRates']    = firmDist['passDividendRates']
            Market['passTaxIncRates']      = firmDist['passTaxIncRates']
            
            Market['bondPrice0']           = 1
            Market['bondPrices']           = np.ones((1,T_model))
            if scenario.UseStaticDebt:
                Market['bondDividendRates']    = budget['debtrates'] #rem: dividendrate is per $ of assets
            
            # Save previous iteration capital for closed economy klRatio
            prev_iter['caps'] = Dynamic['caps']
            
            # Generate dynamic aggregates
            (Dynamic, LABs, savings, DIST, OPTs, DIST_trans) = generate_aggregates(Market, DIST0, [], [], [])
            
            # Calculate and record additional dynamic aggregates
            # (Note that open economy requires capital calculation before debt calculation 
            # while closed economy requires the reverse)
            
            Dynamic['corpDividends']      = firmDist['corpDividends']
            Dynamic['corpTaxs']           = firmDist['corpTaxs']
            Dynamic['corpDebts']          = firmDist['corpDebts']
            Dynamic['corpTaxbase']        = firmDist['corpTaxbase']
            Dynamic['corpCits']           = firmDist['corpCits']
            Dynamic['corpDebtInterest']   = firmDist['corpDebtInterest']
            Dynamic['corpDepreciation']   = firmDist['corpDepreciation']
            Dynamic['corpExpensing']      = firmDist['corpExpensing']
            Dynamic['corpDebtTaxBenefit'] = firmDist['corpDebtTaxBenefit']
            Dynamic['corpDebtCost']       = firmDist['corpDebtCost']
            Dynamic['passDividends']      = firmDist['passDividends']
            Dynamic['passDebts']          = firmDist['passDebts']
            Dynamic['passTaxableIncome']  = firmDist['passTaxableIncome']
            Dynamic['passDebtInterest']   = firmDist['passDebtInterest']
            Dynamic['passDepreciation']   = firmDist['passDepreciation']
            Dynamic['passExpensing']      = firmDist['passExpensing']
            Dynamic['passDebtTaxBenefit'] = firmDist['passDebtTaxBenefit']
            
            Dynamic['corpDividendsForeign']   = firmDist['corpDividendsForeign']
            Dynamic['passDividendsForeign']   = firmDist['passDividendsForeign']
            Dynamic['corpForeignWithholding'] = firmDist['corpForeignWithholding']
            Dynamic['passForeignWithholding'] = firmDist['passForeignWithholding']
            

            Dynamic['revs'] = (Dynamic['pits'] + Dynamic['ssts'] + Dynamic['corpTaxs'] + Dynamic['constax'] +
                           Dynamic['corpForeignWithholding'] + Dynamic['passForeignWithholding'])

            if isSteadyEconomy:
                # Calculate debt, capital, and output
                # (Numerical solver used due to absence of closed form solution)
                f_debts = lambda outs: budget['debttoout'] * outs
                f_caps  = lambda debts: (Dynamic['assetsPM'] - debts) / Market['equityPrices']
                f_outs  = lambda caps: theFirm.output(caps, Dynamic.labeffs )
                x_ = scipy.optimize.fsolve(lambda x: x - np.vstack((f_debts(x[2]), f_caps(x[0]), f_outs(x[1]))), np.zeros((3,1)))
                Dynamic['debts'] = x_[0]
                Dynamic['caps']  = x_[1]
                Dynamic['outs']  = x_[2]

                Dynamic['debts_domestic'] = Dynamic['debts']
                Dynamic['debts_foreign']  = np.zeros((1,T_model))
                Dynamic['deficits']       = np.zeros((1,T_model))
                Dynamic['caps_domestic']  = Dynamic['caps']
                Dynamic['caps_foreign']   = np.zeros((1,T_model))
                Dynamic['invest_foreign'] = np.zeros((1,T_model))
                Dynamic['tot_assetsAM']   = Dynamic['assetsAM']
                Dynamic['tot_assetsPM']   = Dynamic['assetsPM']

                # Calculate income
                Dynamic['labincs'] = Dynamic['labeffs'] * Market['wages']
                Dynamic['capincs'] = Market['MPKs'] * Dynamic['caps']
                Dynamic['GNI']     = Dynamic['outs'] - (Dynamic['corpDividendsForeign'] + Dynamic['passDividendsForeign'])

                # Gross investment in physical capital
                #    and resulting K' for "next" period
                DIST_gs            = np.reshape(np.sum(DIST, axis=4), (nz,nk,nb,T_life,T_model))
                assets_tomorrow    = np.sum(np.sum(np.reshape(DIST_gs * OPTs['SAVINGS'], (-1, T_model)), axis = 0), axis = 2)
                Dynamic['caps_next']  = (Market['capsharesPM'] * (assets_tomorrow - Dynamic['bequests'])) / Market['equityPrices']
                Dynamic['investment'] = Dynamic['caps_next'] - (1 - depreciation) * Dynamic['caps']

                Dynamic['caps_domestic_next'] = Dynamic['caps_next']
                Dynamic['caps_foreign_next'] = 0

                # "Next" period debt. 
                #   TBD: Revise to get real deficit (rem: none in
                #   steady state)
                Dynamic['debts_next'] = (Dynamic['debts'][0] * (1 + Market['bondDividendRates'][0])
                                     + 0.035 * Dynamic['outs'][0]) # TBD: hardcoded deficit for now
                Dynamic['debts_foreign_next'] = 0    # steady-state has no foreigners

                # Update guesses
                rhos      = Dynamic['caps'] / Dynamic['labeffs']
                beqs      = Dynamic['bequests'] / sum(DIST_trans[:])          # Note: capgains is zero in steady state, so bequests don't need to be changed
                invtocaps = Dynamic['investment'] / Dynamic['caps']
                capsharesAM = (Dynamic['assetsAM'] - Dynamic['debts']) / Dynamic['assetsAM']
                capsharesPM = (Dynamic['assetsPM'] - Dynamic['debts']) / Dynamic['assetsPM']
            
            else:
                # Capital(1), Debts(1) come from inital state (K', D')
                Dynamic['caps'][0] = Dynamic0['caps_next']
                Market['effectiveRatesByMaturity'][0:T_model,0:GovtDebt['maxDuration']] = None
                Market['debtDistributionByMaturity'][0:T_model,0:GovtDebt['maxDuration']] = None
                if not scenario['UseStaticDebt']:
                    Market['effectiveRatesByMaturity'][0,:] = Market0['effectiveRatesByMaturity_next']
                    Market['debtDistributionByMaturity'][0,:] = Market0['debtDistributionByMaturity_next']
                    Market['bondDividendRates'][0] = Market0['bondDividendRates_next']

                # Calculate debt and take-up
                #  NOTE: Foreign debt take-up is for next period debt
                #        We cumulate the take-ups to generate foreign
                #        debt holdings
                if scenario['UseStaticDebt']:
                    (Dynamic['debts'], Dynamic['deficits']) = GovtDebt.calculateStaticDebt( Dynamic, Gtilde, Ctilde, Ttilde, Dynamic0['debts_next'], budget, T_model)
                else:
                    (Dynamic['debts'], Dynamic['deficits'], Market['bondDividendRates'], Market['effectiveRatesByMaturity'],Market['debtDistributionByMaturity']) = (
                    GovtDebt.calculateDynamicDebts( Dynamic, Gtilde, Ctilde, Ttilde, Dynamic0['debts_next'], T_model, budget, Market['effectiveRatesByMaturity'][0,:], Market['debtDistributionByMaturity'][0,:], Market['bondDividendRates'][0]))

                debt_foreign_1  = Dynamic0['debts_foreign_next']
                new_debt_issued = Dynamic['deficits'] + (Dynamic['debts'] * Market['bondDividendRates'])
                Dynamic['debts_foreign']   = np.cumsum(np.hstack((debt_foreign_1,                                             
                                            new_debt_issued[0:T_model-1] * international['debtTakeUp'][0:T_model-1])))   
                Dynamic['debts_domestic']  = Dynamic['debts'] - Dynamic['debts_foreign']

                # Calculate capital and output
                # Note: Dynamic.assetsPM represents current assets at new prices.
                klRatio = theFirm.calculateKLRatio(Market['worldAfterTaxReturn'], prev_iter['caps'], Dynamic['labeffs'])
                open_econ_caps          = klRatio * Dynamic['labeffs']
                Dynamic['caps_domestic']   = (Dynamic['assetsPM'] - Dynamic['debts_domestic']) / Market['equityPrices']
                Dynamic['caps_foreign']    = (open_econ_caps - Dynamic['caps_domestic']) * international['capitalTakeUp']
                Dynamic['caps_foreign'][0] = Dynamic0['caps_foreign_next']    # Investment was set t=0 -- cannot be changed
                Dynamic['caps']            = Dynamic['caps_domestic'] + Dynamic['caps_foreign']
                Dynamic['outs']            = theFirm.output(Dynamic['caps'], Dynamic['labeffs'])

                # Converge to find Ctilde which closes the D/Y ratio
                Ctilde_error    = float('Inf')
                Ctilde          = np.zeros((1, T_model))
                while Ctilde_error > 1e-13 :

                    prev_Ctilde = Ctilde

                    # Calculate Ctilde to close debt growth
                    # Rem: This effects outs, so need to converge
                    closure_debttoout   = Dynamic['debts']['closure_year']/Dynamic['outs']['closure_year']
                    cont_Ctilde         = GovtDebt.calculateFixedDebts(closure_debttoout,   
                                                                       Dynamic['deficits'][closure_year-2:T_model],         
                                                                       Dynamic['outs'][closure_year-2:T_model],             
                                                                       Dynamic['debts'][closure_year-2],
                                                                       Market['bondDividendRates'][closure_year-2:T_model])
                    # Update Ctilde
                    Ctilde = np.hstack((np.zeros((1, closure_year-2)), cont_Ctilde))

                    # Recalculate debt and check if D/Y has been fixed
                    #   Note: D/Y for t=ClosureYear is unchanged by Ctilde
                    if scenario['UseStaticDebt']:
                        (Dynamic['debts'], Dynamic['deficits']) = GovtDebt.calculateStaticDebt( Dynamic, Gtilde, Ctilde, Ttilde, Dynamic0.debts_next, budget, T_model)
                    else:
                        (Dynamic['debts'], Dynamic['deficits'], Market['bondDividendRates'], Market['effectiveRatesByMaturity'], Market['debtDistributionByMaturity']) = (
                            GovtDebt.calculateDynamicDebts( Dynamic, Gtilde, Ctilde, Ttilde, Dynamic0['debts_next'], T_model, budget, Market['effectiveRatesByMaturity'][0,:], Market['debtDistributionByMaturity'][0,:], Market['bondDividendRates'][0]))

                    # Re-calculate capital and output
                    new_debt_issued         = Dynamic['deficits'] - Ctilde + (Dynamic['debts'] * Market['bondDividendRates'])
                    Dynamic['debts_foreign']   = np.cumsum(np.hstack((debt_foreign_1,                                             
                                            new_debt_issued[0:T_model-1] * international['debtTakeUp'][0:T_model-1]))) 
                    Dynamic['debts_domestic']  = Dynamic['debts'] - Dynamic['debts_foreign']

                    # We do not recalculate the KL ratio -- so
                    # open_econ_caps stays the same
                    Dynamic['caps_domestic']   = (Dynamic['assetsPM'] - Dynamic['debts_domestic']) / Market['equityPrices']
                    Dynamic['caps_foreign']    = (open_econ_caps - Dynamic['caps_domestic']) * international['capitalTakeUp']
                    Dynamic['caps_foreign'][0] = Dynamic0['caps_foreign_next']    # Investment was set t=0 -- cannot be changed
                    Dynamic['caps']            = Dynamic['caps_domestic'] + Dynamic['caps_foreign']

                    # Check for capital going wacky
                    too_low_caps = np.where(Dynamic['caps'] <= 0)
                    if len(too_low_caps) != 0 :
                        # Ctilde did not fix debt explosion in time
                        raise Exception('Capital becomes negative at t=%u.', too_low_caps[0])
                    too_high_caps = np.where(Dynamic['caps'] > 1e6)
                    if len(too_high_caps) != 0 :
                        # Ctilde did not fix debt explosion in time
                        raise Exception( 'Capital becomes >1e6 at t=%u.', too_high_caps[0])

                    Dynamic['outs'] = theFirm.output(Dynamic['caps'], Dynamic['labeffs'])

                    # Check that Ctilde has converged
                    Ctilde_error = max(abs(Ctilde - prev_Ctilde))

                Dynamic['tot_assetsAM']   = Dynamic['assetsAM'] + Dynamic['debts_foreign']
                Dynamic['tot_assetsPM']   = Dynamic['assetsPM'] + Dynamic['debts_foreign']

                # Calculate income
                Dynamic['labincs'] = Dynamic['labeffs'] * Market['wages']
                Dynamic['capincs'] = Market['MPKs'] * Dynamic['caps']
                Dynamic['GNI']     = Dynamic['outs'] - (Dynamic['corpDividendsForeign'] + Dynamic['passDividendsForeign'])

                # Gross investment in physical capital.
                #   (see above for comment in 'open' economy)
                #   Also, note, that closed econ capital does not
                #   converge to final steady state quickly, so I/K
                #   may have big jump at T_model with this approach.
                #   So, we use last I/K for smoothness.
                # TBD: Using I/K from init steady state
                #   but this should be final steady state (e.g. if
                #   pop_growth changes.)
                Dynamic['investment'] = np.hstack((Dynamic['caps'][1:T_model], 0)) - (1 - depreciation) * np.hstack((Dynamic['caps'][0:T_model-1], 0))
                Dynamic['investment'][T_model] = Market['invtocaps'][T_model-1] * Dynamic['caps'][T_model]

                # Calculate foreign investment series (for reporting)
                Dynamic['invest_foreign'] = (np.hstack((Dynamic['caps_foreign'][1:T_model], Dynamic['caps_foreign'][T_model])) 
                                  - (1 - depreciation) * np.hstack((Dynamic['caps_foreign'][0:T_model-1], Dynamic['caps_foreign'][T_model-1])))

                # Update guesses
                # Note: Dynamic.assets represents current assets at new prices.
                #       Bequests should also be priced according to the new policy.
                #       So we apply today's prices to yesterday's bequests and capshares.
                rhos      = Dynamic['caps'] / Dynamic['labeffs']
                beqs      = np.hstack((Dynamic0['bequests'] * (1 + Market0['capsharesPM'] * Market['capgains'][0]),
                             Dynamic['bequests'][0:T_model-1] * (1 + Market['capsharesPM'][0:T_model-1] * Market['capgains'][1:T_model]))) / Dynamic['pops']
                invtocaps = Dynamic['investment'] / Dynamic['caps']

                # NOTE: capshares/capital calculation potentially
                #   DESTROYS domestic capital. This is because debt
                #   crowds out capital in the HH assets and the g'vt
                #   can create an arbitrarily large increase in debt
                #   from period to period.
                capsharesAM = (Dynamic['assetsAM'] - Dynamic['debts_domestic']) / Dynamic['assetsAM']
                capsharesPM = (Dynamic['assetsPM'] - Dynamic['debts_domestic']) / Dynamic['assetsPM']

            # Store g'vt budget residuals to struct
            Dynamic['Gtilde'] = Gtilde
            Dynamic['Ttilde'] = Ttilde
            Dynamic['Ctilde'] = Ctilde
            
            # Save the out of model adjustment variables
            Dynamic['infraSpending'] = budget['infraSpending']
            Dynamic['lumpSumTaxes'] = budget['lumpSumTaxes']
                   
            # Calculate market clearing series
            clearing['rhos']   = max(abs((Market['rhos']      - rhos)      / rhos     ))
            clearing['beqs']   = max(abs((Market['beqs']      - beqs)      / beqs     ))
            clearing['invtocaps'] = max(abs((Market['invtocaps'] - invtocaps) / invtocaps))
            clearing['capshares'] = max(abs((Market['capsharesPM'] - capsharesPM) / capsharesPM))
                    
            # Check convergence
            isConverged = ((clearing['rhos']      < tolerance['rhos']     ) and (clearing['beqs']      < tolerance['beqs']     ) and 
                          (clearing['invtocaps'] < tolerance['invtocaps']))
            
            # Erase 'RUNNING' text and then print errors
            print('\b\b\b\b\b\b\b')
            print('Errors: K/L = %7.6f beqs = %7.6f I/K = %7.6f (K)/A = %7.6f\n' % (clearing['rhos'], clearing['beqs'], clearing['invtocaps'], clearing['capshares']))
        
        print( '\nFinished at: %s\n' % datetime)

        Dynamic['is_converged'] = isConverged
        # Issue warning if did not converge
        if not Dynamic['is_converged']:
            warnings.warn('Model did not converge.')
        
        # Save market conditions, HH policies, and dynamic aggregates
        scipy.io.savemat(os.join.path(save_dir, 'market.mat'), Market)
        scipy.io.savemat(os.join.path(save_dir, 'dynamics.mat'), Dynamic)
        decisions = {'OPTs': OPTs, 'LABs': LABs, 'savings': savings}
        scipy.io.savemat(os.join.path(save_dir, 'decisions.mat'), decisions)
        if isSteadyEconomy: 
            DIST = {'DIST': DIST, 'DIST_trans': DIST_trans}
            scipy.io.savemat(os.join.path(save_dir, 'distribution.mat'), DIST)
        else:
            scipy.io.savemat(os.join.path(save_dir, 'distribution.mat'), DIST)
        
        ## Elasticity calculation
        if isSteadyEconomy:

            # Calculate capital to output ratio
            captoout = (Dynamic['assetsPM'] - Dynamic['debts']) / Dynamic['outs']

            # Calculate labor elasticity
            LAB_  = LABs[0]
            DIST_ = np.sum(DIST['DIST'][:,:,:,:,:,0], axis = 4)

            workind = (LAB_ > 0.01)

            workmass = np.sum(DIST_[workind])
            frisch   = np.sum(DIST_[workind] * (1 - LAB_[workind]) / LAB_[workind]) * (1 - gamma*(1-sigma))/sigma

            labelas = frisch / workmass

            # Calculate savings elasticity
            ratedev = 0.01
            Market_dev = Market

            # Increase rates of return to HH by ratedev
            #   Note: Cap gains is zero in steady state, so 
            #         return to HH is only on equity + debt.
            Market_dev['equityDividendRates'] = Market['equityDividendRates'] * (1 + ratedev)
            Market_dev['bondDividendRates']   = Market['bondDividendRates']   * (1 + ratedev)

            Dynamic_dev = generate_aggregates(Market_dev, [], [], [], [])

            savelas = (Dynamic_dev['assetsPM'] - Dynamic['assetsPM']) / (Dynamic['assetsPM'] * ratedev)

            # Calculate $GDP/HH
            outperHH = (Dynamic['outs']/Dynamic['pops'])/scenario['modelunit_dollar']

            # Calculate gini
            GiniTable = MomentsGenerator(scenario,DIST['DIST'],Market,OPTs).giniTable
            gini      = GiniTable.model[GiniTable.Gini=='wealth']

            # Save and display elasticities
            to_save = {'captoout': captoout, 'labelas': labelas, 'savelas': savelas, 'outperHH': outperHH, 'gini': gini,
                'beta': beta, 'gamma': gamma, 'sigma': sigma, 'modelunit_dollar': modelunit_dollar, 'bequest_phi_1': bequest_phi_1}
            scipy.io.savemat(os.join.path(save_dir, 'paramsTargets.mat'), to_save)

            print( '\n' )
            for label in [ ['Beta'          , beta              ] ,
                           ['Gamma'         , gamma             ] ,
                           ['Sigma'         , sigma             ] ,
                           ['Model$'        , modelunit_dollar  ] ,
                           ['phi_1'         , bequest_phi_1     ] ]:
                print('\t%-25s= % 7.8f\n' % label)
            
            print( '--------------\n' )
            for label in [['Capital/Output'        , captoout   ] , 
                          ['Labor elasticity'      , labelas    ] ,
                          ['Savings elasticity'    , savelas    ] ,
                          ['Output/HH'             , outperHH   ] ,
                          ['Wealth Gini'           , gini       ]]:
                print('\t%-25s= % 7.4f\n' % label)
            
            print('\n')

        ##        
        # Create completion indicator file
        f = open(os.join.path(save_dir, 'solved'), 'w+')
        f.close()
        
        # Attempt to move scenario solution to its final resting place.
        # (Errors are typically due to multiple simultaneous copy attempts, 
        # which should result in at least one successful copy)
        try: 
            os.rename(save_dir, PathFinder.getCacheDir( scenario ))
        except Exception as e:
            print( 'WARNING! Moving temporary working dir to final cached dir: %s \n' % str(e))
        
        # Release MEX file to avoid locks
        #clear mex; %#ok<CLMEX>
        
    ## 
    # Solve unanticipated shock on transition path
    #
    def solvePolicyShock( scenario, callertag ):
        
        if scenario.isSteady() :
            raise Exception('Scenario cannot be steady-state.' )
        
        # Identify a temporary directory for output which will be copied
        # to the scenario's WorkingDir upon completion.
        (save_dir, callingtag) = PathFinder.getTempWorkingDir(scenario, callertag )

        # Solving 'closed' transition path requires 'open'
        # ModelSolver.solve() does not have the right dependences for shocked paths
        # (need to pass initial_state from 'open' when hardyload runs from
        # 'closed'. Since we don't have that, solve explicitly here.
        if not scenario.isSteady() and not scenario.isOpen():
           ModelSolver.solvePolicyShock( scenario.open(), '' ) 

        # Create baseline path and shocked path 
        scenario1 = scenario.currentPolicy()
        scenario2 = scenario.postShock()
        
        scenario1_generator = lambda: ModelSolver.solve(scenario1, callingtag)
        scenario1_dir       = PathFinder.getCacheDir(scenario1)
        Market1             = hardyload('market.mat'        , scenario1_generator, scenario1_dir)
        Dynamic1            = hardyload('dynamics.mat'      , scenario1_generator, scenario1_dir)
        Distribution1       = hardyload('distribution.mat'  , scenario1_generator, scenario1_dir)
        Decisions1          = hardyload('decisions.mat'     , scenario1_generator, scenario1_dir)

        T = scenario2.TransitionFirstYear - scenario1.TransitionFirstYear
        
        # This is a list of the GovtDebt variables that have a nonstandard
        # size in Market.  They are usually a matrix of (periods, 30) as
        # opposed to a vector or a singleton. Therefore, these have to be
        # separately loaded into Markets.  
        debtVarList = ['effectiveRatesByMaturity', 'debtDistributionByMaturity', 'effectiveRatesByMaturity_next', 'debtDistributionByMaturity_next']
        
        # Find model year for scenario2 t=0 and fetch Market and
        # Dynamic from that year
        
        for p in Market1.keys():
            valuename = p
            # For non-vectors, just copy the whole thing -- these should be
            # single values (or the priceindices struct)
            if valuename in debtVarList:
                vectorLoc = min(Market1[valuename].shape[0],T)
                Market[valuename] = Market1[valuename][vectorLoc,:]
            else:
                if Market1[valuename].shape[1] < T :
                    Market[valuename] = Market1[valuename]
                else:
                    Market[valuename] = Market1[valuename][T]

        for p in Dynamic1.keys():
            valuename = p
            # For non-vectors, just copy the whole thing 
            if Dynamic1[valuename].shape[1] < T :
                Dynamic[valuename] = Dynamic1[valuename]
            else:
                Dynamic[valuename] = Dynamic1[valuename][T]
        
        # Create caps_next and debts_next variables
        Dynamic['caps_next']           = Dynamic1['caps'][ min(T+1, len(Dynamic1['caps'])) - 1 ]
        Dynamic['caps_domestic_next']  = Dynamic1['caps_domestic'][ min(T+1, len(Dynamic1['caps_domestic'])) - 1 ]
        Dynamic['caps_foreign_next']   = Dynamic1['caps_foreign'][ min(T+1, len(Dynamic1['caps_foreign'])) - 1 ]
        Dynamic['debts_next']          = Dynamic1['debts'][ min(T+1, len(Dynamic['debts_next'])) - 1 ]
        Dynamic['debts_domestic_next'] = Dynamic1['debts_domestic'][ min(T+1, len(Dynamic1['debts_domestic'])) - 1 ]
        Dynamic['debts_foreign_next']  = Dynamic1['debts_foreign'][ min(T+1, len(Dynamic1['debts_foreign'])) - 1 ]
        Market['bondDividendRates_next'] = Market1['bondDividendRates'][ min(T+1, len(Market1['bondDividendRates'])) - 1 ]
        Market['effectiveRatesByMaturity_next'] = Market1['effectiveRatesByMaturity'][ min(T+1, len(Market1['effectiveRatesByMaturity'])) - 1, : ]
        Market['debtDistributionByMaturity_next'] = Market1['debtDistributionByMaturity'][ min(T+1, len(Market1['debtDistributionByMaturity'])) - 1, : ]
        
        # Set priceCapital0 so that capgains are correctly calculated
        Market['priceCapital0']        = Market1['equityPrices'][ T ]
        
        # Set average wage history for wage index calculations
        #   Include all history from steady-state, inclusive
        Market['averageWagesHistory']  = np.hstack((Market1['averageWagesHistory'][0], Market1['averageWages']))
        
        # Pull correct time period from DIST, rem: expecting t=1 population
        DIST = Distribution1['DIST'][:,:,:,:,:,T]            
        
        # Package initial_state, solve scenario2
        init_state['Market']       = Market 
        init_state['Dynamic']      = Dynamic 
        init_state['DIST']         = DIST
        init_state['yearSteady']   = scenario1['TransitionFirstYear'] - 1  # year of steady state

        scenario2_generator = lambda: ModelSolver.solve(scenario2, callingtag, init_state)
        scenario2_dir       = PathFinder.getCacheDir(scenario2)
        Market2             = hardyload('market.mat'        , scenario2_generator, scenario2_dir)
        Dynamic2            = hardyload('dynamics.mat'      , scenario2_generator, scenario2_dir)
        Distribution2       = hardyload('distribution.mat'  , scenario2_generator, scenario2_dir)
        Decisions2          = hardyload('decisions.mat'     , scenario2_generator, scenario2_dir)       
        
        Static2             = []
        if not scenario2.isCurrentPolicy():
            Static2         = hardyload('statics.mat', scenario2_generator, scenario2_dir)
        
        # helper function to combine scenario1 files w/ scenario2
        def join_series( J1, J2 ):
            
            # This is a list of the GovtDebt variables that have a nonstandard
            # size in Market.  They are usually a matrix of (periods, 30) as
            # opposed to a vector or a singleton. Therefore, these have to be
            # separately loaded into Markets. Removing the "next" variables
            # for the debt too, as we do not need those anymore.
            debtVarList = ['effectiveRatesByMaturity', 'debtDistributionByMaturity','effectiveRatesByMaturity_next','debtDistributionByMaturity_next']
            
            for p in J1.keys():
                valuename = p
                # Note: Market (though not Dynamic,Static) has some non-vectors, 
                #    also implictly checks length
                if valuename in debtVarList:
                    # because those are stored as 2D matrices and need to
                    # be concatenated differently.
                    #if strfind(valuename, '_next')
                    #    J.(valuename) = J1.(valuename);
                    #else
                    if not '_next' in valuename:
                        J[valuename] = J1[valuename][0:T,:]
                        J2[valuename]
                
                else:
                    if J1[valuename].shape(1) >= T :
                        J[valuename] = np.hstack((J1[valuename][0:T], J2[valuename]))
                    else:
                        if J1[valuename].shape[1] == 1 : 
                            J[valuename] = J1[valuename]
            return J
        
        def join_dist( Dist1, Dist2 ):
            Dist = Dist1['DIST']
            Dist[:,:,:,:,:,T:-1] = Dist2['DIST']
            return Dist
        
        def join_decisions( D1, D2 ):           
            # Get OPTs variables
            for p in D1['OPTs'].keys():
                valuename = p
                decisions['OPTs'][valuename] = D1['OPTs'][valuename]
                decisions['OPTs'][valuename][:,:,:,:,T:-1] = D2['OPTs'][valuename]
            
            # Auxiliary variables for concatenating cohort-based LABs and savings
            T_model1 = D1['OPTs']['LABOR'].shape[4]
            T_model2 = D2['OPTs']['LABOR'].shape[4]
            T        = T_model1 - T_model2
            # Initialize LABs & savings for all cohorts using scenario1
            for i in range(D1['LABs'].shape[0]):
                decisions['LABs'][i]    = D1['LABs'][i]
                decisions['savings'][i] = D1['savings'][i]
            
            # Update LABs & savings from T+1 on for cohorts alive T_model1 periods
            for i in range(T, D1['LABs'].shape[0] - T_model1 + 1):
                decisions['LABs'][i][:,:,:,T:-1]    = D2['LABs'][i-T]
                decisions['savings'][i][:,:,:,T:-1] = D2['savings'][i-T]
            
            # Update LABs & savings for cohorts alive for less than T_model1
            # periods and more than T_model2 periods (born during scenario1)
            gap = T
            for i in range(D1['LABs'].shape[0] - T_model1 + 1, D1['LABs'].shape[0] - T_model2):
                decisions['LABs'][i][:,:,:,(gap-1):-1] = D2['LABs'][i-T]
                decisions['savings'][i][:,:,:,(gap-1):-1] = D2['savings'][i-T]
                gap = gap - 1
            
            # Get LABs & savings for cohorts born during scenario2
            for i in range(D1['LABs'].shape[0] - T_model2, D1['LABs'].shape[0]):
                decisions['LABs'][i]    = D2['LABs'][i-T]
                decisions['savings'][i] = D2['savings'][i-T]                
            
            return decisions
        
        # Combine scenario1 and scenario 2 into single file
        Market    = join_series( Market1 , Market2  )
        Dynamic   = join_series( Dynamic1, Dynamic2 )
        DIST      = join_dist( Distribution1, Distribution2 )
        decisions = join_decisions( Decisions1, Decisions2 )
        Static  = []
        if len(Static2) != 0:
            Static = join_series( Dynamic1, Static2 )
        
        # Save the combined Market, Dynamic
        #   Rem: First build into save_dir, then move
        if os.path.isfile(save_dir):
            shutil.rmtree(save_dir)
        os.mkdir(save_dir)

        scipy.io.savemat(os.join.path(save_dir, 'market.mat'), Market)
        scipy.io.savemat(os.join.path(save_dir, 'dynamics.mat'), Dynamic)
        scipy.io.savemat(os.join.path(save_dir, 'distribution.mat'), DIST)
        scipy.io.savemat(os.join.path(save_dir, 'decisions.mat'), decisions)
        if len(Static) != 0:
            scipy.io.savemat(os.join.path(save_dir, 'statics.mat'), Static)
        
        # Create completion indicator file
        f = open(os.join.path(save_dir, 'solved'), 'w+')
        f.close()
        
        # Attempt to move scenario solution to its final resting place.
        # (Errors are typically due to multiple simultaneous copy attempts, 
        # which should result in at least one successful copy)
        try: 
            os.rename(save_dir, PathFinder.getCacheDir( scenario )) 
        except Exception as e:
            print( 'WARNING! Moving temporary working dir to final cached dir: %s \n' % str(e) )

        return save_dir
    

    ##
    # Create indexes for benefits and taxmax calculations and import CPI index
    #
    #   Inputs:
    #       Market.wages    = T_model-dimension vector, 
    #       Dynamic.labs    = aggregate hours worked
    #       Dynamic.labeffs = aggregate effective labor
    #       nstartyears     = number of cohorts, 
    #       realage_entry   = real age at entry, 
    #       T_model         = number of periods in model run, 
    #       T_life          = maximum life spam,
    def generateIndices( Market, budget, nstartyears, realage_entry, T_steady, T_model, T_life):

        # TBD: FIX!  The timing on nominals is wrong when T_steady <> 0
        index['nominals']        = 1 / budget['deflator']                       # Time-varying reciprocal CPI indexes from CBO
        
        # Build long wage index
        # Index runs from before steady state starts to after model ends
        # and last household there has died
        T                   = T_life + T_steady + T_model - 1 + T_life
        t_1                 = T_life + T_steady                    # t=1 in model
        t_0                 = T_life                               # t of steady state
        wage_index          = np.ones((1,T))
        
        # Overwrite index with historical values, calculated values for model periods
        #  For t<t_steady, wages are what they are in steady
        #  For t>T_model, wages are what they are in T_model (terminal steady state)
        average_wage0           = Market['averageWagesHistory'][0]     # Normalization from steady
        average_wageTmodel      = Market['averageWages'][T_model - 1]
        average_wages           = np.hstack(average_wage0 * np.ones((1,T_life)),          
                                        Market['averageWagesHistory'][1:T_steady],
                                        Market['averageWages'],                     
                                        average_wageTmodel * np.ones((1,T_life)))
        wage_index              = average_wages / average_wage0  

        # Store long index and size/location info
        index['T_steadyStart']     = t_0
        index['T_modelStart']      = t_1
        index['T_modelEnd']        = t_1 + T_model - 1
        index['T_index']           = T                            
        index['wage_inflations']   = wage_index
        
        index['cohort_wages']    = np.zeros((T_model, nstartyears))    # Time- and cohort-varying indexes
        for i in range(nstartyears):
            # Birthyear of i=1   --> T_steady
            # So Birthyear of i  --> T_steady + i - 1
            # So year age60 of i --> T_steady + i - 1 + (60 - realage_entry)
            year_turn60 = T_steady + i - 1 + (60 - realage_entry)
            index['cohort_wages'][:,i] = average_wages[year_turn60] / average_wages[t_1-1:t_1+T_model-1]
        
        return index


##
# Check if file exists and generate before loading if necessary, handling parallel write and read conflicts.
# 
def hardyload(filename, generator, save_dir):
    
    # Check if file exists and generate if necessary
    filepath = os.join.path(save_dir, filename)
    if not os.path.isfile(filepath):
        tagged_dir = generator()
        print( 'INFO: Generated needed dependent scenario to %s \n' % tagged_dir)
    
    # Set maximum number of attempts and maximum pause time in seconds
    maxtries = 200
    maxpause = 2.0
    
    for itry in range(maxtries):
        try:
            # Attempt load
            s = scipy.io.loadmat(filepath)
            break
        except Exception as e:
            if itry == maxtries:
                raise Exception('Failed to load ''%s'' after %d attempts.\n%s' % (filepath, maxtries, str(e)))
            
            # Take a breather
            time.sleep(math.random() * maxpause)
    
    return s