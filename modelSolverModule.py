##
# Dynamic model solver.

from pathFinderModule import PathFinder
from initialGuessModule import InitialGuess
from paramGeneratorModule import ParamGenerator
from firmModule import Firm
from govtDebtModule import GovtDebt
from momentsGeneratorModule import MomentsGenerator
from socialSecurityModule import SocialSecurity
from convergenceModule import Convergence
from helpersModule import makeIterable, checkOneColumn

import os
import shutil
import numpy as np
import scipy
import scipy.optimize as sio
from scipy.sparse import linalg
import math
import time
import warnings
import datetime
import shelve
import copy
import pickle
import random
from numba import jit

class ModelSolver:
    
    # Remove cached Scenario files
    @staticmethod
    def removeCached(scenario):
        
        cacheDir = PathFinder.getCacheDir(scenario)
        if os.path.exists(cacheDir):
            shutil.rmtree(cacheDir)
            
        # For policy shock scenario, remove postShock scenario
        # We leave the baseline scenario
        if scenario.isPolicyShock():
            scenario2 = scenario.postShock()
            cacheDir = PathFinder.getCacheDir(scenario)
            if os.path.exists(cacheDir):
                shutil.rmtree(cacheDir)
                
    ##
    # Solve dynamic model
    #   Inputs:     scenario        -- Scenario to solve
    #               callertag       -- string tag of calling process (if any)
    #               initial_state   -- pre-model state (if any)
    @staticmethod
    def solve(scenario, callertag = '', initial_state = None): ##ok<*FXUP>
        
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
        #@staticmethod
        
        def generate_aggregates(Market, Aggregate, DIST_steady, LABs_static, savings_static, DIST_static):
            # Define dynamic aggregate generation flag
            isdynamic = (len(DIST_static) == 0) or (len(LABs_static) == 0) or (len(savings_static) == 0)
            
            # Set static optimal decision values to empty values for dynamic aggregate generation
            if isdynamic:
                LABs_static = np.array([[] for i in range(nstartyears)])
                savings_static = np.array([[] for i in range(nstartyears)])
                
            # Initialize optimal decision value arrays
            OPTs = {}
            os = ['V', 'LUMP_SUM_TAX', 'LABOR', 'SAVINGS', 'CONSUMPTION', 'CONSUMPTION_TAX', 'AVG_EARNINGS', 'TAXABLE_INC', 'OASI_BENEFITS', 'ORD_LIABILITY', 'PAYROLL_LIABILITY', 'PREF_LIABILITY']
            for o in os:
                OPTs[o] = np.zeros((nz,nk,nb,T_life,T_model))
            
            # Initialize array of cohort optimal labor values
            LABs    = nstartyears * [None] 
            savings = nstartyears * [None]
            
            # Initialize population distribution array
            DIST = np.zeros((nz,nk,nb,T_life,ng,T_model))
            newDIST_trans = np.zeros((nz,nk,nb,T_life,ng,T_model))
            if isSteadyEconomy:
                DIST_trans = np.zeros((nz,nk,nb,T_life,ng,T_model))
            else:
                DIST_trans = np.empty((nz,nk,nb,T_life,ng,T_model))
            
            theSocialSecurity = SocialSecurity(socialsecurity, bv, Market['priceindices'], startyears, realage_entry, T_model)
            
            # Calculate indexed policy variables, use only model periodindices
            ssincmins_indexed   = theSocialSecurity.incomeFloor
            ssincmaxs_indexed   = theSocialSecurity.incomeCeiling
            
            sstax_brackets_indexed  = theSocialSecurity.payrollTaxBrackets
            sstax_rates_indexed     = theSocialSecurity.payrollTaxRates
            sstax_burdens_indexed   = theSocialSecurity.payrollTaxBurdens

            # Package fixed dynamic optimization arguments into anonymous function
            solve_cohort_ = lambda V0, LAB_static, saving_static, T_past, T_shift, T_active, T_works, ssbenefits, cohort_wageindexes: ModelSolver.solve_cohort(
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
                 Market['equityPrices'], np.hstack((makeIterable(Market['equityPrice0']), Market['equityPrices'][0:T_model-1]))  # Equity prices (pm and am)
                 )
                
            # Initialize series of terminal utility values
            V0s = np.zeros((nz,nk,nb,T_life))
            
            # TBD: Solve steady state / post-transition path cohort
            if isdynamic:
                
                ssbenefits = theSocialSecurity.getBenefitsForCohort(-1)

                # Solve dynamic optimization
                # (Note that active time is set to full lifetime)
                OPT = solve_cohort_(V0s[:,:,:,T_life-1], np.array([[[[]]]]), np.array([[[[]]]]), T_pasts[-1], T_shifts[-1], T_life, T_works[-1], ssbenefits, Market['priceindices']['cohort_wages'][:,-1])
                
                # Define series of terminal utility values
                V0s[:,:,:,0:T_life-1] = np.copy(OPT['V'][:,:,:,1:T_life])
            
            if isSteadyEconomy:
                # Store optimal decision values
                for o in os:
                    OPTs[o][:,:,:,:,0] = np.copy(OPT[o])
                LABs[0]    = np.copy(OPT['LABOR'])
                savings[0] = np.copy(OPT['SAVINGS'])
            else:
                # Solve transition path cohorts
                OPTs_cohort = nstartyears * [None]
                
         
                # TBD parallelize for loop
                for i in range(nstartyears): 
                    
                    # Extract terminal utility values
                    V0 = np.copy(V0s[:,:,:,T_ends[i]-1]) ##ok<PFBNS>
                    
                    # Calculate cohort-based year-varying benefits policy
                    ssbenefits = theSocialSecurity.getBenefitsForCohort(i) 
                    # Solve dynamic optimization
                    try:
                        OPTs_cohort[i] = solve_cohort_(V0, LABs_static[i, None, None, None], savings_static[i, None, None, None], T_pasts[i], T_shifts[i], T_actives[i], T_works[i], ssbenefits, Market['priceindices']['cohort_wages'][:,i])
                    except:
                        OPTs_cohort[i] = solve_cohort_(V0, LABs_static[i], savings_static[i], T_pasts[i], T_shifts[i], T_actives[i], T_works[i], ssbenefits, Market['priceindices']['cohort_wages'][:,i])
                        
                    LABs[i]    = np.copy(OPTs_cohort[i]['LABOR'])
                    savings[i] = np.copy(OPTs_cohort[i]['SAVINGS'])
                    
                # Construct optimal decision value arrays
                for i in range(nstartyears):
                    for t in range(T_actives[i]):
                        age  = t + T_pasts[i]
                        year = t + T_shifts[i]
                        for o in os:
                            OPTs[o][:,:,:,age,year] = np.copy(OPTs_cohort[i][o][:,:,:,t])
                  
            if isdynamic:
                if isSteadyEconomy:
                    # Determine steady state age distribution without immigration
                    (_, v) = linalg.eigs(np.vstack((birth_rate[0:T_model]*np.ones(T_life), np.hstack((np.diag(np.squeeze(surv[0:T_model,0:T_life-1])), np.zeros((T_life-1,1)))))), k=1, which='LR')
                    v=v.real
                    DISTage0 = v / np.sum(v)
                    
                    # Define initial population distribution using
                    # steady state distributions over age and
                    # productivity for each of the three groups
                    # (citizens, legals, and illegals).
                    
                    newDIST_next = np.zeros((nz,nk,nb,T_life,ng))
                    newDIST_next[:,:,:,:,g['citizen']] = np.tile(np.reshape(np.tile(np.reshape(sNew['citizen']['initPop'], (1,T_life)), [nz,1]) * newDISTz[:,:,g['citizen']], (nz,1,1,T_life)), [1,nk,nb,1]) / (nk*nb)
                    newDIST_next[:,:,:,:,g['legal']] = np.tile(np.reshape(np.tile(np.reshape(sNew['legal']['initPop'], (1,T_life)), [nz,1]) * newDISTz[:,:,g['legal']], (nz,1,1,T_life)), [1,nk,nb,1]) / (nk*nb)
                    newDIST_next[:,:,:,:,g['illegal']] = np.tile(np.reshape(np.tile(np.reshape(sNew['illegal']['initPop'], (1,T_life)), [nz,1]) * newDISTz[:,:,g['illegal']], (nz,1,1,T_life)), [1,nk,nb,1]) / (nk*nb)
                    
                    # Define initial population distribution using steady state distributions over age and productivity without immigration
                    DIST_next = np.zeros((nz,nk,nb,T_life,ng))
                    DIST_next[:,:,:,:,g['citizen']] = np.tile(np.reshape(np.tile(np.reshape(DISTage0, (1,T_life)), [nz,1]) * DISTz[:,:,g['citizen']], (nz,1,1,T_life)), [1,nk,nb,1]) / (nk*nb)

                    # Specify number of distribution generation years as number of years required to achieve steady state
                    nyears = 153
                else:
                    # Define initial population distribution as steady state distribution
                    DIST_next = np.copy(DIST_steady)

                    # Specify number of distribution generation years as number of model years
                    nyears = T_model
                    
                    # Define initial population distribution as steady state distribution
                    newDIST_next = np.copy(DIST_steady)
                    
                # We are initializing our entry and exit matrices.  These
                # matrices will keep track of how many people enter and
                # leave explicitly.  These are in rates, not in levels.                
                entry          = np.zeros((nz,nk,nb,T_life,ng,nyears))
                exit           = np.zeros((nz,nk,nb,T_life,ng,nyears))
 
                newDIST = np.zeros(newDIST_next.shape + (min(nyears, T_model),))
                
                for year in range(1,nyears+1):
                    
                    # Extract optimal decision values for current year
                    K = np.copy(OPTs['SAVINGS'][:,:,:,:,min(year, T_model)-1])
                    B = np.copy(OPTs['AVG_EARNINGS'][:,:,:,:,min(year, T_model)-1])
                    f = lambda x: np.reshape(x, (nz,1,1,x.shape[1],1))
                    
                    # Initializing the current period distribution and the
                    # mass of agents
                    newDIST_year = np.copy(newDIST_next)
                    newDIST[:,:,:,:,:,min(year, T_model)-1] = np.copy(newDIST_year)
                    newP = np.sum(newDIST_next)
                    
                    # If we're in the first period, we need next period's
                    # birth rates and immigration rates to get the right
                    # population
                    max_Model = T_model
                    if isSteadyEconomy:
                        max_Model = 2

                    # The ENTRY matrix is additive.  generate_distribution
                    # will take this matrix and ADD it to the existing
                    # matrix.
                    # The EXIT matrix is multiplicative.  It works like a
                    # survival rate; we take the existing distribution, and
                    # like survival rates, multiply through by the EXIT
                    # matrix.
                    # The reason for this is that if we ever want to adjust
                    # discount factors for non-death exit, we'll be able to
                    # easily multiply through by the exit matrix to get the
                    # new time value of money.
                    entry[:,0,0,0,g['citizen'], year-1] = np.squeeze(f(newDISTz[:,0,g['citizen'],np.newaxis]) * newP *  sNew['citizen']['birthRate'][min(year+1, max_Model)-1, 0])
                    entry[:,0,0,:,g['legal'], year-1] = np.squeeze(f(newDISTz[:,:,g['legal'],np.newaxis]) * newP * np.tile(np.reshape( sNew['legal']['birthRate'][min(year+1, max_Model)-1,:] + sNew['legal']['immigrationRate'][min(year + 1, max_Model)-1,:], (1,1,1,T_life,1)), [nz,1,1,1,1]))
                    entry[:,0,0,:,g['illegal'], year-1] = np.squeeze(f(newDISTz[:,:,g['illegal'],np.newaxis]) * newP * np.tile(np.reshape( sNew['illegal']['birthRate'][min(year+1, max_Model)-1,:] + sNew['illegal']['immigrationRate'][min(year+1, max_Model)-1,:], (1,1,1,T_life,1)), [nz,1,1,1,1]))
                    entry[:,:,:,42:T_life,:,:] = 0                          # Capping entry at age 62, so people have a chance to earn a work history and accululate assets (else could end up with negative consumption in retirement)
                    exit[:,:,:,:,g['legal'], year-1]    = np.squeeze(np.tile(np.reshape( sNew['legal']['emigrationRate'][min(year, max_Model)-1,:], (1,1,1,T_life,1)), [nz,nk, nb,1,1]))
                    exit[:,:,:,:,g['illegal'], year-1]  = np.squeeze(np.tile(np.reshape( sNew['illegal']['emigrationRate'][min(year, max_Model)-1,:], (1,1,1,T_life,1)), [nz,nk, nb,1,1]))
                    
                    newDIST_next = ModelSolver.generate_distribution_new(newDIST_year, entry[:,:,:,:,:,min(year, T_model)-1], K, B, nz, nk, nb, T_life, ng, np.squeeze(transz[min(year, T_model)-1,:,:,:]), kv, bv, sNew['surv'][min(year, T_model)-1,:],  np.ones((nz,nk,nb,T_life,ng)) - exit[:,:,:,:,:,min(year, T_model)-1])
                    assert (newDIST_next>=0).all(), 'Negative mass of people at DIST_next.'
                    
                    # TBD: Amnesty for new demographic approach. We MAY
                    # want to do this through the microsim model in the
                    # future when we integrate movement among legal
                    # statuses.

                    # Store population distribution for current year
                    DIST_year = np.copy(DIST_next)
                    DIST[:,:,:,:,:,min(year, T_model)-1] = np.copy(DIST_year)

                    # Define population growth distribution
                    DIST_grow     = np.zeros((nz,nk,nb,T_life,ng))
                    P             = np.sum(DIST_year)
                    legal_entry   = np.tile(np.reshape(legal_rate_age[min(year, T_model)-1,:], (1,1,1,T_life,1)), [nz,1,1,1,1])
                    illegal_entry = np.tile(np.reshape(illegal_rate_age[min(year, T_model)-1,:], (1,1,1,T_life,1)), [nz,1,1,1,1])

                    DIST_grow[:,0,0,0,g['citizen']] = np.squeeze(f(DISTz[:,0,g['citizen'],np.newaxis]) * P * birth_rate[min(year, T_model)-1])
                    DIST_grow[:,0,0,:,g['legal']] = np.squeeze(f(DISTz[:,:,g['legal']]) * P * legal_entry)
                    DIST_grow[:,0,0,:,g['illegal']] = np.squeeze(f(DISTz[:,:,g['illegal']]) * P * illegal_entry)

                    # Generate population distribution for next year
                    DIST_next = ModelSolver.generate_distribution(DIST_year, DIST_grow, K, B, nz, nk, nb, T_life, ng,
                                        np.squeeze(transz[min(year, T_model)-1,:,:,:]), kv, bv, surv[min(year, T_model)-1,:])
                    assert (DIST_next>=0).all(), 'Negative mass of people at DIST_next.'

                    # Increase legal immigrant population for amnesty, maintaining distributions over productivity
                    DISTz_legal = DIST_next[:,:,:,:,g['legal']] / np.tile(np.sum(DIST_next[:,:,:,:,g['legal']], axis = 0), [nz,1,1,1])
                    DISTz_legal[np.isnan(DISTz_legal)] = 1/nz
                    
                    DIST_next[:,:,:,:,g['legal']] = DIST_next[:,:,:,:,g['legal']] + np.tile(np.sum(amnesty*DIST_next[:,:,:,:,g['illegal']], axis = 0), [nz,1,1,1]) * DISTz_legal

                    # Reduce illegal immigrant population for amnesty
                    DIST_next[:,:,:,:,g['illegal']] = (1-amnesty) * DIST_next[:,:,:,:,g['illegal']]

                    if isSteadyEconomy:
                        DIST_trans[:,:,:,:,:,0] = np.copy(DIST_next)
                        newDIST_trans[:,:,:,:,:,0] = np.copy(newDIST_next)

            else:
                DIST = np.copy(DIST_static)
                newDIST = np.copy(DIST_static)
            
            newDIST = np.squeeze(newDIST)
            
            # Normalize steady state population distribution
            if isSteadyEconomy:
                DIST_trans = DIST_trans / np.sum(DIST)
                DIST = DIST / np.sum(DIST)
                newDIST = newDIST / np.sum(newDIST)
                
                # We are creating targets for the population by immigration
                # status. Then we will renormalize each cohort of the
                # steady-state population to reach those targets,
                # preserving the distribution of capital, productivity,
                # work history, etc.
                targets = np.transpose(np.vstack((sNew['citizen']['initPop'], sNew['legal']['initPop'], sNew['illegal']['initPop']))/(sum(sNew['citizen']['initPop']) + sum(sNew['legal']['initPop']) + sum(sNew['illegal']['initPop'])))

                for numG in range(ng):
                    for numAge in range(T_life):
                        newDIST[:, :, :, numAge, numG] = newDIST[:, :, :, numAge, numG] * targets[numAge,numG] / np.sum(newDIST[:,:,:,numAge,numG])
                
                # We have to recompute the transition because we've
                # renornalized the age / status distribution to the actual
                # microsim, so next period needs to be updated to reflect
                # these adjustments. We are also recomputing the "entry"
                # matrices b/c they use the "next-period" value of birth
                # and immigration rates.  Exit, by contrast, uses the same
                # value, so we can reuse those matrices.
                f = lambda x: np.reshape(x, (nz,1,1,x.shape[1],1))
                entry[:,0,0,0,g['citizen'], 0] = np.squeeze(f(newDISTz[:,0,g['citizen'],np.newaxis]) *  sNew['citizen']['birthRate'][1, 0])
                entry[:,0,0,:,g['legal'], 0]   = np.squeeze(f(newDISTz[:,:,g['legal']]  ) * np.tile(np.reshape( sNew['legal']['birthRate'][1,:] + sNew['legal']['immigrationRate'][1,:], (1,1,1,T_life,1)), [nz,1,1,1,1]))
                entry[:,0,0,:,g['illegal'], 1] = np.squeeze(f(newDISTz[:,:,g['illegal']]) * np.tile(np.reshape( sNew['illegal']['birthRate'][min(year, T_model-1),:] + sNew['illegal']['immigrationRate'][min(year, T_model-1),:], (1,1,1,T_life,1)), [nz,1,1,1,1]))
                entry[:,:,:,42:T_life,:,:] = 0      # Capping entry at age 62, so people have a chance to earn a work history and accululate assets (else could end up with negative consumption in retirement)
                exit[:,:,:,:,g['legal'], 0]    = np.squeeze(np.tile(np.reshape( sNew['legal']['emigrationRate'][0,:], (1,1,1,T_life,1)), [nz,nk, nb,1,1]))
                exit[:,:,:,:,g['illegal'], 0]  = np.squeeze(np.tile(np.reshape( sNew['illegal']['emigrationRate'][0,:], (1,1,1,T_life,1)), [nz,nk, nb,1,1]))
                
                newDIST_trans = ModelSolver.generate_distribution_new(newDIST, entry[:,:,:,:,:,0], K, B, nz, nk, nb, T_life, ng, np.squeeze(transz[0,:,:,:]), kv, bv, sNew['surv'][0,:],  np.ones((nz,nk,nb,T_life,ng)) - exit[:,:,:,:,:,0])
                assert (newDIST_trans>=0).all(), 'Negative mass of people at DIST_next.'

                # If we're using the new demographic approach, replace the
                # next period distribution from steady state with the new demographics.  
                if scenario.UseNewDemographics:
                    DIST_trans = np.copy(newDIST_trans)
            
            # If we're using the new demographic approach, replace the
            # distribution with the new demographics.  
            if scenario.UseNewDemographics:
                DIST = np.copy(newDIST)
            
            # Generate aggregates
            assert (DIST>=0).all(), 'WARNING! Negative mass of people at DIST.'
            DIST_gs = np.reshape(np.sum(DIST, axis = 4), [nz,nk,nb,T_life,T_model], order = 'F')
                       
            f = lambda F: np.sum(np.reshape(DIST_gs * F, (-1, T_model), order = 'F'), axis=0)
         
            Aggregate['pops']      = f(1)                                                                                        # Population
            Aggregate['lumpSumTaxes'] = f(OPTs['LUMP_SUM_TAX'])
            Aggregate['bequests']  = f(OPTs['SAVINGS'] * np.tile(np.reshape(1-np.transpose(surv), (1,1,1,T_life,T_model)), [nz,nk,nb,1,1]))    # Bequests
            Aggregate['labs']      = f(OPTs['LABOR'])                                                                            # Labor
            Aggregate['labeffs']   = f(OPTs['LABOR'] * np.tile(np.reshape(np.transpose(zs, [2,1,0]), (nz,1,1,T_life,T_model)), [1,nk,nb,1,1]))         # Effective labor
            Aggregate['lfprs']     = f(OPTs['LABOR'] > 0.01) / f(1)                                                           # Labor force participation rate
            Aggregate['incs']      = f(OPTs['TAXABLE_INC'])                                                                  # Income
            Aggregate['pits']      = f(OPTs['ORD_LIABILITY'] + OPTs['PREF_LIABILITY'])                                        # Personal income tax
            Aggregate['ssts']      = f(OPTs['PAYROLL_LIABILITY'])                                                             # Capital income tax
            Aggregate['bens']      = f(OPTs['OASI_BENEFITS'])                                                                 # Social Security benefits
            Aggregate['cons']      = f(OPTs['CONSUMPTION'])                                                                   # Consumption
            Aggregate['constax']   = f(OPTs['CONSUMPTION_TAX'])                                                               # Consumption tax
            Aggregate['assetsAM']  = f(np.tile(np.reshape(kv, (1,nk,1,1,1)), [nz, 1,nb,T_life,T_model]))                      # Assets before re-pricing
            Aggregate['assetsPM']  = (Aggregate['assetsAM'] * (np.ones(T_model) + Market['capgains'])                      # Assets after re-pricing            
                                    * (Market['capsharesAM']/Market['capsharesPM']))                                           # Note: The definition of assetsPM corresponds to beginning of period assets at new policy prices, that is, accounting for eventual capital gains.
            Aggregate['laborIncomes'] = f(                                                                                  # Total labor income
                OPTs['LABOR']                                           
                    * np.reshape(np.transpose(zs, axes=[2,1,0]), [nz,1,1,T_life,T_model])                
                    * np.reshape(Market['wages'], [1,1,1,1,T_model]))
        
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
        
        steadyBaseScenario = np.array([])
        steady_generator = np.array([])
        steady_dir = np.array([])
        DIST0 = np.array([])
        baselineScenario = np.array([])
        base_generator = np.array([])
        base_dir = np.array([])
        openScenario = np.array([])
        open_generator = np.array([])
        open_dir = np.array([])
        
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
        if initial_state != None:
            Market0     = initial_state['Market']
            Dynamic0    = initial_state['Dynamic']
            DIST0       = initial_state['DIST']
            yearSteady  = initial_state['yearSteady']
        else:
            yearSteady      = scenario.TransitionFirstYear - 1
            if not isSteadyEconomy:
                Market0   = ModelSolver.hardyload('market.pkl'      , steady_generator, steady_dir)
                Dynamic0  = ModelSolver.hardyload('dynamics.pkl'    , steady_generator, steady_dir)
                s         = ModelSolver.hardyload('distribution.pkl', steady_generator, steady_dir)
                DIST0     = s['DIST_trans']
            else: # for steady, load initial guess
                guess       = InitialGuess(scenario)
                Market0     = guess.Market
                Dynamic0    = guess.Dynamic
                DIST0       = np.array([])

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
        Dynamic_base = np.array([])
        Market_base = np.array([])
        Dynamic_open = np.array([])
        Market_open = np.array([])
        
        if not isBaseline:
            # Load baseline market and dynamics conditions
            Market_base     = ModelSolver.hardyload('market.pkl'  , base_generator, base_dir)
            Dynamic_base    = ModelSolver.hardyload('dynamics.pkl', base_generator, base_dir)
        if not isSteadyEconomy and not scenario.isOpen():
            # WARNING: This does not work with solvePolicyShock,
            #   ensure these files already exist from solving shocked 'open'
            Market_open     = ModelSolver.hardyload('market.pkl'  , open_generator, open_dir)
            Dynamic_open    = ModelSolver.hardyload('dynamics.pkl', open_generator, open_dir)
        
        
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
        
        T_pasts   = np.array([x if x > 0 else 0 for x in -startyears])                            # Life years before first model year
        T_shifts  = np.array([x if x > 0 else 0 for x in +startyears])                           # Model years before first life year
        T_actives = np.array([x if x < T_model else T_model for x in startyears+T_life]) - T_shifts     # Life years within modeling period
        T_ends    = np.array([x if x < T_life else T_life for x in T_model-startyears])                # Maximum age within modeling period
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
        sNew     = ParamGenerator.demographicsNew( scenario )
        newDISTz = sNew['newDISTz']
        
        # Load Social Security parameters
        socialsecurity  = ParamGenerator.social_security( scenario )
        T_works         = socialsecurity['T_works']      # retirement time per cohort 

        #  Budget: interest rates, expenditures, and debt
        #       Rem -- budget.debttoout only used in steady economy
        budget = ParamGenerator.budget( scenario )

        # Load out-of-model parameters
        s = ParamGenerator.outofmodel( scenario, budget )
        production['TFP'] = s['TFP']
        budget['infraSpending']    = s['infraSpending']
        budget['lumpSumTaxes']     = s['lumpSumTaxes']
        
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
        if os.path.exists(save_dir):
            shutil.rmtree(save_dir)
        os.makedirs(save_dir)

        
        ## Static aggregate generation
        if not isBaseline and not isSteadyEconomy:
            
            # Load baseline optimal decisions and population distribution
            s      = ModelSolver.hardyload('decisions.pkl', base_generator, base_dir)
            LABs_static    = s['LABs']
            savings_static = s['savings']
            
            s      = ModelSolver.hardyload('distribution.pkl', base_generator, base_dir)
            try:
                DIST_static = s['DIST']
            except:
                DIST_static = s
            
            # Generate static aggregates
            # (Intermediary structure used to filter out extraneous fields)
            Static = {}
            
            (Static, _, _, Static_DIST, Static_OPTs, _) = generate_aggregates(
                    Market_base, Static, [], LABs_static, savings_static, DIST_static)

            
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
            
            if scenario.UseStaticDebt:
                (Static['debts'], Static['deficits']) = GovtDebt.calculateStaticDebt( Static, Dynamic_base['debts'][0], budget, T_model)
            else:
                Market0['effectiveRatesByMaturity']   = Market0['effectiveRatesByMaturity_next']
                Market0['debtDistributionByMaturity'] = Market0['debtDistributionByMaturity_next']
                (Static['debts'], Static['deficits'], Market0['bondDividendRates'], Market0['effectiveRatesByMaturity'], Market0['debtDistributionByMaturity']) = GovtDebt.calculateDynamicDebts( Static, Dynamic_base['debts'][0], T_model, budget, Market0['effectiveRatesByMaturity'][0,:], Market0['debtDistributionByMaturity'][0,:], Market0['bondDividendRates'][0]) 
            
            # Total assets
            # Note: tot_assets is a sum of choice variables, those are constant at baseline values
            price               = Market_base['equityPrices']
            Static['tot_assetsAM'] = np.hstack((Market_base['equityPrice0'], price[0:T_model-1])) * Static['caps'] + Static['debts_domestic'] + Static['debts_foreign']
            Static['tot_assetsPM'] = price * Static['caps'] + Static['debts_domestic'] + Static['debts_foreign']
                    
            # Ad field for symmetry w/ Dynamic
            Static['is_converged'] = False
                                
            # Save static aggregates
            with open(os.path.join(save_dir, 'statics.pkl'), 'wb') as handle:
                pickle.dump(Static, handle, protocol=pickle.HIGHEST_PROTOCOL)
            with open(os.path.join(save_dir, 'Static_decisions.pkl'), 'wb') as handle:
                pickle.dump(Static_OPTs, handle, protocol=pickle.HIGHEST_PROTOCOL)
            with open(os.path.join(save_dir, 'Static_distribut.pkl'), 'wb') as handle:
                pickle.dump(Static_DIST, handle, protocol=pickle.HIGHEST_PROTOCOL)        
                
        ## Dynamic aggregate generation
        
        # Set initial guesses (Market_init, Dynamic_init)
        # TBD: Revise how we get initial guesses
        guessSource = 'steady'
        if not isSteadyEconomy:
            guessSource = scenario.OpennessPath
        
        if guessSource == 'steady':
            # rem: steady guess is from guess file
            Market_init     = copy.deepcopy(Market0)
            Dynamic_init    = copy.deepcopy(Dynamic0)
            Market_init['bondDividendRates'] = budget['debtrates'] # Use the long-term average real interest coupon
            
            Market_init['effectiveRatesByMaturity_next'] = np.full((T_model,GovtDebt.maxDuration), np.nan)
            Market_init['debtDistributionByMaturity_next'] = np.full((T_model,GovtDebt.maxDuration), np.nan)
            if not scenario.UseStaticDebt:
                Market_init['bondDividendRates_next'] = (1+sum(GovtDebt.initEffectiveRatesByMaturity * GovtDebt.initDebtDistributionByMaturity)) / budget['steadyStatePriceGrowth'] - 1
                Market_init['effectiveRatesByMaturity_next'][0,:] = GovtDebt.initEffectiveRatesByMaturity
                Market_init['debtDistributionByMaturity_next'][0,:] =  GovtDebt.initDebtDistributionByMaturity
                    
        elif guessSource == 'open':    
            if isBaseline:
                Market_init = {}
                for p in Market0.keys():
                    value0         = Market0[p]
                    Market_init[p] = value0
                    if checkOneColumn(value0):
                        Market_init[p] = value0 * np.ones(T_model)
                
                Dynamic_init = {}
                for p in Dynamic0.keys():
                    value0          = Dynamic0[p]
                    Dynamic_init[p] = value0
                    if checkOneColumn(value0):
                        Dynamic_init[p] = value0 * np.ones(T_model)
            else:
                # Guess is baseline
                Market_init  = copy.deepcopy(Market_base)
                Dynamic_init = copy.deepcopy(Dynamic_base)
               
        else:                
            if isBaseline:
                # Guess is from open_base
                Market_init  = copy.deepcopy(Market_open)
                Dynamic_init = copy.deepcopy(Dynamic_open)
            else:
                # Guess from closed_base
                Market_init  = copy.deepcopy(Market_base)
                Dynamic_init = copy.deepcopy(Dynamic_base)
                
               
        ##       
        # Instantiate the Convergence manager
        
        theConvergence = Convergence( scenario )
        
        # Display header
        print( '%s\n' % scenario.longDescription() )
        print( 'Started at: %s \n' % str(datetime.datetime.now()))
        print( '%s\n' % scenario.shortDescription() )

        # Initialize Market, Dynamic to guesses 
        #   rem: guesses should be right size.
        Market  = copy.deepcopy(Market_init)
        Dynamic = copy.deepcopy(Dynamic_init)
                
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
        
        # Track previous iteration
        prevMarket  = copy.deepcopy(Market)
        prevDynamic = copy.deepcopy(Dynamic)
        
        # Iterate over economy until convergence
        while ( theConvergence.stillConverging() ):
            
            print(' Iteration %2d  ...  RUNNING' % theConvergence.iteration)
                
            # Calculate g'vt budget residuals
            # TBD: Revise where these come from
            prevDynamic['Gtilde']      = 0
            prevDynamic['Ttilde']      = 0
            prevDynamic['Ctilde']      = 0
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

                    prevDynamic['Gtilde'] = budget['outlays_by_GDP'] * baseline_outs - baseline_bens
                    prevDynamic['Ttilde'] = budget['tax_revenue_by_GDP'] * baseline_outs - baseline_revs 
                    prevDynamic['Ctilde'] = np.zeros((1,T_model))

                else:

                    # Set g'vt budget residuals from open economy
                    prevDynamic['Gtilde']      = Dynamic_open['Gtilde']
                    prevDynamic['Ttilde']      = Dynamic_open['Ttilde']
                    prevDynamic['Ctilde']      = Dynamic_open['Ctilde']                
               
            ##
            # Firms sector
            theFirm  = Firm( prevDynamic, prevMarket, taxBusiness, production )
            
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
                Market['priceCapital0']        = theFirm.priceCapital0
                Market['userCostCapital0']     = theFirm.userCostCapital0
                Market['corpLeverageCost']     = theFirm.corpLeverageCost
                Market['passLeverageCost']     = theFirm.passLeverageCost
                
                if not isinstance(Market['equityDividendRates'], np.ndarray):
                    Market['equityDividendRates'] = np.array([Market['equityDividendRates']])
                
                Market['averageWagesHistory']  = Market['averageWages'].item(0)
                Market['worldAfterTaxReturn']  = Market['equityDividendRates'].item(0) * (1 - taxBusiness['rateForeignBusinessIncome'].item(0)) + Market['capgains'].item(0) # Capgains=0 in steady
                Market['investmentToCapital0'] = Dynamic['investment'] / Dynamic['caps']
            
            # Compute price indices
            Market['priceindices']     = ModelSolver.generateIndices(Market,
                  budget, nstartyears, realage_entry, T_steady, T_model, T_life)
            
            # Payouts from firm
            firmDist  = theFirm.distributions()
            
            # 'Price' of assets -- HH own equal shares of both bond & equity funds
            Market['userCostCapital0']     = theFirm.userCostCapital0
            Market['priceCapital0']        = theFirm.priceCapital0
            Market['equityPrice0']         = theFirm.priceCapital0
            Market['equityPrices']         = theFirm.priceCapital()
            Market['equityDividendRates']  = firmDist['equityDividendRates']
            
            Market['corpDividendRates']    = firmDist['corpDividendRates']
            Market['passDividendRates']    = firmDist['passDividendRates']
            Market['passTaxIncRates']      = firmDist['passTaxIncRates']
            
            Market['bondPrice0']           = 1
            Market['bondPrices']           = np.ones(T_model)
            if scenario.UseStaticDebt:
                Market['bondDividendRates']    = np.copy(budget['debtrates']) #rem: dividendrate is per $ of assets
            
            # Save previous iteration capital for closed economy klRatio
            #prev_iter['caps'] = Dynamic['caps']
            
            # Generate dynamic aggregates
            
            (Dynamic, LABs, savings, DIST, OPTs, DIST_trans) = generate_aggregates(Market, prevDynamic, DIST0, [], [], [])
            
            
            
            #load generate_aggregate outputs "manually" from localsave for now instead of recomputing
            
            #LOAD ONLY FOR FIRST PASS THROUGH GENERATE AGGREGATES
            #IF WORKING FILES WERE CREATED THEN comment THIS out
            '''
            shelf = shelve.open('localsave.out')
            
            Dynamic = shelf['Dynamic']
            LABs = shelf['LABs']
            savings = shelf['savings']
            DIST = shelf['DIST']
            OPTs = shelf['OPTs']
            DIST_trans = shelf['DIST_trans']
            
            shelf.close()
            '''
            
            #LOAD ONLY FOR SECOND PASS THROUGH GENERATE AGGREGATES
            #IF WORKING FILES WERE CREATED THEN UNCOMMENT THIS OUT
            '''
            shelf = shelve.open('localsave2.out')
            
            Dynamic = shelf['Dynamic']
            LABs = shelf['LABs']
            savings = shelf['savings']
            DIST = shelf['DIST']
            OPTs = shelf['OPTs']
            DIST_trans = shelf['DIST_trans']
            
            shelf.close()
            '''
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
            # TBD: Add foreigner payments of interest on govdebt
            if isSteadyEconomy:
                # Calculate debt, capital, and output
                # (Numerical solver used due to absence of closed form solution)
                f_debts = lambda outs: budget['debttoout'] * outs
                f_caps  = lambda debts: (Dynamic['assetsPM'] - debts) / Market['equityPrices']
                f_outs  = lambda caps: theFirm.output(caps, Dynamic['labeffs'] )
                f = lambda x: x - np.hstack((f_debts(x[2]), f_caps(x[0]), f_outs(x[1])))
                x_ = scipy.optimize.root(f, np.zeros(3),method='lm')
                x_ = x_.x
                Dynamic['debts'] = x_[0]
                Dynamic['caps']  = x_[1]
                Dynamic['outs']  = x_[2]

                Dynamic['debts_domestic'] = Dynamic['debts']
                Dynamic['debts_foreign']  = np.zeros(T_model)
                Dynamic['deficits']       = np.zeros(T_model)
                Dynamic['caps_domestic']  = Dynamic['caps']
                Dynamic['caps_foreign']   = np.zeros(T_model)
                Dynamic['invest_foreign'] = np.zeros(T_model)
                Dynamic['tot_assetsAM']   = Dynamic['assetsAM']
                Dynamic['tot_assetsPM']   = Dynamic['assetsPM']

                # Calculate income
                Dynamic['labincs'] = Dynamic['labeffs'] * Market['wages']
                Dynamic['capincs'] = Market['MPKs'] * Dynamic['caps']
                Dynamic['GNI']     = Dynamic['outs'] - (Dynamic['corpDividendsForeign'] + Dynamic['passDividendsForeign'])

                # Gross investment in physical capital
                #    and resulting K' for "next" period
                #global DIST
                #global OPTs
                #global DIST_trans
                DIST_gs            = np.reshape(np.sum(DIST, axis=4), (nz,nk,nb,T_life,T_model))
                
                assets_tomorrow    = np.sum(np.reshape(DIST_gs * OPTs['SAVINGS'], (-1, T_model)), axis = 0)
                Dynamic['caps_next']  = (Market['capsharesPM'] * (assets_tomorrow - Dynamic['bequests'])) / Market['equityPrices']
                Dynamic['investment'] = Dynamic['caps_next'] - (1 - depreciation) * Dynamic['caps']

                Dynamic['caps_domestic_next'] = Dynamic['caps_next']
                Dynamic['caps_foreign_next'] = 0

                # "Next" period debt. 
                #   TBD: Revise to get real deficit (rem: none in
                #   steady state)
                Dynamic['debts_next'] = (Dynamic['debts'] * (1 + Market['bondDividendRates'][0])
                                     + 0.035 * Dynamic['outs']) # TBD: hardcoded deficit for now
                Dynamic['debts_foreign_next'] = 0    # steady-state has no foreigners

                # Update guesses for next iteration
                # Note: Dynamic.assets represents current assets at new prices.
                #       Bequests should also be priced according to the new policy.
                #       So we apply today's prices to yesterday's bequests and capshares.
                Market['rhos']      = Dynamic['caps'] / Dynamic['labeffs']
                Market['beqs']      = makeIterable(Dynamic['bequests'] / np.sum(DIST_trans) )  
                Market['invtocaps'] = Dynamic['investment'] / Dynamic['caps']

                # NOTE: capshares/capital calculation potentially
                #   DESTROYS domestic capital. This is because debt
                #   crowds out capital in the HH assets and the g'vt
                #   can create an arbitrarily large increase in debt
                #   from period to period.
                Market['capsharesAM'] = makeIterable((Dynamic['assetsAM'] - Dynamic['debts_domestic']) / Dynamic['assetsAM'])
                Market['capsharesPM'] = (Dynamic['assetsPM'] - Dynamic['debts_domestic']) / Dynamic['assetsPM']
            
            else:
                # Capital(1), Debts(1) come from inital state (K', D')
                Dynamic['caps'][0] = Dynamic0['caps_next']
                Market['effectiveRatesByMaturity'] = np.full((T_model, GovtDebt.maxDuration), float('nan'))
                Market['debtDistributionByMaturity'] = np.full((T_model, GovtDebt.maxDuration), float('nan'))
                if not scenario.UseStaticDebt:
                    Market['effectiveRatesByMaturity'][0,:] = Market0['effectiveRatesByMaturity_next']
                    Market['debtDistributionByMaturity'][0,:] = Market0['debtDistributionByMaturity_next']
                    Market['bondDividendRates'][0] = Market0['bondDividendRates_next']
                
                # Calculate debt and take-up
                #  NOTE: Foreign debt take-up is for next period debt
                #        We cumulate the take-ups to generate foreign
                #        debt holdings
                
                if scenario.UseStaticDebt:
                    (Dynamic['debts'], Dynamic['deficits']) = GovtDebt.calculateStaticDebt( Dynamic, Dynamic0['debts_next'], budget, T_model)
                else:
                    (Dynamic['debts'], Dynamic['deficits'], Market['bondDividendRates'], Market['effectiveRatesByMaturity'],Market['debtDistributionByMaturity']) = (
                    GovtDebt.calculateDynamicDebts( Dynamic, Dynamic0['debts_next'], T_model, budget, Market['effectiveRatesByMaturity'][0], Market['debtDistributionByMaturity'][0], Market['bondDividendRates'][0]))

                debt_foreign_1  = Dynamic0['debts_foreign_next']
                new_debt_issued = Dynamic['deficits'] + (Dynamic['debts'] * Market['bondDividendRates'])
                Dynamic['debts_foreign']   = np.cumsum(np.hstack((debt_foreign_1,                                             
                                            new_debt_issued[0:T_model-1] * international['debtTakeUp'][0:T_model-1])))   
                Dynamic['debts_domestic']  = Dynamic['debts'] - Dynamic['debts_foreign']

                # Calculate capital and output
                # Note: Dynamic.assetsPM represents current assets at new prices.
                (klRatio, _) = theFirm.calculateKLRatio(Market['worldAfterTaxReturn'], prevDynamic['caps'], Dynamic['labeffs'])
                open_econ_caps          = klRatio * Dynamic['labeffs']
                Dynamic['caps_domestic']   = (Dynamic['assetsPM'] - Dynamic['debts_domestic']) / Market['equityPrices']
                Dynamic['caps_foreign']    = (open_econ_caps - Dynamic['caps_domestic']) * international['capitalTakeUp']
                Dynamic['caps_foreign'][0] = Dynamic0['caps_foreign_next']    # Investment was set t=0 -- cannot be changed
                Dynamic['caps']            = Dynamic['caps_domestic'] + Dynamic['caps_foreign']
                Dynamic['outs']            = theFirm.output(Dynamic['caps'], Dynamic['labeffs'])

                # Converge to find Ctilde which closes the D/Y ratio
                Ctilde_error    = float('Inf')
                Ctilde          = np.zeros(T_model)
                while Ctilde_error > 1e-13 :

                    prev_Ctilde = Ctilde

                    # Calculate Ctilde to close debt growth
                    # Rem: This effects outs, so need to converge
                    closure_debttoout   = Dynamic['debts'][closure_year-1]/Dynamic['outs'][closure_year-1]
                    
                    cont_Ctilde         = GovtDebt.calculateFixedDebts(closure_debttoout,   
                                                                       Dynamic['deficits'][closure_year-2:T_model],         
                                                                       Dynamic['outs'][closure_year-2:T_model],             
                                                                       Dynamic['debts'][closure_year-2],
                                                                       Market['bondDividendRates'][closure_year-2:T_model])
                    
                    # Update Ctilde
                    Ctilde = np.append(np.zeros(closure_year-2), cont_Ctilde)
                    Dynamic['Ctilde'] = Ctilde
                  
                    # Recalculate debt and check if D/Y has been fixed
                    #   Note: D/Y for t=ClosureYear is unchanged by Ctilde
                    if scenario.UseStaticDebt:
                        (Dynamic['debts'], Dynamic['deficits']) = GovtDebt.calculateStaticDebt( Dynamic, Dynamic0['debts_next'], budget, T_model)
                    else:
                        (Dynamic['debts'], Dynamic['deficits'], Market['bondDividendRates'], Market['effectiveRatesByMaturity'], Market['debtDistributionByMaturity']) = (
                            GovtDebt.calculateDynamicDebts( Dynamic, Dynamic0['debts_next'], T_model, budget, Market['effectiveRatesByMaturity'][0,:], Market['debtDistributionByMaturity'][0,:], Market['bondDividendRates'][0]))

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
                    too_low_caps = np.where(Dynamic['caps'] <= 0)[0]
                    if len(too_low_caps) != 0 :
                        # Ctilde did not fix debt explosion in time
                        raise Exception('Capital becomes negative at t=%u.', too_low_caps[0])
                    too_high_caps = np.where(Dynamic['caps'] > 1e6)[0]
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
                Dynamic['investment'][T_model-1] = Market['invtocaps'][T_model-2] * Dynamic['caps'][T_model-1]

                # Calculate foreign investment series (for reporting)
                Dynamic['invest_foreign'] = (np.hstack((Dynamic['caps_foreign'][1:T_model], Dynamic['caps_foreign'][T_model-1])) 
                                  - (1 - depreciation) * np.hstack((Dynamic['caps_foreign'][0:T_model-1], Dynamic['caps_foreign'][T_model-2])))

                # Update guesses
                # Note: Dynamic.assets represents current assets at new prices.
                #       Bequests should also be priced according to the new policy.
                #       So we apply today's prices to yesterday's bequests and capshares.
                Market['rhos']      = Dynamic['caps'] / Dynamic['labeffs']
                Market['beqs']      = np.hstack((Dynamic0['bequests'] * (1 + Market0['capsharesPM'] * Market['capgains'][0]),
                             Dynamic['bequests'][0:T_model-1] * (1 + Market['capsharesPM'][0:T_model-1] * Market['capgains'][1:T_model]))) / Dynamic['pops']
                Market['invtocaps'] = Dynamic['investment'] / Dynamic['caps']

                # NOTE: capshares/capital calculation potentially
                #   DESTROYS domestic capital. This is because debt
                #   crowds out capital in the HH assets and the g'vt
                #   can create an arbitrarily large increase in debt
                #   from period to period.
                Market['capsharesAM'] = makeIterable((Dynamic['assetsAM'] - Dynamic['debts_domestic']) / Dynamic['assetsAM'])
                Market['capsharesPM'] = (Dynamic['assetsPM'] - Dynamic['debts_domestic']) / Dynamic['assetsPM']

            # Update Market guesses and check for convergence
            Market      = theConvergence.iterate( Market, prevMarket )
            prevMarket  = copy.deepcopy(Market)
            prevDynamic = copy.deepcopy(Dynamic)
            
            # Erase 'RUNNING' text and then print convergence status
            print('\b\b\b\b\b\b\b')
            print('%s\n' % theConvergence.status())
        
        print( '\nFinished at: %s\n' % str(datetime.datetime.now()))

        Dynamic['is_converged'] = theConvergence.isConverged
        
        # Issue warning if did not converge
        if not Dynamic['is_converged']:
            warnings.warn('Model did not converge.')
        
        # Save market conditions, HH policies, and dynamic aggregates
        with open(os.path.join(save_dir, 'market.pkl'), 'wb') as handle:
            pickle.dump(Market, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open(os.path.join(save_dir, 'dynamics.pkl'), 'wb') as handle:
            pickle.dump(Dynamic, handle, protocol=pickle.HIGHEST_PROTOCOL)
        decisions = {'OPTs': OPTs, 'LABs': LABs, 'savings': savings}
        with open(os.path.join(save_dir, 'decisions.pkl'), 'wb') as handle:
            pickle.dump(decisions, handle, protocol=pickle.HIGHEST_PROTOCOL)
        if isSteadyEconomy: 
            DIST = {'DIST': DIST, 'DIST_trans': DIST_trans}
            with open(os.path.join(save_dir, 'distribution.pkl'), 'wb') as handle:
                pickle.dump(DIST, handle, protocol=pickle.HIGHEST_PROTOCOL)
        else:
            with open(os.path.join(save_dir, 'distribution.pkl'), 'wb') as handle:
                pickle.dump(DIST, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
        targets = {}
        
        ## Elasticity calculation
        if isSteadyEconomy:

            # Calculate capital to output ratio
            captoout = (Dynamic['assetsPM'] - Dynamic['debts']) / Dynamic['outs']
            
            # Calculate labor elasticity
            LAB_  = LABs[0]
            DIST_ = np.sum(DIST['DIST'], axis = 4)

            workind = (LAB_ > 0.01)

            workmass = np.sum(DIST_[workind])
            frisch   = np.sum(DIST_[workind] * (1 - LAB_[workind]) / LAB_[workind]) * (1 - gamma*(1-sigma))/sigma

            laborElasticity = frisch / workmass

            # Calculate savings elasticity
            ratedev = 0.01
            Market_dev = copy.deepcopy(Market)
            # Increase rates of return to HH by ratedev
            #   Note: Cap gains is zero in steady state, so 
            #         return to HH is only on equity + debt.
            Market_dev['equityDividendRates'] = Market['equityDividendRates'] * (1 + ratedev)
            Market_dev['bondDividendRates']   = Market['bondDividendRates']   * (1 + ratedev)
            
            (Dynamic_dev, _, _, _, _, _) = generate_aggregates(Market_dev, Dynamic, [], [], [], [])
            
            '''
            shelf = shelve.open('dd_save.ou')
            globals()['Dynamic_dev'] = shelf['Dynamic_dev']
            shelf.close()
            '''
            
            savingElastcitity = (Dynamic_dev['assetsPM'] - Dynamic['assetsPM']) / (Dynamic['assetsPM'] * ratedev)

            # Calculate $GDP/HH
            outperHH = (Dynamic['outs']/Dynamic['pops'])/scenario.modelunit_dollar

            # Calculate gini
            GiniTable = MomentsGenerator(scenario,DIST['DIST'],Market,OPTs).giniTable()
            wealthGini      = GiniTable.model[GiniTable.Gini=='wealth'].values

            # Save and display elasticities
            targets = {'CapitalToOutput': captoout, 'LaborElasticity': laborElasticity, 'SavingElasticity': savingElastcitity, 'OutputPerHH': outperHH, 'WealthGini': wealthGini}

            print( '\n' )
            for label in [ ('Beta'          , beta              ) ,
                           ('Gamma'         , gamma             ) ,
                           ('Sigma'         , sigma             ) ,
                           ('Model$'        , modelunit_dollar  ) ,
                           ('phi_1'         , bequest_phi_1     ) ]:
                print('\t%-25s= % 7.8f\n' % label)
            
            print( '--------------\n' )
            for label in [('Capital/Output'        , captoout   ) , 
                          ('Labor elasticity'      , laborElasticity    ) ,
                          ('Savings elasticity'    , savingElastcitity    ) ,
                          ('Output/HH'             , outperHH   ) ,
                          ('Wealth Gini'           , wealthGini       )]:
                print('\t%-25s= % 7.4f\n' % label)
            
            print('\n')
        
        
        ##
        # Save scenario info
        with open(os.path.join(save_dir, 'scenario.pkl'), 'wb') as handle:
            to_save = {'targets': targets, 'scenario': vars(scenario)}
            pickle.dump(to_save, handle, protocol=pickle.HIGHEST_PROTOCOL)
  
        # Create completion indicator file
        f = open(os.path.join(save_dir, 'solved'), 'w+')
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
        Market1             = ModelSolver.hardyload('market.pkl'        , scenario1_generator, scenario1_dir)
        Dynamic1            = ModelSolver.hardyload('dynamics.pkl'      , scenario1_generator, scenario1_dir)
        Distribution1       = ModelSolver.hardyload('distribution.pkl'  , scenario1_generator, scenario1_dir)
        Decisions1          = ModelSolver.hardyload('decisions.pkl'     , scenario1_generator, scenario1_dir)

        T = scenario2.TransitionFirstYear - scenario1.TransitionFirstYear
        
        # This is a list of the GovtDebt variables that have a nonstandard
        # size in Market.  They are usually a matrix of (periods, 30) as
        # opposed to a vector or a singleton. Therefore, these have to be
        # separately loaded into Markets.  
        debtVarList = ['effectiveRatesByMaturity', 'debtDistributionByMaturity', 'effectiveRatesByMaturity_next', 'debtDistributionByMaturity_next']
        
        # Find model year for scenario2 t=0 and fetch Market and
        # Dynamic from that year
        Market = {}
        Dynamic = {}
        
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
        init_state = {}
        init_state['Market']       = Market 
        init_state['Dynamic']      = Dynamic 
        init_state['DIST']         = DIST
        init_state['yearSteady']   = scenario1['TransitionFirstYear'] - 1  # year of steady state

        scenario2_generator = lambda: ModelSolver.solve(scenario2, callingtag, init_state)
        scenario2_dir       = PathFinder.getCacheDir(scenario2)
        Market2             = ModelSolver.hardyload('market.pkl'        , scenario2_generator, scenario2_dir)
        Dynamic2            = ModelSolver.hardyload('dynamics.pkl'      , scenario2_generator, scenario2_dir)
        Distribution2       = ModelSolver.hardyload('distribution.pkl'  , scenario2_generator, scenario2_dir)
        Decisions2          = ModelSolver.hardyload('decisions.pkl'     , scenario2_generator, scenario2_dir)       
        
        Static2             = []
        if not scenario2.isCurrentPolicy():
            Static2         = ModelSolver.hardyload('statics.pkl', scenario2_generator, scenario2_dir)
        
        # helper function to combine scenario1 files w/ scenario2
        def join_series( J1, J2 ):
            
            # This is a list of the GovtDebt variables that have a nonstandard
            # size in Market.  They are usually a matrix of (periods, 30) as
            # opposed to a vector or a singleton. Therefore, these have to be
            # separately loaded into Markets. Removing the "next" variables
            # for the debt too, as we do not need those anymore.
            debtVarList = ['effectiveRatesByMaturity', 'debtDistributionByMaturity','effectiveRatesByMaturity_next','debtDistributionByMaturity_next']
            J = {}
            
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
            Dist[:,:,:,:,:,T:] = Dist2['DIST']
            return Dist
        
        def join_decisions( D1, D2 ):  
            global T
            
            # Get OPTs variables
            for p in D1['OPTs'].keys():
                valuename = p
                decisions['OPTs'][valuename] = D1['OPTs'][valuename]
                decisions['OPTs'][valuename][:,:,:,:,T:] = D2['OPTs'][valuename]
            
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
                decisions['LABs'][i][:,:,:,T:]    = D2['LABs'][i-T]
                decisions['savings'][i][:,:,:,T:] = D2['savings'][i-T]
            
            # Update LABs & savings for cohorts alive for less than T_model1
            # periods and more than T_model2 periods (born during scenario1)
            gap = T
            for i in range(D1['LABs'].shape[0] - T_model1 + 1, D1['LABs'].shape[0] - T_model2):
                decisions['LABs'][i][:,:,:,(gap-1):] = D2['LABs'][i-T]
                decisions['savings'][i][:,:,:,(gap-1):] = D2['savings'][i-T]
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
        
        with open(os.path.join(save_dir, 'market.pkl'), 'wb') as handle:
            pickle.dump(Market, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open(os.path.join(save_dir, 'dynamics.pkl'), 'wb') as handle:
            pickle.dump(Dynamic, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open(os.path.join(save_dir, 'distribution.pkl'), 'wb') as handle:
            pickle.dump(DIST, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open(os.path.join(save_dir, 'decisions.pkl'), 'wb') as handle:
            pickle.dump(decisions, handle, protocol=pickle.HIGHEST_PROTOCOL)
        if len(Static) != 0:
            with open(os.path.join(save_dir, 'statics.pkl'), 'wb') as handle:
                pickle.dump(Static, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
        # Create completion indicator file
        f = open(os.path.join(save_dir, 'solved'), 'w+')
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
        index = {}
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
        if not isinstance(Market['averageWagesHistory'], np.ndarray):
            Market['averageWagesHistory'] = np.array([Market['averageWagesHistory']])
        average_wage0           = Market['averageWagesHistory'][0]     # Normalization from steady
        
        average_wageTmodel      = Market['averageWages'].item(T_model - 1)
        average_wages           = np.hstack((average_wage0 * np.ones(T_life),          
                                        Market['averageWagesHistory'][1:T_steady],
                                        Market['averageWages'],                     
                                        average_wageTmodel * np.ones(T_life)))
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
        filepath = os.path.join(save_dir, filename)
        if not os.path.isfile(filepath):
            tagged_dir = generator()
            print( 'INFO: Generated needed dependent scenario to %s \n' % tagged_dir)
    
        # Set maximum number of attempts and maximum pause time in seconds
        maxtries = 200
        maxpause = 2.0
    
        for itry in range(maxtries):
            try:
                # Attempt load
                with open(filepath, 'rb') as handle:
                    s = pickle.load(handle)
                break
            except Exception as e:
                if itry == maxtries:
                    raise Exception('Failed to load ''%s'' after %d attempts.\n%s' % (filepath, maxtries, str(e)))
                    
                # Take a breather
                time.sleep(random.random() * maxpause)
    
        return s
    
    ##
    # Solve dynamic optimization problem for a cohort.
    # 
    
    def solve_cohort(
        V0, LAB_static, saving_static, isdynamic,
        nz, ns, nb, T_past, T_shift, T_active, T_work, T_model,
        zs, transz, sv, bv, beta, gamma, sigma, surv_rates,
        bequest_phi_1, bequest_phi_2, bequest_phi_3,
        ssbenefits, ssincmins, ssincmaxs, sswageindexes,
        sstax_brackets, sstax_burdens, sstax_rates,
        pittax_sscredit, pittax_brackets, pittax_burdens, pittax_rates,
        constax_brackets, constax_burdens, constax_rates,
        lump_sum_taxes,
        preferredCorp_shares, ordinaryCorp_shares, ordinaryPass_shares,
        preftax_brackets, preftax_burdens, preftax_rates,
        capgain_rates, capgain_taxrates,
        beqs,
        wages,
        portfolio_equityshares, corpincome_shares, passincome_shares,
        equity_dividendrates, bond_dividendrates,
        corp_taxableincrates, pass_taxableincrates,
        equity_pricesPM, equity_pricesAM):


        ## Argument verification

        nz_max          = 50
        ns_max          = 50
        nb_max          = 50
        T_max           = 100
        nbrackets_max   = 20
        
        assert V0.dtype == 'float64' and (V0.shape[0] <= nz_max) and (V0.shape[1] <= ns_max) and (V0.shape[2] <= nb_max)
        assert LAB_static.dtype == 'float64' and (LAB_static.shape[0] <= nz_max) and (LAB_static.shape[1] <= ns_max) and (LAB_static.shape[2] <= nb_max) and (LAB_static.shape[3] <= T_max)
        assert saving_static.dtype == 'float64' and (saving_static.shape[0] <= nz_max) and (saving_static.shape[1] <= ns_max) and (saving_static.shape[2] <= nb_max) and (saving_static.shape[3] <= T_max)
        assert isinstance(isdynamic, bool ) #and (isdynamic.shape[0] == 1) and (isdynamic.shape[1] == 1)
        assert isinstance(nz, (float, int)) #and (nz.shape[0] == 1) and (nz.shape[1] == 1)
        assert isinstance(ns, (float, int)) #and (ns.shape[0] == 1) and (ns.shape[1] == 1)
        assert isinstance(nb, (float, int)) #and (nb.shape[0] == 1) and (nb.shape[1] == 1)
        assert isinstance(T_past, (float, int, np.int32)) #and (T_past.shape[0] == 1) and (T_past.shape[1] == 1)
        assert isinstance(T_shift, (float, int, np.int32)) #and (T_shift.shape[0] == 1) and (T_shift.shape[1] == 1)
        assert isinstance(T_active, (float, int, np.int32)) #and (T_active.shape[0] == 1) and (T_active.shape[1] == 1)
        assert isinstance(T_work, (float, int, np.int32)) #and (T_work.shape[0] == 1) and (T_work.shape[1] == 1)
        assert isinstance(T_model, (float, int, np.int32)) #and (T_model.shape[0] == 1) and (T_model.shape[1] == 1)
        assert zs.dtype == 'float64' and (zs.shape[0] <= T_max) and (zs.shape[1] <= T_max) and (zs.shape[2] <= nz_max)
        assert transz.dtype == 'float64' and (transz.shape[0] <= T_max) and (transz.shape[1] <= T_max) and (transz.shape[2] <= nz_max) and (transz.shape[3] <= nz_max)
        assert sv.dtype == 'float64' and (sv.shape[0] <= ns_max) and len(sv.shape) == 1 #and (sv.shape[1] == 1)
        assert bv.dtype == 'float64' and (bv.shape[0] <= nb_max) and len (bv.shape) == 1 #and (bv.shape[1] == 1)
        assert isinstance(beta, float) #and (beta.shape[0] == 1) and (beta.shape[1] == 1)
        assert isinstance(gamma, float) #and (gamma.shape[0] == 1) and (gamma.shape[1] == 1)
        assert isinstance(sigma, float) #and (sigma.shape[0] == 1) and (sigma.shape[1] == 1)
        assert surv_rates.dtype == 'float64' and (surv_rates.shape[0] <= T_max) and (surv_rates.shape[1] <= T_max)

        assert isinstance(bequest_phi_1, (float, int)) #and (bequest_phi_1.shape[0] == 1 ) and (bequest_phi_1.shape[1] == 1 )
        assert isinstance(bequest_phi_2, (float, int)) #and (bequest_phi_2.shape[0] == 1 ) and (bequest_phi_2.shape[1] == 1 )
        assert isinstance(bequest_phi_3, (float, int)) #and (bequest_phi_3.shape[0] == 1 ) and (bequest_phi_3.shape[1] == 1 )
        
        assert ssincmins.dtype == 'float64' and (ssincmins.shape[0] <= T_max) and len(ssincmins.shape) == 1 #and (ssincmins.shape[1] == 1)
        assert ssincmaxs.dtype == 'float64' and (ssincmaxs.shape[0] <= T_max) and len(ssincmaxs.shape) == 1 #and (ssincmaxs[1] == 1)
        assert sswageindexes.dtype == 'float64' and (sswageindexes.shape[0] <= T_max) and len(sswageindexes.shape) == 1 #and (sswageindexes.shape[1] == 1)
        assert ssbenefits.dtype == 'float64' and (ssbenefits.shape[0] <= T_max) and (ssbenefits.shape[1] <= nb_max)

        assert sstax_brackets.dtype == 'float64' and (sstax_brackets.shape[0] <= T_max) and (sstax_brackets.shape[1] <= nbrackets_max)
        assert sstax_burdens.dtype == 'float64' and (sstax_burdens.shape[0] <= T_max) and (sstax_burdens.shape[1] <= nbrackets_max)
        assert sstax_rates.dtype == 'float64' and (sstax_rates.shape[0] <= T_max) and (sstax_rates.shape[1] <= nbrackets_max)

        assert pittax_sscredit.dtype == 'float64' and (pittax_sscredit.shape[0] <= T_max) and len(pittax_sscredit.shape) == 1 #and (pittax_sscredit.shape[1] == 1)
        assert pittax_brackets.dtype == 'float64' and (pittax_brackets.shape[0] <= T_max) and (pittax_brackets.shape[1] <= nbrackets_max)
        assert pittax_burdens.dtype == 'float64' and (pittax_burdens.shape[0] <= T_max) and (pittax_burdens.shape[1] <= nbrackets_max)
        assert pittax_rates.dtype == 'float64' and (pittax_rates.shape[0] <= T_max) and (pittax_rates.shape[1] <= nbrackets_max)

        assert constax_brackets.dtype == 'float64' and (constax_brackets.shape[0] <= T_max) and (constax_brackets.shape[1] <= nbrackets_max)
        assert constax_burdens.dtype == 'float64' and (constax_burdens.shape[0] <= T_max) and (constax_burdens.shape[1] <= nbrackets_max)
        assert constax_rates.dtype in ('float64', 'int64') and (constax_rates.shape[0] <= T_max) and (constax_rates.shape[1] <= nbrackets_max)
        assert lump_sum_taxes.dtype == 'float64' and (lump_sum_taxes.shape[0] <= T_max) and len(lump_sum_taxes.shape) == 1 #and (lump_sum_taxes.shape[1] == 1)

        assert preferredCorp_shares.dtype == 'float64' and (preferredCorp_shares.shape[0] <= T_max) and len(preferredCorp_shares.shape) == 1 #and (preferredCorp_shares.shape[1] == 1)
        assert ordinaryCorp_shares.dtype == 'float64' and (ordinaryCorp_shares.shape[0] <= T_max) and len(ordinaryCorp_shares.shape) == 1 #and (ordinaryCorp_shares.shape[1] == 1)
        assert ordinaryPass_shares.dtype == 'float64' and (ordinaryPass_shares.shape[0] <= T_max) and len(ordinaryPass_shares.shape) == 1 #and (ordinaryPass_shares.shape[1] == 1)

        assert preftax_brackets.dtype == 'float64' and (preftax_brackets.shape[0] <= T_max) and (preftax_brackets.shape[1] <= nbrackets_max)
        assert preftax_burdens.dtype == 'float64' and (preftax_burdens.shape[0] <= T_max) and (preftax_burdens.shape[1] <= nbrackets_max)
        assert preftax_rates.dtype == 'float64' and (preftax_rates.shape[0] <= T_max) and (preftax_rates.shape[1] <= nbrackets_max)

        assert capgain_rates.dtype == 'float64' and (capgain_rates.shape[0] <= T_max) and len(capgain_rates.shape) == 1 
        assert capgain_taxrates.dtype == 'float64' and (capgain_taxrates.shape[0] <= T_max) and len(capgain_taxrates.shape) == 1
        
        assert beqs.dtype == 'float64' and (beqs.shape[0] <= T_max) and len(beqs.shape) == 1
        assert wages.dtype == 'float64' and (wages.shape[0] <= T_max) and len(wages.shape) == 1 #and (wages.shape[1] == 1)
        
        assert portfolio_equityshares.dtype == 'float64' and (portfolio_equityshares.shape[0] <= T_max) and len(portfolio_equityshares.shape) == 1
        assert corpincome_shares.dtype == 'float64' and (corpincome_shares.shape[0] <= T_max) and len(corpincome_shares.shape) == 1
        assert passincome_shares.dtype == 'float64' and (passincome_shares.shape[0] <= T_max) and len(passincome_shares.shape) == 1
        assert equity_dividendrates.dtype == 'float64' and (equity_dividendrates.shape[0] <= T_max) and len(equity_dividendrates.shape) == 1 
        assert bond_dividendrates.dtype == 'float64' and (bond_dividendrates.shape[0] <= T_max) and len(bond_dividendrates.shape) == 1
        assert corp_taxableincrates.dtype == 'float64' and (corp_taxableincrates.shape[0] <= T_max) and len(corp_taxableincrates.shape) == 1
        assert pass_taxableincrates.dtype == 'float64' and (pass_taxableincrates.shape[0] <= T_max) and len(pass_taxableincrates.shape) == 1 
        assert equity_pricesPM.dtype == 'float64' and (equity_pricesPM.shape[0]<= T_max) and len(equity_pricesPM.shape) == 1 
        assert equity_pricesAM.dtype == 'float64' and (equity_pricesAM.shape[0]<= T_max) and len(equity_pricesAM.shape) == 1

        ## Dynamic optimization
        OPT = {}

        # Initialize utility and optimal decision value arrays
        OPT['V']   = np.zeros((nz,ns,nb,T_active))                 # Utility

        OPT['LABOR']                   = np.zeros((nz,ns,nb,T_active))       # Labor level
        OPT['SAVINGS']                 = np.zeros((nz,ns,nb,T_active))       # Savings
        OPT['CONSUMPTION']             = np.zeros((nz,ns,nb,T_active))       # Consumption
        OPT['CONSUMPTION_TAX']         = np.zeros((nz,ns,nb,T_active))       # Consumption tax
        OPT['AVG_EARNINGS']            = np.zeros((nz,ns,nb,T_active))       # Average earnings
        OPT['TAXABLE_INC']             = np.zeros((nz,ns,nb,T_active))       # Taxable income
        OPT['OASI_BENEFITS']           = np.zeros((nz,ns,nb,T_active))       # Social Security benefits

        OPT['ORD_LIABILITY']           = np.zeros((nz,ns,nb,T_active))   # Personal income liability at ordinary rates
        OPT['PREF_LIABILITY']          = np.zeros((nz,ns,nb,T_active))   # Personal income liability at preferred rates
        OPT['PAYROLL_LIABILITY']       = np.zeros((nz,ns,nb,T_active))   # Payroll tax liability
        OPT['LUMP_SUM_TAX']            = np.zeros((nz,ns,nb,T_active))   # Lump sum tax
        
        # Initialize forward-looking utility values
        V_step = np.copy(V0)

        # Specify settings for dynamic optimization subproblems
        # optim_options = optimset('Display', 'off', 'TolFun', 1e-4, 'TolX', 1e-4);
        

        # Solve dynamic optimization problem through backward induction
        for t in range(T_active,0,-1):
    
            # Determine age and year, bounded by modeling period
            age  = int(t + T_past )
            year = int(min(t + T_shift, T_model) )

            # Survival probability by year
            survival_probability = surv_rates[year-1, age-1]
    
            # Wages and Bequests for current year
            beq         = beqs    [year-1]
            wage        = wages   [year-1]
    
            # Social security for current year 
            ssincmax    = ssincmaxs    [year-1]
            ssincmin    = ssincmins    [year-1]
            sswageindex = sswageindexes[year-1]
            ssbenefit   = ssbenefits   [year-1, :]
    
            # Asset params for current year
            portfolio_equityshare = portfolio_equityshares[year-1]
            portfolio_bondshare   = 1 - portfolio_equityshare
            corpincome_share      = corpincome_shares     [year-1]
            passincome_share      = passincome_shares     [year-1]
    
            equity_dividendrate   = equity_dividendrates  [year-1]
            bond_dividendrate     = bond_dividendrates    [year-1]
            corp_taxableincrate   = corp_taxableincrates  [year-1]
            pass_taxableincrate   = pass_taxableincrates  [year-1]
    
            equity_pricePM        = equity_pricesPM       [year-1]
            equity_priceAM        = equity_pricesAM       [year-1]
    
            # Capital gains taxes for this year
            capgain_rate        = capgain_rates         [year-1]
            capgain_taxrate     = capgain_taxrates      [year-1]
    
            # Capital taxation for this year
            preferredCorp_share = preferredCorp_shares  [year-1]
            ordinaryCorp_share  = ordinaryCorp_shares   [year-1]
            ordinaryPass_share  = ordinaryPass_shares   [year-1]
            cappref_brackets    = preftax_brackets      [year-1, :]
            cappref_burdens     = preftax_burdens       [year-1, :]
            cappref_rates       = preftax_rates         [year-1, :]

            # Payroll taxes for this year
            sst_brackets        = sstax_brackets        [year-1, :]
            sst_burdens         = sstax_burdens         [year-1, :]
            sst_rates           = sstax_rates           [year-1, :]
    
            # Personal income taxes for this year
            pit_brackets        = pittax_brackets       [year-1, :]
            pit_burdens         = pittax_burdens        [year-1, :]
            pit_rates           = pittax_rates          [year-1, :]
            pit_sscredit        = pittax_sscredit       [year-1]
            
            # Consumption tax rates for this year
            cons_brackets        = constax_brackets     [year-1, :]
            cons_burdens         = constax_burdens      [year-1, :]
            cons_rates           = constax_rates        [year-1, :]  
    
            # Lump sum tax for this year
            lump_sum_tax         = lump_sum_taxes  [year-1]
            
            # Pre-calculate for speed and conciseness
            bequest_p_1   = beta * (1-survival_probability)* bequest_phi_1;
            bequest_p_2   = bequest_phi_2;
            bequest_p_3   = bequest_phi_3;
            
            
            for ib in range(nb):
                for ic in range(ns):
            
                    # Porfolio holdings in dollars
                    bond_value = sv[ic] * portfolio_bondshare
                    equity_value = sv[ic] * portfolio_equityshare
                    
                    # Actual dollar payments to household as dividend income
                    bond_dividend   = bond_value   * bond_dividendrate
                    equity_dividend = equity_value * equity_dividendrate * (equity_pricePM / equity_priceAM)
                    
                    # Dividend liabilities in dollars
                    corp_taxableinc = equity_value * corp_taxableincrate * (equity_pricePM / equity_priceAM) * corpincome_share
                    pass_taxableinc = equity_value * pass_taxableincrate * (equity_pricePM / equity_priceAM) * passincome_share
            
                    # Retired person
                    if (age > T_work):
                
                        # Labor and income
                        labor   = 0
                        labinc  = 0
                        ssinc   = ssbenefit[ib]
                                
                        # Calculate available resources and tax terms
                        (resources, taxable_inc, ord_liability, payroll_liability, pref_liability) = ModelSolver.calculate_resources(
                            labinc, ssinc,
                            equity_value, bond_value,
                            equity_dividend, bond_dividend, corp_taxableinc, pass_taxableinc,
                            sst_brackets, sst_burdens, sst_rates,
                            pit_sscredit, pit_brackets, pit_burdens, pit_rates,
                            lump_sum_tax,
                            preferredCorp_share, ordinaryCorp_share, ordinaryPass_share,
                            cappref_brackets, cappref_burdens, cappref_rates,
                            capgain_rate,
                            beq )
                        
                
                        if isdynamic:
                    
                            # Calculate expected value conditional on living using forward-looking 
                            #   utility values. Pre-multiply by prob. survival and
                            #   beta to save on computation.
                            EV = survival_probability*beta*np.reshape(V_step[0,:,:], (ns,nb))
                            
                            # Solve dynamic optimization subproblem
                            f = lambda s: ModelSolver.value_retirement(
                                    float(s), resources, EV[:,ib], 
                                    sv,
                                    cons_brackets, cons_burdens, cons_rates,
                                    bequest_p_1, bequest_p_2, bequest_p_3,
                                    sigma, gamma)
                            
                            (s, v, _, _, _) = sio.fmin(f, sv[ic], xtol = 1e-4, ftol = 1e-4, full_output = True, disp = False)
                            s = s[0]
                            
                            # Checks -> only work in the absence of mex file!
                            assert not math.isinf(v), 'v is inf'
                            assert s <= sv[-1], 's (s_next) is too big!'
                            
                            # Record utility and optimal decision values
                            #   s (savings) is set by the optimizer
                            #   v is also from optimizer
                            v = -v  # Rem: flipped for minimization
                            
                        else: # STATIC
                    
                            s     = saving_static[0,ic,ib,t-1]  # First dimension is productivity shock (doesn't matter for retirees)
                            labor = 0
                            v     = float('nan')                      # Utility is not properly defined since consumption can be negative
                    
                        # Record utility, decisions, and other values
                        OPT['V']                [:,ic,ib,t-1] = v
                        OPT['LABOR']            [:,ic,ib,t-1] = labor
                        OPT['SAVINGS']          [:,ic,ib,t-1] = s
                        OPT['CONSUMPTION']      [:,ic,ib,t-1] = resources - s - ModelSolver.find_consumption_tax( resources - s, cons_brackets, cons_burdens, cons_rates )
                        OPT['CONSUMPTION_TAX']  [:,ic,ib,t-1] = ModelSolver.find_consumption_tax( resources - s, cons_brackets, cons_burdens, cons_rates )
                        OPT['AVG_EARNINGS']     [:,ic,ib,t-1] = bv[ib]
                        OPT['TAXABLE_INC']      [:,ic,ib,t-1] = taxable_inc
                        OPT['OASI_BENEFITS']    [:,ic,ib,t-1] = ssinc
                        OPT['ORD_LIABILITY']    [:,ic,ib,t-1] = ord_liability
                        OPT['PREF_LIABILITY']   [:,ic,ib,t-1] = pref_liability
                        OPT['PAYROLL_LIABILITY'][:,ic,ib,t-1] = payroll_liability
                        OPT['LUMP_SUM_TAX']     [:,ic,ib,t-1] = lump_sum_tax
                
                    else:
                        # Working age person
                        
                        # Create local instance of average earnings calculation function with fixed parameters
                        calculate_b_ = lambda labinc: ModelSolver.calculate_b(
                                        labinc, age, bv[ib],
                                        ssincmin, ssincmax, sswageindex)
                
                        # Social Security income
                        ssinc = 0

                        # Create local instance of resource calculation function with fixed parameters
                        calculate_resources_ = lambda labinc: ModelSolver.calculate_resources( 
                                                labinc, ssinc, 
                                                equity_value, bond_value,
                                                equity_dividend, bond_dividend, corp_taxableinc, pass_taxableinc, 
                                                sst_brackets, sst_burdens, sst_rates, 
                                                pit_sscredit, pit_brackets, pit_burdens, pit_rates, 
                                                lump_sum_tax,
                                                preferredCorp_share, ordinaryCorp_share, ordinaryPass_share, 
                                                cappref_brackets, cappref_burdens, cappref_rates, 
                                                capgain_rate, 
                                                beq)
                      
                            
                        for iz in range(nz):
                    
                            # Calculate effective wage
                            wage_eff = wage * zs[year-1,age-1,iz]
                    
                            if isdynamic:
                        
                                # Calculate expected value conditional on living using forward-looking 
                                #   utility values. Pre-multiply by prob. survival and
                                #   beta to save on computation.
                                squeezed = np.squeeze(transz[year-1,age-1,iz,:])[:, np.newaxis, np.newaxis]
                                EV = survival_probability * beta * np.reshape(np.sum(np.tile(squeezed, [1,ns,nb]) * V_step, axis=0), (ns,nb))
                                
                                # Solve dynamic optimization subproblem
                                labor0 = 0.5
                                s0   = max(sv[ic], min(sv[-1], 0.1 * wage_eff * labor0))   # Assumes taxation will not exceed 90% of labor income and at the same time forces k to be in the grid
                        
                                f = lambda x: ModelSolver.value_working(
                                         x, EV,
                                         sv, 
                                         bv, wage_eff,
                                         cons_brackets, cons_burdens, cons_rates,
                                         bequest_p_1, bequest_p_2, bequest_p_3,
                                         sigma, gamma, 
                                         calculate_b_, calculate_resources_ )
                                (x, v, _, _, _) = sio.fmin(f, [s0, labor0], xtol = 1e-4, ftol = 1e-4, full_output = True, disp = False)
                        
                                s     = x[0]
                                labor = x[1]       
                        
                                # Checks -> only work in the absence of mex file!
                                assert not math.isinf(v)   , 'v is inf'
                                assert s <= sv[-1], 's (s_next) is too big!'
                        
                                # Record utility and optimal decision values
                                #   s (savings) is set by the optimizer
                                #   v is also from optimizer
                                v = -v  # Rem: flipped for minimization
                        
                            else:   # STATIC
                        
                                s     = saving_static[iz,ic,ib,t-1]
                                labor = LAB_static[iz,ic,ib,t-1]
                                v     = float('nan')                       # Utility is not properly defined since consumption can be negative
                    
                            # Calculate resources
                            labinc      = wage_eff * labor
                            (resources, taxable_inc, ord_liability, payroll_liability, pref_liability) = calculate_resources_(labinc)

                            # Record utility, decisions, and other values
                            OPT['V']                [iz,ic,ib,t-1] = v
                            OPT['LABOR']            [iz,ic,ib,t-1] = labor
                            OPT['SAVINGS']          [iz,ic,ib,t-1] = s
                            OPT['CONSUMPTION']      [iz,ic,ib,t-1] = resources - s - ModelSolver.find_consumption_tax( resources - s, cons_brackets, cons_burdens, cons_rates )
                            OPT['CONSUMPTION_TAX']  [iz,ic,ib,t-1] = ModelSolver.find_consumption_tax( resources - s, cons_brackets, cons_burdens, cons_rates )
                            OPT['AVG_EARNINGS']     [iz,ic,ib,t-1] = calculate_b_(labinc)
                            OPT['TAXABLE_INC']      [iz,ic,ib,t-1] = taxable_inc
                            OPT['OASI_BENEFITS']    [iz,ic,ib,t-1] = 0
                            OPT['ORD_LIABILITY']    [iz,ic,ib,t-1] = ord_liability
                            OPT['PREF_LIABILITY']   [iz,ic,ib,t-1] = pref_liability
                            OPT['PAYROLL_LIABILITY'][iz,ic,ib,t-1] = payroll_liability
                            OPT['LUMP_SUM_TAX']     [iz,ic,ib,t-1] = lump_sum_tax
                    
            # Update forward-looking utility values
            V_step = OPT['V'][:,:,:,t-1]
         
        return OPT
    
    # Retirement age value function
    def value_retirement(
        s, resources, EV_ib,
        sv,
        cons_brackets, cons_burdens, cons_rates,
        bequest_p_1, bequest_p_2, bequest_p_3,
        sigma, gamma):

        # TBD Enforce function inlining for C code generation
        # coder.inline('always');

        # Calculate consumption net of CONSUMPTION tax
        #   Since con_lump = con_rate = 0, we don't return con_liability
        #       where con_liability = con_lump + con_rate * consumption;
        #   consumption + consumption_taxes = total money available after savings
    
        consumption = resources - s - ModelSolver.find_consumption_tax( resources - s, cons_brackets, cons_burdens, cons_rates )
        
        # Perform bound checks
        if sv[0] <= s and s <= sv[-1] and 0 <= consumption:

            # Residual value of bequest.
            # NOTE: (1) bequest is assets chosen for next period,
            #       (2) bequest_p_1 is beta*prob_death*bequest_phi_1
            value_bequest = bequest_p_1 * (1 + s/bequest_p_2)**(1-bequest_p_3)

            # Calculate utility
            v = ((consumption**(gamma*(1-sigma)))*(1/(1-sigma))     # flow utility
                + scipy.interpolate.interp1d(sv, EV_ib)(s)          # continuation value of life
                + value_bequest)                                    # value of bequest

            # Negate utility for minimization and force to scalar for C code generation
            v = -v

        else:
            v = float('Inf')

        return v
    
    # Working age value function
    def value_working(
        x, EV,
        sv,
        bv, wage_eff, 
        cons_brackets, cons_burdens, cons_rates,
        bequest_p_1, bequest_p_2, bequest_p_3,
        sigma, gamma,
        calculate_b_, calculate_resources_ ):

        # TBD Enforce function inlining for C code generation
        # coder.inline('always');

        # Define decision variables and perform bound checks
        s   = x[0]
        lab = x[1]
        
        labinc = wage_eff * lab

        if not ((0 <= lab) and (lab <= 1) and (sv[0] <= s) and (s <= sv[-1])):
            v = float('inf')
            return v

        b = calculate_b_(labinc)

        # Calculate available resources
        (resources, _, _, _, _) = calculate_resources_(labinc)

        # Calculate consumption and perform bound check
        # Since con_lump = con_rate = 0, we don't return con_liability
        #       where con_liability = con_lump + con_rate * consumption;
        #   consumption + consumption_taxes = total money available after savings
        #consumption = (resources - con_lump - s) / (1 + con_rate);
    
        consumption = resources - s - ModelSolver.find_consumption_tax( resources - s, cons_brackets, cons_burdens, cons_rates )
        
        if not (0 <= consumption):
            v = float('inf')
            return v

        # Residual value of bequest.
        # NOTE: (1) bequest is assets chosen for next period,
        #       (2) bequest_p_1 is beta*prob_death*bequest_phi_1
        value_bequest = bequest_p_1 * (1 + s/bequest_p_2) ** (1-bequest_p_3)
        
        # Calculate utility
        v = ((((consumption**gamma)*((1-lab)**(1-gamma)))**(1-sigma))*(1/(1-sigma))           # flow utility
            + scipy.interpolate.interp2d(sv, bv, np.transpose(EV), fill_value = float('nan'))(s, b)       # continuation value of life
            + value_bequest )                                                                 # value of bequest
        
        # Negate utility for minimization and force to scalar for C code generation
        v = -v[0]
        
        return v
    
    # Average earnings calculation function
    def calculate_b(
        labinc, age, bv_ib,
        ssincmin, ssincmax, sswageindex):

        # Enforce function inlining for C code generation
        # coder.inline('always')

        # Calculate average earnings, cap them to a policy ceiling, and deflate
        # them at each period using Market.priceindices.cohort_wages(:,i) = sswageindex
        # for household of cohort i. 
        if labinc > (ssincmin + 10*np.spacing(1)):
            b = (bv_ib*(age-1) + sswageindex*min(labinc, ssincmax)) / age - 10*np.spacing(1)
        else:
            b = bv_ib
        
        return b
    
    # Resource and tax calculation function
    def calculate_resources(labinc, ssinc,
                            equity_value, bond_value,
                            equity_dividend, bond_dividend, corp_dividend, pass_taxableinc,
                            sst_brackets, sst_burdens, sst_rates,
                            pit_sscredit, pit_brackets, pit_burdens, pit_rates,
                            lump_sum_tax,
                            preferredCorp_share, ordinaryCorp_share, ordinaryPass_share,
                            cappref_brackets, cappref_burdens, cappref_rates,
                            capgain_rate,
                            beq):

        # Enforce function inlining for C code generation
        # coder.inline('always');

        # Cap gains are ONLY for equityfund for now
        #   NOTE that capgain_rate is in form (p_t - p_(t-1))/p_(t-1)
        equitycapgain   = equity_value * capgain_rate 
    
        # Calculate PIT taxable income
        #   We do not allow negative incomes
        pit_inc = max( 0,                            
                      ordinaryPass_share * pass_taxableinc 
                      + ordinaryCorp_share * corp_dividend   
                      + bond_dividend                        
                      + (1 - pit_sscredit)*ssinc             
                      + labinc)
        ord_liability = ModelSolver.find_tax_liability( pit_inc, pit_brackets, pit_burdens, pit_rates )

        # Calculate Social Security tax from wage income
        payroll_liability = ModelSolver.find_tax_liability( labinc, sst_brackets, sst_burdens, sst_rates );

        # Calculate preferred rates tax 
        prefinc         = max(0, preferredCorp_share * corp_dividend)
        pref_liability  = ModelSolver.find_tax_liability( prefinc, cappref_brackets, cappref_burdens, cappref_rates )
    
        # Since capgain_taxrate = 0, we don't add capgain_liability to any liability series
        # capgain_liability = equitycapgain * capgain_taxrate
    
        # Calculate available resources
        resources = (labinc + ssinc + beq                              # Labor income + Social Security benefits + bequests
                    + equity_value + bond_value                        # Asset holdings
                    + equity_dividend + bond_dividend                  # Income from assets
                    + equitycapgain                                    # Capital gains due to price changes
                    - (ord_liability + payroll_liability + pref_liability + lump_sum_tax)) # Tax liability without consumption taxes

        return (resources, pit_inc, ord_liability, payroll_liability, pref_liability)
    
    ##
    #  Helper function to find tax liability from brackets & rates
    #       NOTE:   income, burdens, & brackets are in US dollars; 
    #           calculated liability is in also in US dollars.
    #           rates apply for income between brackets(i-1) and brackets(i)
    #           burdens(i) are pre-calculated total tax liability at brackets(i)
    #       IMPORTANT:  Expect 
    #                   (1) equal-size vectors with brackets(1)=0
    #                   (2) brackets are in ascending order
    #                   (3) rates, burdens match brackets
    def find_tax_liability( income, brackets, burdens, rates ):

        # Enforce function inlining for C code generation
        # coder.inline('always');

        # Linear search since assume relatively small size vectors
        numbrackets = len(brackets)
        thebracket  = 0
        while (thebracket < numbrackets) and (brackets[thebracket] <= income) :
            thebracket = thebracket + 1
        
        thebracket = thebracket - 1
    
        tax = burdens[thebracket] + rates[thebracket]*(income - brackets[thebracket])

        return tax

    def find_consumption_tax( consumptionInclTax, brackets, burdens, rates ):
    #     INPUTS: consumptionInclTax = resources minus savings (consumption + taxes)
    #             rates apply for consumptionInclTax between brackets(i-1) and brackets(i)
    #             burdens(i) are pre-calculated total tax liability at brackets(i)


        # Linear search since assume relatively small size vectors
        numbrackets = len(brackets)
        # Initialization
        thebracket  = 0
        tmpTax = 0
        actualConsumption = 10000000000
    
        # Find the tax associated with each bracket
        #     Loop through each bracket, computes tax & consumption consistent
        #     with that bracket, checks if consumption is inside that bracket.
        #     If not, moves to the next bracket
        while (thebracket < numbrackets) and (actualConsumption > brackets[thebracket]):
            tmpTax = (burdens[thebracket] + rates[thebracket]*(consumptionInclTax - brackets[thebracket])) / (1 + rates[thebracket]) 
            actualConsumption = consumptionInclTax - tmpTax
            thebracket = thebracket + 1
    
        tax = tmpTax

        return tax

    ##
    # Generate population distribution for next year from population distribution for current year.
    # 
    ##
    @jit
    def generate_distribution(DIST_year, DIST_grow, K, B, nz, nk, nb, T_life, ng, transz, kv, bv, surv): #codegen

        ## Argument verification

        nz_max =  50
        nk_max =  50
        nb_max =  50
        T_max  = 100
        ng_max =  10
        
        assert DIST_year.dtype == 'float64' and (DIST_year.shape[0]    <= nz_max  ) and (DIST_year.shape[1]    <= nk_max  ) and (DIST_year.shape[2]     <= nb_max  ) and (DIST_year.shape[3]     <= T_max   ) and (DIST_year.shape[4]    <= ng_max  ) 
        assert DIST_grow.dtype == 'float64' and (DIST_grow.shape[0]    <= nz_max  ) and (DIST_grow.shape[1]    <= nk_max  ) and (DIST_grow.shape[2]     <= nb_max  ) and (DIST_grow.shape[3]     <= T_max   ) and (DIST_grow.shape[4]    <= ng_max  ) 
        assert K.dtype == 'float64' and (K.shape[0]            <= nz_max  ) and (K.shape[1]            <= nk_max  ) and (K.shape[2]             <= nb_max  ) and (K.shape[3]             <= T_max   ) 
        assert B.dtype == 'float64' and (B.shape[0]            <= nz_max  ) and (B.shape[1]            <= nk_max  ) and (B.shape[2]             <= nb_max  ) and (B.shape[3]             <= T_max   ) 
        assert isinstance(nz, (float, int)) # and (nz.shape[0]           == 1       ) and (nz.shape[1]            == 1      )
        assert isinstance(nk, (float, int)) # and (nk.shape[0]           == 1       ) and (nk.shape[1]            == 1      )
        assert isinstance(nb, (float, int)) # and (nb.shape[0]           == 1       ) and (nb.shape[1]            == 1      )
        assert isinstance(T_life, (float, int)) # and (T_life.shape[0]       == 1       ) and (T_life.shape[1]        == 1      )
        assert isinstance(ng, (float, int)) # and (ng.shape[0]           == 1       ) and (ng.shape[1]            == 1      )
        assert transz.dtype == 'float64' and (transz.shape[0]       <= T_max   ) and (transz.shape[1]       <= nz_max  ) and (transz.shape[2]        <= nz_max  )
        assert kv.dtype == 'float64' and (kv.shape[0]           <= nk_max  ) and (len(kv.shape)            == 1      )
        assert bv.dtype == 'float64' and (bv.shape[0]           <= nb_max  ) and (len(bv.shape)           == 1      )
        assert surv.dtype == 'float64' and (surv.shape[0]         <= T_max       ) and (len(surv.shape)         <= T_max   )



        ## Distribution generation

        # Initialize population distribution for next year using population growth distribution
        DIST_next = np.copy(DIST_grow)

        for age in range(1,T_life):
            
            # Extract optimal decision values
            k_t = np.copy(K[:,:,:,age-1])
            b_t = np.copy(B[:,:,:,age-1])
            k_shape = k_t.shape
            b_shape = b_t.shape
            k_t_flat = k_t.flatten(order = 'F')
            b_t_flat = b_t.flatten(order = 'F')
    
            # Find indices of nearest values in decision value discretization vectors
            jk_lt = np.ones(len(k_t_flat))
            first = kv[0:-1]
            for elem in range(len(k_t_flat)):
                second = k_t_flat[elem]
                jk_lt[elem] = np.where(first <= second)[0][-1]
                
            jk_lt = np.reshape(jk_lt, k_shape, order = 'F').astype('int')
            jk_gt = jk_lt + 1
            
            jb_lt = np.ones(len(b_t_flat))
            first = bv[0:-1]
            for elem in range(len(b_t_flat)):
                second = b_t_flat[elem]
                jb_lt[elem] = np.where(first <= second)[0][-1]
            
            jb_lt = np.reshape(jb_lt, b_shape, order = 'F').astype('int')
            jb_gt = jb_lt + 1
            
            # Calculate linear weights for nearest values
            wk_lt = (kv[jk_gt] - k_t) / (kv[jk_gt] - kv[jk_lt])
            wk_gt = 1 - wk_lt
    
            wb_lt = (bv[jb_gt] - b_t) / (bv[jb_gt] - bv[jb_lt])
            wb_gt = 1 - wb_lt
            
            # Checks -> only work in the absence of mex file!
            assert (wk_lt>=0).all() and (wk_gt>=0).all() and (wb_lt>=0).all() and (wb_gt>=0).all(), 'Negative weights to compute households distribution.'     
    
            for jz in range(nz):
        
                # Apply survival and productivity transformations to population distribution for current year
                DIST_transz = DIST_year[:,:,:,age-1,:] * surv[age-1] * np.tile(np.reshape(transz[age,:,jz], (nz,1,1,1)), [1,nk,nb,ng])
                assert (DIST_transz>=0).all(), 'Negative mass of people at DIST_transz.'
        
                # Redistribute population distribution from current year to next year according to target indices and weights
                for ib in range(nb):
                    for ik in range(nk):
                        for iz in range(nz): #ok<ALIGN>
  
                          DIST_next[jz, jk_lt[iz,ik,ib], jb_lt[iz,ik,ib], age, :] = DIST_next[jz, jk_lt[iz,ik,ib], jb_lt[iz,ik,ib], age, :] + wk_lt[iz,ik,ib]*wb_lt[iz,ik,ib]*DIST_transz[iz,ik,ib,:]
                          DIST_next[jz, jk_gt[iz,ik,ib], jb_lt[iz,ik,ib], age, :] = DIST_next[jz, jk_gt[iz,ik,ib], jb_lt[iz,ik,ib], age, :] + wk_gt[iz,ik,ib]*wb_lt[iz,ik,ib]*DIST_transz[iz,ik,ib,:]
                          DIST_next[jz, jk_lt[iz,ik,ib], jb_gt[iz,ik,ib], age, :] = DIST_next[jz, jk_lt[iz,ik,ib], jb_gt[iz,ik,ib], age, :] + wk_lt[iz,ik,ib]*wb_gt[iz,ik,ib]*DIST_transz[iz,ik,ib,:]
                          DIST_next[jz, jk_gt[iz,ik,ib], jb_gt[iz,ik,ib], age, :] = DIST_next[jz, jk_gt[iz,ik,ib], jb_gt[iz,ik,ib], age, :] + wk_gt[iz,ik,ib]*wb_gt[iz,ik,ib]*DIST_transz[iz,ik,ib,:]

        return DIST_next
    
    
    ## Generate population distribution for next year from population distribution for current year.
    #
    @jit
    def generate_distribution_new(DIST_year, DIST_grow, K, B, nz, nk, nb, T_life, ng, transz, kv, bv, surv, stay): #codegen
        
        ## Argument verification
        
        nz_max =  50
        nk_max =  50
        nb_max =  50
        T_max  = 100
        ng_max =  10
        
        assert DIST_year.dtype == 'float64' and (DIST_year.shape[0]     <= nz_max  ) and (DIST_year.shape[1]     <= nk_max  ) and (DIST_year.shape[2]     <= nb_max  ) and (DIST_year.shape[3]     <= T_max   ) and (DIST_year.shape[4]     <= ng_max  ) 
        assert DIST_grow.dtype == 'float64' and (DIST_grow.shape[0]     <= nz_max  ) and (DIST_grow.shape[1]     <= nk_max  ) and (DIST_grow.shape[2]     <= nb_max  ) and (DIST_grow.shape[3]     <= T_max   ) and (DIST_grow.shape[4]     <= ng_max  ) 
        assert K.dtype == 'float64' and (K.shape[0]             <= nz_max  ) and (K.shape[1]             <= nk_max  ) and (K.shape[2]             <= nb_max  ) and (K.shape[3]             <= T_max   ) 
        assert B.dtype == 'float64' and (B.shape[0]             <= nz_max  ) and (B.shape[1]             <= nk_max  ) and (B.shape[2]             <= nb_max  ) and (B.shape[3]             <= T_max   ) 
        
        assert isinstance(nz          , (int, float) ) #and (nz.shape[0]           == 1       ) and (nz.shape[1]            == 1       ) 
        assert isinstance(nk          , (int, float) ) #and (nk.shape[0]           == 1       ) and (nk.shape[1]            == 1       ) 
        assert isinstance(nb          , (int, float) ) #and (nb.shape[0]           == 1       ) and (nb.shape[1]            == 1       ) 
        assert isinstance(T_life      , (int, float) ) #and (T_life.shape[0]       == 1       ) and (T_life.shape[1]        == 1       ) 
        assert isinstance(ng          , (int, float) ) #and (ng.shape[0]           == 1       ) and (ng.shape[1]            == 1       ) 
        assert transz.dtype == 'float64' and (transz.shape[0]       <= T_max   ) and (transz.shape[1]        <= nz_max  ) and (transz.shape[2]      <= nz_max  )
        assert kv.dtype == 'float64' and (kv.shape[0]           <= nk_max  ) and (len(kv.shape)            == 1       ) 
        assert bv.dtype == 'float64' and (bv.shape[0]           <= nb_max  ) and (len(bv.shape)            == 1       )
        assert surv.dtype == 'float64' and (len(surv.shape)         == 1       ) and (surv.shape[0]          <= T_max   ) 
        assert stay.dtype == 'float64' and (stay.shape[0]         <= nz_max  ) and (stay.shape[1]          <= nk_max  ) and (stay.shape[2]     <= nb_max  ) and (stay.shape[3]     <= T_max   ) and (stay.shape[4]    <= ng_max  )

        ## Distribution generation

        # Initialize population distribution for next year using population growth distribution
        DIST_next = np.copy(DIST_grow)
        
        for age in range(1,T_life):
            
            # Extract optimal decision values
            k_t = np.copy(K[:,:,:,age-1])
            b_t = np.copy(B[:,:,:,age-1])
            k_shape = k_t.shape
            b_shape = b_t.shape
            k_t_flat = k_t.flatten(order = 'F')
            b_t_flat = b_t.flatten(order = 'F')
    
            # Find indices of nearest values in decision value discretization vectors
            jk_lt = np.ones(len(k_t_flat))
            first = kv[0:-1]
            for elem in range(len(k_t_flat)):
                second = k_t_flat[elem]
                jk_lt[elem] = np.where(first <= second)[0][-1]
                
            jk_lt = np.reshape(jk_lt, k_shape, order = 'F').astype('int')
            jk_gt = jk_lt + 1
            
            jb_lt = np.ones(len(b_t_flat))
            first = bv[0:-1]
            for elem in range(len(b_t_flat)):
                
                second = b_t_flat[elem]
                jb_lt[elem] = np.where(first <= second)[0][-1]
            
            jb_lt = np.reshape(jb_lt, b_shape, order = 'F').astype('int')
            jb_gt = jb_lt + 1
            
            # Calculate linear weights for nearest values
            wk_lt = (kv[jk_gt] - k_t) / (kv[jk_gt] - kv[jk_lt])
            wk_gt = 1 - wk_lt
    
            wb_lt = (bv[jb_gt] - b_t) / (bv[jb_gt] - bv[jb_lt])
            wb_gt = 1 - wb_lt
            
            # Checks -> only work in the absence of mex file!
            assert (wk_lt>=0).all() and (wk_gt>=0).all() and (wb_lt>=0).all() and (wb_gt>=0).all(), 'Negative weights to compute households distribution.'     
            
            for jz in range(nz):
                
                # Apply survival and productivity transformations to population distribution for current year
                DIST_transz = DIST_year[:,:,:,age-1,:] * surv[age-1] * stay[:,:,:,age-1,:] * np.tile(np.reshape(transz[age,:,jz], (nz,1,1,1)), [1,nk,nb,ng])
                assert (DIST_transz>=0).all(), 'Negative mass of people at DIST_transz.'
        
                # Redistribute population distribution from current year to next year according to target indices and weights
                for ib in range(nb):
                    for ik in range(nk):
                        for iz in range(nz): #ok<ALIGN>
  
                          DIST_next[jz, jk_lt[iz,ik,ib], jb_lt[iz,ik,ib], age, :] = DIST_next[jz, jk_lt[iz,ik,ib], jb_lt[iz,ik,ib], age, :] + wk_lt[iz,ik,ib]*wb_lt[iz,ik,ib]*DIST_transz[iz,ik,ib,:]
                          DIST_next[jz, jk_gt[iz,ik,ib], jb_lt[iz,ik,ib], age, :] = DIST_next[jz, jk_gt[iz,ik,ib], jb_lt[iz,ik,ib], age, :] + wk_gt[iz,ik,ib]*wb_lt[iz,ik,ib]*DIST_transz[iz,ik,ib,:]
                          DIST_next[jz, jk_lt[iz,ik,ib], jb_gt[iz,ik,ib], age, :] = DIST_next[jz, jk_lt[iz,ik,ib], jb_gt[iz,ik,ib], age, :] + wk_lt[iz,ik,ib]*wb_gt[iz,ik,ib]*DIST_transz[iz,ik,ib,:]
                          DIST_next[jz, jk_gt[iz,ik,ib], jb_gt[iz,ik,ib], age, :] = DIST_next[jz, jk_gt[iz,ik,ib], jb_gt[iz,ik,ib], age, :] + wk_gt[iz,ik,ib]*wb_gt[iz,ik,ib]*DIST_transz[iz,ik,ib,:]

        return DIST_next