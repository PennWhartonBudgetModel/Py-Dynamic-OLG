##
# Generate moments distribution for the transition path of a counterfactual economy whose working folder is complete, that is, after solving for it

import pathFinderModule as PathFinder
import paramGeneratorModule as ParamGenerator

import pickle
import os
import numpy as np
import openpyxl
import scipy
import pandas as pd

class TransitionMoments:
    
    @staticmethod
    def showDistribution(scenario):
        
        # Scenarios
        sc_steady  = scenario.currentPolicy().steady()
        if scenario == sc_steady:
            raise Exception('TransitionMoments scenario economy cannot be steady state.')
        if scenario.isOpen():
            sc_base = scenario.currentPolicy().open()
        else:
            sc_base = scenario.currentPolicy().closed()
        if scenario == sc_base:
            raise Exception('Scenario must use a counterfactual policy plan.')
        
        # Directories
        steady_dir = PathFinder.getCacheDir(sc_steady)
        sc_dir     = PathFinder.getCacheDir(scenario)
        base_dir   = PathFinder.getCacheDir(sc_base)
        
        ##PARAMETERS
        
        # Define time constants
        s = ParamGenerator.timing(scenario)
        T_life  = s['T_life']    # Total life years
        T_model = s['T_model']   # Transition path model years

        # Define grids
        s = ParamGenerator.grids(scenario)
        ng     = s['ng']         # num groups
        nz     = s['nz']         # num labor productivity shocks
        zs     = s['zs']         # shocks grid (by demographic type and age)
        nk     = s['nk']         # num asset points
        nb     = s['nb']         # num avg. earnings points
        kv     = s['kv']         # asset grid
        
        # Import total income without Social Security transfers in the baseline economy and define percentiles
        
        (totinc_base, _) = TransitionMoments.append_decisions('TOTINCbase', nz, nk, nb, T_life, ng, T_model, kv, zs, steady_dir, base_dir, [])
        (dist_base, _) = TransitionMoments.append_dist(nz, nk, nb, T_life, ng, T_model, steady_dir, base_dir, 1)
        (_, sort_base, index_base) = TransitionMoments.get_percentiles(totinc_base, dist_base, T_model)

        ## Households distribution

        (dist, dist_static) = TransitionMoments.append_dist(nz, nk, nb, T_life, ng, T_model, steady_dir, sc_dir, 0)

        ## Generate groups for other variables based on the distribution of total income without SS in baseline economy
        groups = {}

        for var_name in ['LABOR', 'SAVINGS', 'CONSUMPTION', 'ASSETS', 'TAXABLE_INC', 'PAYROLL_LIABILITY', 'OASI_BENEFITS', 'ORD_LIABILITY', 'PREF_LIABILITY', 'LABINC', 'TAX', 'TOTINC', 'TOTINCwSS', 'AINC']:

            # append steady state values
            (var, var_static) = TransitionMoments.append_decisions(var_name, nz, nk, nb, T_life, ng, T_model, kv, zs, steady_dir, sc_dir, base_dir)

            # generate percentile-like groups (dynamic and static)
            groups[var_name] = TransitionMoments.generate_groups(var, dist, sort_base, index_base, T_model)
            groups[var_name + '_static'] = TransitionMoments.generate_groups(var_static, dist_static, sort_base, index_base, T_model)

            # generate deltas
            groups[var_name + '_delta']  = groups[var_name] / groups[var_name + '_static']

        with open(os.path.join(sc_dir, 'groups.pkl')) as f:
            pickle.dump(groups, f)
            
        # Save deltas in a spreadsheet
        header = np.array(['year'])
        delta_table = np.arange(2016,2017 + T_model).reshape((-1,1))
        for v in groups.keys():
            if '_delta' in v:
                for i in range(7):
                    header_name = '%s_%i' % (v, i)
                    header.append(header_name)
                    delta_table = np.hstack(groups[v][:,i], delta_table)

        delta_table = pd.DataFrame(delta_table, columns = header)

        # Save counter variables in a spreadsheet
        header = np.array(['year'])
        counter_table = np.arange(2017,2017 + T_model).reshape((-1,1))
        for var_name in ['LABOR', 'SAVINGS', 'CONSUMPTION', 'ASSETS', 'TAXABLE_INC', 'PAYROLL_LIABILITY', 'OASI_BENEFITS', 'ORD_LIABILITY', 'PREF_LIABILITY', 'LABINC', 'TAX', 'TOTINC', 'TOTINCwSS', 'AINC']:
            for i in range(7):
                header_name = '%s_%i' % (var_name, i)
                header.append(header_name)
                counter_table = np.hstack(counter_table, groups[var_name][:,i])

        counter_table = pd.DataFrame(counter_table, columns = header)
        
        #use opepyxl to create workbook
        filepath = os.path.join(sc_dir, 'moments.xlsx')
        wb = openpyxl.Workbook()
        wb.save(filepath)
        
        #use pandas to save dataframes to workbook
        writer = pd.Excelwriter(filepath)
        delta_table.to_excel(writer, 'delta')
        counter_table.to_excel(writer, 'counter')

        ## FUNCTIONS
        # Pre-appends steady state distribution
        @staticmethod
        def append_dist(nz, nk, nb, T_life, ng, T_model, dir_ss, dir_sc, base):

        # Inputs:  (nz, nk, nb, T_life, ng, T_model) = dimensions of DIST (number of states for each state variable)
        #          dir_ss = steady state directory
        #          dir_sc = directory of the scenario of interest
        # Outputs: x        = dynamic distribution array with pre-appended steady state distribution
        #          x_static = static distribution array with pre-appended steady state distribution

            x        = np.zeros((nz,nk,nb,T_life,ng,T_model+1))

            with open(os.path.join(dir_ss, 'distribution.pkl'), 'rb') as handle:
                s = pickle.load(handle)
            x[:,:,:,:,:,0,:] = s['DIST']

            if not base:
                x_static = np.zeros((nz,nk,nb,T_life,ng,T_model+1))
                x_static[:,:,:,:,:,0,:] = s['DIST']
                with open(os.path.join(dir_sc, 'Static_distribut.pkl'), 'rb') as handle:
                    s = pickle.load(handle)
                x_static[:,:,:,:,:,1:,:] = s['Static_DIST']
            else:
                x_static = []

            with open(os.path.join(dir_sc, 'distribution.pkl'), 'rb') as handle:
                s = pickle.load(handle)
            x[:,:,:,:,:,1:,:] = s['DIST']

            return (x, x_static)

        # Pre-appends steady state variables from decisions mat file
        @staticmethod
        def append_decisions(x_name, nz, nk, nb, T_life, ng, T_model, kv, zs, dir_ss, dir_sc, dir_bs):

        # Inputs:  (nz, nk, nb, T_life, ng, T_model) = dimensions of DIST (number of states for each state variable)
        #          x_name = name of variable of interest
        #          kv     = capital vector
        #          zs     = productivity shocks matrix
        #          dir_ss = steady state directory
        #          dir_sc = directory of the scenario of interest
        #          dir_bs = directory of the baseline case of the scenario of interest
        # Outputs: x        = dynamic variable array with pre-appended steady state value
        #          x_static = static variable array with pre-appended steady state value

            # Assets case
            if x_name == 'ASSETS':

                x        = np.tile(np.reshape(kv, (1,nk,1,1,1,1)),[nz,1,nb,T_life,ng,T_model+1])
                x_static = np.tile(np.reshape(kv, (1,nk,1,1,1,1)),[nz,1,nb,T_life,ng,T_model+1])

            # Labor income case
            elif x_name == 'LABINC':

                x        = np.zeros((nz,nk,nb,T_life,ng,T_model+1))
                x_static = np.zeros((nz,nk,nb,T_life,ng,T_model+1))

                with open(os.path.join(dir_ss, 'market.pkl'), 'rb') as handle:
                    s        = pickle.load(handle)
                wages_ss = s['wages']
                f = lambda X: np.tile(np.reshape(X, (nz,nk,nb,T_life,1,1)), [1,1,1,1,ng,1])
                with open(os.path.join(dir_ss, 'decisions.pkl'), 'rb') as handle:
                    s = pickle.load(handle)
                x[:,:,:,:,:,0,:]        = wages_ss*f(s['OPTs']['LABOR'] * np.tile(np.reshape(np.transpose(zs, [2,1,0]), (nz,1,1,T_life,1)), [1,nk,nb,1,1]))
                x_static[:,:,:,:,:,0,:] = wages_ss*f(s['OPTs']['LABOR'] * np.tile(np.reshape(np.transpose(zs, [2,1,0]), (nz,1,1,T_life,1)), [1,nk,nb,1,1]))

                with open(os.path.join(dir_sc, 'market.pkl'), 'rb') as handle:
                    s     = pickle.load(handle)
                wages = s['wages']
                f = lambda X: np.tile(np.reshape(X, (nz,nk,nb,T_life,1,T_model)), [1,1,1,1,ng,1])
                with open(os.path.join(dir_sc, 'decisions.pkl'), 'rb') as handle:
                    s = pickle.load(handle)
                x[:,:,:,:,:,1:,:] = f(np.tile(np.reshape(wages, (1,1,1,1,T_model)), [nz,nk,nb,T_life,1]) * s['OPTs']['LABOR'] * np.tile(np.reshape(np.transpose(zs, [2,1,0]), (nz,1,1,T_life,T_model)), [1,nk,nb,1,1]))
                with open(os.path.join(dir_bs, 'market.pkl'), 'rb') as handle:
                    s = pickle.load(handle)
                wages_static = s['wages']
                with open(os.path.join(dir_sc, 'Static_decisions.pkl'), 'rb') as handle:
                    s = pickle.load(handle)
                x_static[:,:,:,:,:,1:,:] = f(np.tile(np.reshape(wages_static, [1,1,1,1,T_model]), (nz,nk,nb,T_life,1)) * s['LABOR'] * np.tile(np.reshape(np.transpose(zs, [3,2,1]), (nz,1,1,T_life,T_model)), [1,nk,nb,1,1]))

            # Assets income case
            elif x_name == 'AINC':

                x        = np.zeros((nz,nk,nb,T_life,ng,T_model+1))
                x_static = np.zeros((nz,nk,nb,T_life,ng,T_model+1))
                
                with open(os.path.join(dir_ss, 'market.pkl'), 'rb') as handle:
                    s        = pickle.load(handle)
                totrates_ss = s['capsharesPM'] * s['equityDividendRates'] + (1 - s['capsharesPM']) * s['bondDividendRates']
                x[:,:,:,:,:,0,:]        = totrates_ss * np.tile(np.reshape(kv, (1,nk,1,1,1,1)),[nz,1,nb,T_life,ng,1])
                x_static[:,:,:,:,:,0,:] = totrates_ss * np.tile(np.reshape(kv, (1,nk,1,1,1,1)),[nz,1,nb,T_life,ng,1])

                f = lambda X: np.tile(np.reshape(X, (1, 1,1,1,1,T_model)), [nz,nk,nb,T_life,ng,1])
                g = lambda X: np.tile(np.reshape(X, (1,nk,1,1,1,1)), [nz, 1,nb,T_life,ng,T_model])
                
                with open(os.path.join(dir_sc, 'market.pkl'), 'rb') as handle:
                    s = pickle.load(handle)
                totrates = s['capsharesPM'] * s['equityDividendRates'] + (1 - s['capsharesPM']) * s['bondDividendRates']
                x[:,:,:,:,:,1:,:] = f(totrates) * g(kv)
                with open(os.path.join(dir_bs, 'market.pkl'), 'rb') as handle:
                    s = pickle.load(handle)
                totrates_static = s['capsharesPM'] * s['equityDividendRates'] + (1 - s['capsharesPM']) * s['bondDividendRates']
                x_static[:,:,:,:,:,1:,:] = f(totrates_static) * g(kv)

            # Total taxes paid case
            elif x_name == 'TAX':

                x        = np.zeros((nz,nk,nb,T_life,ng,T_model+1))
                x_static = np.zeros((nz,nk,nb,T_life,ng,T_model+1))

                f = lambda X: np.tile(np.reshape(X, (nz,nk,nb,T_life,1,1)), [1,1,1,1,ng,1])
                with open(os.path.join(dir_ss, 'decisions.pkl'), 'rb') as handle:
                    s = pickle.load(handle)
                x[:,:,:,:,:,0,:]        = f(s['OPTs']['PAYROLL_LIABILITY'] + s['OPTs']['ORD_LIABILITY'] + s['OPTs']['PREF_LIABILITY'])
                x_static[:,:,:,:,:,0,:] = f(s['OPTs']['PAYROLL_LIABILITY'] + s['OPTs']['ORD_LIABILITY'] + s['OPTs']['PREF_LIABILITY'])

                f = lambda X: np.tile(np.reshape(X, (nz,nk,nb,T_life,1,T_model)), [1,1,1,1,ng,1])
                with open(os.path.join(dir_sc, 'decisions.pkl'), 'rb') as handle:
                    s = pickle.load(handle)
                x[:,:,:,:,:,1:,:] = f(s['OPTs']['PAYROLL_LIABILITY'] + s['OPTs']['ORD_LIABILITY'] + s['OPTs']['PREF_LIABILITY'])
                with open(os.path.join(dir_sc, 'Static_decisions.pkl'), 'rb') as handle:
                    s = pickle.load(handle)
                x_static[:,:,:,:,:,1:,:] = f(s['PAYROLL_LIABILITY'] + s['ORD_LIABILITY'] + s['PREF_LIABILITY'])

            # Total income case
            elif x_name == 'TOTINC':

                x        = np.zeros((nz,nk,nb,T_life,ng,T_model+1))
                x_static = np.zeros((nz,nk,nb,T_life,ng,T_model+1))

                with open(os.path.join(dir_ss, 'market.pkl'), 'rb') as handle:
                    s = pickle.load(handle)
                wages_ss = s['wages']
                totrates_ss = s['capsharesPM'] * s['equityDividendRates'] + (1 - s['capsharesPM']) * s['bondDividendRates']
                f = lambda X: np.tile(np.reshape(X, (nz,nk,nb,T_life,1,1)), [1,1,1,1,ng,1])
                with open(os.path.join(dir_ss, 'decisions.pkl'), 'rb') as handle:
                    s = pickle.load(handle)
                x[:,:,:,:,:,0,:]        = f(wages_ss*s['OPTs']['LABOR'] * np.tile(np.reshape(np.transpose(zs, [2,1,0]), (nz,1,1,T_life,1)), [1,nk,nb,1,1]) +
                                          totrates_ss*np.tile(np.reshape(kv, (1,nk,1,1,1)),[nz,1,nb,T_life,1]))
                x_static[:,:,:,:,:,0,:] = f(wages_ss*s['OPTs']['LABOR'] * np.tile(np.reshape(np.transpose(zs, [2,1,0]), (nz,1,1,T_life,1)), [1,nk,nb,1,1]) +
                                          totrates_ss*np.tile(np.reshape(kv, (1,nk,1,1,1)),[nz,1,nb,T_life,1]))

                with open(os.path.join(dir_sc, 'market.pkl'), 'rb') as handle:
                    s     = pickle.load(handle)
                wages = s['wages']
                totrates = s['capsharesPM'] * s['equityDividendRates'] + (1 - s['capsharesPM']) * s['bondDividendRates']
                f = lambda X: np.tile(np.reshape(X, (nz,nk,nb,T_life,1,T_model)), [1,1,1,1,ng,1])
                with open(os.path.join(dir_sc, 'decisions.pkl'), 'rb') as handle:
                    s = pickle.load(handle)
                x[:,:,:,:,:,1:,:] = f(np.tile(np.reshape(wages, (1,1,1,1,T_model)), [nz,nk,nb,T_life,1]) * s['OPTs']['LABOR'] * np.tile(np.reshape(np.transpose(zs, [2,1,0]), (nz,1,1,T_life,T_model)), [1,nk,nb,1,1]) +
                                         np.tile(np.reshape(totrates, (1,1,1,1,T_model)), [nz,nk,nb,T_life,1]) * np.tile(np.reshape(kv, (1,nk,1,1,1)),[nz,1,nb,T_life,T_model]))
                with open(os.path.join(dir_bs, 'market.pkl'), 'rb') as handle:
                    s     = pickle.load(handle)
                wages_static = s['wages']
                totrates_static = s['capsharesPM'] * s['equityDividendRates'] + (1 - s['capsharesPM']) * s['bondDividendRates']
                with open(os.path.join(dir_sc, 'Static_decisions.pkl'), 'rb') as handle:
                    s = pickle.load(handle)
                x_static[:,:,:,:,:,1:,:] = f(np.tile(np.reshape(wages_static, (1,1,1,1,T_model)), [nz,nk,nb,T_life,1]) * s['LABOR'] * np.tile(np.reshape(np.tile(zs, [2,1,0]), (nz,1,1,T_life,T_model)), [1,nk,nb,1,1]) +
                                                np.tile(np.reshape(totrates_static, (1,1,1,1,T_model)), [nz,nk,nb,T_life,1]) * np.tile(np.reshape(kv, [1,nk,1,1,1]),[nz,1,nb,T_life,T_model]))

            # Total income with Social Security transfers case
            elif x_name == 'TOTINCwSS':

                x        = np.zeros((nz,nk,nb,T_life,ng,T_model+1))
                x_static = np.zeros((nz,nk,nb,T_life,ng,T_model+1))

                with open(os.path.join(dir_ss, 'market.pkl'), 'rb') as handle:
                    s        = pickle.load(handle)
                wages_ss = s['wages']
                totrates_ss = s['capsharesPM']*s['equityDividendRates'] + (1-s['capsharesPM'])*s['bondDividendRates']
                f = lambda X: np.tile(np.reshape(X, (nz,nk,nb,T_life,1,1)), [1,1,1,1,ng,1])
                with open(os.path.join(dir_ss, 'decisions.pkl'), 'rb') as handle:
                    s = pickle.load(handle)
                x[:,:,:,:,:,0,:]        = f(wages_ss*s['OPTs']['LABOR'] * np.tile(np.reshape(np.transpose(zs, [2,1,0]), (nz,1,1,T_life,1)), [1,nk,nb,1,1]) + 
                                          totrates_ss*np.tile(np.reshape(kv, (1,nk,1,1,1)),[nz,1,nb,T_life,1]) + s['OPTs']['OASI_BENEFITS'])
                x_static[:,:,:,:,:,0,:] = f(wages_ss*s['OPTs']['LABOR'] * np.tile(np.reshape(np.transpose(zs, [2,1,0]), (nz,1,1,T_life,1)), [1,nk,nb,1,1]) +
                                          totrates_ss*np.tile(np.reshape(kv, (1,nk,1,1,1)),[nz,1,nb,T_life,1]) + s['OPTs']['OASI_BENEFITS'])

                with open(os.path.join(dir_sc, 'market.pkl'), 'rb') as handle:
                    s     = pickle.load(handle)
                wages = s['wages']
                totrates = s['capsharesPM'] * s['equityDividendRates'] + (1 - s['capsharesPM']) * s['bondDividendRates']
                f = lambda X: np.tile(np.reshape(X, (nz,nk,nb,T_life,1,T_model)), [1,1,1,1,ng,1])
                with open(os.path.join(dir_sc, 'decisions.pkl'), 'rb') as handle:    
                    s = pickle.load(handle)
                x[:,:,:,:,:,1:,:] = f(np.tile(np.reshape(wages, (1,1,1,1,T_model)), [nz,nk,nb,T_life,1]) * s['OPTs']['LABOR'] * np.tile(np.reshape(np.transpose(zs, [2,1,0]), (nz,1,1,T_life,T_model)), [1,nk,nb,1,1]) + 
                                         np.tile(np.reshape(totrates, (1,1,1,1,T_model)), [nz,nk,nb,T_life,1]) * np.tile(np.reshape(kv, (1,nk,1,1,1)),[nz,1,nb,T_life,T_model])
                                          + s['OPTs']['OASI_BENEFITS'])
                with os.path.join((dir_bs, 'market.pkl'), 'rb') as handle: 
                    s     = pickle.load(handle)
                wages_static = s['wages']
                totrates_static = s['capsharesPM']*s['equityDividendRates'] + (1-s['capsharesPM']) * s['bondDividendRates']
                with os.path.join((dir_sc, 'Static_decisions.pkl'), 'rb') as handle: 
                    s = pickle.load(handle)
                x_static[:,:,:,:,:,1:,:] = f(np.tile(np.reshape(wages_static, (1,1,1,1,T_model)), [nz,nk,nb,T_life,1]) * s['LABOR'] * np.tile(np.reshape(np.transpose(zs, [2,1,0]), (nz,1,1,T_life,T_model)), [1,nk,nb,1,1]) +
                                                np.tile(np.reshape(totrates_static, (1,1,1,1,T_model)), [nz,nk,nb,T_life,1]) * np.tile(np.reshape(kv, (1,nk,1,1,1)),[nz,1,nb,T_life,T_model])
                                                + s['OASI_BENEFITS'])

            # Total income without Social Security transfers case for baseline economy
            elif x_name == 'TOTINCbase':

                x        = np.zeros((nz,nk,nb,T_life,ng,T_model+1))
                x_static = np.array([])

                with os.path.join((dir_ss, 'market.pkl'), 'rb') as handle: 
                    s        = pickle.load(handle)
                wages_ss = s['wages']
                totrates_ss = s['capsharesPM']*s['equityDividendRates'] + (1-s['capsharesPM'])*s['bondDividendRates']
                f = lambda X: np.tile(np.reshape(X, (nz,nk,nb,T_life,1,1)), [1,1,1,1,ng,1])
                with os.path.join((dir_ss, 'decisions.pkl'), 'rb') as handle: 
                    s = pickle.load(handle)
                x[:,:,:,:,:,0,:]        = f(wages_ss*s['OPTs']['LABOR'] * np.tile(np.reshape(np.transpose(zs, [2,1,0]), (nz,1,1,T_life,1)), [1,nk,nb,1,1]) +
                                          totrates_ss*np.tile(np.reshape(kv, (1,nk,1,1,1)),[nz,1,nb,T_life,1]))

                with os.path.join((dir_sc, 'market.pkl'), 'rb') as handle: 
                    s     = pickle.load(handle)
                wages = s['wages']
                totrates = s['capsharesPM']*s['equityDividendRates'] + (1-s['capsharesPM']) * s['bondDividendRates']
                f = lambda X: np.tile(np.reshape(X, (nz,nk,nb,T_life,1,T_model)), [1,1,1,1,ng,1])
                with os.path.join((dir_sc, 'decisions.pkl'), 'rb') as handle: 
                    s = pickle.load(handle)
                x[:,:,:,:,:,1:,:] = f(np.tile(np.reshape(wages, (1,1,1,1,T_model)), [nz,nk,nb,T_life,1]) * s['OPTs']['LABOR'] * np.tile(np.reshape(np.transpose(zs, [2,1,0]), (nz,1,1,T_life,T_model)), [1,nk,nb,1,1]) + 
                                         np.tile(np.reshape(totrates, (1,1,1,1,T_model)), [nz,nk,nb,T_life,1]) * np.tile(np.reshape(kv, (1,nk,1,1,1)),[nz,1,nb,T_life,T_model]))

            # All other variables (already stored in decisions.mat file)    
            else:

                x        = np.zeros((nz,nk,nb,T_life,ng,T_model+1))
                x_static = np.zeros((nz,nk,nb,T_life,ng,T_model+1))

                f = lambda X: np.tile(np.reshape(X, (nz,nk,nb,T_life,1,1)), [1,1,1,1,ng,1])
                with os.path.join((dir_ss, 'decisions.pkl'), 'rb') as handle: 
                    s = pickle.load(handle)
                x[:,:,:,:,:,0,:]        = f(s['OPTs'][x_name])
                x_static[:,:,:,:,:,0,:] = f(s['OPTs'][x_name])

                f = lambda X: np.tile(np.reshape(X, (nz,nk,nb,T_life,1,T_model)), [1,1,1,1,ng,1])
                with os.path.join((dir_sc, 'decisions.pkl'), 'rb') as handle: 
                    s = pickle.load(handle)
                x[:,:,:,:,:,1:,:] = f(s['OPTs'][x_name])
                with os.path.join((dir_sc, 'Static_decisions.pkl'), 'rb') as handle: 
                    s = pickle.load(handle)
                x_static[:,:,:,:,:,1:,:] = f(s[x_name])

            return (x, x_static)

        # Get percentiles function
        @staticmethod
        def get_percentiles(x, dist, T_model):

        # Inputs:  x       = array with the variable of interest
        #          dist    = measure of households at each state of x
        #          T_model = number of transition periods
        # Outputs: percentiles = total amount of x held by each percentiles
        #          sort_x    = array on how to sort x in ascending order at each period
        #          index_x   = array with the cutoffs delimiting each quintile

            percentiles = np.zeros((T_model+1, 7))
            index_x   = np.zeros((T_model+1, 7))
            sort_x = np.zeros((x.shape[0],T_model+1))

            for t in range(T_model+1):

                # Vectorize variables
                x_t    = x[:,:,:,:,:,t,:]
                x_t    = x_t[:]
                dist_t = dist[:,:,:,:,:,t,:]
                dist_t = dist_t[:]
                dist_t = dist_t/sum(dist_t[:])

                # Sort variables (TBD not sure if I picked the right axis to sort on)
                x_t = np.sort(x_t, axis = 1)
                sort_x[:,t] = np.argsort(x_t, axis = 1)
                dist_t = dist_t[sort_x[:,t]]

                # Find percentiles
                x_cum = np.zeros((1,7))
                q = 0
                for percentile in [0.2, 0.4, 0.6, 0.8, 0.9, 0.95]:

                    # Total income with SS distribution
                    i = np.where(np.cumsum(dist_t, axis = 1) >= percentile)[0]
                    x_cum[0,q]   = np.sum(x_t[0:i] * dist_t[0:i])
                    index_x[t,q] = i

                    # Counter
                    q = q + 1

                # Top percentile
                x_cum[0,q]   = np.sum(x_t * dist_t)
                index_x[t,q] = len(x_t)

                # Groups
                percentiles[t,:] = np.diff(np.hstack((0, x_cum[0,:])))

            return (percentiles, sort_x, index_x)

        # Generate percentile-like groups according to how households are distributed wrt total income with SS percentiles in the baseline economy
        def generate_groups(x, dist, sortx, q_index, T_model):

        # Inputs:  x        = array with the variable of interest
        #          dist     = measure of households at each state of x
        #          sort_x   = array on how to sort x in total income without SS ascending order at each period
        #          index_x  = array with the cutoffs delimiting each total income without SS percentile
        #          T_model  = number of transition periods
        # Outputs: x_groups = total amount of x held by each 'percentile-like group'

            x_groups = np.zeros((T_model+1,7))

            for t in range(T_model+1):

                # Vectorize variables
                x_t    = x[:,:,:,:,:,t,:]
                dist_t = dist[:,:,:,:,:,t,:]
                x_t    = x_t[:]
                dist_t = dist_t[:]
                dist_t = dist_t/sum(dist_t[:])

                # Sort variables
                x_t    = x_t[sortx[:,1]]
                dist_t = dist_t[sortx[:,1]]

                # Find percentiles
                x_cum = np.zeros((1,7))
                q = 0
                for percentile in [0.2, 0.4, 0.6, 0.8, 0.9, 0.95]:

                    # Asset distributed according to total income without SS
                    i = q_index[1,q]
                    x_cum[0,q] = sum(x_t[0:i] * dist_t[0:i])
                    # Counter
                    q = q + 1

                # Top quintile
                x_cum[0,q] = sum(x_t[0:] * dist_t[0:])

                # Groups
                x_groups[t,:] = np.diff(np.hstack((0, x_cum[0,:])))

            return x_groups