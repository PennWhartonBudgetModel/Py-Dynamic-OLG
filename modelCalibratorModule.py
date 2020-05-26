 ##
# Dynamic model baseline parameter calibrator.
# 

from pathFinderModule import PathFinder
from scenarioModule import Scenario
from modelSolverModule import ModelSolver
from paramGeneratorModule import ParamGenerator

import os
import shutil
import numpy as np
import pandas as pd
import pickle
import scipy
import matplotlib.pyplot as plt
import datetime

class ModelCalibrator:
    
    # Define list of parameters which define the steady state
    paramlist = ['beta', 'gamma', 'sigma', 'modelunit_dollar']
    
    # Define list of targets
    targetlist = ['captoout', 'labelas', 'savelas', 'outperHH']
    ntarget    = len(targetlist)
    
    # Define number of discretization points for each dimension of the calibration grid
    ngrid = 15
    
    # Determine total number of calibration points
    #   There are 3 dimensions for the calibration grid -- beta, sigma, gamma
    npoint = ngrid ** 3
    
    # Define calibration point directory and calibration point file path
    pointdir  = os.path.join(PathFinder.getSourceDir(), 'CalibrationPoints')
    pointfile = lambda ipoint: os.path.join(ModelCalibrator.pointdir, 'point%05d.pkl' % ipoint)
    
    # Define the moment targets for the reports on how we did
    #   Cell array: { Variable Name, Value, Description }
    moment_targets   = [['r',        0.05,   'Return on capital'],
                     ['PIT',      0.08,   'PIT/GDP'],
                     ['SSTax',    0.05,   'SSTax/GDP'],
                     ['KbyY',     3.0,    'Capital/GDP'],
                     ['outperHH', 7.98e4, 'GDP$/adult']]
    
    # Define calibration points
    @staticmethod
    def definePoints():
        
        assert ModelCalibrator.npoint <= 75000, 'Number of calibration points exceeds HPCC task array job size limit.'
        
        # Clear or create calibration point directory
        if os.path.exists(ModelCalibrator.pointdir):
            shutil.rmtree(ModelCalibrator.pointdir)
        os.mkdir(ModelCalibrator.pointdir)
        
        # Specify parameter lower and upper bounds
        lb = {}
        ub = {}
        lb['beta'] = 0.920
        lb['gamma'] = 0.150
        lb['sigma'] = 1.20
        ub['beta'] = 1.050
        ub['gamma'] = 0.900
        ub['sigma'] = 9.00
        
        # Construct vectors of parameter values
        v = {}
        v['beta']  = np.linspace(lb['beta'], ub['beta'], num=ModelCalibrator.ngrid)
        v['gamma'] = np.linspace(lb['gamma'], ub['gamma'], num=ModelCalibrator.ngrid)
        v['sigma'] = np.logspace(np.log10(lb['sigma']), np.log10(ub['sigma']), num=ModelCalibrator.ngrid)
        
        # Generate calibration points as unique combinations of parameter values
        grid = {}
        (grid['beta'], grid['gamma'], grid['sigma']) = np.meshgrid(v['beta'], v['gamma'], v['sigma'])
        for ipoint in range(ModelCalibrator.npoint):
            params = {}
            for p in ['beta', 'gamma', 'sigma']:
                params[p] = grid[p][ipoint] #ok<STRNU>
            with open(ModelCalibrator.pointfile(ipoint)) as f:
                pickle.dump(params, f)
    
    # Solve calibration point
    @staticmethod
    def calibratePoint(ipoint):
        
        # Load parameter values for calibration point
        with open(ModelCalibrator.pointfile(ipoint)) as f:
            s = pickle.load(f)
        params = s['params']
        
        try:
            # Calibrate steady state on modelunit_dollar
            (targets, modelunit_dollar, solved) = ModelCalibrator.calibrate_dollar(params) #ok<ASGLU>   
        except:
            print('Error encountered calibrating point %u:\n\t \nSaving placeholder solution values.\n' % ipoint)
            
            for o in ModelCalibrator.targetlist:
                targets[o] = None
            modelunit_dollar = None
            solved = 0 #ok<NASGU>
            
        # Extend parameters structure
        params['modelunit_dollar'] = modelunit_dollar
        
        # Save parameters, targets, and solution condition to calibration point file
        with open(ModelCalibrator.pointfile(ipoint)) as f:
            pickle.dump(params)
            pickle.dump(targets)
            pickle.dump(solved)

    # Consolidate solved calibration points
    @staticmethod
    def consolidatePoints():
        
        # Clear or create calibration output directory
        outputdir = PathFinder.getCalibrationOutputDir()
        if os.path.exists(outputdir):
            shutil.rmtree(outputdir)
        os.mkdir(outputdir)
        
        paramv = {}
        targetv = {}
        
        # Initialize vectors of parameters, targets, and solution conditions
        for o in ModelCalibrator.paramlist:
            paramv[o] = np.empty(shape=(1, ModelCalibrator.npoint))
        for o in ModelCalibrator.targetlist:
            targetv[o] = np.empty(shape=(1, ModelCalibrator.npoint))
        solved = np.zeros(shape=(1, ModelCalibrator.npoint))
        
        # Load and consolidate calibration points
        for i in range(ModelCalibrator.npoint):
            
            print('Reading calibration point %5d of %5d\n' % (i, ModelCalibrator.npoint))
            
            s = {}
            
            with open(ModelCalibrator.pointfile[i]) as f:
                s['params'] = pickle.load(f)
                s['targets'] = pickle.load(f)
                s['solved'] = pickle.load(f)
                
            for o in ModelCalibrator.paramlist:
                paramv[o][i] = s['params'][o]
            for o in ModelCalibrator.targetlist:
                targetv[o][i] = s['targets'][o]
            solved[i] = s['solved']
        
        # Save consolidated points to calibration output directory
        with open(os.path.join(outputdir, 'calibration.pkl')) as f:
            pickle.dump(paramv)
            pickle.dump(targetv)
            pickle.dump(solved)
            
        # Initialize plot of calibration point solution conditions
        fig = plt.figure()
        ax = fig.add_subplot(111, projection = '3d')
        
        # Determine colors
        cv = np.zeros((ModelCalibrator.npoint, 3))
        devs = min(abs(targetv['captoout'][solved] - 3) ** 0.5, 1)
        cv[solved, :] = np.hstack((devs, np.ones(devs.shape), devs)) * 180 / 256        # Gray to green
        cv[np.logical_not(solved), :] = np.tile([200/256, 0, 0], (sum(np.logical_not(solved)), 1))    # Red
        
        # Plot solution conditions
        ax.scatter(paramv['beta'], paramv['gamma'], paramv['sigma'], s = 40, c = cv, marker = 'o')
        
        # Format axes
        plt.axis('tight')
        ax.set_frame_on(True)
        ax.grid(b = True, which = 'minor')
        ax.set_aspect(num = 1)
        ax.set_xlabel('beta')
        ax.set_xscale('linear')
        ax.set_xticks(np.linspace(ax.get_xlim[0], ax.get_ylim[1], num = 3))
        ax.set_ylabel('gamma')
        ax.set_yscale('linear')
        ax.set_yticks(np.linspace(ax.get_ylim[0], ax.get_ylim[1], num = 3))
        ax.set_zlabel('sigma')
        ax.set_zscale('log')
        ax.set_zticks(np.logspace(np.log10(ax.get_zlim[0]), np.log10(ax.get_zlim[1]), num = 3))
        ax.minorticks_off()
        ax.grid(which = 'minor')
        
        # Save plot to calibration output directory
        plt.savefig(fig, os.path.join(outputdir, 'conditions.fig'))
    
    ##
    #  Single loop to calibrate on modelunit_dollar targets
    def calibrate_dollar(gridpoint):

        # Set target = $gdp/adult
        #     from Alex $79.8k for 2016
        #     REM: In moment_targets, 
        #        col 1 = varname, col 2 = value, col 3 = description 
        target_outperHH_index = np.where(ModelCalibrator.moment_targets[:, 0] == 'outperHH')[0]
        target_outperHH       = np.array([ModelCalibrator.moment_targets[target_outperHH_index, 1]])
        
        # Set initial modelunit_dollar.
        # In the future, we could apply a heuristic better initial guess.
        modelunit_dollar    = 4.0e-05  

        tolerance           = 0.01    # as ratio 
        err_size            = 1
        iter_num            = 1
        iter_max            = 8       # iterations for modelunit_dollar

        while err_size > tolerance and iter_num <= iter_max:

            # Create Scenario to run
            scenario    = Scenario( { 'economy'           : 'steady'          ,
                                      'beta'              : gridpoint.beta    ,
                                      'gamma'             : gridpoint.gamma   ,
                                      'sigma'             : gridpoint.sigma   ,
                                      'modelunit_dollar'  : modelunit_dollar  ,
                                      'bequest_phi_1'     : 0                 } )
            save_dir    = ModelSolver.solve( scenario )

            # find target -- $gdp/pop
            with open(os.path.join(save_dir, 'paramsTargets.pkl'), 'rb') as handle:
                s_paramsTargets = pickle.load(handle)
            run_outperHH    = s_paramsTargets['outperHH']
            
            err_size        = abs( run_outperHH/target_outperHH - 1 )
            print( '...MODELUNIT_DOLLAR iteration %u   error=%f\n ' % (iter_num, err_size) )

            # package up answer
            targets         = {   'savelas':      s_paramsTargets['savelas'],  
                                  'labelas':      s_paramsTargets['labelas'],  
                                  'captoout':     s_paramsTargets['captoout'],  
                                  'outperHH':     run_outperHH             }                       

            # Update by percent shift, reduced a bit as number of 
            # iterations increases. This approach slows the update rate
            # in case of slow convergence -- we're usually bouncing around then.
            exp_reduce        = max( 0.5, 1.0 - iter_num *0.07 )
            modelunit_dollar = modelunit_dollar*((run_outperHH/target_outperHH)**exp_reduce)

            # Find if converged
            #    This only needs to be done after the loop, but
            #    we're about to wipe out the run's files.
            with open(os.path.join(save_dir, 'dynamics.pkl' ), 'rb') as handle:
                s_dynamics     = pickle.load(handle)
            is_converged   = s_dynamics['is_converged']

            # Delete save directory along with parent directories
            shutil.rmtree(os.path.join(save_dir, '..', '..'))
            
            iter_num = iter_num + 1
   
        # Keep last successful run with modelunit_dollar
        modelunit_dollar   = scenario.modelunit_dollar
        
        # Check solution condition.
        # Stable solution identified as:
        #  1. Robust solver convergence rate
        #  2. modelunit_dollar convergence
        is_solved = is_converged and ( err_size <= tolerance )
        if iter_num > iter_max :
           print( '...MODELUNIT_DOLLAR -- max iterations (%u) reached.\n' % iter_max )

        return (targets, modelunit_dollar, is_solved)

    ##
    #  Print moments info on a particular steady state
    def report_moments( save_dir, targets = None ):

        delimiter   = [chr(13), chr(10)]  # end-of-line 
        
        filepath    = os.path.join(save_dir % 'iterations.csv')
        T           = pd.read_csv(filepath)
        iters       = T.iloc[:,0].values
        iterations  = iters[-1]

        with open(os.path.join(save_dir, 'dynamics.pkl' ), 'rb') as handle:
            s_dynamics      = pickle.load(handle)
        with open(os.path.join(save_dir, 'paramsTargets.pkl'), 'rb') as handle:
            s_paramsTargets = pickle.load(handle)
        with open(os.path.join(save_dir, 'market.pkl'), 'rb') as handle:
            s_markets       = pickle.load(handle)
        
        # Define some helper vars for clarity
        pop             = s_dynamics['pops']    
        gdp             = s_dynamics['outs']    
        dollar          = 1/s_paramsTargets['modelunit_dollar'] 

        if targets == None:
            targets = ModelCalibrator.moment_targets
            targets = np.vstack((targets, ['labelas', 1, 'Labor elasticity']))
            targets = np.vstack((targets, ['savelas', 1, 'Savings elasticity']))

        # helper function to format results
        myTargetPrint = (lambda lbl, modelResult, targetResult:
                    '   %20s = %f (%f) error = %0.1f%%' % (lbl, modelResult, targetResult, (modelResult/targetResult - 1)*100.0 ))
        myParamPrint  = lambda lbl, modelInput : '   %20s = %f' % (lbl, modelInput )
        
        # Make PARAMS section
        params  = {'beta'   : s_paramsTargets['beta'],
                   'sigma'  : s_paramsTargets['sigma'],
                   'gamma'  : s_paramsTargets['gamma'],
                   'model$' : s_paramsTargets['modelunit_dollar'] }
        
        param_part = '%s   PARAMS%s' % (delimiter, delimiter)
        for i in range(len(params)):             
            result  = params[i, 1]
            lbl     = params[i, 0]
            line    = myParamPrint(lbl, result)
            param_part = '%s%s%s' % (param_part, line, delimiter)
        
        # Make structure for results 
        #   targets has been passed in (or set to default)
        model_results = {'r'          : s_markets['MPKs'],
                         'PIT'        : s_dynamics['pits']/gdp,
                         'SSTax'      : s_dynamics['ssts']/gdp,
                         'KbyY'       : s_paramsTargets['captoout'],
                         'outperHH'   : gdp * dollar/pop,
                         'labelas'    : s_paramsTargets['labelas'],
                         'savelas'    : s_paramsTargets['savelas']}
        
        # Make TARGETS section
        target_part = '%s   TARGETS%s' % (delimiter, delimiter)
        for i in range(len(model_results[:, 0])):  
            m_index = np.where( targets[:, 0] == model_results[i,0])[0]
            target  = targets[ m_index, 1]
            lbl     = targets[ m_index, 2]
            result  = model_results[i, 1]
            
            line        = myTargetPrint(lbl, result, target)
            target_part = '%s%s%s' % (target_part, line, delimiter)
        
        # Make convergence part
        if s_dynamics['is_converged'] :
            s_iter = 'Converged in %u iterations' % iterations 
        else :
            s_iter = 'DID NOT converge in %u iterations.' % iterations 
        
        converge_part = '%s CONVERGENCE: %s %s' % (delimiter, s_iter, delimiter)      

        # Concatentate for full report
        outstr = '%s%s%s' % ( param_part, target_part, converge_part )
        return outstr

    ##
    #   Make a report of various moments for the 16 baselines
    def report_baseline_moments():
        
        outputfilename      = os.path.join(PathFinder.getSourceDir(), 'BaselineMoments.txt')
        f             = open(outputfilename,'w+')
        
        f.write('-------------BASELINE MOMENTS-------------')
        f.write('%s \r\n' % str(datetime.datetime.now()))
        
        # load the matrix and get inverter function
        (_, f_invert) = ParamGenerator.invert();
        
        for labelas in np.arange(0.25,1.0,0.25):
            for savelas in np.arange(0.25,1.0,0.25):
                target = {'labelas': labelas, 'savelas': savelas}
                f.write('\r\nBASELINE labor elas = %0.2f  savings elas = %0.2f \r\n' % (labelas, savelas) ) 
                inverse = f_invert(target)
                
                scenario = Scenario({'economy': 'steady',
                                     'beta': inverse['beta'],
                                     'gamma': inverse['gamma'],
                                     'sigma': inverse['sigma'],
                                     'modelunit_dollar': inverse['modelunit_dollar'],
                                     'bequest_phi_1': 0})
                
                save_dir = ModelSolver.solve(scenario)
                
                targets = ModelCalibrator.moment_targets
                targets = np.vstack((targets, ['labelas', labelas, 'Labor elasticity']))
                targets = np.vstack((targets, ['savelas', savelas, 'Savings elasticity']))
                outstr  = ModelCalibrator.report_moments( save_dir, targets )
                f.write( '%s \r\n' % outstr )
                f.write( '-------------------------------------\r\n' )                
            
        f.write( ' ==== DONE ===== \r\n' )    
        f.close()

    ## 
    #   Adjust calibration grid boundaries
def adjust_grid():
    
    epsilon = 1e-4
    (_, f_invert) = ParamGenerator.invert()
    green = [0, 180/256, 0]
    cv    = np.tile(np.reshape(green, (3,1)), [1,16])  #zeros(16,3)
    grid_beta  = [0.950, 1.100]
    grid_gamma = [0.150, 0.900]
    grid_sigma = [1.200, 9.000]
    delta_beta  = np.zeros((16,2))
    delta_gamma = np.zeros((16,2))
    delta_sigma = np.zeros((16,2))
    labelasv = np.zeros((16,1))
    savelasv = np.zeros((16,1))
    iter = 0
    
    for labelas in np.arange(0.25,1,0.25):
        for savelas in np.arange(0.25,1,0.25):
            
            target = {'labelas': labelas, 'savelas': savelas}
            inverse = f_invert(target)
            delta_beta[iter,:]  = np.hstack(((inverse['beta']  - grid_beta[0])/inverse['beta'] , (grid_beta[1] - inverse['beta'])/inverse['beta']))
            delta_gamma[iter,:] = np.hstack(((inverse['gamma'] - grid_gamma[0])/inverse['gamma'], (grid_gamma[1] - inverse['gamma'])/inverse['gamma']))
            delta_sigma[iter,:] = np.hstack(((inverse['sigma'] - grid_sigma[0])/inverse['sigma'], (grid_sigma[1] - inverse['sigma'])/inverse['sigma']))
            delta = min(np.minimum(np.minimum(delta_beta[iter,:], delta_gamma[iter,:]), delta_sigma[iter,:]))
            if delta <= epsilon :
                cv[iter,:] = [1, delta, 0] * 200/256
            
            labelasv[iter,0] = labelas
            savelasv[iter,0] = savelas
            iter = iter + 1
            
    # Plot
    plt.figure()
    plt.scatter(labelasv, savelasv, s = 40, c = cv, marker = 'o')
    plt.xlabel('labor elasticity', fontsize = 13)
    plt.xticks(ticks = np.arange(0,1.00,0.25)) 
    plt.ylabel('savings elasticity', fontsize = 13)
    plt.yticks(ticks = np.arange(0,1.00,0.25))
    plt.grid(b = True)

    # Adjust grids
    if min(delta_beta[:,0])  <= epsilon :
        grid_beta[0]  = 0.9*grid_beta[0]
    if min(delta_beta[:,1])  <= epsilon :
        grid_beta[1]  = 1.1*grid_beta[1]
    if min(delta_gamma[:,0]) <= epsilon :
        grid_gamma[0] = 0.9*grid_gamma[0] 
    if min(delta_gamma[:,1]) <= epsilon :
        grid_gamma[1] = 1.1*grid_gamma[1]
    if min(delta_sigma[:,0]) <= epsilon :
        grid_sigma[0] = max(1.01, 0.9*grid_sigma[0])
    if min(delta_sigma[:,1]) <= epsilon :
        grid_sigma[1] = 1.1*grid_sigma[1]
  
    return (grid_beta, grid_gamma, grid_sigma)

