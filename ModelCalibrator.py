##
# Dynamic model baseline parameter calibrator.
# 

import os
import PathFinder
import shutil
import numpy as np
import pickle

class ModelCalibrator:
    
    # Define list of parameters which define the steady state
    paramlist = ['beta', 'gamma', 'sigma', 'modelunit_dollar']
    
    # Define list of targets
    targetlist = ['captoout', 'labelas', 'savelas', 'outperHH']
    ntarget    = len(ModelCalibrator.targetlist)
    
    # Define number of discretization points for each dimension of the calibration grid
    ngrid = 15
    
    # Determine total number of calibration points
    #   There are 3 dimensions for the calibration grid -- beta, sigma, gamma
    npoint = ModelCalibrator.ngrid ^ 3
    
    # Define calibration point directory and calibration point file path
    pointdir  = os.path.join(PathFinder.getSourceDir(), 'CalibrationPoints')
    pointfile = lambda ipoint: os.path.join(ModelCalibrator.pointdir, 'point%05d.mat' % ipoint)
    
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
        lb['beta'] = 0.920
        lb['gamma'] = 0.150
        lb['sigma'] = 1.20
        ub['beta'] = 1.050
        ub['gamma'] = 0.900
        ub['sigma'] = 9.00
        
        # Construct vectors of parameter values
        v['beta']  = np.linspace(lb['beta'], ub['beta'], num=ModelCalibrator.ngrid)
        v['gamma'] = np.linspace(lb['gamma'], ub['gamma'], num=ModelCalibrator.ngrid)
        v['sigma'] = np.logspace(np.log10(lb['sigma']), np.log10(ub['sigma']), num=ModelCalibrator.ngrid)
        
        # Generate calibration points as unique combinations of parameter values
        # TBD: figure out what this grid thing is about
        (grid['beta'], grid['gamma'], grid['sigma']) = np.meshgrid(v['beta'], v['gamma'], v['sigma'])
        for ipoint in range(ModelCalibrator.npoint):
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
        
        # Initialize vectors of parameters, targets, and solution conditions
        for o in ModelCalibrator.paramlist:
            paramv[o] = np.empty(shape=(1, ModelCalibrator.npoint))
        for o in ModelCalibrator.targetlist:
            targetv[o] = np.empty(shape=(1, ModelCalibrator.npoint))
        solved = np.zeros(shape=(1, ModelCalibrator.npoint))
        
        # Load and consolidate calibration points
        for i in range(ModelCalibrator.npoint):
            
            print('Reading calibration point %5d of %5d\n' % (i, ModelCalibrator.npoint))
            
            with open(ModelCalibrator.pointfile(i)) as f:
                s['params'] = pickle.load(f)
                s['targets'] = pickle.load(f)
                s['solved'] = pickle.load(f)
                
            for o in ModelCalibrator.paramlist:
                paramv[o][i] = s['params'][o]
            for o in ModelCalibrator.targetlist:
                targetv[o][i] = s['targets'][o]
            solved(i) = s['solved']
        
        # Save consolidated points to calibration output directory
        with open(os.path.join(outputdir, 'calibration.mat')) as f:
            pickle.dump(paramv)
            pickle.dump(targetv)
            pickle.dump(solved)
            
        # Initialize plot of calibration point solution conditions
        fig = figure; ax = axes;
        
        # Determine colors
        cv = np.zeros((ModelCalibrator.npoint, 3))
        devs = min(abs(targetv.captoout(solved) - 3).^0.5, 1);
        cv( solved, :) = [devs, ones(size(devs)), devs]*180/256;        % Gray to green
        cv(~solved, :) = repmat([200/256, 0, 0], [sum(~solved), 1]);    % Red
        
        % Plot solution conditions
        scatter3(paramv.beta, paramv.gamma, paramv.sigma, 40, cv, 'filled');
        
        % Format axes
        axis(ax, 'tight'), box(ax, 'on'), grid(ax, 'minor'), view(3), pbaspect([1,1,1])
        xlabel('beta' ), ax.XScale = 'linear'; ax.XTickMode = 'manual'; ax.XTick = linspace(ax.XLim(1)       , ax.XLim(2)       , 3); ax.XMinorTick = 'off';
        ylabel('gamma'), ax.YScale = 'linear'; ax.YTickMode = 'manual'; ax.YTick = linspace(ax.YLim(1)       , ax.YLim(2)       , 3); ax.YMinorTick = 'off';
        zlabel('sigma'), ax.ZScale = 'log'   ; ax.ZTickMode = 'manual'; ax.ZTick = logspace(log10(ax.ZLim(1)), log10(ax.ZLim(2)), 3); ax.ZMinorTick = 'off';
        grid(ax, 'minor')
        
        % Save plot to calibration output directory
        savefig(fig, fullfile(outputdir, 'conditions.fig'));
        
    end
    
    
    
    %%
    %   Single loop to calibrate on modelunit_dollar targets
    function [ targets, modelunit_dollar, is_solved ] = calibrate_dollar( gridpoint )

        % Set target = $gdp/adult
        %     from Alex $79.8k for 2016
        %     REM: In moment_targets, 
        %        col 1 = varname, col 2 = value, col 3 = description 
        target_outperHH_index = find( strcmp( ModelCalibrator.moment_targets(:, 1), 'outperHH' ), 1 );
        target_outperHH       = cell2mat( ModelCalibrator.moment_targets( target_outperHH_index, 2 ) );
        
        % Set initial modelunit_dollar.
        % In the future, we could apply a heuristic better initial guess.
        modelunit_dollar    = 4.0e-05;  

        tolerance           = 0.01;    % as ratio 
        err_size            = 1;
        iter_num            = 1;
        iter_max            = 8;       % iterations for modelunit_dollar

        while (( err_size > tolerance ) && (iter_num <= iter_max) )

            % Create Scenario to run
            scenario    = Scenario( struct( 'economy'           , 'steady'          ...
                                         ,  'beta'              , gridpoint.beta    ...
                                         ,  'gamma'             , gridpoint.gamma   ...
                                         ,  'sigma'             , gridpoint.sigma   ...
                                         ,  'modelunit_dollar'  , modelunit_dollar  ...
                                         ,  'bequest_phi_1'     , 0                 ...
                                   ));
            save_dir    = ModelSolver.solve( scenario );

            % find target -- $gdp/pop
            s_paramsTargets = load( fullfile(save_dir, 'paramsTargets.mat' ) );
            run_outperHH    = s_paramsTargets.outperHH;
            
            err_size        = abs( run_outperHH/target_outperHH - 1 );
            fprintf( '...MODELUNIT_DOLLAR iteration %u   error=%f\n ', iter_num, err_size );

            % package up answer
            targets         = struct(   'savelas',      s_paramsTargets.savelas  ...
                                     ,  'labelas',      s_paramsTargets.labelas  ...
                                     ,  'captoout',     s_paramsTargets.captoout ... 
                                     ,  'outperHH',     run_outperHH             );                       

            % Update by percent shift, reduced a bit as number of 
            % iterations increases. This approach slows the update rate
            % in case of slow convergence -- we're usually bouncing around then.
            exp_reduce        = max( 0.5, 1.0 - iter_num *0.07 );
            modelunit_dollar = modelunit_dollar*((run_outperHH/target_outperHH)^exp_reduce);

            % Find if converged
            %    This only needs to be done after the loop, but
            %    we're about to wipe out the run's files.
            s_dynamics     = load( fullfile(save_dir, 'dynamics.mat' ) );
            is_converged   = s_dynamics.is_converged;

            % Delete save directory along with parent directories
            rmdir(fullfile(save_dir, '..', '..'), 's')
            clear( 's_ParamsTargets' ); clear( 's_dynamics' );
            
            iter_num = iter_num + 1;
        end % while
   
        % Keep last successful run with modelunit_dollar
        modelunit_dollar   = scenario.modelunit_dollar;
        
        % Check solution condition.
        % Stable solution identified as:
        %  1. Robust solver convergence rate
        %  2. modelunit_dollar convergence
        is_solved = is_converged && ( err_size <= tolerance );
        if( iter_num > iter_max )
           fprintf( '...MODELUNIT_DOLLAR -- max iterations (%u) reached.\n', iter_max );
        end 

    end % calibrate_dollar

    %%
    %  Print moments info on a particular steady state
    function outstr = report_moments( save_dir, targets )

        delimiter   = [char(13) char(10)];  % end-of-line 
        
        filepath    = fullfile(save_dir, 'iterations.csv');
        T           = readtable(filepath, 'Format', '%u%f');
        iters       = table2array(T(:,1));  
        iterations  = iters(end);

        s_dynamics      = load( fullfile(save_dir, 'dynamics.mat' ) );
        s_paramsTargets = load( fullfile(save_dir, 'paramsTargets.mat' ) );
        s_markets       = load( fullfile(save_dir, 'market.mat' ) );
        
        % Define some helper vars for clarity
        pop             = s_dynamics.pops;    
        gdp             = s_dynamics.outs;    
        dollar          = 1/s_paramsTargets.modelunit_dollar; 

        if( ~exist('targets', 'var') || isempty(targets) )
            targets = ModelCalibrator.moment_targets;
            targets(end+1,:) = {'labelas', 1, 'Labor elasticity'};
            targets(end+1,:) = {'savelas', 1, 'Savings elasticity'};
        end

        % helper function to format results
        myTargetPrint = @( lbl, modelResult, targetResult ) ...
                    sprintf( '   %20s = %f (%f) error = %0.1f%%', lbl, modelResult, targetResult, (modelResult/targetResult - 1)*100.0 );
        myParamPrint  = @( lbl, modelInput ) ...
                    sprintf( '   %20s = %f', lbl, modelInput );
        
        % Make PARAMS section
        params  = {'beta'   , s_paramsTargets.beta;
                   'sigma'  , s_paramsTargets.sigma;
                   'gamma'  , s_paramsTargets.gamma;
                   'model$' , s_paramsTargets.modelunit_dollar; };
        
        param_part = sprintf('%s   PARAMS%s', delimiter, delimiter );
        for i = 1:length(params)             
            result  = cell2mat( params(i, 2) );
            lbl     = cell2mat( params(i, 1) );
            
            line    = myParamPrint( lbl, result );
            param_part = sprintf( '%s%s%s', param_part, line, delimiter );
        end        
        
        % Make structure for results 
        %   targets has been passed in (or set to default)
        model_results = {'r'          , s_markets.MPKs;
                         'PIT'        , s_dynamics.pits/gdp;
                         'SSTax'      , s_dynamics.ssts/gdp;
                         'KbyY'       , s_paramsTargets.captoout;
                         'outperHH'   , gdp*dollar/pop;
                         'labelas'    , s_paramsTargets.labelas;
                         'savelas'    , s_paramsTargets.savelas;};
        
        % Make TARGETS section
        target_part = sprintf('%s   TARGETS%s', delimiter, delimiter );
        for i = 1:length(model_results(:, 1))  
            m_index = find( strcmp( targets(:, 1), model_results(i,1) ), 1 );
            target  = cell2mat( targets( m_index, 2) );
            lbl     = cell2mat( targets( m_index, 3) );
            result  = cell2mat( model_results(i, 2) );
            
            line        = myTargetPrint( lbl, result , target );
            target_part = sprintf( '%s%s%s', target_part, line, delimiter );
        end
        
        % Make convergence part
        if( s_dynamics.is_converged )
            s_iter = sprintf('Converged in %u iterations', iterations );
        else
            s_iter = sprintf('DID NOT converge in %u iterations.', iterations );
        end
        converge_part = sprintf( '%s CONVERGENCE: %s %s', delimiter, s_iter, delimiter );      

        % Concatentate for full report
        outstr = sprintf( '%s%s%s', param_part, target_part, converge_part );
    end % report_moments

    %%
    %   Make a report of various moments for the 16 baselines
    function [] = report_baseline_moments()
        
        outputfilename      = fullfile(PathFinder.getSourceDir(), 'BaselineMoments.txt');
        fileID              = fopen(outputfilename,'w');
        
        fprintf( fileID, '-------------BASELINE MOMENTS-------------' );
        fprintf( fileID, '%s \r\n', datestr(now));
        
        % load the matrix and get inverter function
        [~, f_invert] = ParamGenerator.invert();
        
        for labelas = 0.25:0.25:1.0
            for savelas = 0.25:0.25:1.0
                target = struct('labelas', labelas, 'savelas', savelas);
                fprintf( fileID, '\r\nBASELINE labor elas = %0.2f  savings elas = %0.2f \r\n', labelas, savelas ); 
                inverse = f_invert(target);
                
                scenario = Scenario(struct('economy','steady','beta',inverse.beta,'gamma',inverse.gamma,...
                                           'sigma',inverse.sigma,'modelunit_dollar',inverse.modelunit_dollar,...
                                           'bequest_phi_1',0));
                
                save_dir = ModelSolver.solve(scenario);
                
                targets  = ModelCalibrator.moment_targets;
                targets(end+1,:) = {'labelas', labelas, 'Labor elasticity'};
                targets(end+1,:) = {'savelas', savelas, 'Savings elasticity'};
                outstr   = ModelCalibrator.report_moments( save_dir, targets );
                fprintf( fileID, '%s \r\n', outstr );
                fprintf( fileID, '-------------------------------------\r\n' );                
            end % for 
        end % for 
        fprintf( fileID, ' ==== DONE ===== \r\n' );    
        fclose( fileID );
    end % function report_baseline_moments

    %% 
    %   Adjust calibration grid boundaries
function [grid_beta, grid_gamma, grid_sigma] = adjust_grid()
    
    epsilon = 1e-4;
    [~, f_invert] = ParamGenerator.invert();
    green = [0 180/256 0];
    cv    = repmat(reshape(green, [3,1]), [1,16])';%zeros(16,3);
    grid_beta  = [0.950 1.100];
    grid_gamma = [0.150 0.900];
    grid_sigma = [1.200 9.000];
    delta_beta  = zeros(16,2);
    delta_gamma = zeros(16,2);
    delta_sigma = zeros(16,2);
    labelasv = zeros(16,1);
    savelasv = zeros(16,1);
    iter = 1;
    
    for labelas = 0.25:0.25:1.00
        for savelas = 0.25:0.25:1.00
            
            target  = struct('labelas', labelas, 'savelas', savelas);
            inverse = f_invert(target);
            delta_beta(iter,:)  = [(inverse.beta  - grid_beta(1))/inverse.beta , (grid_beta(2) - inverse.beta)/inverse.beta  ];
            delta_gamma(iter,:) = [(inverse.gamma - grid_gamma(1))/inverse.gamma, (grid_gamma(2) - inverse.gamma)/inverse.gamma];
            delta_sigma(iter,:) = [(inverse.sigma - grid_sigma(1))/inverse.sigma, (grid_sigma(2) - inverse.sigma)/inverse.sigma];
            delta = min(min(min(delta_beta(iter,:), delta_gamma(iter,:)), delta_sigma(iter,:)));
            if (delta <= epsilon)
                cv(iter,:) = [1, delta, 0]*200/256;
            end
            labelasv(iter,1) = labelas;
            savelasv(iter,1) = savelas;
            iter = iter + 1;
        end
    end
            
    % Plot
    figure
    scatter(labelasv, savelasv, 40, cv, 'filled');
    xlabel('labor elasticity'  ,'FontSize',13); set(gca,'XTick',0:0.25:1.00) 
    ylabel('savings elasticity','FontSize',13); set(gca,'YTick',0:0.25:1.00)
    grid on

    % Adjust grids
    if (min(delta_beta(:,1))  <= epsilon), grid_beta(1)  = 0.9*grid_beta(1) ; end
    if (min(delta_beta(:,2))  <= epsilon), grid_beta(2)  = 1.1*grid_beta(2) ; end
    if (min(delta_gamma(:,1)) <= epsilon), grid_gamma(1) = 0.9*grid_gamma(1); end
    if (min(delta_gamma(:,2)) <= epsilon), grid_gamma(2) = 1.1*grid_gamma(2); end
    if (min(delta_sigma(:,1)) <= epsilon), grid_sigma(1) = max(1.01, 0.9*grid_sigma(1)); end
    if (min(delta_sigma(:,2)) <= epsilon), grid_sigma(2) = 1.1*grid_sigma(2); end
    
    
end

end % methods


end