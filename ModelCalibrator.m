%%
% Dynamic model baseline parameter calibrator.
% 
%%
classdef ModelCalibrator

properties (Constant)
    
    % Define list of parameters which define the steady state
    paramlist = {'beta', 'gamma', 'sigma', 'modelunit_dollar'};
    
    % Define list of targets
    targetlist  = {'captoout', 'labelas', 'savelas', 'outperHH'};
    ntarget     = length(ModelCalibrator.targetlist);
    
    % Define number of discretization points for each parameter
    npoint = 10;
    
    % Define number of parameter sets per batch
    batchsize = 1;
    
    % Determine number of parameter sets and number of batches
    %   REM: There are 3 dimensions for the calibration grid:
    %              beta, sigma, gamma
    nset   = ModelCalibrator.npoint ^ 3;
    nbatch = ceil(ModelCalibrator.nset / ModelCalibrator.batchsize);
    
    % Define batch directory and batch file path
    batch_dir  = fullfile(PathFinder.getSourceDir(), 'Batches');
    batch_file = @(ibatch) fullfile(ModelCalibrator.batch_dir, sprintf('batch%05d.mat', ibatch));
    
    % Define the moment targets for the reports on how we did
    %    Cell array: Col1=varname, Col2=value, Col3=description
    moment_targets   = {'r',        0.05,   'Return on capital';
                        'PIT',      0.08,   'PIT/GDP';
                        'SSTax',    0.05,   'SSTax/GDP';
                        'KbyY',     3.0,    'Capital/GDP';
                        'outperHH', 7.98e4, 'GDP$/adult';};
    
 end

methods (Static)
    
    % Define batches of parameter sets
    function [] = define_batches()
        
        % Specify parameter lower and upper bounds
        lb.beta = 0.950; lb.gamma = 0.150; lb.sigma = 1.20;
        ub.beta = 1.100; ub.gamma = 0.900; ub.sigma = 9.00;
        
        % Construct vectors of parameter values
        v.beta  = linspace(lb.beta        , ub.beta        , ModelCalibrator.npoint);
        v.gamma = linspace(lb.gamma       , ub.gamma       , ModelCalibrator.npoint);
        v.sigma = logspace(log10(lb.sigma), log10(ub.sigma), ModelCalibrator.npoint);
        
        % Generate grid parameter sets as unique combinations of parameter values
        [grid.beta, grid.gamma, grid.sigma] = ndgrid(v.beta, v.gamma, v.sigma);
        for i = 1:ModelCalibrator.nset
            for p = {'beta', 'gamma', 'sigma'}
                paramsets(i).(p{1}) = grid.(p{1})(i); %#ok<AGROW>
            end 
        end 
 
        % Clear or create batch directory
        if exist(ModelCalibrator.batch_dir, 'dir'), rmdir(ModelCalibrator.batch_dir, 's'), end, mkdir(ModelCalibrator.batch_dir)
        
        % Extract and save batches of parameter sets
        for ibatch = 1:ModelCalibrator.nbatch
            
            shift = (ibatch-1)*ModelCalibrator.batchsize;
            params = paramsets(shift+1 : min(shift+ModelCalibrator.batchsize, end)); %#ok<NASGU>
            save(ModelCalibrator.batch_file(ibatch), 'params')
            
        end
        
    end
    
    
    % Solve dynamic model steady states for all parameter sets in a batch
    function [] = solve_batch(ibatch)
        
        % Load batch of parameter sets
        s = load(ModelCalibrator.batch_file(ibatch)); 
        params = s.params; clear('s');
        nparamsets = length(params);
        
        parfor i = 1:nparamsets  
            % Calibrate steady state on modelunit_dollar
            [ targets(i), modelunit_dollar(i), solved(i) ] = ModelCalibrator.calibrate_dollar( params(i) ); %#ok<NASGU,PFOUS,ASGLU>
        end
        
        % Add modelunit_dollar to the params
        for i = 1:nparamsets
            params(i).modelunit_dollar = modelunit_dollar(i);
        end
        
        % Save elasticity sets and solution conditions to batch file
        % Rem: Overwrite parambatch since modelunit_dollar will have been modified
        save(ModelCalibrator.batch_file(ibatch), '-append' ...
                            , 'targets', 'solved', 'params' );
        
    end
    
    
    % Consolidate solutions from batch runs
    function [] = consolidate_batches(clean)
        
        % Initialize vectors of parameters, elasticities, and solution conditions
        for p = ModelCalibrator.paramlist , paramv.(p{1})  = []; end
        for e = ModelCalibrator.targetlist, targetv.(e{1}) = []; end
        solved = [];
        
        % Load and consolidate solutions from batch files
        for ibatch = 1:ModelCalibrator.nbatch
            
            fprintf('Reading batch %5d of %5d\n', ibatch, ModelCalibrator.nbatch)
            
            s = load(fullfile(ModelCalibrator.batch_dir, sprintf('batch%05d.mat', ibatch)));
            
            for p = ModelCalibrator.paramlist , paramv.(p{1})   = [paramv.(p{1}) , s.params.(p{1}) ]; end
            for e = ModelCalibrator.targetlist, targetv.( e{1}) = [targetv.(e{1}), s.targets.(e{1})]; end
            solved = [solved, s.solved]; %#ok<AGROW>
            
        end
        solved = boolean(solved); %#ok<NASGU>
        
        % Save solutions to new calibration file in new input version directory
        outdir = PathFinder.getCalibrationOutputDir();
        if exist(outdir, 'dir'), rmdir(outdir, 's'), end, mkdir(outdir)
        save(fullfile(outdir, 'calibration.mat'), 'paramv', 'targetv', 'solved');
        
        % Delete batch directory
        if (exist('clean', 'var') && clean), rmdir(ModelCalibrator.batch_dir, 's'), end
        
    end
    
    
    % Plot calibration solution conditions
    function [] = plot_conditions()
        
        % Load calibration solutions
        cal_dir = PathFinder.getCalibrationInputDir();
        s       = load(fullfile(cal_dir, 'calibration.mat'));
        paramv  = s.paramv;
        targetv = s.targetv ;
        solved  = s.solved;
        
        % Initialize figure
        figure(1); ax = axes;
        
        % Determine colors
        cv = zeros(ModelCalibrator.nset, 3);
        devs = min(abs(targetv.captoout(solved)' - 3).^0.5, 1);
        cv( solved, :) = [devs, ones(size(devs)), devs]*180/256;        % Gray to green
        cv(~solved, :) = repmat([200/256, 0, 0], [sum(~solved), 1]);    % Red
        
        % Plot parameter sets
        scatter3(paramv.beta, paramv.gamma, paramv.sigma, 40, cv, 'filled');
        
        % Format axes
        axis(ax, 'tight'), box(ax, 'on'), grid(ax, 'minor'), view(3), pbaspect([1,1,1])
        xlabel('beta' ), ax.XScale = 'linear'; ax.XTickMode = 'manual'; ax.XTick = linspace(ax.XLim(1)       , ax.XLim(2)       , 3); ax.XMinorTick = 'off';
        ylabel('gamma'), ax.YScale = 'linear'; ax.YTickMode = 'manual'; ax.YTick = linspace(ax.YLim(1)       , ax.YLim(2)       , 3); ax.YMinorTick = 'off';
        zlabel('sigma'), ax.ZScale = 'log'   ; ax.ZTickMode = 'manual'; ax.ZTick = logspace(log10(ax.ZLim(1)), log10(ax.ZLim(2)), 3); ax.ZMinorTick = 'off';
        grid(ax, 'minor')
        
    end
    
    
    % Invert target elasticities, constructing a reusable elasticity inverter in the process
    function [inverse, f] = invert(target)
        
        % Load calibration solutions
        cal_dir = PathFinder.getCalibrationInputDir();
        s       = load(fullfile(cal_dir, 'calibration.mat'));
        paramv  = s.paramv;
        targetv = s.targetv ;
        solved  = s.solved;
        
        % Construct inverse interpolants that map individual targets to parameters
        for p = ModelCalibrator.paramlist
            interp.(p{1}) = scatteredInterpolant(targetv.captoout(solved)', targetv.labelas(solved)', targetv.savelas(solved)', paramv.(p{1})(solved)', 'nearest');
        end
        
        % Construct elasticity inverter by consolidating inverse interpolants
        function [inverse] = f_(target)
            captoout = 3;
            for p_ = ModelCalibrator.paramlist, inverse.(p_{1}) = interp.(p_{1})(captoout, target.labelas, target.savelas); end
        end
        f = @f_;
        
        % Invert target
        if exist('target', 'var'), inverse = f(target); else, inverse = struct(); end
        
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
        model_results = {'r'          , s_markets.caprates;
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
        [~, f_invert] = ModelCalibrator.invert();               
        
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
    [~,f_invert] = ModelCalibrator.invert();
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