%%
% Dynamic model baseline parameter calibrator.
% 
%%


classdef modelCalibrator

properties (Constant)
    
    % Define list of parameters to be calibrated
    %    NOTE: modelunit_dollar is added as a param during solve_batch
    %          since paramlist defines the grid (and modelunit_dollars  
    %          is not in grid).
    paramlist = {'beta', 'gamma', 'sigma'};
    nparam    = length(modelCalibrator.paramlist);
    
    % Define list of parameters which define the steady state
    %    This is with modelunit_dollars.
    fullparamlist = {'beta', 'gamma', 'sigma', 'modelunit_dollars'};
    nfullparam    = length(modelCalibrator.fullparamlist);
    
    % Define list of elasticities to be targeted
    targetlist  = {'captoout', 'labelas', 'savelas', 'outperHH'};
    ntarget     = length(modelCalibrator.targetlist);
    
    % Define number of discretization points for each parameter
    npoint = 30;
    
    % Define number of parameter sets per batch
    batchsize = 4;
    
    % Determine number of parameter sets and number of batches
    nset   = modelCalibrator.npoint ^ modelCalibrator.nparam;
    nbatch = ceil(modelCalibrator.nset / modelCalibrator.batchsize);
    
    % Define batch directory and batch file path
    batch_dir  = fullfile(dirFinder.source(), 'Batches');
    batch_file = @(ibatch) fullfile(modelCalibrator.batch_dir, sprintf('batch%05d.mat', ibatch));
    
    % Define the calibration targets:
    %   Currently, just output per adult
    target_outperHH  = 7.98e4;
    
    % Define the moment targets for the reports on how we did
    moment_targets   = struct('r', 0.05, 'PIT', 0.08, 'SSTax', 0.05 ...
                            , 'KbyY', 3.0, 'GDPperHH', 7.98e4);
    
 end

methods (Static)
    
    % Define batches of parameter sets
    function [] = define_batches()
        
        % Specify parameter lower and upper bounds
        lb.beta = 0.990; lb.gamma = 0.150; lb.sigma =  1.50;
        ub.beta = 1.170; ub.gamma = 0.900; ub.sigma = 30.00;
        
        % Construct vectors of parameter values
        v.beta  = linspace(lb.beta        , ub.beta        , modelCalibrator.npoint);
        v.gamma = linspace(lb.gamma       , ub.gamma       , modelCalibrator.npoint);
        v.sigma = logspace(log10(lb.sigma), log10(ub.sigma), modelCalibrator.npoint);
        
        % Generate parameter sets as unique combinations of parameter values
        [grid.beta, grid.gamma, grid.sigma] = ndgrid(v.beta, v.gamma, v.sigma);
        for i = 1:modelCalibrator.nset, for p = modelCalibrator.paramlist, paramsets(i).(p{1}) = grid.(p{1})(i); end, end %#ok<AGROW>
 
        % Add modelunit_dollars to the params
        % and set to a default initial value.
        % In the future, we could apply a heuristic better initial guess.
        for i = 1:modelCalibrator.nset
            paramsets(i).modelunit_dollars = 4.199e-05;
        end

        % Clear or create batch directory
        if exist(modelCalibrator.batch_dir, 'dir'), rmdir(modelCalibrator.batch_dir, 's'), end, mkdir(modelCalibrator.batch_dir)
        
        % Extract and save batches of parameter sets
        for ibatch = 1:modelCalibrator.nbatch
            
            shift = (ibatch-1)*modelCalibrator.batchsize;
            params = paramsets(shift+1 : min(shift+modelCalibrator.batchsize, end)); %#ok<NASGU>
            save(modelCalibrator.batch_file(ibatch), 'params')
            
        end
        
    end
    
    
    % Solve dynamic model steady states for all parameter sets in a batch
    function [] = solve_batch(ibatch)
        
        % Load batch of parameter sets
        s = load(modelCalibrator.batch_file(ibatch)); params = s.params; clear('s')
        
        parfor i = 1:length(params)
            % Calibrate steady state on modelunit_dollars
            [ targets(i), modelunit_dollars, solved(i) ] = modelCalibrator.calibrate_dollar( params(i) );
            % overwrite parambatch w/ new modelunit_dollars
            params(i).modelunit_dollars = modelunit_dollars;
        end
        
        
        % Save elasticity sets and solution conditions to batch file
        % Rem: Overwrite parambatch since modelunit_dollar will have been modified
        save(modelCalibrator.batch_file(ibatch), '-append' ...
                            , 'targets', 'solved', 'params' );
        
    end
    
    
    % Consolidate solutions from batch runs
    function [] = consolidate_batches(clean)
        
        % Initialize vectors of parameters, elasticities, and solution conditions
        for p = modelCalibrator.fullparamlist, paramv.(p{1})  = []; end
        for e = modelCalibrator.targetlist   , targetv.(e{1}) = []; end
        solved = [];
        
        % Load and consolidate solutions from batch files
        for ibatch = 1:modelCalibrator.nbatch
            
            fprintf('Reading batch %5d of %5d\n', ibatch, modelCalibrator.nbatch)
            
            s = load(fullfile(modelCalibrator.batch_dir, sprintf('batch%05d.mat', ibatch)));
            
            for p = modelCalibrator.fullparamlist, paramv.(p{1})   = [paramv.(p{1}) , s.params.(p{1}) ]; end
            for e = modelCalibrator.targetlist   , targetv.( e{1}) = [targetv.(e{1}), s.targets.(e{1})]; end
            solved = [solved, s.solved]; %#ok<AGROW>
            
        end
        solved = boolean(solved); %#ok<NASGU>
        
        % Save solutions into calibration file
        % (Note that this file should be copied to the Parameter directory and committed anew following a completed calibration)
        save(fullfile(dirFinder.saveroot(), 'calibration.mat'), 'paramv', 'targetv', 'solved');
        
        % Delete batch directory
        if (exist('clean', 'var') && clean), rmdir(modelCalibrator.batch_dir, 's'), end
        
    end
    
    
    % Plot calibration solution conditions
    function [] = plot_conditions()
        
        % Load calibration solutions
        s = load(fullfile(dirFinder.param(), 'calibration.mat'));
        paramv  = s.paramv;
        targetv = s.targetv ;
        solved  = s.solved;
        
        % Initialize figure
        figure(1); ax = axes;
        
        % Determine colors
        cv = zeros(modelCalibrator.nset, 3);
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
        s = load(fullfile(dirFinder.param(), 'calibration.mat'));
        paramv  = s.paramv;
        targetv = s.targetv ;
        solved  = s.solved;
        
        % Construct inverse interpolants that map individual elasticities to parameters
        for p = modelCalibrator.fullparamlist
            interp.(p{1}) = scatteredInterpolant(targetv.captoout(solved)', targetv.labelas(solved)', targetv.savelas(solved)', paramv.(p{1})(solved)', 'nearest');
        end
        
        % Construct elasticity inverter by consolidating inverse interpolants
        function [inverse] = f_(target)
            captoout = 3;
            for p_ = modelCalibrator.fullparamlist, inverse.(p_{1}) = interp.(p_{1})(captoout, target.labelas, target.savelas); end
        end
        f = @f_;
        
        % Invert target
        if exist('target', 'var'), inverse = f(target); else, inverse = struct(); end
        
    end


%%
%   Single loop to calibrate on modelunit_dollar targets
function [ targets, modelunit_final, is_solved ] = calibrate_dollar( basedef )

    % Set target = $gdp/adult
    %     from Alex $79.8k for 2016
    target              = modelCalibrator.target_outperHH;
    modelunit_dollars   = basedef.modelunit_dollars;
    
    tolerance           = 0.01;    % as ratio 
    err_size            = 1;
    iter_num            = 1;
    iter_max            = 8;   % iterations for modelunit_dollars
 
    while ( err_size > tolerance ) 

        modelunit_final             = modelunit_dollars;  % for func output
        basedef.modelunit_dollars   = modelunit_dollars;
        save_dir                    = dynamicSolver.steady( basedef );

        % find target -- $gdp/pop
        s_elas          = load( fullfile(save_dir, 'elasticities.mat' ) );
        actual_value    = s_elas.outperHH;
        
        err_size        = abs( actual_value/target - 1 );
        fprintf( '...MODELUNIT_DOLLAR iteration %u   error=%f\n ', iter_num, err_size );
        
        % package up answer
        targets         = struct(   'savelas',      s_elas.savelas ...
                                 ,  'labelas',      s_elas.labelas ...
                                 ,  'outperHH',     actual_value ...                       
                                 ,  'captoout',     s_elas.captoout );
        
        % Check solution condition
        % (Stable solution identified as reasonable capital-to-output ratio and robust solver convergence rate)
        is_solved = (targets.captoout > 0.5) && (size(csvread(fullfile(save_dir, 'iterations.csv')), 1) < 25); %#ok<NASGU,PFOUS>
        if( iter_num >= iter_max )
            is_solved = false;
            fprintf( '...MODELUNIT_DOLLAR -- max iterations (%u) reached.', iter_max );
            break;
        end 
        
        % Update by percent shift, reduced a bit as number of 
        % iterations increases. This approach slows the update rate
        % in case of slow convergence -- we're usually bouncing around then.
        exp_reduce        = max( 0.5, 1.0 - iter_num *0.07 );
        modelunit_dollars = modelunit_dollars*((actual_value/target)^exp_reduce);
        
        % Delete save directory along with parent directories
        rmdir(fullfile(save_dir, '..', '..'), 's')

        iter_num = iter_num + 1;

    end % while

end % calibrate_dollar

%%
%  Print moments info on a particular steady state
function outstr = report_moments( save_dir, moment_targets )
    
    filepath    = fullfile(save_dir, 'iterations.csv');
    T           = readtable(filepath, 'Format', '%u%f');
    iters       = table2array(T(:,1));  
    iterations  = iters(end);
 
    s_dynamics = load( fullfile(save_dir, 'dynamics.mat' ) );
    s_elas     = load( fullfile(save_dir, 'elasticities.mat' ) );
    s_markets  = load( fullfile(save_dir, 'market.mat' ) );
    pop        = s_dynamics.pops;    % helper variable
    gdp        = s_dynamics.outs;    % helper variable
    dollar     = 1/s_elas.modelunit_dollars; % helper variable
    
    if( isempty(moment_targets) )
        moment_targets = modelCalibrator.moments_targets;
        moment_targets.labelas = 1; moment_targets.savelas = 1; % dummy vals
    end
    
    % helper function to format results
    myPrint = @( lbl, modelResult, targetResult ) ...
                sprintf( '%s %f (%f) error = %0.1f%%', lbl, modelResult, targetResult, (modelResult/targetResult - 1)*100.0 );
    
    % make vector of outputs then concatenate 
    ss     = string(17);
    ss(1)  =          '   PARAMS ';
    ss(2)  = sprintf( '   beta               = %f', s_elas.beta );
    ss(3)  = sprintf( '   sigma              = %f', s_elas.sigma );
    ss(4)  = sprintf( '   gamma              = %f', s_elas.gamma );
    ss(5)  = sprintf( '   model$             = %f', s_elas.modelunit_dollars );
    ss(6)  = ' ';
    ss(7)  =          '   TARGETS';
    ss(8)  = myPrint( '   lab elas           =', s_elas.labelas, moment_targets.labelas );
    ss(9)  = myPrint( '   sav elas           =', s_elas.savelas, moment_targets.savelas );
    ss(10) = myPrint( '   K/Y                =', s_elas.captoout, moment_targets.KbyY );
    ss(11) = myPrint( '   SStax/Y            =', s_dynamics.ssts/gdp, moment_targets.SSTax );
    ss(12) = myPrint( '   PIT/Y              =', s_dynamics.pits/gdp, moment_targets.PIT );
    ss(13) = myPrint( '   r                  =', s_markets.caprates, moment_targets.r );
    ss(14) = myPrint( '   $GDP/adult         =', gdp*dollar/pop, moment_targets.GDPperHH );
    ss(15) = ' ';
    if( iterations == 25 )
        s_iter = 'DID NOT converge in 25 iterations.';
    else
        s_iter = sprintf('Converged in %u iterations', iterations );
    end
    ss(16) = sprintf( ' CONVERENCE: %s', s_iter );      
 
    outstr = strjoin( ss , '\r\n' );
end % report_moments

%%
%   Make a report of various moments for the 16 baselines
function [] = report_baseline_moments()
    
    outputfilename      = fullfile(dirFinder.saveroot(), 'BaselineMoments.txt');
    fileID              = fopen(outputfilename,'w');

    fprintf( fileID, '-------------BASELINE MOMENTS-------------' );
    fprintf( fileID, '%s \r\n', datestr(now));
    
    loaded = false;
    for labelas = 0.25:0.25:1.0
        for savelas = 0.25:0.25:1.0
            target = struct('labelas', labelas, 'savelas', savelas);
            fprintf( fileID, '\r\nBASELINE labor elas = %0.2f  savings elas = %0.2f \r\n', labelas, savelas ); 
            if( loaded == false )
                [inverse, f_invert] = modelCalibrator.invert(target);
            else
                inverse = f_invert(target);
            end
            
            save_dir = dynamicSolver.steady(inverse);
            
            moment_targets = modelCalibrator.moment_targets;
            moment_targets.labelas = labelas; moment_targets.savelas = savelas;
            outstr   = modelCalibrator.report_moments( save_dir, moment_targets );
            fprintf( fileID, '%s \r\n', outstr );
            fprintf( fileID, '-------------------------------------\r\n' );
        end % for 
    end % for 
    fprintf( fileID, ' ==== DONE ===== \r\n' );    
    fclose( fileID );
end % function report_baseline_moments

end % methods


end