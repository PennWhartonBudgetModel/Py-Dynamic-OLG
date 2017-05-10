%%
% Dynamic model baseline parameter calibrator.
% 
%%


classdef modelCalibrator

properties (Constant)
    
    % Define list of parameters to be calibrated
    %    NOTE: modelunit_dollar is added as a param during solve_batch
    paramlist = {'beta', 'gamma', 'sigma'};
    nparam = length(modelCalibrator.paramlist);
    
    % Define list of elasticities to be targeted
    elaslist = {'captoout', 'labelas', 'savelas', 'outperHH'};
    nelas = length(modelCalibrator.elaslist);
    
    % Define number of discretization points for each parameter
    npoint = 10; % 30;
    
    % Define number of parameter sets per batch
    batchsize = 10;
    
    % Determine number of parameter sets and number of batches
    nset   = modelCalibrator.npoint ^ modelCalibrator.nparam;
    nbatch = ceil(modelCalibrator.nset / modelCalibrator.batchsize);
    
    % Define batch directory and batch file path
    batch_dir  = fullfile(dirFinder.source(), 'Batches');
    batch_file = @(ibatch) fullfile(modelCalibrator.batch_dir, sprintf('batch%05d.mat', ibatch));
    
end

methods (Static)
    
    % Define batches of parameter sets
    function [] = define_batches()
        
        % Specify parameter lower and upper bounds
        lb.beta = 0.990; lb.gamma = 0.150; lb.sigma =  1.50;
        ub.beta = 1.180; ub.gamma = 0.900; ub.sigma = 32.00;
        
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
            paramsets(i).modelunit_dollars = 2.88e-05;
        end

        % Clear or create batch directory
        if exist(modelCalibrator.batch_dir, 'dir'), rmdir(modelCalibrator.batch_dir, 's'), end, mkdir(modelCalibrator.batch_dir)
        
        % Extract and save batches of parameter sets
        for ibatch = 1:modelCalibrator.nbatch
            
            shift = (ibatch-1)*modelCalibrator.batchsize;
            parambatch = paramsets(shift+1 : min(shift+modelCalibrator.batchsize, end)); %#ok<NASGU>
            
            save(modelCalibrator.batch_file(ibatch), 'parambatch')
            
        end
        
    end
    
    
    % Solve dynamic model steady states for all parameter sets in a batch
    function [] = solve_batch(ibatch)
        
        % Load batch of parameter sets
        s = load(modelCalibrator.batch_file(ibatch)); params = s.parambatch; clear('s')
        
        parfor i = 1:length(params)
            % Calibrate steady state on modelunit_dollar
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
        for p = modelCalibrator.paramlist , paramv.(p{1})  = []; end
        for e = modelCalibrator.targetlist, targetv.(e{1}) = []; end
        solved = [];
        
        % Load and consolidate solutions from batch files
        for ibatch = 1:modelCalibrator.nbatch
            
            fprintf('Reading batch %5d of %5d\n', ibatch, modelCalibrator.nbatch)
            
            s = load(fullfile(modelCalibrator.batch_dir, sprintf('batch%05d.mat', ibatch)));
            
            for p = modelCalibrator.paramlist , paramv.(p{1})   = [paramv.(p{1}) , s.params.(p{1}) ]; end
            for e = modelCalibrator.targetlist, targetv.( e{1}) = [targetv.(e{1}), s.targets.(e{1})]; end
            solved = [solved, s.solved]; %#ok<AGROW>
            
        end
        solved = boolean(solved); %#ok<NASGU>
        
        % Save solutions into calibration file
        % (Note that this file should be copied to the Parameter directory and committed anew following a completed calibration)
        save(fullfile(dirFinder.saveroot(), 'calibration.mat'), 'paramv', 'rargetv', 'solved');
        
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
        for p = modelCalibrator.paramlist
            interp.(p{1}) = scatteredInterpolant(targetv.captoout(solved)', targetv.labelas(solved)', targetv.savelas(solved)', paramv.(p{1})(solved)', 'nearest');
        end
        
        % Construct elasticity inverter by consolidating inverse interpolants
        function [inverse] = f_(target)
            captoout = 3;
            for p_ = modelCalibrator.paramlist, inverse.(p_{1}) = interp.(p_{1})(captoout, target.labelas, target.savelas); end
        end
        f = @f_;
        
        % Invert target
        if exist('target', 'var'), inverse = f(target); else, inverse = struct(); end
        
    end

 
    
    
%%
%   Make a report of various moments for the 16 baselines
function [] = view_baseline_moments()
    
    moment_targets      = struct('r', 0.07, 'PIT', 0.08, 'SSTax', 0.05 ...
                            , 'KbyY', 3.0, 'GDPperHH', 1e5);
    
    outputfilename      = fullfile(dirFinder.saveroot(), 'BaselineMoments.txt');
    fileID              = fopen(outputfilename,'w');

    % helper function to print results
    myPrint = @( modelResult, targetResult ) ...
                fprintf( fileID, ' %f (%f) error = %0.1f%% \r\n', modelResult, targetResult, (modelResult/targetResult - 1)*100.0 );
    
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
            
            % TODO: this is temp: inverse doesn't have modelunit_dollars yet
            inverse(:).modelunit_dollars = 2.88e-05;
            
            save_dir = dynamicSolver.steady(inverse);
            % find iterations
            filepath    = fullfile(save_dir, 'iterations.csv');
            T           = readtable(filepath, 'Format', '%u%f');
            iterations  = table2array(T(:,1));            
            fprintf( fileID, '   Computed in %u iterations \r\n', iterations(end) );
            
            s_dynamics = load( fullfile(save_dir, 'dynamics.mat' ) );
            s_elas     = load( fullfile(save_dir, 'elasticities.mat' ) );
            s_markets  = load( fullfile(save_dir, 'market.mat' ) );
            pop        = s_dynamics.pops;    % helper variable
            gdp        = s_dynamics.outs;    % helper variable
            dollar     = 1/s_elas.modelunit_dollars; % helper variable
            
            % print reports
            fprintf( fileID, '   beta               = %f \r\n', inverse.beta );
            fprintf( fileID, '   sigma              = %f \r\n', inverse.sigma );
            fprintf( fileID, '   gamma              = %f \r\n', inverse.gamma );
            
            %  this will come from inverse struct
            fprintf( fileID, '   model$             = %f \r\n', s_elas.modelunit_dollars );
            
            fprintf( fileID, '   lab elas           = ' ); myPrint( s_elas.labelas, labelas );
            fprintf( fileID, '   sav elas           = ' ); myPrint( s_elas.savelas, savelas );
            fprintf( fileID, '   K/Y                = ' ); myPrint( s_elas.captoout, moment_targets.KbyY );
            fprintf( fileID, '   SStax/Y            = ' ); myPrint( s_dynamics.ssts/gdp, moment_targets.SSTax );
            fprintf( fileID, '   PIT/Y              = ' ); myPrint( s_dynamics.pits/gdp, moment_targets.PIT );
            fprintf( fileID, '   r                  = ' ); myPrint( s_markets.caprates, moment_targets.r );
            fprintf( fileID, '   $GDP/worker        = ' ); myPrint( gdp*dollar/pop, moment_targets.GDPperHH );
            fprintf( fileID, '-------------------------------------\r\n' );
        end % for 
    end % for 
    fprintf( fileID, ' ==== DONE ===== \r\n' );    
    fclose( fileID );
end % function view_baseline_moments

%%
%   Single loop to calibrate on modelunit_dollar targets
function [ targets, modelunit_final, is_solved ] = calibrate_dollar( basedef )

    % Set target = $gdp/worker
    target              = 1.14e05;
    modelunit_dollars   = basedef.modelunit_dollars;
    
    tolerance           = 0.01;    % as ratio 
    err_size            = 1;
    iter_num            = 1;
 
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
            
        % Update by percent shift, reduced a bit as it takes 
        % longer to converge. This approach updates faster for 
        % "normal" points, which converge quickly, and updates
        % slower for "odd" points, which converge slower. 
        exp_reduce        = max( 0.5, 1.0 - iter_num *0.07 );
        modelunit_dollars = modelunit_dollars*((actual_value/target)^exp_reduce);
        
        % Delete save directory along with parent directories
        rmdir(fullfile(save_dir, '..', '..'), 's')

        iter_num = iter_num + 1;

        if( iter_num > 8 )
            is_solved = false;
            fprintf( '...MODELUNIT_DOLLAR -- max iterations reached.' );
            break;
        end 
        
    end % while

end % calibrate_dollar


end % methods


end