%%
% Dynamic model baseline parameter calibrator.
% 
%%


classdef modelCalibrator

properties (Constant)
    
    % Define list of parameters to be calibrated
    paramlist = {'beta', 'gamma', 'sigma'};
    nparam = length(modelCalibrator.paramlist);
    
    % Define list of elasticities to be targeted
    elaslist = {'captoout', 'labelas', 'savelas'};
    nelas = length(modelCalibrator.elaslist);
    
    % Define number of discretization points for each parameter
    npoint = 30;
    
    % Define number of parameter sets per batch
    batchsize = 20;
    
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
        s = load(modelCalibrator.batch_file(ibatch)); parambatch = s.parambatch; clear('s')
        
        parfor i = 1:length(parambatch)
            
            % Solve steady state
            save_dir = dynamicSolver.steady(parambatch(i));
            
            % Extract elasticity set
            elasbatch(i) = load(fullfile(save_dir, 'elasticities.mat')); %#ok<PFOUS>
            
            % Check solution condition
            % (Stable solution identified as reasonable capital-to-output ratio and robust solver convergence rate)
            solvedbatch(i) = (elasbatch(i).captoout > 0.5) && (size(csvread(fullfile(save_dir, 'iterations.csv')), 1) < 25); %#ok<NASGU,PFOUS>
            
            % Delete save directory along with parent directories
            rmdir(fullfile(save_dir, '..', '..'), 's')
            
        end
        
        % Save elasticity sets and solution conditions to batch file
        save(modelCalibrator.batch_file(ibatch), '-append', 'elasbatch', 'solvedbatch')
        
    end
    
    
    % Construct elasticity inverter
    function [] = construct_inverter(clean)
        
        % Initialize vectors of parameters, elasticities, and solution conditions
        for p = modelCalibrator.paramlist, paramv.(p{1}) = []; end
        for e = modelCalibrator.elaslist , elasv.( e{1}) = []; end
        solved = [];
        
        % Load parameter and elasticity sets from batch files
        for ibatch = 1:modelCalibrator.nbatch
            
            fprintf('Reading batch %5d of %5d\n', ibatch, modelCalibrator.nbatch)
            
            s = load(fullfile(modelCalibrator.batch_dir, sprintf('batch%05d.mat', ibatch)));
            
            for p = modelCalibrator.paramlist, paramv.(p{1}) = [paramv.(p{1}), s.parambatch.(p{1})]; end
            for e = modelCalibrator.elaslist , elasv.( e{1}) = [elasv.( e{1}), s.elasbatch.( e{1})]; end
            solved = [solved, s.solvedbatch]; %#ok<AGROW>
            
        end
        solved = boolean(solved);
        
        % Construct inverse interpolants that map elasticities to parameters
        for p = modelCalibrator.paramlist, interp.(p{1}) = scatteredInterpolant(elasv.captoout(solved)', elasv.labelas(solved)', elasv.savelas(solved)', paramv.(p{1})(solved)', 'nearest'); end
        
        % Construct elasticity inverter by consolidating inverse interpolants
        function [inverse] = invert_(target)
            captoout = 3;
            for p_ = modelCalibrator.paramlist, inverse.(p_{1}) = interp.(p_{1})(captoout, target.labelas, target.savelas); end
        end
        invert = @invert_; %#ok<NASGU>
        
        % Save elasticity inverter
        save(fullfile(dirFinder.saveroot(), 'invert.mat'), 'invert', 'paramv', 'elasv', 'solved');
        
        % Delete batch directory
        if (exist('clean', 'var') && clean), rmdir(modelCalibrator.batch_dir, 's'), end
        
    end
    
    
    % Plot calibration solution conditions
    function [] = plot_conditions()
        
        % Load vectors of parameters, elasticities, and solution conditions
        s = load(fullfile(dirFinder.saveroot(), 'invert.mat'));
        paramv = s.paramv;
        elasv  = s.elasv ;
        solved = s.solved;
        
        % Initialize figure
        fig = figure; ax = axes;
        
        % Determine colors
        cv = zeros(modelCalibrator.nset, 3);
        devs = min(abs(elasv.captoout(solved)' - 3).^0.5, 1);
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
        
        % Save figure
        savefig(fig, fullfile(dirFinder.saveroot(), 'invert.fig'))
        
    end
    
end

end