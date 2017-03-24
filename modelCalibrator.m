%%
% Dynamic model calibrator.
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
    npoint = 40;
    
    % Define number of parameter sets per batch
    batchsize = 5;
    
    % Determine number of parameter sets and number of batches
    nset   = modelCalibrator.npoint^modelCalibrator.nparam;
    nbatch = ceil(modelCalibrator.nset / modelCalibrator.batchsize);
    
    % Define batch directory
    batch_dir = fullfile(dirFinder.source, 'Batches');
    
end

methods (Static)
    
    % Initialize batch files
    function [] = initialize_batches()
        
        % Specify parameter lower and upper bounds
        lb.beta = 0.970; lb.gamma = 0.150; lb.sigma =  1.50;
        ub.beta = 1.110; ub.gamma = 0.900; ub.sigma = 22.00;
        
        % Construct vectors of parameter values
        for p = modelCalibrator.paramlist, v.(p{1}) = linspace(lb.(p{1}), ub.(p{1}), modelCalibrator.npoint); end
        
        % Generate parameter sets as unique combinations of parameter values
        [grid.beta, grid.gamma, grid.sigma] = ndgrid(v.beta, v.gamma, v.sigma);
        for i = 1:modelCalibrator.nset, for p = modelCalibrator.paramlist, paramsets(i).(p{1}) = grid.(p{1})(i); end, end %#ok<AGROW>
        
        % Clear or create batch directory
        if exist(modelCalibrator.batch_dir, 'dir'), rmdir(modelCalibrator.batch_dir, 's'), end, mkdir(modelCalibrator.batch_dir)
        
        % Extract and save batches of parameter sets
        for ibatch = 1:modelCalibrator.nbatch
            shift = (ibatch-1)*modelCalibrator.batchsize;
            parambatch = paramsets(shift+1 : min(shift+modelCalibrator.batchsize, end)); %#ok<NASGU>
            save(fullfile(modelCalibrator.batch_dir, sprintf('batch%05d.mat', ibatch)), 'parambatch')
        end
        
    end
    
    
    % Solve dynamic model steady states for all parameter sets in a batch
    function [] = solve_batch(ibatch)
        
        % Identify batch file
        batchfile = fullfile(modelCalibrator.batch_dir, sprintf('batch%05d.mat', ibatch));
        
        % Load batch of parameter sets
        s = load(batchfile);
        parambatch = s.parambatch;
        
        parfor i = 1:length(parambatch)
            
            % Solve steady state
            save_dir = dynamicSolver.steady(parambatch(i));
            
            % Extract elasticity set
            elasbatch(i) = load(fullfile(save_dir, 'elasticities.mat')); %#ok<PFOUS>
            
            % Check if stable solution found
            solved(i) = (elasbatch(i).captoout > 0.5) && (size(csvread(fullfile(save_dir, 'iterations.csv')), 1) < 25); %#ok<NASGU,PFOUS>
            
            % Delete save directory
            rmdir(save_dir, 's')
            
            % Delete parent directories if empty
            [~] = rmdir(fullfile(save_dir, '..'));
            [~] = rmdir(fullfile(save_dir, '..', '..'));
            [~] = rmdir(fullfile(save_dir, '..', '..', '..'));
            
        end
        
        % Save elasticity sets and stable solution checks to batch file
        save(batchfile, '-append', 'elasbatch', 'solved')
        
    end
    
    
    % Construct elasticity inverter
    function [] = construct_inverter(clean)
        
        if ~exist('clean', 'var'), clean = false; end
        
        % Load parameter and elasticity sets from batch files
        for p = modelCalibrator.paramlist, paramsets.(p{1}) = NaN; end
        for e = modelCalibrator.elaslist , elassets.( e{1}) = NaN; end
        
        for ibatch = 1:modelCalibrator.nbatch
            s = load(fullfile(modelCalibrator.batch_dir, sprintf('batch%05d.mat', ibatch)));
            paramsets = [paramsets, s.parambatch(s.solved)]; %#ok<AGROW>
            elassets  = [elassets , s.elasbatch( s.solved)]; %#ok<AGROW>
        end
        paramsets(1) = [];
        elassets (1) = [];
        
        % Construct inverse interpolants that map elasticities to parameters
        for p = modelCalibrator.paramlist, interp.(p{1}) = scatteredInterpolant([elassets.captoout]', [elassets.labelas]', [elassets.savelas]', [paramsets.(p{1})]', 'nearest'); end
        
        % Construct elasticity inverter by consolidating inverse interpolants
        function [inverse] = invert_(target)
            for p_ = modelCalibrator.paramlist, inverse.(p_{1}) = interp.(p_{1})(target.captoout, target.labelas, target.savelas); end
        end
        invert = @invert_; %#ok<NASGU>
        
        % Save elasticity inverter
        save(fullfile(dirFinder.source, 'invert.mat'), 'invert');
        
        % Delete batch directory
        if clean, rmdir(modelCalibrator.batch_dir, 's'), end
        
    end
    
end

end