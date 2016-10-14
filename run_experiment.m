%%
% Run a singleton policy experiment.
% 
% Arguments:
% 
%   deep_params
%       beta, gamma, and sigma preference parameters collected in a 1 x 3 vector [beta, gamma, sigma].
% 
%   plan
%       Tax plan string identifier {'base', 'clinton', 'trump', 'ryan'}.
% 
%   gcut
%       Percentage government expenditure reduction, with positive values corresponding to reductions.
%       Specify as [] for open economy.
% 
%   exp_param_set
%       Parameters to adjust for experiment specified as a cell array of strings -- e.g. {'avg_deduc', 'coefs'}.
% 
%%


function [] = run_experiment(deep_params, plan, gcut, exp_param_set)

% Set default deep parameters if none provided
if ~exist('deep_params', 'var') || isempty(deep_params)
    beta  = 1.005;
    gamma = 0.390;
    sigma = 06.00;
    deep_params = [beta, gamma, sigma];
end

% Set base plan as default
if ~exist('plan', 'var') || isempty(plan)
    plan = 'base';
end

% Identify open economy by absence of gcut
isopen = ~exist('gcut', 'var') || isempty(gcut);

% Set null experiment as default
if ~exist('exp_param_set', 'var') || isempty(exp_param_set)
    exp_param_set = {};
end


% Identify parameter directory
param_dir = dirFinder.param;


% Define experiment plan
% (Will be equivalent to plan if no experiment parameters specified)
exp_plan = sprintf('%s%s%s', plan, repmat('-', ~isempty(exp_param_set)), strjoin(exp_param_set, ','));


if ~isempty(exp_param_set)
    
    % Load baseline tax parameters
    s_inctax = load(fullfile(param_dir, 'param_inctax_base.mat'));
    s_bustax = load(fullfile(param_dir, 'param_bustax_base.mat'));

    % Load experiment tax parameters to override baseline parameters
    s_inctax_plan = load(fullfile(param_dir, sprintf('param_inctax_%s.mat', plan)));
    s_bustax_plan = load(fullfile(param_dir, sprintf('param_bustax_%s.mat', plan)));

    for i = 1:length(exp_param_set)
        param = exp_param_set{i};
        if     isfield(s_inctax_plan, param), s_inctax.(param) = s_inctax_plan.(param);
        elseif isfield(s_bustax_plan, param), s_bustax.(param) = s_bustax_plan.(param);
        else   error('''%s'' not found amongst tax parameters.', param)
        end
    end
    
    % Create experiment parameter files
    exp_inctax_file = fullfile(param_dir, sprintf('param_inctax_%s.mat', exp_plan));
    exp_bustax_file = fullfile(param_dir, sprintf('param_bustax_%s.mat', exp_plan));
    
    save(exp_inctax_file, '-struct', 's_inctax')
    save(exp_bustax_file, '-struct', 's_bustax')
    
end


% Call solver
if isopen
    solve_open(  deep_params, exp_plan, false)
else
    solve_closed(deep_params, exp_plan, gcut, false)
end


% Clean up experiment parameter files
if ~isempty(exp_param_set)
    delete(exp_inctax_file)
    delete(exp_bustax_file)
end


end