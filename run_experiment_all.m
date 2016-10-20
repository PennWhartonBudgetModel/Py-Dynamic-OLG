%%
% Run all singleton policy experiments.
% 
%%


function [] = run_experiment_all()

% Specify experiment parameter sets
exp_param_sets = { {'limit', 'X'}         ;
                   {'cap_tax_share'}      ; 
                   {'avg_deduc', 'coefs'} ;
                   {'tau_cap'}            ;
                   {'exp_share'}          };
n_sets = length(exp_param_sets);

% Specify plans
plans = {'trump', 'clinton', 'ryan'};
n_plans = length(plans);

% Specify deep parameter set
deep_params = inddeep_to_params(6);

% Specify government expenditure reduction
gcut = +0.05;


parfor i = 1:n_sets*n_plans

    i_set  = mod(i-1, n_sets)+1;
    i_plan = ceil(i/n_sets);

	% Run open economy experiment
    run_experiment(deep_params, plans{i_plan}, [],   exp_param_sets{i_set}) %#ok<PFBNS>
		
	% Run closed economy experiment
    run_experiment(deep_params, plans{i_plan}, gcut, exp_param_sets{i_set})

end


end