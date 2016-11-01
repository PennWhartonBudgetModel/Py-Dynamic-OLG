function [] = generate_firm_parameters()


% Determine productivity shock process
n_prodshocks = 5;
ar1_coef     = .85;
var_innov    = .1;
addpath('..')
[prod_shocks, prod_transprob, ~] = cooper_mex(n_prodshocks,0,ar1_coef,sqrt(var_innov));
prod_shocks   = exp(prod_shocks);


% Other parameters
depreciation    = .08; %#ok<*NASGU>
adj_cost_param  = .15;
prod_func_param = .7;
nk              = 10;
kgrid           = linspace(.1,50,nk);


% Create structure
firm_params.n_prodshocks    = n_prodshocks;
firm_params.prod_transprob  = prod_transprob;
firm_params.prod_shocks     = prod_shocks;
firm_params.depreciation    = depreciation;
firm_params.adj_cost_param  = adj_cost_param;
firm_params.prod_func_param = prod_func_param; 
firm_params.nk              = nk; 
firm_params.kgrid           = kgrid; %#ok<*STRNU>

save('firm_parameters.mat','firm_params')

end