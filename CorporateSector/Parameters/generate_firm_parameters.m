function [] = generate_firm_parameters()


% Determine productivity shock process
n_prodshocks = 5;
ar1_coef     = .85;
var_innov    = .1;
addpath('..')
[prod_shocks, prod_transprob, ~] = cooper_mex(n_prodshocks,0,ar1_coef,sqrt(var_innov));
prod_shocks   = exp(prod_shocks);


% Other parameters
depreciation    = .03; %#ok<*NASGU>
adj_cost_param  = .15;
prod_func_param = .8;
discount_factor = .95;
nk              = 15;
kgrid           = linspace(.1,1500,nk);


% Create structure
firm_params.n_prodshocks    = n_prodshocks;
firm_params.prod_transprob  = prod_transprob;
firm_params.prod_shocks     = prod_shocks;
firm_params.depreciation    = depreciation;
firm_params.adj_cost_param  = adj_cost_param;
firm_params.prod_func_param = prod_func_param; 
firm_params.nk              = nk; 
firm_params.kgrid           = kgrid; 
firm_params.discount_factor = discount_factor; %#ok<*STRNU>

save('firm_parameters.mat','firm_params')

end