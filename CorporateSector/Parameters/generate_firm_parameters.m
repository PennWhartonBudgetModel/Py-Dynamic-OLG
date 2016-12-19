function [] = generate_firm_parameters()


% Determine productivity shock process
n_prodshocks = 5;
ar1_coef     = .9;
var_innov    = .118;
addpath('..')
[prod_shocks, prod_transprob, ~] = cooper_mex(n_prodshocks,0,ar1_coef,sqrt(var_innov));
prod_shocks   = exp(prod_shocks);


% Other parameters
depreciation                = .085; %#ok<*NASGU>
adj_cost_param              = .05;
returns_to_scale_adjustment = .85;    % Want: f(k,n) = [(k^a)*(n^(1-a))]^b = (k^(a*b))*(n^((1-a)*b)).
capital_share               = .33;
capital_share               = capital_share*returns_to_scale_adjustment;
labor_share                 = (1-capital_share)*returns_to_scale_adjustment;
discount_factor             = .95;
nk                          = 15;
% kgrid                       = linspace(.1,2000,nk);
kgrid                       = logspace(log(.1)/log(10),log(2000)/log(10),nk);


% Create structure
firm_params.n_prodshocks    = n_prodshocks;
firm_params.prod_transprob  = prod_transprob;
firm_params.prod_shocks     = prod_shocks;
firm_params.depreciation    = depreciation;
firm_params.adj_cost_param  = adj_cost_param;
firm_params.capital_share   = capital_share; 
firm_params.labor_share     = labor_share; 
firm_params.nk              = nk; 
firm_params.kgrid           = kgrid; 
firm_params.discount_factor = discount_factor; %#ok<*STRNU>

save('firm_parameters.mat','firm_params')

end