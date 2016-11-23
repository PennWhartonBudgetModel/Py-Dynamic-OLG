function [] = generate_hh_parameters()


% Determine productivity shock process
n_prodshocks = 2;
ar1_coef     = .75;
var_innov    = .20;
addpath('..')
[prod_shocks, prod_transprob, ~] = cooper_mex(n_prodshocks,0,ar1_coef,sqrt(var_innov));
prod_shocks   = exp(prod_shocks);



% Other parameters
crra = 2.5;
discount_factor = .96;
na = 100;
% asset_grid = [0,logspace(-1,log(1000000)/log(10),na-1)];
asset_grid = linspace(0,100000,na);


% Create structure
hh_params.n_prodshocks     = n_prodshocks;
hh_params.prod_transprob   = prod_transprob;
hh_params.prod_shocks      = prod_shocks;
hh_params.crra             = crra;
hh_params.discount_factor  = discount_factor; 
hh_params.na               = na; 
hh_params.asset_grid        = asset_grid; %#ok<*STRNU>

save('hh_parameters.mat','hh_params')

end