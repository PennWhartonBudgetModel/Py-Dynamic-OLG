function [] = generate_hh_parameters(tauchen, heathcote)


% Determine productivity shock process

if tauchen
    n_prodshocks = 2;
    ar1_coef     = .75;
    var_innov    = .20;
    addpath('..')
    [prod_shocks, prod_transprob, ~] = cooper_mex(n_prodshocks,0,ar1_coef,sqrt(var_innov));
    prod_shocks   = exp(prod_shocks)';
end
if heathcote
    n_prodshocks = 3;
    prod_shocks = [0.167, 0.839, 5.087];
    prod_transprob = zeros(n_prodshocks);
    diagonal_probabilities = [.9, .99, .9];
    prod_transprob(1,1) = diagonal_probabilities(1);
    prod_transprob(1,2) = 1-diagonal_probabilities(1);
    prod_transprob(2,1) = (1-diagonal_probabilities(2))/2;
    prod_transprob(2,2) = diagonal_probabilities(2);
    prod_transprob(2,3) = (1-diagonal_probabilities(2))/2;
    prod_transprob(3,2) = 1-diagonal_probabilities(3);
    prod_transprob(3,3) = diagonal_probabilities(3);
end
    


% Other parameters
crra = 2;
discount_factor = .96;
ns = 50;
shares_grid = [0,logspace(-1,log(10000)/log(10),ns-1)];
% shares_grid = linspace(0,10000,ns);


% Create structure
hh_params.n_prodshocks     = n_prodshocks;
hh_params.prod_transprob   = prod_transprob;
hh_params.prod_shocks      = prod_shocks;
hh_params.crra             = crra;
hh_params.discount_factor  = discount_factor; 
hh_params.ns               = ns; 
hh_params.shares_grid      = shares_grid; %#ok<*STRNU>

save('hh_parameters.mat','hh_params')

end