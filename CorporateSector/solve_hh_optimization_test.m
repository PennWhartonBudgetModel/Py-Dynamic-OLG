% solve_hh_optimization_test
s = load(fullfile('Parameters','hh_parameters.mat'));
s_tax = load(fullfile('Parameters','tax_parameters.mat'));

% Prices
prices.fund        = 10;
prices.wage        = 1;
prices.rate        = 1.03;

% Other inputs
dividend  = .1;
tolerance = 1;

[V, sopt, shares_total, consumption_total, tax_total, V_total, dist] = solve_hh_optimization_mex(s.hh_params, prices, s_tax.taxes,  dividend, tolerance);
