% solve_hh_optimization_test
s = load(fullfile('Parameters','hh_parameters.mat'));

% Prices
prices.fund        = 10;
prices.wage        = 1;
prices.rate        = 1.03;

% Other inputs
dividend  = .1;
tolerance = 1;

[V, sopt, shares_total, consumption_total, dist] = solve_hh_optimization(s.hh_params, prices, dividend, tolerance);
