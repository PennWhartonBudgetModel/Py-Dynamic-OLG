
% solve_firm_optimization_test
s = load(fullfile('Parameters','firm_parameters.mat'));

prices.consumption = 1;

[capital_total, eq_total, V_total, dist] = solve_firm_optimization_mex(prices, s.firm_params);
