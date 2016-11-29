
% solve_firm_optimization_test
s = load(fullfile('Parameters','firm_parameters.mat'));

prices.consumption = 1;
prices.rate        = 1.05;
prices.wage        = 1;

s.firm_params.discount_factor = 1/prices.rate;

[capital_total, eq_total, V_total, output_total, dist] = solve_firm_optimization_mex(prices, s.firm_params);
