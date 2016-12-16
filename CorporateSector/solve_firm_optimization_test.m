
% solve_firm_optimization_test
s = load(fullfile('Parameters','firm_parameters.mat'));

prices.fund        = 0.98;
prices.wage        = 1;
prices.rate        = .98;

s.firm_params.discount_factor = prices.rate;

[capital_total, labor_total, eq_total, inv_total, adjcost_total, V_total, output_total, dist] = solve_firm_optimization(prices, s.firm_params);
