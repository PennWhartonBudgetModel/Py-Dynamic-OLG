% solve_hh_optimization_test
s = load(fullfile('Parameters','hh_parameters.mat'));

prices.consumption = 1;
prices.rate        = 1.03;
prices.wage        = 1;
tolerance = 0.0001;

[V, aopt, assets_total, consumption_total, dist] = solve_hh_optimization_mex(s.hh_params,prices, tolerance);
