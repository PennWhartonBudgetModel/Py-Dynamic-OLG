% solve_hh_optimization_test
s = load(fullfile('Parameters','hh_parameters.mat'));

rate = 1.04;

[V, aopt, assets_total, dist] = solve_hh_optimization(s.hh_params,rate);
