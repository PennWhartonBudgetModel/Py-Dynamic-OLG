function [] = solve_equilibrium_()


% Generate grid of interest rates
rate_lb = 1.02;
rate_ub = 1.0415;
n_rates = 5;
rates = linspace(rate_lb,rate_ub,n_rates);

% Test bounds of algorithm
excess_demand_lb = excess_demand(rate_lb, rate_lb, rate_ub);    % Should be positive
excess_demand_ub = excess_demand(rate_ub, rate_lb, rate_ub);    % Should be negative
if excess_demand_lb<0 
    error('Excess demand negative at lower bound.  Decrease lower bound and set upper bound to lower bound value.')
elseif excess_demand_ub>0
    error('Excess demand positive at lower bound.  Increase upper bound and set lower bound to upper bound value.')
end

% Solve equilibrium rate of return
eqm_r = fsolve(@(x) excess_demand(x, rate_lb, rate_ub),(rate_lb+rate_ub)/2)

[~,eqm_discount] = asset_demand(eqm_r);
eqm_discount

beta_times_r = eqm_r*eqm_discount


% Solve household optimization at equilibrium solution
s_hh = load(fullfile('Parameters','hh_parameters.mat'));
[V, aopt, assets_total, consumption_total, dist] = solve_hh_optimization_mex(s_hh.hh_params,eqm_r);

% Solve firm optimization at equilibrium solution
s_firm = load(fullfile('Parameters','firm_parameters.mat'));
s_firm.firm_params.discount_factor = eqm_discount;
prices.consumption = 1;
[capital_total, eq_total, V_total, output_total, dist] = solve_firm_optimization_mex( prices, s_firm.firm_params); %#ok<*ASGLU>

if s_hh.hh_params.n_prodshocks==1
    [assets_total, consumption_total] = asset_supply(eqm_r);
    aggregate_consumption = consumption_total %#ok<*NOPRT,*NASGU>
    aggregate_savings     = assets_total
elseif s_hh.hh_params.n_prodshocks~=1
    aggregate_consumption = consumption_total
    aggregate_savings     = assets_total
end
mutual_fund_value     = V_total
aggregate_capital     = capital_total
aggregate_output      = output_total

end


function [ex_dem] = excess_demand(x, rate_lb, rate_ub)

if (x<rate_lb)||(x>rate_ub)
    ex_dem = Inf;
    return
end

ex_dem = asset_demand(x) - asset_supply(x);

end
