function [] = solve_equilibrium()


% Generate grid of interest rates
s_hh = load(fullfile('Parameters','hh_parameters.mat'));
rate_ub = (s_hh.hh_params.discount_factor)^(-1)-0.001;
rate_lb = 1.0055;
n_rates = 20;
rates = linspace(rate_lb,rate_ub,n_rates);

assets_demanded = zeros(size(rates));
assets_supplied = zeros(size(rates));

for ir = 1:n_rates
%     percent_done = ir/n_rates
    assets_demanded(ir) = asset_demand(rates(ir));
    assets_supplied(ir) = asset_supply(rates(ir));

end


plot(rates,assets_demanded,rates,assets_supplied,'LineWidth',4)
legend('Asset Demand','Asset Supply')

eqm_r = fsolve(@(x) excess_demand(x,assets_demanded,assets_supplied,rates),(rate_lb+rate_ub)/2)

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


function [ex_dem] = excess_demand(x,assets_demanded,assets_supplied,rates)

ex_dem = interp1(rates,assets_demanded,x) - interp1(rates,assets_supplied,x);

end
