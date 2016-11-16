function [] = solve_equilibrium()


% Generate grid of interest rates
rate_lb = 1.0055;
rate_ub = 1.05;
n_rates = 10;
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

eqm_r = fsolve(@(x) excess_demand(x,assets_demanded,assets_supplied,rates),1.02)

[~,eqm_discount] = asset_demand(eqm_r);
eqm_discount

beta_times_r = eqm_r*eqm_discount

end


function [ex_dem] = excess_demand(x,assets_demanded,assets_supplied,rates)

ex_dem = interp1(rates,assets_demanded,x) - interp1(rates,assets_supplied,x);

end
