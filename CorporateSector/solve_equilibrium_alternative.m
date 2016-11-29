function [prices, quantities] = solve_equilibrium_alternative()


% Generate grid of interest rates
s_hh = load(fullfile('Parameters','hh_parameters.mat'));
rate_ub = (s_hh.hh_params.discount_factor)^(-1);
rate_lb = 1;
optim_options = optimset('Display', 'off', 'TolFun', 1e-6, 'TolX', 1e-6);


% Solve equilibrium price of consumption
pc_ub = .75;
pc_lb = .65;
while true
    prices.consumption = (pc_ub + pc_lb)/2;
    prices.rate = 1.02;    % This value will be replaced in the solver.
    prices.wage        = 1;
    % Solve equilibrium rate of return
    eqm_prices = fsolve( @(x) excess_demand(x, s_hh.hh_params,prices) , (rate_lb+rate_ub)/2, optim_options);
    prices.rate        = eqm_prices(1);
    
    [asset_demand, output_total] = firm_quantities(prices);
    quantities.firm.value_total  = asset_demand;
    quantities.firm.output_total = output_total;
    
    [~, ~, asset_supply, consumption_total, ~] = solve_hh_optimization_mex(s_hh.hh_params,prices);
    quantities.hh.assets_total     = asset_supply;
    quantities.hh.consumtion_total = consumption_total;
    
    fprintf('\nInterest Rate = %0.4f\n',prices.rate)
    fprintf('\nAsset Demand = %0.2f\n',asset_demand)
    fprintf('\nAsset Supply = %0.2f\n',asset_supply)
    fprintf('\nConsumption Price = %0.4f\n',prices.consumption)
    fprintf('\nTotal Output = %0.2f\n',output_total)
    fprintf('\nTotal Consumption = %0.2f\n',consumption_total)
    
    error = abs(consumption_total - output_total)
    
    if error<.1, break, end
    
    if consumption_total>output_total
        pc_lb = prices.consumption;
    elseif output_total>consumption_total
        pc_ub = prices.consumption;
    end
    
end
    


% if s_hh.hh_params.n_prodshocks==1
%     [assets_total, consumption_total] = asset_supply(eqm_r);
%     aggregate_consumption = consumption_total 
%     aggregate_savings     = assets_total
% elseif s_hh.hh_params.n_prodshocks~=1
%     aggregate_consumption = consumption_total
%     aggregate_savings     = assets_total
% end
% mutual_fund_value     = V_total
% aggregate_capital     = capital_total
% aggregate_output      = output_total

end


function [ex_dem] = excess_demand(x, hh_params, prices)%, rate_lb, rate_ub)


% if (x(2)<rate_lb)||(x(2)>rate_ub)
%     ex_dem = Inf;
%     return
% end

prices.rate        = x;

[asset_demand, output_total] = firm_quantities(prices); %#ok<ASGLU>
[~, ~, asset_supply, consumption_total, ~] = solve_hh_optimization_mex(hh_params,prices); %#ok<ASGLU>


% ex_dem(1) = consumption_total - output_total;
ex_dem = asset_demand - asset_supply;

end
