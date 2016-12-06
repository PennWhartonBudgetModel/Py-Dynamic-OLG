function [prices, quantities] = solve_equilibrium_alternative()


% Prepare starting points and settings for optimization routine.
s_hh = load(fullfile('Parameters','hh_parameters.mat'));
rate_ub = (s_hh.hh_params.discount_factor)^(-1);
rate_lb = 1;
optim_options = optimset('Display', 'off', 'TolFun', 1e-6, 'TolX', 1e-6);


% Initiating price structure
prices.consumption = 1;
prices.rate        = 1.02;
prices.wage        = 1;


% Solve equilibrim wage
w_ub = 10;
w_lb = 0;
labor_supply = 1;    % Update if labor supply becomes elastically demanded.
wage_iter = 0;
while true
    % Update wage.
    prices.wage = (w_ub+w_lb)/2;

    % Solve equilibrium price of consumption.
    pc_ub = 6;
    pc_lb = 0;
    wage_iter = wage_iter+1;
    goods_iter = 0;
    while true
        goods_iter = goods_iter+1;
        % Update price.
        prices.consumption = (pc_ub + pc_lb)/2;
        
        % Solve equilibrium rate of return.
        eqm_prices = fsolve( @(x) excess_demand(x, s_hh.hh_params,prices) , (rate_lb+rate_ub)/2, optim_options);
        prices.rate        = eqm_prices(1);
        
        % Give values at optimal solution.
        [asset_demand, output_total, labor_demand] = firm_quantities(prices);
        quantities.firm.value_total  = asset_demand;
        quantities.firm.output_total = output_total;

        [~, ~, asset_supply, consumption_total, ~] = solve_hh_optimization_mex(s_hh.hh_params,prices);
        quantities.hh.assets_total     = asset_supply;
        quantities.hh.consumtion_total = consumption_total;
        
        % Calculate error.
        goods_market_error = abs(consumption_total - output_total);
        labor_market_error = abs(labor_demand - labor_supply);
        
        % Print values to screen.
        fprintf('\nInterest Rate = %0.4f\n',prices.rate)
        fprintf('\nAsset Demand = %0.2f\n',asset_demand)
        fprintf('\nAsset Supply = %0.2f\n',asset_supply)
        fprintf('\nConsumption Price = %0.4f\n',prices.consumption)
        fprintf('\nTotal Output = %0.2f\n',output_total)
        fprintf('\nTotal Consumption = %0.2f\n',consumption_total)
        fprintf('\nLabor Demand = %0.2f\n',labor_demand)
        fprintf('\nLabor Supply = %0.2f\n',labor_supply)
        fprintf('\nGoods Market Error = %0.4f\n',goods_market_error)
        fprintf('\nLabor Market Error = %0.4f\n',labor_market_error)
        fprintf('\nWage Iteration = %0.0f\n',wage_iter)
        fprintf('\nGoods Iteration = %0.0f\n\n',goods_iter)
        

        if goods_market_error<.01, break, end

        if consumption_total>output_total
            pc_lb = prices.consumption;
        elseif output_total>consumption_total
            pc_ub = prices.consumption;
        end

    end
    
    
    if labor_market_error<.01, break, end
    
    if labor_demand>labor_supply
        w_lb = prices.wage;
    elseif labor_supply>labor_demand
        w_ub = prices.wage;
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

[asset_demand, output_total, labor_demand] = firm_quantities(prices);
[~, ~, asset_supply, consumption_total, ~] = solve_hh_optimization_mex(hh_params,prices);


% ex_dem(1) = consumption_total - output_total;
ex_dem = asset_demand - asset_supply;

end





end
