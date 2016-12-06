function [prices, quantities] = solve_ss_equilibrium(print_info, eqm_prices, global_opt) %#ok<INUSD,INUSL>


% Prepare starting points and settings for optimization routine.
s_hh = load(fullfile('Parameters','hh_parameters.mat'));
rate_ub = (s_hh.hh_params.discount_factor)^(-1);
rate_lb = 1;
optim_options1 = optimset('Display', 'off', 'TolFun', 1e-6, 'TolX', 1e-6); %#ok<NASGU>
optim_options2 = optimset('TolFun', 1e-6, 'TolX', 1e-6,'MaxFunEvals',10000); %#ok<NASGU>
optim_options3 = optimoptions(@fmincon,'UseParallel',true);


% Initiating price structure and providing first guess
prices.consumption = .5;
prices.rate        = (rate_lb+rate_ub)/2;
prices.wage        = 1;


% Solve equilibrim wage
labor_supply = 1;    % Update if labor supply becomes elastically demanded.

% Solve equilibrium rate of return.
if ~exist('eqm_prices','var')
%     if ~exist('global_opt','var')
%         eqm_prices = fminsearch( @(x) excess_demand(x, s_hh.hh_params,prices, labor_supply) , [prices.rate, prices.wage], optim_options2);
        eqm_prices = fmincon( @(x) excess_demand(x, s_hh.hh_params,prices, labor_supply) , [prices.rate, prices.wage],...
            [],[],[],[],[rate_lb,.01],[rate_ub,10],[],optim_options3);
% %         eqm_prices = simulannealbnd( @(x) excess_demand(x, s_hh.hh_params,prices, labor_supply) , [prices.rate, prices.wage],[rate_lb,.1],[rate_ub,10], optim_options2);
%     else
%         % Employ multi-start global optimization routine.
%         fun_to_minimize = @(x) excess_demand(x, s_hh.hh_params,prices, labor_supply);
%         lb1 = [];
%         ub1 = [];
%         initial_guess = [prices.rate, prices.wage];
%         
        
        
    
end

% prices.consumption = eqm_prices(1);
prices.rate        = eqm_prices(1);
prices.wage        = eqm_prices(2);

% Give values at optimal solution.
[asset_demand, output_total, labor_demand] = firm_quantities(prices);
quantities.firm.value_total  = asset_demand;
quantities.firm.output_total = output_total;

[~, ~, asset_supply, consumption_total, ~] = solve_hh_optimization_mex(s_hh.hh_params,prices);
quantities.hh.assets_total     = asset_supply;
quantities.hh.consumtion_total = consumption_total;

% Calculate error.
goods_market_error = abs(consumption_total - output_total);
asset_market_error = abs(asset_demand - asset_supply);
labor_market_error = abs(labor_demand - labor_supply);

if ~exist('print_info','var')
    return
else
    print_information
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





function [ex_dem] = excess_demand(x, hh_params, prices, labor_supply)

% prices.consumption = x(1);
prices.rate        = x(1);
prices.wage        = x(2);

[asset_demand, output_total, labor_demand] = firm_quantities(prices);
[~, ~, asset_supply, consumption_total, ~] = solve_hh_optimization_mex(hh_params,prices);


% ex_dem(1) = consumption_total - output_total;
ex_dem(1) = asset_demand - asset_supply;
ex_dem(2) = labor_demand - labor_supply;

ex_dem = norm(ex_dem,2);    % implemented for fminsearch ONLY

end


function [] = print_information()

% Print values to screen.
fprintf('\nInterest Rate = %0.4f\n'     ,prices.rate)
fprintf('\nAsset Demand = %0.2f\n'      ,asset_demand)
fprintf('\nAsset Supply = %0.2f\n'      ,asset_supply)
fprintf('\nConsumption Price = %0.4f\n' ,prices.consumption)
fprintf('\nTotal Output = %0.2f\n'      ,output_total)
fprintf('\nTotal Consumption = %0.2f\n' ,consumption_total)
fprintf('\nWage = %0.4f\n'              ,prices.wage)
fprintf('\nLabor Demand = %0.2f\n'      ,labor_demand)
fprintf('\nLabor Supply = %0.2f\n'      ,labor_supply)
fprintf('\nGoods Market Error = %0.4f\n',goods_market_error)
fprintf('\nAsset Market Error = %0.4f\n',asset_market_error)
fprintf('\nLabor Market Error = %0.4f\n',labor_market_error)

end
        


end
