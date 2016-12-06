function [prices, quantities] = solve_ss_equilibrium_global()


% Prepare starting points and settings for optimization routine.
s_hh = load(fullfile('Parameters','hh_parameters.mat'));
rate_ub = (s_hh.hh_params.discount_factor)^(-1);
rate_lb = 1;
% optim_options1 = optimset('Display', 'off', 'TolFun', 1e-6, 'TolX', 1e-6); %#ok<NASGU>
% optim_options2 = optimset('TolFun', 1e-6, 'TolX', 1e-6,'MaxFunEvals',10000); %#ok<NASGU>
% optim_options3 = optimoptions(@fmincon,'UseParallel',true);


% Initiating price structure
prices.consumption = 1;
prices.rate        = 1.02;
prices.wage        = 1;

% Determine invariant distribution of productivity shocks
[~, ~, ~, ~, dist] = solve_hh_optimization_mex(s_hh.hh_params,prices);

% Solve equilibrim wage
labor_supply = sum(sum(dist).*s_hh.hh_params.prod_shocks');    % Update if labor supply becomes elastically demanded.

% Solve equilibrium rate of return.
if ~exist('eqm_prices','var')
    
    % Solve equilibrium price of consumption.
    pc_ub = 2.5;
    pc_lb = 2.25;
    pc_iter = 0;
    while true
        pc_iter = pc_iter + 1;
        prices.consumption = (pc_ub + pc_lb)/2;
        
        pw_ub = 8;
        pw_lb = 6;
        pw_iter = 0;
        while true
            pw_iter = pw_iter + 1;
            prices.wage = (pw_ub + pw_lb)/2;
            
            pr_ub = rate_ub;
            pr_lb = rate_lb;
            pr_iter = 0;
            while true
                pr_iter = pr_iter + 1;
                prices.rate = (pr_ub + pr_lb)/2;
                
                % Solve agent problems
                [asset_demand, output_total, labor_demand] = firm_quantities(prices);
                quantities.firm.value_total  = asset_demand;
                quantities.firm.output_total = output_total;
                
                [~, ~, asset_supply, consumption_total, ~] = solve_hh_optimization_mex(s_hh.hh_params,prices);
                quantities.hh.assets_total     = asset_supply;
                quantities.hh.consumtion_total = consumption_total;
                
                goods_market_error = abs(consumption_total - output_total);
                labor_market_error = abs(labor_demand - labor_supply);
                asset_market_error = abs(asset_demand-asset_supply);
                
                print_information
                
                if (asset_market_error<.001)||(pr_iter>50), break, end
                
                if asset_demand>asset_supply
                    pr_lb = prices.rate;
                elseif asset_demand<asset_supply
                    pr_ub = prices.rate;
                end
                
            end
            
            if (labor_market_error<.01)||(pw_iter>50), break, end
    
            if labor_demand>labor_supply
                pw_lb = prices.wage;
            elseif labor_supply>labor_demand
                pw_ub = prices.wage;
            end
            
        end
        
        if (goods_market_error<.01)||(pc_iter>20), break, end

        if consumption_total>output_total
            pc_lb = prices.consumption;
        elseif output_total>consumption_total
            pc_ub = prices.consumption;
        end
        
    end
            
            

        
    
end









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
fprintf('\nInterest Rate = %0.4f\n'         ,prices.rate)
fprintf('\nAsset Demand = %0.2f\n'          ,asset_demand)
fprintf('\nAsset Supply = %0.2f\n'          ,asset_supply)
fprintf('\nWage = %0.4f\n'                  ,prices.wage)
fprintf('\nLabor Demand = %0.2f\n'          ,labor_demand)
fprintf('\nLabor Supply = %0.2f\n'          ,labor_supply)
fprintf('\nConsumption Price = %0.4f\n'     ,prices.consumption)
fprintf('\nTotal Output = %0.2f\n'          ,output_total)
fprintf('\nTotal Consumption = %0.2f\n'     ,consumption_total)
fprintf('\nGoods Market Error = %0.4f\n'    ,goods_market_error)
fprintf('\nLabor Market Error = %0.4f\n'    ,labor_market_error)
fprintf('\nAsset Market Error = %0.4f\n'    ,asset_market_error)
fprintf('\nGood Market Iteration = %0.0f\n' ,pc_iter)
fprintf('\nLabor Market Iteration = %0.0f\n',pw_iter)
fprintf('\nAsset Market Iteration = %0.0f\n',pr_iter)

end
        


end
