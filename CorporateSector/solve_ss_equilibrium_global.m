function [prices, quantities] = solve_ss_equilibrium_global(prices,print_info) %#ok<INUSD>


% Prepare starting points and settings for optimization routine.
s_hh = load(fullfile('Parameters','hh_parameters.mat'));
rate_lb = s_hh.hh_params.discount_factor;
rate_ub = 1;
hh_tolerance = 0.0001;
% optim_options1 = optimset('Display', 'off', 'TolFun', 1e-6, 'TolX', 1e-6); %#ok<NASGU>
% optim_options2 = optimset('TolFun', 1e-6, 'TolX', 1e-6,'MaxFunEvals',10000); %#ok<NASGU>
% optim_options3 = optimoptions(@fmincon,'UseParallel',true);

% Total shares: should add to one.  Only change if shares are explicitly distributed for equity financing.
shares_outstanding = 1;

% Initiating price structure
if ~exist('prices','var')
    prices.fund        = 10;
    prices.wage        = 1;
    prices.rate        = 0.98;
    dividend           = .2;

    % Determine invariant distribution of productivity shocks
    [~, ~, ~, ~, dist] = solve_hh_optimization_mex(s_hh.hh_params, prices, dividend, .1);

    % Determine labor supply when inelastic
    labor_supply = sum(sum(dist).*s_hh.hh_params.prod_shocks');    % Update if labor supply becomes elastically demanded.

    % Solve equilibrium rate of return.

    pw_ub = 2.75;
    pw_lb = 2.25;
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
            [asset_value, output_total, capital_total, labor_demand, dividend, inv_total, adjcost_total] = firm_quantities(prices);
            quantities.firm.ex_div_value_total  = asset_value;
            quantities.firm.output_total        = output_total;
            quantities.firm.capital_total       = capital_total;
            quantities.firm.inv_total           = inv_total;
            quantities.firm.labor_demand        = labor_demand;
            quantities.firm.dividend            = dividend;
            quantities.firm.adjcost_total       = adjcost_total;

            prices.fund = asset_value;

            [~, ~, shares_total, consumption_total, ~] = solve_hh_optimization_mex(s_hh.hh_params, prices, dividend, hh_tolerance);
            quantities.hh.assets_total     = shares_total;
            quantities.hh.consumtion_total = consumption_total;
            quantities.hh.labor_supply     = labor_supply;

            goods_market_error = abs(consumption_total + inv_total + adjcost_total - output_total);
            labor_market_error = abs(labor_demand - labor_supply);
            asset_market_error = abs(shares_total - shares_outstanding);

            print_information

            if (asset_market_error<.0001)||(pr_iter>30), break, end

            if shares_outstanding > shares_total
                pr_ub = prices.rate;
            elseif shares_outstanding < shares_total
                pr_lb = prices.rate;
            end

        end

        if (labor_market_error<.001)||(pw_iter>50), break, end

        if labor_demand>labor_supply
            pw_lb = prices.wage;
        elseif labor_supply>labor_demand
            pw_ub = prices.wage;
        end

    end
        

    
    save(fullfile('Results','results.mat'),'prices','quantities')
            
elseif exist('prices','var')
    
    % Give values at provided prices.
    [asset_value, output_total, capital_total, labor_demand, dividend, inv_total] = firm_quantities(prices);
    quantities.firm.value_total   = asset_value;
    quantities.firm.output_total  = output_total;
    quantities.firm.capital_total = capital_total;
    quantities.firm.inv_total     = inv_total;
    quantities.firm.labor_demand  = labor_demand;
    quantities.firm.dividend      = dividend;
    quantities.firm.adjcost_total = adjcost_total;
    
    prices.fund = asset_value;

    [~, ~, shares_total, consumption_total, ~] = solve_hh_optimization_mex(s_hh.hh_params, prices, dividend, hh_tolerance);
    quantities.hh.assets_total      = shares_total;
    quantities.hh.consumption_total = consumption_total;
    quantities.hh.labor_supply      = labor_supply;

    % Calculate error.
    goods_market_error = abs(consumption_total + inv_total + adjcost_total - output_total);
    asset_market_error = abs(shares_total - shares_outstanding);
    labor_market_error = abs(labor_demand - labor_supply);

    pw_iter = [];
    pr_iter = [];


end










if ~exist('print_info','var')
    return
else
    print_information
end





function [] = print_information()

% Print values to screen.
fprintf('\nInterest Rate = %0.4f\n'         , prices.rate)
fprintf('\nFund Value = %0.2f\n'            , prices.fund)
fprintf('\nShares Demanded = %0.2f\n'       , shares_total)
fprintf('\nWage = %0.4f\n'                  , prices.wage)
fprintf('\nLabor Demand = %0.2f\n'          , labor_demand)
% fprintf('\nLabor Supply = %0.2f\n'          , labor_supply)
fprintf('\nTotal Capital = %0.2f\n'         , capital_total)
fprintf('\nTotal Output = %0.2f\n'          , output_total)
fprintf('\nTotal Adjustment Costs = %0.4f\n', adjcost_total)
fprintf('\nTotal Investment = %0.2f\n'      , inv_total)
fprintf('\nTotal Consumption = %0.2f\n'     , consumption_total)
fprintf('\nGoods Market Error = %0.4f\n'    , goods_market_error)
fprintf('\nLabor Market Error = %0.4f\n'    , labor_market_error)
fprintf('\nAsset Market Error = %0.4f\n'    , asset_market_error)
fprintf('\nLabor Market Iteration = %0.0f\n', pw_iter)
fprintf('\nAsset Market Iteration = %0.0f\n', pr_iter)

end
        


end
