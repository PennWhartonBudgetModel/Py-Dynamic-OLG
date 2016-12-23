function [prices, quantities] = solve_ss_equilibrium_global(prices, taxes, print_info) %#ok<INUSD>


%% Prepare starting points and settings for optimization routine.
s_hh = load(fullfile('Parameters','hh_parameters.mat'));

if ~exist('taxes','var')||isempty(taxes)
    s_tax = load(fullfile('Parameters','tax_parameters.mat'));
    taxes = s_tax.taxes;
end

rate_lb = s_hh.hh_params.discount_factor;
rate_ub = 1;
hh_tolerance = 0.0001;

% Total shares: should add to one.  Only change if shares are explicitly distributed for equity financing.
shares_outstanding = 1;

%% Solve equilibrium
if ~exist('prices','var')||isempty(prices)
    prices.fund        = 10;
    prices.wage        = 1;
    prices.rate        = 0.98;
    dividend           = .2;

    % Determine invariant distribution of productivity shocks
    if s_hh.hh_params.n_prodshocks>1
        [~, ~, ~, ~, ~, ~, ~, ~, dist] = solve_hh_optimization_mex(s_hh.hh_params, prices, taxes, dividend, .1);
    elseif s_hh.hh_params.n_prodshocks==1
        dist = 1;
    end

    % Determine labor supply when inelastic
    labor_supply = sum(sum(dist).*s_hh.hh_params.prod_shocks);    % Update if labor supply becomes elastically demanded.

    % Solve equilibrium rate of return.

    pw_ub = 2.25;
    pw_lb = 2.075;
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
            [asset_value, output_total, capital_total, labor_demand, dividend, inv_total, adjcost_total, firm_tax_total] = firm_quantities(prices, taxes);
            quantities.firm.ex_div_value_total  = asset_value;
            quantities.firm.output_total        = output_total;
            quantities.firm.capital_total       = capital_total;
            quantities.firm.inv_total           = inv_total;
            quantities.firm.labor_demand        = labor_demand;
            quantities.firm.dividend            = dividend;
            quantities.firm.adjcost_total       = adjcost_total;
            quantities.firm.tax_total           = firm_tax_total;

            prices.fund = asset_value;
            if s_hh.hh_params.n_prodshocks>1
                [~, ~, ~, ~, shares_total, consumption_total, hh_tax_total, V_total, ~]    = solve_hh_optimization_mex(s_hh.hh_params, prices, taxes, dividend, hh_tolerance);
            elseif s_hh.hh_params.n_prodshocks==1
                [V, sopt, copt, taxopt, ~, ~, ~, ~, ~] = solve_hh_optimization_mex(s_hh.hh_params, prices, taxes, dividend, hh_tolerance);
                shares_total = fsolve(@(x) find_shares_fixedpoint(x,sopt',s_hh.hh_params.shares_grid),1);
                V_total           = interp1(s_hh.hh_params.shares_grid, V     , shares_total, 'linear', 'extrap');
                consumption_total = interp1(s_hh.hh_params.shares_grid, copt  , shares_total, 'linear', 'extrap');
                hh_tax_total      = interp1(s_hh.hh_params.shares_grid, taxopt, shares_total, 'linear', 'extrap');
            end
            
            quantities.hh.assets_total     = shares_total;
            quantities.hh.consumtion_total = consumption_total;
            quantities.hh.labor_supply     = labor_supply;
            quantities.hh.tax_total        = hh_tax_total;
            quantities.hh.welfare          = V_total;

            tax_revenue = firm_tax_total + hh_tax_total;
            government_expenditures = tax_revenue;
            welfare     = V_total;
            
            quantities.government.tax_revenue = tax_revenue;
            quantities.government.government_expenditures = government_expenditures;
            
            goods_market_error = abs(consumption_total + inv_total + adjcost_total + government_expenditures - output_total);
            labor_market_error = abs(labor_demand - labor_supply);
            asset_market_error = abs(shares_total - shares_outstanding);
            
            quantities.equilibrium.goods_market_error = goods_market_error;
            quantities.equilibrium.labor_market_error = labor_market_error;
            quantities.equilibrium.asset_market_error = asset_market_error;
            
            if exist('print_info','var')
                print_information
            end

            if (asset_market_error<.00001)||(pr_iter>30), break, end

            if shares_outstanding > shares_total
                pr_ub = prices.rate;
            elseif shares_outstanding < shares_total
                pr_lb = prices.rate;
            end

        end

        if (labor_market_error<.0001)||(pw_iter>50), break, end

        if labor_demand>labor_supply
            pw_lb = prices.wage;
        elseif labor_supply>labor_demand
            pw_ub = prices.wage;
        end

    end
        

    
    save(fullfile('Results','results.mat'),'prices','quantities')
            
elseif exist('prices','var')&&~isempty('prices')
    labor_supply = [];
    
    % Give values at provided prices.
    [asset_value, output_total, capital_total, labor_demand, dividend, inv_total, adjcost_total, firm_tax_total] = firm_quantities(prices, taxes);
    quantities.firm.value_total   = asset_value;
    quantities.firm.output_total  = output_total;
    quantities.firm.capital_total = capital_total;
    quantities.firm.inv_total     = inv_total;
    quantities.firm.labor_demand  = labor_demand;
    quantities.firm.dividend      = dividend;
    quantities.firm.adjcost_total = adjcost_total;
    quantities.firm.tax_total     = firm_tax_total;
    
    prices.fund = asset_value;

    [~, ~, shares_total, consumption_total, hh_tax_total, V_total, ~] = solve_hh_optimization_mex(s_hh.hh_params, prices, taxes, dividend, hh_tolerance);
    quantities.hh.assets_total      = shares_total;
    quantities.hh.consumption_total = consumption_total;
    quantities.hh.labor_supply      = labor_supply;
    quantities.hh.tax_total         = hh_tax_total;
    quantities.hh.welfare           = V_total;
    
    tax_revenue = firm_tax_total + hh_tax_total;
    government_expenditures = tax_revenue;
    welfare     = V_total;

    % Calculate error.
    goods_market_error = abs(consumption_total + inv_total + adjcost_total + government_expenditures - output_total);
    asset_market_error = abs(shares_total - shares_outstanding);
    labor_market_error = abs(labor_demand - labor_supply);

    pw_iter = [];
    pr_iter = [];


end



function error = find_shares_fixedpoint(x,sopt,shares_grid)
    error = interp1(shares_grid, sopt - shares_grid, x, 'linear', 'extrap');
end




function [] = print_information()

% Print values to screen.
fprintf('\nInterest Rate = %0.4f\n'         , prices.rate)
fprintf('\nFund Value = %0.2f\n'            , prices.fund)
fprintf('\nShares Demanded = %0.2f\n'       , shares_total)
fprintf('\nWage = %0.4f\n'                  , prices.wage)
fprintf('\nLabor Demand = %0.2f\n'          , labor_demand)
fprintf('\nTotal Capital = %0.2f\n'         , capital_total)
fprintf('\nTotal Output = %0.2f\n'          , output_total)
fprintf('\nTotal Adjustment Costs = %0.4f\n', adjcost_total)
fprintf('\nTotal Investment = %0.2f\n'      , inv_total)
fprintf('\nTotal Consumption = %0.2f\n'     , consumption_total)
fprintf('\nGoods Market Error = %0.4f\n'    , goods_market_error)
fprintf('\nLabor Market Error = %0.4f\n'    , labor_market_error)
fprintf('\nAsset Market Error = %0.4f\n'    , asset_market_error)
fprintf('\nTax Revenue = %0.2f\n'           , tax_revenue)
fprintf('\nTotal Welfare = %0.2f\n'         , welfare)
fprintf('\nLabor Market Iteration = %0.0f\n', pw_iter)
fprintf('\nAsset Market Iteration = %0.0f\n', pr_iter)

end
        


end
