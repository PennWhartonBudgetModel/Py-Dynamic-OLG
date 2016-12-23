function [] = solve_experiment1()
%  Dividend tax versus dividend income share

%% Find the level of government expenditures for the given taxes.

s_tax = load(fullfile('Parameters','tax_parameters.mat'));
taxes = s_tax.taxes;
[~, s_quantities] = solve_ss_equilibrium_global([], taxes);
benchmark_expenditures = s_quantities.government.government_expenditures;
clear('s_quantities')


% For a fixed policy instrument, find the tax that clears the government budget constraint.

%% Solving optimal income vs dividend tax
div_tax_ub = taxes.hh.dividend;
div_tax_lb = 0;
n_divs = 20;
dividend_taxes = linspace(div_tax_lb,div_tax_ub,n_divs);

maxiter = 40;

dividend_income_shares = zeros(1,n_divs);
welfares = zeros(1,n_divs);

parfor id = 1:n_divs
    
    exp_taxes = taxes;
    exp_taxes.hh.dividend = dividend_taxes(id);
    
    share_ub = 1;
    share_lb = 0;
    iter = 0;
    while true
        iter = iter+1;
        % Determine tax revenue
        exp_taxes.hh.div_inc_share = (share_lb + share_ub)/2;
        [~, quantities] = solve_ss_equilibrium_global([], exp_taxes);
        rev = quantities.government.tax_revenue;
        
        % Determine error and break if complete
        error = abs(rev - benchmark_expenditures);
        if (error<0.001)||(iter>maxiter) , break, end 
        
        % Update criteria
        if rev > benchmark_expenditures
            share_ub = exp_taxes.hh.div_inc_share;
        elseif rev<= benchmark_expenditures
            share_lb = exp_taxes.hh.div_inc_share;
        end
        
    end
    dividend_income_shares(id) = exp_taxes.hh.div_inc_share;
    welfares(id) = quantities.hh.welfare;
    
end

dividend_experiment_results.dividend_taxes         = dividend_taxes;
dividend_experiment_results.dividend_income_shares = dividend_income_shares;
dividend_experiment_results.welfares               = welfares;


% Solve optimal solution (using grid-search)
[~,max_location] = max(welfares);
taxes.hh.dividend = dividend_taxes(max_location);
taxes.hh.div_inc_share = dividend_income_shares(max_location);
[prices, quantities] = solve_ss_equilibrium_global([], taxes);
dividend_experiment_results.prices     = prices;
dividend_experiment_results.quantities = quantities; %#ok<STRNU>

plot(dividend_income_shares,welfares,'LineWidth',4)
    
save(fullfile('Results','experiment_results.mat'),'dividend_experiment_results')      





end