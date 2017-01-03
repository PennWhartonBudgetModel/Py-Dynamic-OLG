function [] = solve_experiment4()
% This experiment determines whether to tax corporate income with partial expensing or non-linear dividend taxation.

% Generates benchmark tax parameter values.
taxes = paramGenerator.tax;
[~, s_quantities] = solve_ss_equilibrium_global([], taxes);
benchmark_expenditures = s_quantities.government.government_expenditures;
benchmark_welfare = s_quantities.hh.welfare;
clear('s_quantities')


%% Solving optimal income vs corporate income tax
cit_ub = taxes.firm.income;
cit_lb = 0;
n_cits = 20;
cits = linspace(cit_lb, cit_ub, n_cits);

max_iter = 40;

dividend_income_shares = zeros(1,n_cits);
welfares = zeros(1, n_cits);

parfor ic = 1:n_cits
    
    experiment_taxes = taxes;
    experiment_taxes.firm.income = cits(ic);
    
    share_ub = 1;
    share_lb = 0;
    
    iter = 0;
    while true
        iter = iter + 1;
        % Determine tax revenue
        experiment_taxes.hh.div_inc_share = (share_lb + share_ub)/2;
        [~, quantities] = solve_ss_equilibrium_global([], experiment_taxes); 
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
    
    dividend_income_shares(ic) = exp_taxes.hh.div_inc_share;
    welfares(ic) = quantities.hh.welfare;
    
end

plot(dividend_income_shares,welfares,'LineWidth',4)

experiment4.cits                  = cits;
experiment4.dividend_income_share = dividend_income_share;
experiment4.welfares              = welfares; %#ok<STRNU>

save(fullfile('Results','experiment4_results.mat'),'experiment4')
    





end
