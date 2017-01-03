function [] = solve_experiment3()
% This experiment solves the optimal dividend and corporate income tax when expensing is different from 1.



%% From the zero-tax structure, set CIT and expense share, solve equilibrium.
taxes = paramGenerator.tax(true);
taxes.firm.income    = .15;
taxes.firm.exp_share = .58;

[~, s_quantities] = solve_ss_equilibrium_global([], taxes);
benchmark_expenditures = s_quantities.government.government_expenditures;
benchmark_welfare = s_quantities.hh.welfare;
clear('s_quantities')




%% Find the dividend tax that clears market
taxes.firm.income = 0;

dub = .25;
dlb = .145;
tolerance = .00001;
while true
    taxes.hh.dividend = (dub+dlb)/2;
    [prices, quantities] = solve_ss_equilibrium_global([], taxes); %#ok<ASGLU>
    rev = quantities.government.tax_revenue;
    error = abs(rev - benchmark_expenditures);
    if error<tolerance, break, end
    fprintf('\nDividend tax upper bound = %0.6f\n', dub)
    fprintf('\nDividend tax lower bound = %0.6f\n', dlb)
    fprintf('\nError = %0.6f\n', error)
    if rev>benchmark_expenditures
        dub = taxes.hh.dividend;
    elseif rev<=benchmark_expenditures
        dlb = taxes.hh.dividend;
    end
end

counterfactual_welfare = quantities.hh.welfare;
fprintf('\nBenchmark Welfare = %0.6f\n', benchmark_welfare)
fprintf('\nCounterfactual Welfare = %0.6f\n', counterfactual_welfare)
fprintf('\nEquilibrium dividend tax = %0.6f\n', taxes.hh.dividend)

save(fullfile('Results','experiment3_results.mat'),'prices','quantities','taxes')




end