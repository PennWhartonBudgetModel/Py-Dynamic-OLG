function [] = solve_experiment2()
%% Dividend versus corporate (with full expensing) when all other are zero.


%% From the zero-tax structure, set CIT and solve equilibrium.
addpath('Parameters')
taxes = generate_tax_parameters(true);
rmpath('Parameters')

taxes.firm.income = .15;
[~, s_quantities] = solve_ss_equilibrium_global([], taxes);
benchmark_expenditures = s_quantities.government.government_expenditures;
benchmark_welfare = s_quantities.hh.welfare;
clear('s_quantities')


%% Find the dividend tax that clears market
taxes.firm.income = 0;

dub = .175;
dlb = .125;
tolerance = .00001;
while true
    taxes.hh.dividend = (dub+dlb)/2;
    [~, quantities] = solve_ss_equilibrium_global([], taxes);
    rev = quantities.government.tax_revenue;
    error = abs(rev - benchmark_expenditures);
    if error<tolerance, break, end
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

save(fullfile('Results','experiment2_results.mat'),'prices','quantities','taxes')





end