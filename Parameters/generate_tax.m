%%
% Generate tax policy parameters according to predefined plans.
% 
%%


function [] = generate_tax()

% Identify parameter directory
param_dir = dirFinder.param();

% Define anonymous function for reading from TPC xlsx file
readTPC = @(sheet, range) xlsread(fullfile(param_dir, 'TPC_Inputs.xlsx'), sheet, range);


function [sum_sq] = gouveiastrauss(phi, x1, y1)
    
    % Calculate taxable income
    ytax = phi(3)*(x1 - (x1.^-phi(1) + phi(2)).^(-1/phi(1)));
    
    % Fit function to effective income
    sum_sq = sum(((ytax - y1)./x1).^2);
    
end


for taxplan_ = {'base', 'trump', 'ryan'}, taxplan = taxplan_{1};
    
    %% Income Tax
    
    % Read taxable income ranges
    income_range_bins = readTPC(1, 'A3:B622');  % (Column A = lower bounds, Column B = upper bounds)
    
    % Read tax rates
    switch taxplan
        case 'base' , col = 'C';
        case 'trump', col = 'E';
        case 'ryan' , col = 'F';
    end
    range = sprintf('%1$c3:%1$c622', col);
    effective_tax_rates = readTPC(1, range);    % Effective tax rates on taxable income
    
    % Estimate income tax function
    n_testpoints = 1000;   % (Determines how finely to discretize income data points for estimation)
    
    % Generate variables for income tax estimation
    income_upper_bound = income_range_bins(end,2);
    discretized_income = linspace(0, income_upper_bound, n_testpoints);
    
    inctax_bills = zeros(1,n_testpoints);
    
    for i = 1:n_testpoints
        income_bin_index = find(discretized_income(i) <= income_range_bins(:,2), 1, 'first');   % (Returns integer corresponding to taxable income range provided by TPC)
        inctax_bills(i) = discretized_income(i)*effective_tax_rates(income_bin_index);
    end
    
    phi = fminsearch(@(phi) gouveiastrauss(phi, discretized_income(2:end), inctax_bills(2:end)), ...
                     [0.8, 0.01, 0.36], optimset('TolFun', 1e-10, 'MaxIter', 1e11, 'MaxFunEvals', 1e11));
    
    s.(taxplan).limit = phi(3);
    s.(taxplan).X     = phi(1:2); 
    
    
    %% Deductions
    
    income_upper_bound = 2e6;
    income_upper_bound = min(income_range_bins(end,2), income_upper_bound);    % (Enforces upper bound provided by TPC)
    a_log = 1;
    b_log = log10(income_upper_bound);
    n_testpoints = 1000;
    discretized_income = logspace(a_log, b_log, n_testpoints);
    n = length(discretized_income);
    
    switch taxplan
        case 'base' , col = 'C';
        case 'trump', col = 'E';
        case 'ryan' , col = 'F';
    end
    range = sprintf('%1$c3:%1$c622', col);
    taxable_income_ratio = readTPC(2, range);
    
    deductions = zeros(1,n);
    for i = 1:n
        income_bin_index = find(discretized_income(i) <= income_range_bins(:,2), 1, 'first');   % (Returns integer corresponding to taxable income range provided by TPC)
        deductions(i) = discretized_income(i)*(1-taxable_income_ratio(income_bin_index));
    end
    
    regressors = [discretized_income.^0; discretized_income.^1; discretized_income.^(1/2)]';    % (Functional form: f(x) = b_1*x^0 + b_2*x^1 + b_3*x^(1/2))
    coefficients = regress(deductions',regressors);
    
    n_regs = min(size(regressors));     % (Assumes that the number of variables < number of observations, which must be true for estimation)
    
    % Limit order of polynomial
    if n_regs < 5
        s.(taxplan).coefs = [coefficients(2:n_regs)', zeros(1,4-(n_regs-1))];
    else
        s.(taxplan).coefs = coefficients(2:5);
    end
    
    s.(taxplan).avg_deduc = coefficients(1);
    
    
    %% Business Tax
    
    % Read business tax parameters
    switch taxplan
        case 'base' , col = 'B';
        case 'trump', col = 'D';
        case 'ryan' , col = 'E';
    end
    range = sprintf('%1$c2:%1$c5', col);
    bus_params = readTPC(3, range);
    
    s.(taxplan).cap_tax_share = bus_params(2); 
    s.(taxplan).exp_share     = bus_params(4); 
    s.(taxplan).tau_cap       = bus_params(1); 
    s.(taxplan).tau_capgain   = 0; 
    
end


% Save parameters
save(fullfile(param_dir, 'tax.mat'), '-struct', 's');


end