%%
% Generate tax policy parameters according to predefined plans.
% 
%%


function [] = generate_tax()

% Identify parameter directory
param_dir = dirFinder.param();

% Define anonymous function for reading from TPC data file
readTPC = @(sheet, range) xlsread(fullfile(param_dir, 'TPC_Inputs.xlsx'), sheet, range);


function [sum_sq] = gouveiastrauss(pit_coefs, x1, y1)
    
    % Calculate taxable income
    ytax = pit_coefs(1)*(x1 - (x1.^-pit_coefs(2) + pit_coefs(3)).^(-1/pit_coefs(2)));
    
    % Fit function to effective income
    sum_sq = sum(((ytax - y1)./x1).^2);
    
end


for taxplan_ = {'base', 'trump', 'ryan'}, taxplan = taxplan_{1};
    
    %% Income Tax
    
    % Read taxable income thresholds
    income_thresholds = readTPC(1, 'B3:B622');
    
    % Read tax rates for plans
    switch taxplan
        case 'base' , col = 'C';
        case 'trump', col = 'E';
        case 'ryan' , col = 'F';
    end
    range = sprintf('%1$c3:%1$c622', col);
    effective_tax_rates = readTPC(1, range);    % Effective tax rates on taxable income
    
    % Estimate income tax function
    discretized_income = linspace(1, income_thresholds(end), 1e3);
    
    inctax_bills = zeros(1,1e3);
    for i = 1:1e3
        income_bin_index = find(discretized_income(i) <= income_thresholds, 1, 'first');   % (Returns integer corresponding to taxable income range provided by TPC)
        inctax_bills(i) = discretized_income(i)*effective_tax_rates(income_bin_index);
    end
    
    pit_coefs = fminsearch(@(pit_coefs) gouveiastrauss(pit_coefs, discretized_income, inctax_bills), ...
                           [0.36, 0.8, 0.01], optimset('TolFun', 1e-10, 'MaxIter', 1e11, 'MaxFunEvals', 1e11));
    
    s.(taxplan).pit_coefs = pit_coefs; 
    
    
    %% Deductions
    
    % (Upper bound provided by TPC)
    discretized_income = logspace(0, log10(2e6), 1e3);
    
    switch taxplan
        case 'base' , col = 'C';
        case 'trump', col = 'E';
        case 'ryan' , col = 'F';
    end
    range = sprintf('%1$c3:%1$c622', col);
    taxable_income_ratio = readTPC(2, range);
    
    deductions = zeros(1,1e3);
    for i = 1:1e3
        income_bin_index = find(discretized_income(i) <= income_thresholds, 1, 'first');   % (Returns integer corresponding to taxable income range provided by TPC)
        deductions(i) = discretized_income(i)*(1 - taxable_income_ratio(income_bin_index));
    end
    
    % Define regression functional form as f(x) = b_1*x^0 + b_2*x^1 + b_3*x^(1/2)
    regressors = [discretized_income.^0; discretized_income.^1; discretized_income.^(1/2)]';
    coefficients = regress(deductions', regressors);
    
    n_regs = min(size(regressors));     % (Assumes that the number of variables < number of observations, which must be true for estimation)
    
    % Limit order of polynomial
    if n_regs < 5
        s.(taxplan).coefs = [coefficients(2:n_regs)', zeros(1,4-(n_regs-1))];
    else
        s.(taxplan).coefs = coefficients(2:5);
    end
    
    s.(taxplan).avg_deduc = coefficients(1);
    
    
    %% Business Tax
    
    % Read business tax parameters from TPC data file
    switch taxplan
        case 'base' , col = 'B';
        case 'trump', col = 'D';
        case 'ryan' , col = 'E';
    end
    range = sprintf('%1$c2:%1$c5', col);
    p = readTPC(3, range);
    
    s.(taxplan).cap_tax_share = p(2); 
    s.(taxplan).exp_share     = p(4); 
    s.(taxplan).tau_cap       = p(1); 
    s.(taxplan).tau_capgain   = 0; 
    
end


% Save parameters
save(fullfile(param_dir, 'tax.mat'), '-struct', 's');


end