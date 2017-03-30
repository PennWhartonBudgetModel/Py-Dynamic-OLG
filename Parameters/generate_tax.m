%%
% Generate tax policy parameters according to predefined plans.
% 
%%


function [] = generate_tax()

% Identify parameter directory
param_dir = dirFinder.param();

% Define anonymous function for reading from TPC xlsx file
readTPC = @(sheet, range) xlsread(fullfile(param_dir, 'TPC_Inputs.xlsx'), sheet, range);


function sum_sq = gouveiastrauss(phi,x1,y1)
    ytax = phi(3).*(x1 - (x1.^(-phi(1)) + phi(2)).^(-1/phi(1)));    % Calculate taxable income
    ybar = ytax./x1;    % Calculate effective income
    sum_sq = sum((ybar - y1./x1).^2);    % Fit function to effective income.
end


for taxplan_ = {'base', 'trump', 'ryan'}, taxplan = taxplan_{1};
    
    %% Income Tax

    % Read taxable income ranges
    income_range_bins = readTPC(1, 'A3:B622');  % Column A = lower bounds, Column B = upper bounds

    % Read tax rates
    if strcmp(taxplan,'base')
        range = 'C3:C622';
    elseif strcmp(taxplan,'trump')    
        range = 'E3:E622';
    elseif strcmp(taxplan,'ryan')
        range = 'F3:F622';
    end
    effective_tax_rates = readTPC(1, range);    % Effective tax rates on taxable income

    % Estimate income tax function
    n_testpoints = 1000;   % Determines how finely to discretize income "data" points for estimation.

    % Setting optimization options.
    opts1 = optimset('TolFun',1e-10,'MaxIter',10e10,'MaxFunEvals',10e10);

    % Generating variables for income tax estimation.
    income_upper_bound = income_range_bins(end,2);
    discretized_income = linspace(0,income_upper_bound,n_testpoints);

    inctax_bills = zeros(1,n_testpoints);

    for i1 = 1:n_testpoints
        income_bin_index = find(discretized_income(i1)<=income_range_bins(:,2),1,'first');   % returns integer corresponding to taxable income range provided by TPC
        inctax_bills(i1) = discretized_income(i1)*effective_tax_rates(income_bin_index);
    end

    X = fminsearch(@(phi) gouveiastrauss(phi,discretized_income(2:end),inctax_bills(2:end)),[0.8, 0.01, 0.36], opts1);

    s.(taxplan).limit = X(3);
    s.(taxplan).X     = X(1:2); 


    %% Deductions

    if strcmp(taxplan,'base')
        range = 'C3:C622';
    elseif strcmp(taxplan,'trump')
        range = 'E3:E622';
    elseif strcmp(taxplan,'ryan')
        range = 'F3:F622';
    end

    income_upper_bound = 2*10e5;
    income_upper_bound = min(income_range_bins(end,2), income_upper_bound);    % Enforcing upper bound provided by TPC.
    a_log = 1;
    b_log = log(income_upper_bound)/log(10);
    n_testpoints = 1000;
    discretized_income = logspace(a_log,b_log,n_testpoints);
    n = length(discretized_income);

    taxable_income_ratio = readTPC(2, range);
    deductions = zeros(1,n);
    for i1 = 1:n
        income_bin_index = find(discretized_income(i1)<=income_range_bins(:,2),1,'first');   % Returns integer corresponding to taxable income range provided by TPC.
        deductions(i1) = discretized_income(i1)*(1-taxable_income_ratio(income_bin_index));
    end


    regressors = [discretized_income.^0; discretized_income.^1; discretized_income.^(1/2)]';   % Note the functional form implied here: f(x) = b_1*x^0 + b_2*x^1 + b_3*x^(1/2)
    coefficients = regress(deductions',regressors);


    n_regs = min(size(regressors));    % Assumes that the number of variables < number of observations.  Must be true for estimation.
    % Below changes the estimates to the order of polynomial currently hard coded
    if n_regs<5
        s.(taxplan).coefs = [coefficients(2:n_regs)', zeros(1,4-(n_regs-1))];
    elseif n_regs>=5
        warning('Deduction polynomial order potential mismatch. Curtailing at 4.')
        s.(taxplan).coefs = coefficients(2:5);
    end

    s.(taxplan).avg_deduc = coefficients(1);


    %% Business Tax

    if strcmp(taxplan,'base')
        range = 'B2:B5';
    elseif strcmp(taxplan,'trump')
        range = 'D2:D5';
    elseif strcmp(taxplan,'ryan')
        range = 'E2:E5';
    end

    bus_params = readTPC(3, range);    % gives business parameters of policy proposal

    s.(taxplan).cap_tax_share = bus_params(2); 
    s.(taxplan).exp_share     = bus_params(4); 
    s.(taxplan).tau_cap       = bus_params(1); 
    s.(taxplan).tau_capgain   = 0; 

end

% Save parameters
save(fullfile(param_dir, 'tax.mat'), '-struct', 's');

end