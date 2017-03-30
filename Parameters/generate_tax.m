%%
% Generate tax policy parameters according to predefined plans.
% 
%%


function [] = generate_tax(plan)

% Identify parameter directory
param_dir = dirFinder.param();

readTPC = @(sheet, range) xlsread(fullfile(param_dir, 'TPC_Inputs.xlsx'), sheet, range);



%% Income Tax

% Read taxable income ranges
income_range_bins = readTPC(1, 'A3:B622');  % Column A = lower bounds, Column B = upper bounds

% Read tax rates
if strcmp(plan,'base')
    range = 'C3:C622';
elseif strcmp(plan,'trump')    
    range = 'E3:E622';
elseif strcmp(plan,'ryan')
    range = 'F3:F622';
end
effective_tax_rates = readTPC(1, range);    % Effective tax rates on taxable income

% Estimate income tax function
n_testpoints = 1000;   % Determines how finely to discretize income "data" points for estimation.
[X, limit] = estimate_inctax_function(income_range_bins, effective_tax_rates, n_testpoints); 


%% Deductions

if strcmp(plan,'base')
    range = 'C3:C622';
elseif strcmp(plan,'trump')
    range = 'E3:E622';
elseif strcmp(plan,'ryan')
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


regressors = [discretized_income.^0; discretized_income.^1; discretized_income.^(1/2)]';   % Note the functional for specification implied here: f(x) = b_1*x^0 + b_2*x^1 + b_3*x^(1/2)
[coefficients,~,~,~,~] = regress(deductions',regressors);


n_regs = min(size(regressors));    % Assumes that the number of variables < number of observations.  Must be true for estimation.
% Below changes the estimates to the order of polynomial currently hard coded
if n_regs<5
    coefs = [coefficients(2:n_regs)', zeros(1,4-(n_regs-1))];
elseif n_regs>=5
    warning('Deduction polynomial order potential mismatch. Curtailing at 4.')
    coefs = coefficients(2:5);
end

avg_deduc = coefficients(1);



%% Business Tax

if strcmp(plan,'base')
    range = 'B2:B5';
elseif strcmp(plan,'trump')
    range = 'D2:D5';
elseif strcmp(plan,'ryan')
    range = 'E2:E5';
end

bus_params    = readTPC(3, range);    % gives business parameters of policy proposal
cap_tax_share = bus_params(2); 
exp_share     = bus_params(4); 
tau_cap       = bus_params(1); 
tau_capgain   = 0; 


% Generating .mat files for future access by non-Windows OS's if running on Windows OS and save_results is true.
save(fullfile(param_dir, ['param_inctax_' plan '.mat']), 'avg_deduc', 'coefs', 'limit', 'X')
save(fullfile(param_dir, ['param_bustax_' plan '.mat']), 'cap_tax_share', 'exp_share', 'tau_cap', 'tau_capgain')

   
end


function [X, limit] = estimate_inctax_function(income_range_bins, effective_tax_rates, n_testpoints)

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

[X,~] = fminsearch(@(phi) gouveiastrauss(phi,discretized_income(2:end),inctax_bills(2:end)),[.8 .01 .36], opts1);
limit = X(3);
X = X(1:2); 

function sum_sq = gouveiastrauss(phi,x1,y1)
    ytax = phi(3).*(x1 - (x1.^(-phi(1)) + phi(2)).^(-1/phi(1)));    % Calculate taxable income
    ybar = ytax./x1;    % Calculate effective income
    sum_sq = sum((ybar - y1./x1).^2);    % Fit function to effective income.
end

end