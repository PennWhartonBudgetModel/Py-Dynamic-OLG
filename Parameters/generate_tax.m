%%
% Generate tax policy parameters according to predefined plans.
% 
%%


function [] = generate_tax()

% Identify parameter directory
param_dir = dirFinder.param();

% Define anonymous function for reading from TPC data file
readTPC = @(sheet, range) xlsread(fullfile(param_dir, 'TPC_Inputs.xlsx'), sheet, range);


function [sqe] = gouveiastrauss(pit_coefs, incv, inct)
    
    % Calculate analytical tax rate
    inct_ = pit_coefs(1)*(incv - (incv.^-pit_coefs(2) + pit_coefs(3)).^(-1/pit_coefs(2))) ./ incv;
    
    % Calculate squared error
    sqe = sum((inct_ - inct).^2);
    
end


for taxplan_ = {'base', 'trump', 'ryan'}, taxplan = taxplan_{1};
        
    % Read taxable income thresholds
    incthresholds = readTPC(1, 'B3:B622');
    
    % Read effective tax rates and taxable income ratios from TPC data file
    switch taxplan
        case 'base' , col = 'C';
        case 'trump', col = 'E';
        case 'ryan' , col = 'F';
    end
    range = sprintf('%1$c3:%1$c622', col);
    
    taxrates  = readTPC(1, range);
    taxratios = readTPC(2, range);
    
    % Define vector of sample incomes and find corresponding tax rates
    incv = linspace(1, incthresholds(end), 1e3)';
    inct = arrayfun(@(inc) taxrates(find(inc <= incthresholds, 1)), incv);
    
    % Fit income tax function using least squares
    pit_coefs = fminsearch(@(pit_coefs) gouveiastrauss(pit_coefs, incv, inct), ...
                           [0.36, 0.8, 0.01], optimset('TolFun', 1e-10, 'MaxIter', 1e11, 'MaxFunEvals', 1e11));
    
    % Define vector of sample incomes and calculate corresponding deductions
    % (Upper bound on deductable income provided by TPC)
    incv = logspace(0, log10(2e6), 1e3)';
    incd = arrayfun(@(inc) inc * (1 - taxratios(find(inc <= incthresholds, 1))), incv);
    
    % Fit deduction function using multilinear regression
    exps = [0, 1, 1/2]; % f(x) = b_1 + b_2*x + b_3*x^(1/2)
    deduc_coefs = regress(incd, repmat(incv, [1,length(exps)]) .^ repmat(exps, [length(incv),1]))';
    
    % Store fitted income tax and deduction function coefficients
    s.(taxplan).pit_coefs   = pit_coefs  ;
    s.(taxplan).deduc_coefs = deduc_coefs;
    
    % Read business tax parameters from TPC data file
    switch taxplan
        case 'base' , col = 'B';
        case 'trump', col = 'D';
        case 'ryan' , col = 'E';
    end
    range = sprintf('%1$c2:%1$c5', col);
    p = readTPC(3, range);
    
    % Store business tax parameters
    s.(taxplan).captaxshare = p(2); 
    s.(taxplan).expshare    = p(4); 
    s.(taxplan).taucap      = p(1); 
    s.(taxplan).taucapgain  = 0   ;
    
end


% Save parameters
save(fullfile(param_dir, 'tax.mat'), '-struct', 's');


end