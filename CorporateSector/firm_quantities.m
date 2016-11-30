function [Available_assets, output_total, labor_total] = firm_quantities(prices)

if ~exist('prices','var')
    prices.consumption = 1;
    prices.rate        = 1.03;
    prices.wage        = 1;
end

% Load firm and household parameter structures
s_firm = load(fullfile('Parameters','firm_parameters.mat'));

% Specify settings for dynamic optimization subproblems
optim_options = optimset('Display', 'off', 'TolFun', 1e-4, 'TolX', 1e-4); %#ok<NASGU>

s_firm.firm_params.discount_factor = 1/prices.rate;
[capital_total, labor_total, eq_total, V_total, output_total, dist] = solve_firm_optimization_mex(prices, s_firm.firm_params); %#ok<ASGLU>
R_mutual = V_total/(V_total - eq_total); %#ok<NASGU>
Available_assets = V_total - eq_total;

end

