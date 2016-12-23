function [Available_assets, output_total, capital_total, labor_total, eq_total, inv_total, adjcost_total, tax_total] = firm_quantities(prices, taxes)

if ~exist('prices','var')
    s_results = load(fullfile('Results','results.mat'));
    prices = s_results.prices;
end
if ~exist('taxes','var')
    s_taxes = load(fullfile('Parameters','tax_parameters.mat'));
    taxes = s_taxes.taxes;
end


% Load firm and household parameter structures
s_firm = load(fullfile('Parameters','firm_parameters.mat'));

% Specify settings for dynamic optimization subproblems
optim_options = optimset('Display', 'off', 'TolFun', 1e-4, 'TolX', 1e-4); %#ok<NASGU>
tolerance = 0.000001;

s_firm.firm_params.discount_factor = prices.rate;
[capital_total, labor_total, eq_total, inv_total, adjcost_total, V_total, output_total, tax_total, dist] = solve_firm_optimization_mex(prices, taxes, s_firm.firm_params, tolerance); %#ok<ASGLU>
R_mutual = V_total/(V_total - eq_total); %#ok<NASGU>
Available_assets = V_total;% - eq_total;

end

