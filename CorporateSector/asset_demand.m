function [Available_assets, discount_opt] = asset_demand(R)

% Load firm and household parameter structures
s_firm = load(fullfile('Parameters','firm_parameters.mat'));

% Generating discount factor grid
discount_lb = .9;
discount_ub = .9999;
n_disc = 10;
discount_factors = linspace(discount_lb,discount_ub,n_disc);
R_m = zeros(size(discount_factors));

% Specify settings for dynamic optimization subproblems
optim_options = optimset('Display', 'off', 'TolFun', 1e-4, 'TolX', 1e-4);


% Solving mutual fund returns over grid of discount factors
for id = 1:length(discount_factors)
    % Guess a firm discount factor
    prices.consumption = 1;
    s_firm.firm_params.discount_factor = discount_factors(id);
    [capital_total, eq_total, V_total, dist] = solve_firm_optimization_mex(prices, s_firm.firm_params); %#ok<ASGLU>
    R_m(id) = V_total/(V_total - eq_total);    % Mutual fund return
end

if R_m(end)>R
    error('Target interest rate below minimum portfolio return.  Recommend increasing upper bound of discount grid.')
end
if R_m(1)<R
    error('Target interest rate above maximum portfolio return.  Recommend reducing lower bound of discount grid.')
end

discount_opt = fsolve(@(x) find_eqm_discount(x,discount_factors,R_m,R), (discount_lb + discount_ub)/2, optim_options);
s_firm.firm_params.discount_factor = discount_opt;
[capital_total, eq_total, V_total, dist] = solve_firm_optimization_mex(prices, s_firm.firm_params); %#ok<ASGLU>
R_mutual = V_total/(V_total - eq_total); %#ok<NASGU>
Available_assets = V_total - eq_total;



end




function root_finder = find_eqm_discount(x,discount_factors,R_m,R)

    root_finder = interp1(discount_factors,R_m,x) - R;

end
