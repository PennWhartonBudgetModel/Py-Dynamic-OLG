function [assets_supplied, consumption_total] = hh_quantities(prices)

s_hh = load(fullfile('Parameters','hh_parameters.mat'));

[V, aopt, assets_supplied, consumption_total, dist] = solve_hh_optimization_mex(s_hh.hh_params,prices); %#ok<ASGLU>
if max(dist(end,:))>1.0e-8
    % Prevent solving constrained optimization problem.
    warning('Distribution endpoint has non-zero mass.  Increase grid upper bound.')
    dist(end,:)
end
if R*s_hh.hh_params.discount_factor>=1
    warning('b*R>=1.  Reduce R and try solving again.')
    display(sprintf('R = %0.4f',R))
end

optim_options = optimset('Display', 'off', 'TolFun', 1e-8, 'TolX', 1e-4);

if s_hh.hh_params.n_prodshocks==1
    
    assets_supplied = fsolve(@(x) find_ss_assets(x,aopt,s_hh.hh_params.asset_grid'),s_hh.hh_params.asset_grid(end)-.1,optim_options);
    consumption_total = s_hh.hh_params.prod_shocks(1) + (R - 1)*assets_supplied;
    
end

end



function ss_deviation = find_ss_assets(x,aopt,asset_grid)

ss_deviation = interp1(asset_grid,aopt - asset_grid,x,'spline');

end


