function [assets_supplied] = asset_supply(R)

s_hh = load(fullfile('Parameters','hh_parameters.mat'));

[V, aopt, assets_supplied, dist] = solve_hh_optimization_mex(s_hh.hh_params,R); %#ok<ASGLU>

optim_options = optimset('Display', 'off', 'TolFun', 1e-4, 'TolX', 1e-4);

if s_hh.hh_params.n_prodshocks==1
    
    assets_supplied = fsolve(@(x) find_ss_assets(x,aopt,s_hh.hh_params.asset_grid'),s_hh.hh_params.asset_grid(end)-.1,optim_options);
    
end

end



function ss_deviation = find_ss_assets(x,aopt,asset_grid)

ss_deviation = interp1(asset_grid,aopt - asset_grid,x,'linear');

end


