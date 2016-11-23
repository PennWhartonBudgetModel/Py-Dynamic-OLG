% function determine_hh_properties()




% Generate grid of interest rates
s_hh = load(fullfile('Parameters','hh_parameters.mat'));
rate_ub = (s_hh.hh_params.discount_factor)^(-1)-0.005;
rate_lb = 1.0055;
n_rates = 20;
rates = linspace(rate_lb,rate_ub,n_rates);

assets_supplied = zeros(size(rates));

for ir = 1:n_rates
%     percent_done = ir/n_rates
    assets_supplied(ir) = asset_supply(rates(ir));

end

plot(assets_supplied,rates)


% end