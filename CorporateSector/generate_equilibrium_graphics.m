function [] = generate_equilibrium_graphics()


% Load equilibrium prices
s = load(fullfile('Results','results.mat'),'prices');
prices = s.prices;
clear('s')

% Set limitations on equilibrium interest rates
s_hh = load(fullfile('Parameters','hh_parameters.mat'));
rate_ub = (s_hh.hh_params.discount_factor)^(-1);
rate_lb = 1;

% Determine invariant distribution of productivity shocks
tolerance = 1;
[~, ~, ~, ~, dist] = solve_hh_optimization_mex(s_hh.hh_params, prices, tolerance);

% Determine labor supply when inelastic
labor_supply = sum(sum(dist).*s_hh.hh_params.prod_shocks');    % Update if labor supply becomes elastically demanded

price_names = fieldnames(prices);
n_markets = length(price_names);
percent_deviation = .1;
n_deviations = 10;

supply = zeros(n_deviations,n_markets);
demand = zeros(n_deviations,n_markets);
x_prices = [];
tolerance = 0.0001;  % Tolerance on hh optimization.
for im = 1:n_markets
    im
    ub = prices.(price_names{im})*(1 + percent_deviation);
    lb = prices.(price_names{im})*(1 - percent_deviation);
    if strcmp(price_names{im},'rate')
        ub = 1.035;
        lb = 1.005;
    end
    
    price_deviation = linspace(lb,ub,n_deviations);
    test_prices = prices;
    for id = 1:n_deviations
        test_prices.(price_names{im}) = price_deviation(id);
        
        [asset_demand, output_total, labor_demand, inv_total] = firm_quantities(test_prices);
        [~, ~, asset_supply, consumption_total, ~] = solve_hh_optimization_mex(s_hh.hh_params, test_prices, tolerance);
        
        if strcmp(price_names{im},'consumption')
            supply(id,im) = output_total;
            demand(id,im) = consumption_total + inv_total;
        elseif strcmp(price_names{im},'rate')
            supply(id,im) = asset_supply;
            demand(id,im) = asset_demand;
        elseif strcmp(price_names{im},'wage')
            supply(id,im) = labor_supply;
            demand(id,im) = labor_demand;
        end
    end
    x_prices = [x_prices; price_deviation]; %#ok<AGROW>
    
end


figure
subplot(3,1,1)
plot(x_prices(1,:),supply(:,1),x_prices(1,:),demand(:,1),'LineWidth',3)
ylabel('Consumption')
xlabel('Output/Goods Price')

subplot(3,1,2)
plot(x_prices(2,:),supply(:,2),x_prices(2,:),demand(:,2),'LineWidth',3)
ylabel('Assets')
xlabel('Interest Rate')

subplot(3,1,3)
plot(x_prices(3,:),supply(:,3),x_prices(3,:),demand(:,3),'LineWidth',3)
ylabel('Labor')
xlabel('Wage')

        
    


desktop













end