function [] = explore_results()
% Access agent optimization files
addpath('..')

% Load taxes and equilibrium results
s_tax = load(fullfile('..','Parameters','tax_parameters.mat'));
taxes = s_tax.taxes;
s = load('results.mat');
prices = s.prices;
quantities = s.quantities;
s_hh = load(fullfile('..','Parameters','hh_parameters.mat'));



k_ratio = quantities.firm.capital_total/quantities.firm.output_total;
i_ratio = quantities.firm.inv_total/quantities.firm.output_total;

t_y = quantities.government.tax_revenue/quantities.firm.output_total;
cit_y = quantities.firm.tax_total/quantities.firm.output_total;
hh_y = quantities.hh.tax_total/quantities.firm.output_total;
pit_y = (quantities.hh.tax_total-(1 - taxes.hh.div_inc_share)*taxes.hh.dividend*quantities.firm.dividend)/quantities.firm.output_total;
effective_pit = (quantities.hh.tax_total-(1 - taxes.hh.div_inc_share)...
                *taxes.hh.dividend*quantities.firm.dividend)/prices.wage;


print_information

income_gini = calculate_gini('income');
fprintf('\nIncome Gini  = %0.4f\n', income_gini)

wealth_gini = calculate_gini('wealth');
fprintf('\nWealth Gini  = %0.4f\n', wealth_gini)





function gini = calculate_gini(gini_variable)
tolerance = 0.0001;
[~, ~, ~, ~, ~, ~, dist] = solve_hh_optimization_mex(s_hh.hh_params, prices, taxes, quantities.firm.dividend, tolerance);
    
if strcmp(gini_variable,'income')    
    inc = repmat(s_hh.hh_params.prod_shocks',[1,s_hh.hh_params.ns])'...
            + repmat(s_hh.hh_params.shares_grid',[1,s_hh.hh_params.n_prodshocks])...
            .*quantities.firm.dividend;
        
    gini_var = inc;
    
elseif strcmp(gini_variable,'wealth')
    wth = repmat(s_hh.hh_params.shares_grid',[1,s_hh.hh_params.n_prodshocks]).*prices.fund;
    
    gini_var = wth;
   
end


% Calculate the gini coefficient
[gini_var, index] = sort(gini_var(:));
dist = dist(index);

if dist(end)<10e-8
    ub_index = find(fliplr(dist')>10e-8,1,'first')-1;    % Starts from the back and looks for first positive
elseif dist(end)>=10e-8
    ub_index = 0;
end
gini_var = gini_var(1:end-ub_index);
dist     = dist    (1:end-ub_index);

n_bins = 25;
var_grid = linspace(gini_var(1),gini_var(end),n_bins);

lorenz = zeros(1,n_bins);
measure = zeros(1,n_bins);
var_total = sum(gini_var.*dist);
for ib = 1:n_bins
    varlim_index = find(gini_var<=var_grid(ib),1,'last');
    lorenz(ib) = sum(gini_var(1:varlim_index).*dist(1:varlim_index))/var_total;
    measure(ib) = sum(dist(1:varlim_index));
end

lorenz  = [0, lorenz];
measure = [0, measure];

figure
plot(measure,lorenz,measure,measure,'LineWidth',4)
title(gini_variable)


% Solve gini


gini = 2*integral(@(x) gini_diff(x,lorenz,measure),0,1);


function diff = gini_diff(x,lorenz,measure)
    diff = interp1(measure-lorenz,measure,x,'linear');
end

end











function [] = print_information()

% General Aggregate Ratios
fprintf('\nCapital to Output = %0.4f\n'    , k_ratio)
fprintf('\nInvestment to Output = %0.4f\n' , i_ratio)

% Tax Ratios
fprintf('\nTotal Tax Revenue to Output = %0.4f\n'    , t_y           )
fprintf('\nTotal HH Tax Revenue to Output = %0.4f\n' , hh_y          )
fprintf('\nCIT to Output = %0.4f\n'                  , cit_y         )
fprintf('\nPIT to Output = %0.4f\n'                  , pit_y         )
fprintf('\nEffective PIT rate = %0.4f\n'             , effective_pit )



end


end