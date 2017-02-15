
% solve_firm_optimization_test
s = paramGenerator.firm;
s_tax = paramGenerator.tax;

prices.fund = 0.98;
prices.wage = [1, 1];
prices.rate = .98;
tolerance   = 50;    

s.discount_factor = prices.rate;

[capital_total, labor_total, eq_total, inv_total, adjcost_total, V_total, output_total, tax_total, dist] =...
    solve_firm_optimization_mex(prices, s_tax, s, tolerance);
