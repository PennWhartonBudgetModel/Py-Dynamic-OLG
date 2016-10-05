% Jorge | 2016-09-28
% 
% Tax function.  To be expanded to include income and consumption calculations.
% 
% Any dependent mex functions should be regenerated after making changes to this function.  These include but may not be limited to:
% 
%       solve_dynamic_optimization
%       calculate_static_taxes
% 
% 


function [fincome, ftax, sstax, fcap] ...
    ...
    = calculate_tax(fincome_lab, kgrid_ik, cap_share, rate_cap, debt_share, rate_gov, cap_tax_share, tau_cap, tau_capgain, exp_subsidy, ...
                    avg_deduc, coefs, limit, X, mpci, rpci, tau_ss, v_ss_max, clinton, year, q_tobin, q_tobin0) %#codegen

% Enforce function inlining for C code generation
coder.inline('always');

fincome     = (rpci/mpci)*(cap_share*(rate_cap-1)*kgrid_ik*(1-cap_tax_share) + debt_share*(rate_gov-1)*kgrid_ik + fincome_lab);
fincome     = max(0,fincome);
deduc       = max(0, avg_deduc + coefs(1)*fincome + coefs(2)*fincome^0.5);
fincome_eff = max(fincome - deduc, 0);
fincome     = (mpci/rpci)*fincome;

ftax = limit*(fincome_eff - (fincome_eff.^(-X(1)) + (X(2))).^(-1/X(1)));
mtr  = limit - (limit*(fincome_eff^(-X(1)) + X(2))^(-1/X(1))) / (X(2)*(fincome_eff^(X(1)+1)) + fincome_eff);
ftax = (mpci/rpci)*(ftax + clinton*max(0, mtr-0.28)*deduc);

sstax = tau_ss*min(fincome_lab, v_ss_max);

fcap = cap_share*kgrid_ik*( tau_cap*((rate_cap - 1) - exp_subsidy)*cap_tax_share + tau_capgain*(year == 1)*(q_tobin - q_tobin0)/q_tobin0 );

end