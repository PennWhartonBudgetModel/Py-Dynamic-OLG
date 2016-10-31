%%
% Resource calculation function.  Also calculates tax terms.  Applies to both working and retirement ages.
% 
% Any dependent mex functions should be regenerated after making changes to this function.  These include but may not be limited to:
% 
%       solve_dynamic_optimization
%       calculate_static_taxes
% 
%%

function [resources, fincome, ftax, sstax, fcap] ...
    ...
    = calculate_resources(fincome_lab, kgrid_ik, ...
                          cap_share, rate_cap, debt_share, rate_gov, cap_tax_share, tau_cap, tau_capgain, exp_subsidy, ...
                          avg_deduc, coefs, limit, X, mpci, rpci, tau_ss, v_ss_max, clinton, year, q_tobin, q_tobin0, ...
                          cap_inc_ik, beq, ss_tax_cred) %#codegen

% Enforce function inlining for C code generation
coder.inline('always');

% Calculate taxable income
fincome     =   cap_share *(rate_cap-1)*kgrid_ik*(1-cap_tax_share) ...
              + debt_share*(rate_gov-1)*kgrid_ik ...
              + (1-ss_tax_cred)*fincome_lab;
fincome     = max(0, (rpci/mpci)*fincome);
deduc       = max(0, avg_deduc + coefs(1)*fincome + coefs(2)*fincome^0.5);
fincome_eff = max(fincome - deduc, 0);
fincome     = (mpci/rpci)*fincome;

% Calculate personal income tax (PIT)
ftax = limit*(fincome_eff - (fincome_eff.^(-X(1)) + (X(2))).^(-1/X(1)));
mtr  = limit - (limit*(fincome_eff^(-X(1)) + X(2))^(-1/X(1))) / (X(2)*(fincome_eff^(X(1)+1)) + fincome_eff);
ftax = (mpci/rpci)*(ftax + clinton*max(0, mtr-0.28)*deduc);

% Calculate Social Security tax
sstax = tau_ss*min(fincome_lab, v_ss_max);

% Calculate capital income tax (CIT)
fcap = cap_share*kgrid_ik*( tau_cap*((rate_cap-1) - exp_subsidy)*cap_tax_share ...
       + tau_capgain*(year == 1)*(q_tobin - q_tobin0)/q_tobin0 );

% Calculate available resources
resources = cap_inc_ik + fincome_lab - (ftax + sstax + fcap) + beq ...
            + (year == 1)*kgrid_ik*cap_share*(q_tobin - q_tobin0)/q_tobin0;

end