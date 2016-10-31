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
    = calculate_resources(fincome_lab, kgrid_ik_, ...
                          cap_share_, rate_cap_, debt_share_, rate_gov_, cap_tax_share_, tau_cap_, tau_capgain_, exp_subsidy_, ...
                          avg_deduc_, coefs_, limit_, X_, mpci_, rpci_, tau_ss_, v_ss_max_, clinton_, year_, q_tobin_, q_tobin0_, ...
                          cap_inc_ik_, beq_, ss_tax_cred_) %#codegen

% Enforce function inlining for C code generation
coder.inline('always');

% Define parameters as persistent variables
persistent initialized
persistent kgrid_ik
persistent cap_share
persistent rate_cap
persistent debt_share
persistent rate_gov
persistent cap_tax_share
persistent tau_cap
persistent tau_capgain
persistent exp_subsidy
persistent avg_deduc
persistent coefs
persistent limit
persistent X
persistent mpci
persistent rpci
persistent tau_ss
persistent v_ss_max
persistent clinton
persistent year
persistent q_tobin
persistent q_tobin0
persistent cap_inc_ik
persistent beq
persistent ss_tax_cred

% Initialize parameters
if isempty(initialized)
    
    kgrid_ik      = 0;
    cap_share     = 0;
    rate_cap      = 0;
    debt_share    = 0;
    rate_gov      = 0;
    cap_tax_share = 0;
    tau_cap       = 0;
    tau_capgain   = 0;
    exp_subsidy   = 0;
    avg_deduc     = 0;
    coefs         = 0;
    limit         = 0;
    X             = 0;
    mpci          = 0;
    rpci          = 0;
    tau_ss        = 0;
    v_ss_max      = 0;
    clinton       = 0;
    year          = 0;
    q_tobin       = 0;
    q_tobin0      = 0;
    cap_inc_ik    = 0;
    beq           = 0;
    ss_tax_cred   = 0;
    
    initialized = true;
    
end

% Set parameters if provided
if (nargin > 1)
    
    kgrid_ik      = kgrid_ik_;
    cap_share     = cap_share_;
    rate_cap      = rate_cap_;
    debt_share    = debt_share_;
    rate_gov      = rate_gov_;
    cap_tax_share = cap_tax_share_;
    tau_cap       = tau_cap_;
    tau_capgain   = tau_capgain_;
    exp_subsidy   = exp_subsidy_;
    avg_deduc     = avg_deduc_;
    coefs         = coefs_;
    limit         = limit_;
    X             = X_;
    mpci          = mpci_;
    rpci          = rpci_;
    tau_ss        = tau_ss_;
    v_ss_max      = v_ss_max_;
    clinton       = clinton_;
    year          = year_;
    q_tobin       = q_tobin_;
    q_tobin0      = q_tobin0_;
    cap_inc_ik    = cap_inc_ik_;
    beq           = beq_;
    ss_tax_cred   = ss_tax_cred_;
    
    if isempty(fincome_lab), return, end
    
end

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