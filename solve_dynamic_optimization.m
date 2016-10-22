%%
% Dynamic optimization solver.  Generalized to support both steady state and transition path solvers.
% 
% For open economy transition paths, the following arguments should not vary with year -- i.e. all columns should be identical:
% 
%       wages
%       cap_shares
%       debt_shares
%       rate_caps
%       rate_govs
%       rate_tots
%       cap_incs
% 
% For steady states, the same is true for the following additional arguments:
% 
%       beqs
%       ss_tax
%       taxmax
%       ss_benefit
% 
% q_tobin0 should equal q_tobin for steady states and transition paths using the base tax plan.
% 
%%


function [V, Vss, kopt, koptss, labopt, bopt, ...
          fedincome, fedincomess, ...
          fitax,     fitaxss,     ...
          fsstax,    fsstaxss,    ...
          fcaptax,   fcaptaxss]   ...
          ...
            = solve_dynamic_optimization(...
                beqs, T, Tr, Tss, birthyear, ...
                bgrid, kgrid, nb, nk, nz, idem, ...
                mpci, rpci, cap_tax_share, ss_tax_cred, ...
                surv, tau_cap, tau_capgain, ss_tax, tr_z, taxmax, z, ...
                wages, cap_shares, debt_shares, rate_caps, rate_govs, rate_tots, cap_incs, exp_subsidys, q_tobin, q_tobin0, Vbeq, ...
                avg_deduc, clinton, coefs, limit, X, ss_benefit, ...
                beta, gamma, sigma) %#codegen


% Define argument properties for C code generation
T_max   = 100;
Tss_max =  75;
nb_max  =  50;
nk_max  =  50;
nz_max  =  50;

assert( isa(beqs,           'double') && (size(beqs,            1) == 1         ) && (size(beqs,            2) <= Tss_max   ) );
assert( isa(T,              'double') && (size(T,               1) == 1         ) && (size(T,               2) == 1         ) );
assert( isa(Tr,             'double') && (size(Tr,              1) == 1         ) && (size(Tr,              2) == 1         ) );
assert( isa(Tss,            'double') && (size(Tss,             1) == 1         ) && (size(Tss,             2) == 1         ) );
assert( isa(birthyear,      'double') && (size(birthyear,       1) == 1         ) && (size(birthyear,       2) == 1         ) );

assert( isa(bgrid,          'double') && (size(bgrid,           1) <= nb_max    ) && (size(bgrid,           2) == 1         ) );
assert( isa(kgrid,          'double') && (size(kgrid,           1) <= nk_max    ) && (size(kgrid,           2) == 1         ) );
assert( isa(nb,             'double') && (size(nb,              1) == 1         ) && (size(nb,              2) == 1         ) );
assert( isa(nk,             'double') && (size(nk,              1) == 1         ) && (size(nk,              2) == 1         ) );
assert( isa(nz,             'double') && (size(nz,              1) == 1         ) && (size(nz,              2) == 1         ) );
assert( isa(idem,           'double') && (size(idem,            1) == 1         ) && (size(idem,            2) == 1         ) );

assert( isa(mpci,           'double') && (size(mpci,            1) == 1         ) && (size(mpci,            2) == 1         ) );
assert( isa(rpci,           'double') && (size(rpci,            1) == 1         ) && (size(rpci,            2) == 1         ) );
assert( isa(cap_tax_share,  'double') && (size(cap_tax_share,   1) == 1         ) && (size(cap_tax_share,   2) == 1         ) );
assert( isa(ss_tax_cred,    'double') && (size(ss_tax_cred,     1) == 1         ) && (size(ss_tax_cred,     2) == 1         ) );

assert( isa(surv,           'double') && (size(surv,            1) == 1         ) && (size(surv,            2) <= T_max     ) );
assert( isa(tau_cap,        'double') && (size(tau_cap,         1) == 1         ) && (size(tau_cap,         2) == 1         ) );
assert( isa(tau_capgain,    'double') && (size(tau_capgain,     1) == 1         ) && (size(tau_capgain,     2) == 1         ) );
assert( isa(ss_tax,         'double') && (size(ss_tax,          1) == 1         ) && (size(ss_tax,          2) <= Tss_max   ) );
assert( isa(tr_z,           'double') && (size(tr_z,            1) <= nz_max    ) && (size(tr_z,            2) <= nz_max    ) );
assert( isa(taxmax,         'double') && (size(taxmax,          1) == 1         ) && (size(taxmax,          2) <= Tss_max   ) );
assert( isa(z,              'double') && (size(z,               1) <= nz_max    ) && (size(z,               2) <= T_max     ) && (size(z, 3) == 2) );

assert( isa(wages,          'double') && (size(wages,           1) == 1         ) && (size(wages,           2) <= Tss_max   ) );
assert( isa(cap_shares,     'double') && (size(cap_shares,      1) == 1         ) && (size(cap_shares,      2) <= Tss_max   ) );
assert( isa(debt_shares,    'double') && (size(debt_shares,     1) == 1         ) && (size(debt_shares,     2) <= Tss_max   ) );
assert( isa(rate_caps,      'double') && (size(rate_caps,       1) == 1         ) && (size(rate_caps,       2) <= Tss_max   ) );
assert( isa(rate_govs,      'double') && (size(rate_govs,       1) == 1         ) && (size(rate_govs,       2) <= Tss_max   ) );
assert( isa(rate_tots,      'double') && (size(rate_tots,       1) == 1         ) && (size(rate_tots,       2) <= Tss_max   ) );
assert( isa(cap_incs,       'double') && (size(cap_incs,        1) <= nk_max    ) && (size(cap_incs,        2) <= Tss_max   ) );
assert( isa(exp_subsidys,   'double') && (size(exp_subsidys,    1) == 1         ) && (size(exp_subsidys,    2) <= Tss_max   ) );
assert( isa(q_tobin,        'double') && (size(q_tobin,         1) == 1         ) && (size(q_tobin,         2) == 1         ) );
assert( isa(q_tobin0,       'double') && (size(q_tobin0,        1) == 1         ) && (size(q_tobin0,        2) == 1         ) );
assert( isa(Vbeq,           'double') && (size(Vbeq,            1) <= nk_max    ) && (size(Vbeq,            2) == 1         ) );

assert( isa(avg_deduc,      'double') && (size(avg_deduc,       1) == 1         ) && (size(avg_deduc,       2) == 1         ) );
assert( isa(clinton,        'double') && (size(clinton,         1) == 1         ) && (size(clinton,         2) == 1         ) );
assert( isa(coefs,          'double') && (size(coefs,           1) == 1         ) && (size(coefs,           2) <= 10        ) );
assert( isa(limit,          'double') && (size(limit,           1) == 1         ) && (size(limit,           2) == 1         ) );
assert( isa(X,              'double') && (size(X,               1) == 1         ) && (size(X,               2) <= 10        ) );
assert( isa(ss_benefit,     'double') && (size(ss_benefit,      1) <= nb_max    ) && (size(ss_benefit,      2) <= Tss_max   ) );

assert( isa(beta,           'double') && (size(beta,            1) == 1         ) && (size(beta,            2) == 1         ) );
assert( isa(gamma,          'double') && (size(gamma,           1) == 1         ) && (size(gamma,           2) == 1         ) );
assert( isa(sigma,          'double') && (size(sigma,           1) == 1         ) && (size(sigma,           2) == 1         ) );


%% Initialization

% Define past years and effective retirement age
% (birthyear is defined relative to the present and can take both positive and negative values)
T_past = -min(0, birthyear);
Tr_eff = max(Tr, T_past);

% Find effective numbers of working and retirement years
N_w = Tr_eff - T_past;
N_r = T      - Tr_eff;

% Initialize arrays
V           = zeros(nk,nz,nb,N_w);

kopt        = zeros(nk,nz,nb,N_w);
labopt      = zeros(nk,nz,nb,N_w);
bopt        = zeros(nk,nz,nb,N_w);

fedincome   = zeros(nk,nz,nb,N_w);
fitax       = zeros(nk,nz,nb,N_w);
fsstax      = zeros(nk,nz,nb,N_w);
fcaptax     = zeros(nk,nz,nb,N_w);

Vss         = zeros(nk,nb,N_r+1);   % (1 extra time slice for initialization of backward induction)

koptss      = zeros(nk,nb,N_r);

fedincomess = zeros(nk,nb,N_r);
fitaxss     = zeros(nk,nb,N_r);
fsstaxss    = zeros(nk,nb,N_r);
fcaptaxss   = zeros(nk,nb,N_r);

% Specify settings for dynamic optimization subproblems
optim_options = optimset('Display', 'off', 'TolFun', 1e-4, 'TolX', 1e-4);


%% Retirement age dynamic optimization

for t = N_r:-1:1
    
    % Determine age and year, bounded by projection period
    age  = t + Tr_eff;
    year = max(1, min(age + birthyear, Tss));
    
    % Extract annual parameters
    ben         = ss_benefit(:,  year);
    beq         = beqs          (year);
    rate_cap    = rate_caps     (year);
    rate_gov    = rate_govs     (year);
    cap_share   = cap_shares    (year);
    debt_share  = debt_shares   (year);
    rate_tot    = rate_tots     (year);
    exp_subsidy = exp_subsidys  (year);
    
    EV = (1-surv(age))*Vbeq*ones(1,nb) + surv(age)*beta*Vss(:,:,t+1);
    
    for ib = 1:nb
        for ik = 1:nk
            
            % Calculate tax terms
            fincome_lab = (1-ss_tax_cred)*ben(ib);
            [fincome, ftax, ~, fcap] ...
                = calculate_tax(fincome_lab, kgrid(ik), ...
                                cap_share, rate_cap, debt_share, rate_gov, cap_tax_share, tau_cap, tau_capgain, exp_subsidy,...
                                avg_deduc, coefs, limit, X, mpci, rpci, 0, 0, clinton, year, q_tobin, q_tobin0);
            
            % Calculate effective income
            income  = rate_tot*kgrid(ik) + ben(ib) + (year == 1)*cap_share*kgrid(ik)*(q_tobin - q_tobin0)/q_tobin0 ;
            eff_inc = income - ftax - fcap + beq;
            
            % Call value function to set parameters
            valfun_retire([], kgrid, eff_inc, EV(:,ib), sigma, gamma);
            
            % Solve dynamic optimization subproblem
            [k_opt, V_min] = fminsearch(@valfun_retire, kgrid(ik), optim_options);
            
            % Record values
            Vss        (ik,ib,t) = -V_min;
            
            koptss     (ik,ib,t) = k_opt;
            
            fedincomess(ik,ib,t) = fincome;
            fitaxss    (ik,ib,t) = ftax;
            fsstaxss   (ik,ib,t) = 0;
            fcaptaxss  (ik,ib,t) = fcap;
            
        end
    end
end


%% Working age dynamic optimization

% (0 iterations if T_past >= Tr -- i.e. if there is no working period)
for t = N_w:-1:1
    
    % Determine age and year, bounded by projection period
    age  = t + T_past;
    year = max(1, min(age + birthyear, Tss));
    
    % Extract annual parameters
    v_ss_max    = taxmax        (year);
    tau_ss      = ss_tax        (year);
    beq         = beqs          (year);
    wage        = wages         (year);
    rate_cap    = rate_caps     (year);
    rate_gov    = rate_govs     (year);
    cap_share   = cap_shares    (year);
    debt_share  = debt_shares   (year);
    cap_inc     = cap_incs(:,    year);
    exp_subsidy = exp_subsidys  (year);
    
    for ib = 1:nb
        for iz = 1:nz
            
            eff_wage = wage * z(iz,age,idem);
            
            % Calculate expected value curve using values for next time step
            if (age == Tr_eff)
                V_step = sum(tr_z(:))*Vss(:,:,1);
            else
                V_step = reshape( sum(bsxfun(@times, tr_z(iz,:), V(:,:,:,t+1)), 2), [nk,nb] );
            end
            EV = (1-surv(age))*repmat(Vbeq, [1,nb]) + surv(age)*beta*V_step;
            
            for ik = 1:nk
                
                % Call value function to set parameters
                valfun_work([], kgrid, bgrid, cap_inc(ik), cap_share, rate_cap, debt_share, rate_gov, cap_tax_share, tau_cap, tau_capgain, exp_subsidy, eff_wage, beq, EV, sigma, gamma, avg_deduc, coefs, limit, X, mpci, rpci, tau_ss, v_ss_max, age, ib, ik, clinton, year, q_tobin, q_tobin0);
                
                % Solve dynamic optimization subproblem
                lab0 = 0.5;
                k0   = max(kgrid(ik), 0.1 * eff_wage * lab0);   % (Assumes taxation will not exceed 90% of labor income)
                
                [x_opt, V_min] = fminsearch(@valfun_work, [k0, lab0], optim_options);
                
                k_opt   = x_opt(1);
                lab_opt = x_opt(2);

                % Calculate tax terms
                fincome_lab = eff_wage * lab_opt;
                [fincome, ftax, sstax, fcap] ...
                    = calculate_tax(fincome_lab, kgrid(ik), ...
                                    cap_share, rate_cap, debt_share, rate_gov, cap_tax_share, tau_cap, tau_capgain, exp_subsidy, ...
                                    avg_deduc, coefs, limit, X, mpci, rpci, tau_ss, v_ss_max, clinton, year, q_tobin, q_tobin0);
                
                % Store values
                V        (ik,iz,ib,t) = -V_min;
                
                kopt     (ik,iz,ib,t) = k_opt;
                labopt   (ik,iz,ib,t) = lab_opt;
                bopt     (ik,iz,ib,t) = lab_to_b(lab_opt, age, bgrid(ib), v_ss_max, eff_wage);
                
                fedincome(ik,iz,ib,t) = fincome;
                fitax    (ik,iz,ib,t) = ftax;
                fsstax   (ik,iz,ib,t) = sstax;
                fcaptax  (ik,iz,ib,t) = fcap;

            end
        end
    end
end

end




function V_tilde = valfun_retire(k_prime, kgrid_, eff_inc_, EV_ib_, sigma_, gamma_)

% Enforce function inlining for C code generation
coder.inline('always');

% Define parameters as persistent variables
persistent initialized
persistent kgrid
persistent eff_inc
persistent EV_ib
persistent sigma
persistent gamma

% Initialize parameters
if isempty(initialized)
    
    kgrid      = 0;
    eff_inc    = 0;
    EV_ib      = 0;
    sigma      = 0;
    gamma      = 0;
    
    initialized = true;
    
end

% Set parameters if provided
if (nargin > 1)
    
    kgrid      = kgrid_;
    eff_inc    = eff_inc_;
    EV_ib      = EV_ib_;
    sigma      = sigma_;
    gamma      = gamma_;
    
    V_tilde = [];
    return
    
end

% Calculate consumption
cons = eff_inc - k_prime;

% Perform bound checks
if (kgrid(1) <= k_prime) && (k_prime <= kgrid(end)) && (0 <= cons)
    
    % Perform linear interpolation
    V_star = interp1(kgrid, EV_ib, k_prime, 'linear');
    
    % Calculate value function
    V_tilde = V_star + (cons^(gamma*(1-sigma)))/(1-sigma);
    
    % Negate for minimization and force to scalar for C code generation
    V_tilde = -V_tilde(1);
    
else
    V_tilde = Inf;
end


end




function V_tilde = valfun_work(x_prime, kgrid_, bgrid_, cap_inc_ik_, cap_share_, rate_cap_, debt_share_, rate_gov_, cap_tax_share_, tau_cap_, tau_capgain_, exp_subsidy_, eff_wage_, beq_, EV_, sigma_, gamma_, avg_deduc_, coefs_, limit_, X_, mpci_, rpci_, tau_ss_, v_ss_max_, age_, ib_, ik_, clinton_, year_, q_tobin_, q_tobin0_)

% Enforce function inlining for C code generation
coder.inline('always');

% Define parameters as persistent variables
persistent initialized
persistent kgrid
persistent bgrid
persistent cap_inc_ik
persistent cap_share
persistent rate_cap
persistent debt_share
persistent rate_gov
persistent cap_tax_share
persistent tau_cap
persistent tau_capgain
persistent exp_subsidy
persistent eff_wage
persistent beq
persistent EV
persistent sigma
persistent gamma
persistent avg_deduc
persistent coefs
persistent limit
persistent X
persistent mpci
persistent rpci
persistent tau_ss
persistent v_ss_max
persistent age
persistent ib
persistent ik
persistent clinton
persistent year
persistent q_tobin
persistent q_tobin0

% Initialize parameters
if isempty(initialized)
    
    kgrid         = 0;
    bgrid         = 0;
    cap_inc_ik    = 0;
    cap_share     = 0;
    rate_cap      = 0;
    debt_share    = 0;
    rate_gov      = 0;
    cap_tax_share = 0;
    tau_cap       = 0;
    tau_capgain   = 0;
    exp_subsidy   = 0;
    eff_wage      = 0;
    beq           = 0;
    EV            = 0;
    sigma         = 0;
    gamma         = 0;
    avg_deduc     = 0;
    coefs         = 0;
    limit         = 0;
    X             = 0;
    mpci          = 0;
    rpci          = 0;
    tau_ss        = 0;
    v_ss_max      = 0;
    age           = 0;
    ib            = 0;
    ik            = 0;
    clinton       = 0;
    year          = 0;
    q_tobin       = 0;
    q_tobin0      = 0;
    
    initialized = true;
    
end

% Set parameters if provided
if (nargin > 1)
    
    kgrid         = kgrid_;
    bgrid         = bgrid_;
    cap_inc_ik    = cap_inc_ik_;
    cap_share     = cap_share_;
    rate_cap      = rate_cap_;
    debt_share    = debt_share_;
    rate_gov      = rate_gov_;
    cap_tax_share = cap_tax_share_;
    tau_cap       = tau_cap_;
    tau_capgain   = tau_capgain_;
    exp_subsidy   = exp_subsidy_;
    eff_wage      = eff_wage_;
    beq           = beq_;
    EV            = EV_;
    sigma         = sigma_;
    gamma         = gamma_;
    avg_deduc     = avg_deduc_;
    coefs         = coefs_;
    limit         = limit_;
    X             = X_;
    mpci          = mpci_;
    rpci          = rpci_;
    tau_ss        = tau_ss_;
    v_ss_max      = v_ss_max_;
    age           = age_;
    ib            = ib_;
    ik            = ik_;
    clinton       = clinton_;
    year          = year_;
    q_tobin       = q_tobin_;
    q_tobin0      = q_tobin0_;
    
    V_tilde = [];
    return
    
end

% Define decision variables and perform bound checks
k_prime   = x_prime(1);
lab_prime = x_prime(2);
b_prime   = lab_to_b(lab_prime, age, bgrid(ib), v_ss_max, eff_wage);

if ~((0 <= lab_prime) && (lab_prime <= 1) ...
     && (kgrid(1) <= k_prime) && (k_prime <= kgrid(end)) ...
     && (bgrid(1) <= b_prime) && (b_prime <= bgrid(end)))
    
    V_tilde = Inf;
    return
    
end

% Calculate tax terms
fincome_lab = eff_wage*lab_prime;
[~, ftax, sstax, fcap] = calculate_tax(fincome_lab, kgrid(ik), cap_share, rate_cap, debt_share, rate_gov, cap_tax_share, tau_cap, tau_capgain, exp_subsidy, avg_deduc, coefs, limit, X, mpci, rpci, tau_ss, v_ss_max, clinton, year, q_tobin, q_tobin0);

% Calculate consumption and perform bound check
income  = cap_inc_ik + eff_wage*lab_prime + (year == 1)*kgrid(ik)*cap_share*( q_tobin - q_tobin0 ) / q_tobin0;
cons    = income - ftax - fcap - sstax - k_prime + beq;

if ~(0 <= cons)
    V_tilde = Inf;
    return
end

% Perform linear interpolation on expected value curve
V_star = interp2(kgrid', bgrid, EV', k_prime, b_prime, 'linear');

% Calculate value function
V_tilde = V_star + (1/(1-sigma))*((cons^gamma)*((1-lab_prime)^(1-gamma)))^(1-sigma);

% Negate for minimization and force to scalar for C code generation
V_tilde = -V_tilde(1);

end




function [b_prime] = lab_to_b(lab_prime, age, bgrid_ib, v_ss_max, eff_wage)

% Enforce function inlining for C code generation
coder.inline('always');

b_prime = (1/age)*(bgrid_ib*(age-1) + min(v_ss_max, eff_wage*lab_prime));

end



