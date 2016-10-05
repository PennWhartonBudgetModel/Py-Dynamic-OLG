% Jorge | 2016-09-26
% 
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
% 


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
lab0 = 1e-2;


%% Retirement age dynamic optimization

for t = N_r:-1:1
    
    % Determine age and year, bounded by projection period
    age  = t + Tr_eff;
    year = max(1, min(age + birthyear, Tss));
    
    % Extract annual parameters
    ben         = ss_benefit(:,year);
    beq         = beqs(year);
    rate_cap    = rate_caps(year);
    rate_gov    = rate_govs(year);
    cap_share   = cap_shares(year);
    debt_share  = debt_shares(year);
    rate_tot    = rate_tots(year);
    exp_subsidy = exp_subsidys(year);
    
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
    v_ss_max    = taxmax(year);
    tau_ss      = ss_tax(year);
    beq         = beqs(year);
    wage        = wages(year);
    rate_cap    = rate_caps(year);
    rate_gov    = rate_govs(year);
    cap_share   = cap_shares(year);
    debt_share  = debt_shares(year);
    cap_inc     = cap_incs(:,year);
    exp_subsidy = exp_subsidys(year);
    
    for ib = 1:nb
        for iz = 1:nz
            
            eff_wage = wage*max(z(iz,age,idem), 0);
            
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
                [x_opt, V_min] = fminsearch(@valfun_work, [kgrid(ik), lab0], optim_options);
                
                k_opt   = x_opt(1);
                lab_opt = x_opt(2);
                
                % Calculate tax terms
                fincome_lab = eff_wage*lab_opt;
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



