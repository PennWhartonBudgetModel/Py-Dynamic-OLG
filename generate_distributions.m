%%
% Solve the dynamic optimization problem for a cohort then generate distributions and cohort aggregates using the optimal decision values.
% 
%%


function [labopt, dist, Kalive, Kdead, ELab, Lab, Lfpr, Fedincome, Fedit, SSrev, Fcaptax, SSexp] ...
          ...
            = generate_distributions(...
                beta, gamma, sigma, T_life, T_work, T_model_opt, T_model_dist, birthyear, ...
                z, tr_z, kgrid, bgrid, nz, nk, nb, idem, ...
                mpci, rpci, cap_tax_share, ss_tax_cred, surv, tau_cap, tau_capgain, ss_tax, taxmax, ...
                beqs, wages, cap_shares, debt_shares, rate_caps, rate_govs, rate_tots, exp_subsidys, q_tobin, q_tobin0, Vbeq, ...
                avg_deduc, clinton, coefs, limit, X, ss_benefit, ...
                dist0, mu2_idem, mu3_idem, ...
                labopt_static, dist_static) %#codegen


% Define argument properties for C code generation
T_life_max  = 100;
T_model_max =  75;
nz_max      =  50;
nk_max      =  50;
nb_max      =  50;

assert( isa(beta,           'double') && (size(beta,            1) == 1         ) && (size(beta,            2) == 1             ) );
assert( isa(gamma,          'double') && (size(gamma,           1) == 1         ) && (size(gamma,           2) == 1             ) );
assert( isa(sigma,          'double') && (size(sigma,           1) == 1         ) && (size(sigma,           2) == 1             ) );
assert( isa(T_life,         'double') && (size(T_life,          1) == 1         ) && (size(T_life,          2) == 1             ) );
assert( isa(T_work,         'double') && (size(T_work,          1) == 1         ) && (size(T_work,          2) == 1             ) );
assert( isa(T_model_dist,   'double') && (size(T_model_dist,    1) == 1         ) && (size(T_model_dist,    2) == 1             ) );
assert( isa(T_model_opt,    'double') && (size(T_model_opt,     1) == 1         ) && (size(T_model_opt,     2) == 1             ) );
assert( isa(birthyear,      'double') && (size(birthyear,       1) == 1         ) && (size(birthyear,       2) == 1             ) );

assert( isa(z,              'double') && (size(z,               1) <= nz_max    ) && (size(z,               2) <= T_life_max    ) && (size(z, 3) == 2) );
assert( isa(tr_z,           'double') && (size(tr_z,            1) <= nz_max    ) && (size(tr_z,            2) <= nz_max        ) );
assert( isa(kgrid,          'double') && (size(kgrid,           1) <= nk_max    ) && (size(kgrid,           2) == 1             ) );
assert( isa(bgrid,          'double') && (size(bgrid,           1) <= nb_max    ) && (size(bgrid,           2) == 1             ) );
assert( isa(nz,             'double') && (size(nz,              1) == 1         ) && (size(nz,              2) == 1             ) );
assert( isa(nk,             'double') && (size(nk,              1) == 1         ) && (size(nk,              2) == 1             ) );
assert( isa(nb,             'double') && (size(nb,              1) == 1         ) && (size(nb,              2) == 1             ) );
assert( isa(idem,           'double') && (size(idem,            1) == 1         ) && (size(idem,            2) == 1             ) );

assert( isa(mpci,           'double') && (size(mpci,            1) == 1         ) && (size(mpci,            2) == 1             ) );
assert( isa(rpci,           'double') && (size(rpci,            1) == 1         ) && (size(rpci,            2) == 1             ) );
assert( isa(cap_tax_share,  'double') && (size(cap_tax_share,   1) == 1         ) && (size(cap_tax_share,   2) == 1             ) );
assert( isa(ss_tax_cred,    'double') && (size(ss_tax_cred,     1) == 1         ) && (size(ss_tax_cred,     2) == 1             ) );
assert( isa(surv,           'double') && (size(surv,            1) == 1         ) && (size(surv,            2) <= T_life_max    ) );
assert( isa(tau_cap,        'double') && (size(tau_cap,         1) == 1         ) && (size(tau_cap,         2) == 1             ) );
assert( isa(tau_capgain,    'double') && (size(tau_capgain,     1) == 1         ) && (size(tau_capgain,     2) == 1             ) );
assert( isa(ss_tax,         'double') && (size(ss_tax,          1) == 1         ) && (size(ss_tax,          2) <= T_model_max   ) );
assert( isa(taxmax,         'double') && (size(taxmax,          1) == 1         ) && (size(taxmax,          2) <= T_model_max   ) );

assert( isa(beqs,           'double') && (size(beqs,            1) == 1         ) && (size(beqs,            2) <= T_model_max   ) );
assert( isa(wages,          'double') && (size(wages,           1) == 1         ) && (size(wages,           2) <= T_model_max   ) );
assert( isa(cap_shares,     'double') && (size(cap_shares,      1) == 1         ) && (size(cap_shares,      2) <= T_model_max   ) );
assert( isa(debt_shares,    'double') && (size(debt_shares,     1) == 1         ) && (size(debt_shares,     2) <= T_model_max   ) );
assert( isa(rate_caps,      'double') && (size(rate_caps,       1) == 1         ) && (size(rate_caps,       2) <= T_model_max   ) );
assert( isa(rate_govs,      'double') && (size(rate_govs,       1) == 1         ) && (size(rate_govs,       2) <= T_model_max   ) );
assert( isa(rate_tots,      'double') && (size(rate_tots,       1) == 1         ) && (size(rate_tots,       2) <= T_model_max   ) );
assert( isa(exp_subsidys,   'double') && (size(exp_subsidys,    1) == 1         ) && (size(exp_subsidys,    2) <= T_model_max   ) );
assert( isa(q_tobin,        'double') && (size(q_tobin,         1) == 1         ) && (size(q_tobin,         2) == 1             ) );
assert( isa(q_tobin0,       'double') && (size(q_tobin0,        1) == 1         ) && (size(q_tobin0,        2) == 1             ) );
assert( isa(Vbeq,           'double') && (size(Vbeq,            1) <= nk_max    ) && (size(Vbeq,            2) == 1             ) );

assert( isa(avg_deduc,      'double') && (size(avg_deduc,       1) == 1         ) && (size(avg_deduc,       2) == 1             ) );
assert( isa(clinton,        'double') && (size(clinton,         1) == 1         ) && (size(clinton,         2) == 1             ) );
assert( isa(coefs,          'double') && (size(coefs,           1) == 1         ) && (size(coefs,           2) <= 10            ) );
assert( isa(limit,          'double') && (size(limit,           1) == 1         ) && (size(limit,           2) == 1             ) );
assert( isa(X,              'double') && (size(X,               1) == 1         ) && (size(X,               2) <= 10            ) );
assert( isa(ss_benefit,     'double') && (size(ss_benefit,      1) <= nb_max    ) && (size(ss_benefit,      2) <= T_model_max   ) );

assert( isa(dist0,          'double') && (size(dist0,           1) <= nz_max    ) && (size(dist0,           2) <= nk_max        ) && (size(dist0,         3) <= nb_max) && (size(dist0,         4) <= T_life_max ) );
assert( isa(mu2_idem,       'double') && (size(mu2_idem,        1) == 1         ) && (size(mu2_idem,        2) <= T_life_max    ) );
assert( isa(mu3_idem,       'double') && (size(mu3_idem,        1) == 1         ) && (size(mu3_idem,        2) <= T_life_max    ) );

assert( isa(labopt_static,  'double') && (size(labopt_static,   1) <= nz_max    ) && (size(labopt_static,   2) <= nk_max        ) && (size(labopt_static, 3) <= nb_max) && (size(labopt_static, 4) <= T_life_max ) );
assert( isa(dist_static,    'double') && (size(dist_static,     1) <= nz_max    ) && (size(dist_static,     2) <= nk_max        ) && (size(dist_static,   3) <= nb_max) && (size(dist_static,   4) <= T_model_max) );



%% Dynamic optimization

% Define dynamic aggregate generation flag
isdynamic = isempty(dist_static);

% Find number of past years, effective living years, and effective working years
T_past     = max(-birthyear,      0);
T_life_eff = max(T_life - T_past, 0);
T_work_eff = max(T_work - T_past, 0);

% Initialize optimal decision value arrays
V           = zeros(nz,nk,nb,T_life_eff+1);  % (1 extra time slice for initialization of backward induction)

kopt        = zeros(nz,nk,nb,T_life_eff);
labopt      = zeros(nz,nk,nb,T_life_eff);
bopt        = zeros(nz,nk,nb,T_life_eff);

fedincome   = zeros(nz,nk,nb,T_life_eff);
fitax       = zeros(nz,nk,nb,T_life_eff);
fsstax      = zeros(nz,nk,nb,T_life_eff);
fcaptax     = zeros(nz,nk,nb,T_life_eff);
benefits    = zeros(nz,nk,nb,T_life_eff);

% Specify settings for dynamic optimization subproblems
optim_options = optimset('Display', 'off', 'TolFun', 1e-4, 'TolX', 1e-4);


% Solve dynamic optimization problem through backward induction
for t = T_life_eff:-1:1
    
    % Determine age and year, bounded by projection period
    age  = t + T_past;
    year = max(1, min(age + birthyear, T_model_opt));
    
    % Extract annual parameters
    ben         = ss_benefit(:,  year);
    v_ss_max    = taxmax        (year);
    tau_ss      = ss_tax        (year);
    beq         = beqs          (year);
    wage        = wages         (year);
    rate_cap    = rate_caps     (year);
    rate_gov    = rate_govs     (year);
    cap_share   = cap_shares    (year);
    debt_share  = debt_shares   (year);
    rate_tot    = rate_tots     (year);
    exp_subsidy = exp_subsidys  (year);
    
    
    for ib = 1:nb
        for ik = 1:nk
            
            if (t > T_work_eff)
                
                % Calculate expected value curve using values for next time step
                EV = (1-surv(age))*repmat(Vbeq, [1,nb]) + surv(age)*beta*reshape(V(1,:,:,t+1), [nk,nb]);
                
                % Calculate available resources and tax terms
                fincome_lab = ben(ib);
                [resources, fincome, ftax, ~, fcap] ...
                    ...
                    = calculate_resources(fincome_lab, kgrid(ik), ...
                                          cap_share, rate_cap, debt_share, rate_gov, cap_tax_share, tau_cap, tau_capgain, exp_subsidy,...
                                          avg_deduc, coefs, limit, X, mpci, rpci, 0, 0, clinton, year, q_tobin, q_tobin0, ...
                                          rate_tot, beq, ss_tax_cred);
                
                if isdynamic
                    
                    % Call retirement age value function to set parameters
                    valfun_retire([], kgrid, resources, EV(:,ib), sigma, gamma);
                    
                    % Solve dynamic optimization subproblem
                    [k_opt, V_min] = fminsearch(@valfun_retire, kgrid(ik), optim_options);
                    
                    % Record values
                    V    (:,ik,ib,t) = -V_min     ;
                    kopt (:,ik,ib,t) = k_opt      ;
                    
                end
                
                labopt   (:,ik,ib,t) = 0          ;
                bopt     (:,ik,ib,t) = bgrid(ib)  ;
                
                fedincome(:,ik,ib,t) = fincome    ;
                fitax    (:,ik,ib,t) = ftax       ;
                fsstax   (:,ik,ib,t) = 0          ;
                fcaptax  (:,ik,ib,t) = fcap       ;
                benefits (:,ik,ib,t) = fincome_lab;
                
            else
                
                for iz = 1:nz
                    
                    % Calculate expected value curve using values for next time step
                    EV = (1-surv(age))*repmat(Vbeq, [1,nb]) + surv(age)*beta*reshape(sum(repmat(tr_z(iz,:)', [1,nk,nb]) .* V(:,:,:,t+1), 1), [nk,nb]);
                    
                    % Calculate effective wage
                    eff_wage = wage * z(iz,age,idem);
                    
                    % Call resource calculation function to set parameters
                    calculate_resources([], kgrid(ik), ...
                                        cap_share, rate_cap, debt_share, rate_gov, cap_tax_share, tau_cap, tau_capgain, exp_subsidy, ...
                                        avg_deduc, coefs, limit, X, mpci, rpci, tau_ss, v_ss_max, clinton, year, q_tobin, q_tobin0, ...
                                        rate_tot, beq, 0);
                    
                    if isdynamic

                        % Call working age value function and average earnings calculation function to set parameters
                        valfun_work([], kgrid, bgrid, eff_wage, EV, sigma, gamma)
                        calculate_b([], age, bgrid(ib), v_ss_max)
                        
                        % Solve dynamic optimization subproblem
                        lab0 = 0.5;
                        k0   = max(kgrid(ik), 0.1 * eff_wage * lab0);   % (Assumes taxation will not exceed 90% of labor income)
                        
                        [x_opt, V_min] = fminsearch(@valfun_work, [k0, lab0], optim_options);
                        
                        k_opt   = x_opt(1);
                        lab_opt = x_opt(2);
                        
                        % Record values
                        V     (iz,ik,ib,t) = -V_min ;
                        kopt  (iz,ik,ib,t) = k_opt  ;
                        
                    else
                        lab_opt = labopt_static(iz,ik,ib,t);
                    end
                    
                    % Calculate tax terms for optimal decision values
                    fincome_lab_opt = eff_wage * lab_opt;
                    [~, fincome, ftax, sstax, fcap] = calculate_resources(fincome_lab_opt);
                    
                    labopt    (iz,ik,ib,t) = lab_opt;
                    bopt      (iz,ik,ib,t) = calculate_b(fincome_lab_opt);
                    
                    fedincome (iz,ik,ib,t) = fincome;
                    fitax     (iz,ik,ib,t) = ftax   ;
                    fsstax    (iz,ik,ib,t) = sstax  ;
                    fcaptax   (iz,ik,ib,t) = fcap   ;
                    benefits  (iz,ik,ib,t) = 0      ;
                    
                end
                
            end
            
        end
    end
end



%% Distribution generation

if isdynamic
    
    % Find number of distribution years
    T_dist = min(birthyear+T_life, T_model_dist) - max(birthyear, 0);
    
    % Initialize distributions
    dist = zeros(nz,nk,nb,T_dist);
    dist(:,:,:,1) = dist0(:,:,:,T_past+1);
    
    % Find distributions through forward propagation
    for t = 1:T_dist-1
        
        % Extract optimal k and b values
        k_t = kopt(:,:,:,t);
        b_t = bopt(:,:,:,t);
        
        % Find indices of nearest values in kgrid and bgrid series
        jk_lt = ones(size(k_t));
        for elem = 1:length(k_t(:))
            jk_lt(elem) = find(kgrid(1:end-1) <= k_t(elem), 1, 'last');
        end
        jk_gt = jk_lt + 1;
        
        jb_lt = ones(size(b_t));
        for elem = 1:length(b_t(:))
            jb_lt(elem) = find(bgrid(1:end-1) <= b_t(elem), 1, 'last');
        end
        jb_gt = jb_lt + 1;
        
        % Calculate linear weights for nearest values
        wk_lt = (kgrid(jk_gt) - k_t) ./ (kgrid(jk_gt) - kgrid(jk_lt));
        wk_gt = 1 - wk_lt;
        
        wb_lt = (bgrid(jb_gt) - b_t) ./ (bgrid(jb_gt) - bgrid(jb_lt));
        wb_gt = 1 - wb_lt;
        
        for jz = 1:nz
            
            % Perform productivity transformation
            dist_step = repmat(tr_z(:,jz), [1,nk,nb]) .* dist(:,:,:,t);
            
            % Calculate distributions for next time step
            for elem = 1:numel(dist_step)
                dist(jz, jk_lt(elem), jb_lt(elem), t+1) = dist(jz, jk_lt(elem), jb_lt(elem), t+1) + wb_lt(elem) * wk_lt(elem) * dist_step(elem);
                dist(jz, jk_gt(elem), jb_lt(elem), t+1) = dist(jz, jk_gt(elem), jb_lt(elem), t+1) + wb_lt(elem) * wk_gt(elem) * dist_step(elem);
                dist(jz, jk_lt(elem), jb_gt(elem), t+1) = dist(jz, jk_lt(elem), jb_gt(elem), t+1) + wb_gt(elem) * wk_lt(elem) * dist_step(elem);
                dist(jz, jk_gt(elem), jb_gt(elem), t+1) = dist(jz, jk_gt(elem), jb_gt(elem), t+1) + wb_gt(elem) * wk_gt(elem) * dist_step(elem);
            end
        
        end
        
    end
    
    % Adjust distributions based on demographics
    dist = repmat(shiftdim(mu2_idem(T_past+(1:T_dist)), -2), [nz,nk,nb]) .* dist;
    
else
    
    dist = dist_static;
    T_dist = size(dist, 4);
    
end



%% Aggregate generation

Kalive    = sum(reshape(kopt            (:,:,:,1:T_dist) .* dist, [], T_dist), 1);

Kdead     = Kalive .* mu3_idem(T_past+(1:T_dist)) ./ mu2_idem(T_past+(1:T_dist));

ELab      = sum(reshape(labopt(:,:,:,1:T_dist) .* repmat(reshape(z(:,T_past+(1:T_dist),idem), [nz,1,1,T_dist]), [1,nk,nb,1]) .* dist, [], T_dist), 1);

Lab       = sum(reshape(labopt          (:,:,:,1:T_dist) .* dist, [], T_dist), 1);

Lfpr      = sum(reshape((labopt(:,:,:,1:T_dist) >= 0.01) .* dist, [], T_dist), 1);

Fedincome = sum(reshape(fedincome       (:,:,:,1:T_dist) .* dist, [], T_dist), 1);

Fedit     = sum(reshape(fitax           (:,:,:,1:T_dist) .* dist, [], T_dist), 1);

SSrev     = sum(reshape(fsstax          (:,:,:,1:T_dist) .* dist, [], T_dist), 1);

Fcaptax   = sum(reshape(fcaptax         (:,:,:,1:T_dist) .* dist, [], T_dist), 1);

SSexp     = sum(reshape(benefits        (:,:,:,1:T_dist) .* dist, [], T_dist), 1);


end




% Resource and tax calculation function
function [resources, fincome, ftax, sstax, fcap] ...
    ...
    = calculate_resources(fincome_lab, kgrid_ik_, ...
                          cap_share_, rate_cap_, debt_share_, rate_gov_, cap_tax_share_, tau_cap_, tau_capgain_, exp_subsidy_, ...
                          avg_deduc_, coefs_, limit_, X_, mpci_, rpci_, tau_ss_, v_ss_max_, clinton_, year_, q_tobin_, q_tobin0_, ...
                          rate_tot_, beq_, ss_tax_cred_) %#codegen

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
persistent rate_tot
persistent beq
persistent ss_tax_cred

% Initialize parameters
if isempty(initialized)
    
    kgrid_ik        = 0;
    cap_share       = 0;
    rate_cap        = 0;
    debt_share      = 0;
    rate_gov        = 0;
    cap_tax_share   = 0;
    tau_cap         = 0;
    tau_capgain     = 0;
    exp_subsidy     = 0;
    avg_deduc       = 0;
    coefs           = 0;
    limit           = 0;
    X               = 0;
    mpci            = 0;
    rpci            = 0;
    tau_ss          = 0;
    v_ss_max        = 0;
    clinton         = 0;
    year            = 0;
    q_tobin         = 0;
    q_tobin0        = 0;
    rate_tot        = 0;
    beq             = 0;
    ss_tax_cred     = 0;
    
    initialized     = true;
    
end

% Set parameters if provided
if (nargin > 1)
    
    kgrid_ik        = kgrid_ik_;
    cap_share       = cap_share_;
    rate_cap        = rate_cap_;
    debt_share      = debt_share_;
    rate_gov        = rate_gov_;
    cap_tax_share   = cap_tax_share_;
    tau_cap         = tau_cap_;
    tau_capgain     = tau_capgain_;
    exp_subsidy     = exp_subsidy_;
    avg_deduc       = avg_deduc_;
    coefs           = coefs_;
    limit           = limit_;
    X               = X_;
    mpci            = mpci_;
    rpci            = rpci_;
    tau_ss          = tau_ss_;
    v_ss_max        = v_ss_max_;
    clinton         = clinton_;
    year            = year_;
    q_tobin         = q_tobin_;
    q_tobin0        = q_tobin0_;
    rate_tot        = rate_tot_;
    beq             = beq_;
    ss_tax_cred     = ss_tax_cred_;
    
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
resources = rate_tot*kgrid_ik + fincome_lab - (ftax + sstax + fcap) + beq ...
            + (year == 1)*kgrid_ik*cap_share*(q_tobin - q_tobin0)/q_tobin0;

end




% Retirement age value function
function V_tilde = valfun_retire(k_prime, kgrid_, resources_, EV_ib_, sigma_, gamma_)

% Enforce function inlining for C code generation
coder.inline('always');

% Define parameters as persistent variables
persistent initialized
persistent kgrid
persistent resources
persistent EV_ib
persistent sigma
persistent gamma

% Initialize parameters
if isempty(initialized)
    
    kgrid           = 0;
    resources       = 0;
    EV_ib           = 0;
    sigma           = 0;
    gamma           = 0;
    
    initialized = true;
    
end

% Set parameters if provided
if (nargin > 1)
    
    kgrid           = kgrid_;
    resources       = resources_;
    EV_ib           = EV_ib_;
    sigma           = sigma_;
    gamma           = gamma_;
    
    return
    
end

% Calculate consumption
consumption = resources - k_prime;

% Perform bound checks
if (kgrid(1) <= k_prime) && (k_prime <= kgrid(end)) && (0 <= consumption)
    
    % Perform linear interpolation
    V_star = interp1(kgrid, EV_ib, k_prime, 'linear');
    
    % Calculate value function
    V_tilde = V_star + (consumption^(gamma*(1-sigma)))/(1-sigma);
    
    % Negate for minimization and force to scalar for C code generation
    V_tilde = -V_tilde(1);
    
else
    V_tilde = Inf;
end


end




% Working age value function
function V_tilde = valfun_work(x_prime, kgrid_, bgrid_, eff_wage_, EV_, sigma_, gamma_)

% Enforce function inlining for C code generation
coder.inline('always');

% Define parameters as persistent variables
persistent initialized
persistent kgrid
persistent bgrid
persistent eff_wage
persistent EV
persistent sigma
persistent gamma

% Initialize parameters
if isempty(initialized)
    
    kgrid           = 0;
    bgrid           = 0;
    eff_wage        = 0;
    EV              = 0;
    sigma           = 0;
    gamma           = 0;
    
    initialized     = true;
    
end

% Set parameters if provided
if (nargin > 1)
    
    kgrid           = kgrid_;
    bgrid           = bgrid_;
    eff_wage        = eff_wage_;
    EV              = EV_;
    sigma           = sigma_;
    gamma           = gamma_;
    
    return
    
end

% Define decision variables and perform bound checks
k_prime   = x_prime(1);
lab_prime = x_prime(2);

fincome_lab = eff_wage * lab_prime;
b_prime     = calculate_b(fincome_lab);

if ~((0 <= lab_prime) && (lab_prime <= 1) ...
     && (kgrid(1) <= k_prime) && (k_prime <= kgrid(end)) ...
     && (bgrid(1) <= b_prime) && (b_prime <= bgrid(end)))
    
    V_tilde = Inf;
    return
    
end

% Calculate available resources
resources = calculate_resources(fincome_lab);

% Calculate consumption and perform bound check
consumption = resources - k_prime;

if ~(0 <= consumption)
    V_tilde = Inf;
    return
end

% Perform linear interpolation on expected value curve
V_star = interp2(kgrid', bgrid, EV', k_prime, b_prime, 'linear');

% Calculate value function
V_tilde = V_star + (1/(1-sigma))*((consumption^gamma)*((1-lab_prime)^(1-gamma)))^(1-sigma);

% Negate for minimization and force to scalar for C code generation
V_tilde = -V_tilde(1);

end




% Average earnings calculation function
function [b_prime] = calculate_b(fincome_lab, age_, bgrid_ib_, v_ss_max_)

% Define parameters as persistent variables
persistent initialized
persistent v_ss_max
persistent age
persistent bgrid_ib

% Initialize parameters
if isempty(initialized)
    
    v_ss_max        = 0;
    age             = 0;
    bgrid_ib        = 0;
    
    initialized     = true;
    
end

% Set parameters if provided
if (nargin > 1)
    
    v_ss_max        = v_ss_max_;
    age             = age_;
    bgrid_ib        = bgrid_ib_;
    
    return
    
end

% Enforce function inlining for C code generation
coder.inline('always');

b_prime = (1/age)*(bgrid_ib*(age-1) + min(fincome_lab, v_ss_max));

end



