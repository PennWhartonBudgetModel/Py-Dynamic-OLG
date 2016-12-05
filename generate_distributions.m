%%
% Solve the dynamic optimization problem for a cohort then generate distributions and cohort aggregates using the optimal decision values.
% 
%%


function [labopt, dist_w, dist_r, N_w, N_r, Kalive, Kdead, ELab, Lab, Lfpr, Fedincome, Fedit, SSrev, Fcaptax, SSexp] ...
          ...
            = generate_distributions(...
                beta, gamma, sigma, T, Tr, Tss, birthyear, ...
                kgrid, z, tr_z, bgrid, nk, nz, nb, idem, ...
                mpci, rpci, cap_tax_share, ss_tax_cred, surv, tau_cap, tau_capgain, ss_tax, taxmax, ...
                beqs, wages, cap_shares, debt_shares, rate_caps, rate_govs, rate_tots, exp_subsidys, q_tobin, q_tobin0, Vbeq, ...
                avg_deduc, clinton, coefs, limit, X, ss_benefit, ...
                dist_w0, dist_r0, mu2_idem, mu3_idem) %#codegen


% Define argument properties for C code generation
T_max   = 100;
Tss_max = 100;  % (Should be greater than or equal to T_max since Tss = T for steady state)
nk_max  =  50;
nz_max  =  50;
nb_max  =  50;

assert( isa(beta,           'double') && (size(beta,            1) == 1         ) && (size(beta,            2) == 1         ) );
assert( isa(gamma,          'double') && (size(gamma,           1) == 1         ) && (size(gamma,           2) == 1         ) );
assert( isa(sigma,          'double') && (size(sigma,           1) == 1         ) && (size(sigma,           2) == 1         ) );
assert( isa(T,              'double') && (size(T,               1) == 1         ) && (size(T,               2) == 1         ) );
assert( isa(Tr,             'double') && (size(Tr,              1) == 1         ) && (size(Tr,              2) == 1         ) );
assert( isa(Tss,            'double') && (size(Tss,             1) == 1         ) && (size(Tss,             2) == 1         ) );
assert( isa(birthyear,      'double') && (size(birthyear,       1) == 1         ) && (size(birthyear,       2) == 1         ) );

assert( isa(kgrid,          'double') && (size(kgrid,           1) <= nk_max    ) && (size(kgrid,           2) == 1         ) );
assert( isa(z,              'double') && (size(z,               1) <= nz_max    ) && (size(z,               2) <= T_max     ) && (size(z, 3) == 2) );
assert( isa(tr_z,           'double') && (size(tr_z,            1) <= nz_max    ) && (size(tr_z,            2) <= nz_max    ) );
assert( isa(bgrid,          'double') && (size(bgrid,           1) <= nb_max    ) && (size(bgrid,           2) == 1         ) );
assert( isa(nk,             'double') && (size(nk,              1) == 1         ) && (size(nk,              2) == 1         ) );
assert( isa(nz,             'double') && (size(nz,              1) == 1         ) && (size(nz,              2) == 1         ) );
assert( isa(nb,             'double') && (size(nb,              1) == 1         ) && (size(nb,              2) == 1         ) );
assert( isa(idem,           'double') && (size(idem,            1) == 1         ) && (size(idem,            2) == 1         ) );

assert( isa(mpci,           'double') && (size(mpci,            1) == 1         ) && (size(mpci,            2) == 1         ) );
assert( isa(rpci,           'double') && (size(rpci,            1) == 1         ) && (size(rpci,            2) == 1         ) );
assert( isa(cap_tax_share,  'double') && (size(cap_tax_share,   1) == 1         ) && (size(cap_tax_share,   2) == 1         ) );
assert( isa(ss_tax_cred,    'double') && (size(ss_tax_cred,     1) == 1         ) && (size(ss_tax_cred,     2) == 1         ) );
assert( isa(surv,           'double') && (size(surv,            1) == 1         ) && (size(surv,            2) <= T_max     ) );
assert( isa(tau_cap,        'double') && (size(tau_cap,         1) == 1         ) && (size(tau_cap,         2) == 1         ) );
assert( isa(tau_capgain,    'double') && (size(tau_capgain,     1) == 1         ) && (size(tau_capgain,     2) == 1         ) );
assert( isa(ss_tax,         'double') && (size(ss_tax,          1) == 1         ) && (size(ss_tax,          2) <= Tss_max   ) );
assert( isa(taxmax,         'double') && (size(taxmax,          1) == 1         ) && (size(taxmax,          2) <= Tss_max   ) );

assert( isa(beqs,           'double') && (size(beqs,            1) == 1         ) && (size(beqs,            2) <= Tss_max   ) );
assert( isa(wages,          'double') && (size(wages,           1) == 1         ) && (size(wages,           2) <= Tss_max   ) );
assert( isa(cap_shares,     'double') && (size(cap_shares,      1) == 1         ) && (size(cap_shares,      2) <= Tss_max   ) );
assert( isa(debt_shares,    'double') && (size(debt_shares,     1) == 1         ) && (size(debt_shares,     2) <= Tss_max   ) );
assert( isa(rate_caps,      'double') && (size(rate_caps,       1) == 1         ) && (size(rate_caps,       2) <= Tss_max   ) );
assert( isa(rate_govs,      'double') && (size(rate_govs,       1) == 1         ) && (size(rate_govs,       2) <= Tss_max   ) );
assert( isa(rate_tots,      'double') && (size(rate_tots,       1) == 1         ) && (size(rate_tots,       2) <= Tss_max   ) );
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

assert( isa(dist_w0,        'double') && (size(dist_w0,         1) <= nk_max    ) && (size(dist_w0,         2) <= nz_max    ) && (size(dist_w0, 3) <= nb_max) && (size(dist_w0, 4) <= T_max) );
assert( isa(dist_r0,        'double') && (size(dist_r0,         1) <= nk_max    ) && (size(dist_r0,         2) <= nb_max    ) && (size(dist_r0, 3) <= T_max) );
assert( isa(mu2_idem,       'double') && (size(mu2_idem,        1) == 1         ) && (size(mu2_idem,        2) <= T_max     ) );
assert( isa(mu3_idem,       'double') && (size(mu3_idem,        1) == 1         ) && (size(mu3_idem,        2) <= T_max     ) );



%% Dynamic optimization

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


% Solve retirement age dynamic optimization problem
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
            
            % Calculate available resources and tax terms
            fincome_lab = ben(ib);
            [resources, fincome, ftax, ~, fcap] ...
                ...
                = calculate_resources(fincome_lab, kgrid(ik), ...
                                      cap_share, rate_cap, debt_share, rate_gov, cap_tax_share, tau_cap, tau_capgain, exp_subsidy,...
                                      avg_deduc, coefs, limit, X, mpci, rpci, 0, 0, clinton, year, q_tobin, q_tobin0, ...
                                      rate_tot, beq, ss_tax_cred);
            
            % Call value function to set parameters
            valfun_retire([], kgrid, resources, EV(:,ib), sigma, gamma);
            
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


% Solve working age dynamic optimization problem
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
    rate_tot    = rate_tots     (year);
    exp_subsidy = exp_subsidys  (year);
    
    for ib = 1:nb
        for iz = 1:nz
            
            eff_wage = wage * z(iz,age,idem);
            
            % Calculate expected value curve using values for next time step
            if (age == Tr_eff)
                V_step = sum(tr_z(:))*Vss(:,:,1);
            else
                V_step = reshape( sum(repmat(tr_z(iz,:), [nk,1,nb]) .* V(:,:,:,t+1), 2), [nk,nb] );
            end
            EV = (1-surv(age))*repmat(Vbeq, [1,nb]) + surv(age)*beta*V_step;
            
            for ik = 1:nk
                
                % Call average earnings calculation function, resource calculation function, and value function to set parameters
                calculate_b([], age, bgrid(ib), v_ss_max)
                
                calculate_resources([], kgrid(ik), ...
                                    cap_share, rate_cap, debt_share, rate_gov, cap_tax_share, tau_cap, tau_capgain, exp_subsidy, ...
                                    avg_deduc, coefs, limit, X, mpci, rpci, tau_ss, v_ss_max, clinton, year, q_tobin, q_tobin0, ...
                                    rate_tot, beq, 0);
                
                valfun_work([], kgrid, bgrid, eff_wage, EV, sigma, gamma)
                
                % Solve dynamic optimization subproblem
                lab0 = 0.5;
                k0   = max(kgrid(ik), 0.1 * eff_wage * lab0);   % (Assumes taxation will not exceed 90% of labor income)
                
                [x_opt, V_min] = fminsearch(@valfun_work, [k0, lab0], optim_options);
                
                k_opt   = x_opt(1);
                lab_opt = x_opt(2);

                % Calculate tax terms for optimal decision values
                fincome_lab_opt = eff_wage * lab_opt;
                [~, fincome, ftax, sstax, fcap] = calculate_resources(fincome_lab_opt);
                
                % Store values
                V        (ik,iz,ib,t) = -V_min;
                
                kopt     (ik,iz,ib,t) = k_opt;
                labopt   (ik,iz,ib,t) = lab_opt;
                bopt     (ik,iz,ib,t) = calculate_b(fincome_lab_opt);
                
                fedincome(ik,iz,ib,t) = fincome;
                fitax    (ik,iz,ib,t) = ftax;
                fsstax   (ik,iz,ib,t) = sstax;
                fcaptax  (ik,iz,ib,t) = fcap;

            end
        end
    end
end



%% Distribution generation

% Define past years, effective lifespan, and effective retirement age
% (Unlike dynamic optimization, projections truncated at year Tss)
T_past = -min(0, birthyear);
T_eff  = min(T, Tss - birthyear);
Tr_eff = max(min(Tr, T_eff), T_past);

% Find effective numbers of working and retirement years
N_w = Tr_eff - T_past;
N_r = T_eff  - Tr_eff;

% Initialize distributions
dist_w = zeros(nk,nz,nb,N_w);
dist_r = zeros(nk,   nb,N_r);

if (N_w > 0)
    dist_w(:,:,:,1) = dist_w0(:,:,:,T_past   +1);
else
    dist_r(:,  :,1) = dist_r0(:,  :,T_past-Tr+1);
end


% Find distributions for working years
for t = 1:N_w
    
    % Extract optimal k and b values
    % (Note that k and b should already bounded by the ranges of kgrid and bgrid respectively)
    k_t = kopt(:,:,:,t);
    b_t = bopt(:,:,:,t);
    
    % Find indices of nearest values in kgrid and bgrid series
    % (Using arrayfun here leads to nontrivial anonymous function overhead)
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
    
    % Calculate linear weights for lower and upper nearest values
    wk_lt = (kgrid(jk_gt) - k_t) ./ (kgrid(jk_gt) - kgrid(jk_lt));
    wk_gt = 1 - wk_lt;
    
    wb_lt = (bgrid(jb_gt) - b_t) ./ (bgrid(jb_gt) - bgrid(jb_lt));
    wb_gt = 1 - wb_lt;
    
    if t ~= N_w
        
        for jz = 1:nz
            
            % Perform productivity transformation
            dist_step = repmat(tr_z(:,jz)', [nk,1,nb]) .* dist_w(:,:,:,t);
            
            % Calculate distributions for next time step
            for elem = 1:numel(dist_step)
                dist_w(jk_lt(elem), jz, jb_lt(elem), t+1) = dist_w(jk_lt(elem), jz, jb_lt(elem), t+1) + wb_lt(elem) * wk_lt(elem) * dist_step(elem);
                dist_w(jk_gt(elem), jz, jb_lt(elem), t+1) = dist_w(jk_gt(elem), jz, jb_lt(elem), t+1) + wb_lt(elem) * wk_gt(elem) * dist_step(elem);
                dist_w(jk_lt(elem), jz, jb_gt(elem), t+1) = dist_w(jk_lt(elem), jz, jb_gt(elem), t+1) + wb_gt(elem) * wk_lt(elem) * dist_step(elem);
                dist_w(jk_gt(elem), jz, jb_gt(elem), t+1) = dist_w(jk_gt(elem), jz, jb_gt(elem), t+1) + wb_gt(elem) * wk_gt(elem) * dist_step(elem);
            end
            
        end
        
    else
        
        % Perform transition from working years to retirement years
        if (N_r > 0)
            
            dist_step = dist_w(:,:,:,N_w);
            
            for elem = 1:numel(dist_step)
                dist_r(jk_lt(elem), jb_lt(elem), 1) = dist_r(jk_lt(elem), jb_lt(elem), 1) + wb_lt(elem) * wk_lt(elem) * dist_step(elem);
                dist_r(jk_gt(elem), jb_lt(elem), 1) = dist_r(jk_gt(elem), jb_lt(elem), 1) + wb_lt(elem) * wk_gt(elem) * dist_step(elem);
                dist_r(jk_lt(elem), jb_gt(elem), 1) = dist_r(jk_lt(elem), jb_gt(elem), 1) + wb_gt(elem) * wk_lt(elem) * dist_step(elem);
                dist_r(jk_gt(elem), jb_gt(elem), 1) = dist_r(jk_gt(elem), jb_gt(elem), 1) + wb_gt(elem) * wk_gt(elem) * dist_step(elem);
            end
            
        end
        
    end

end


% Find distributions for retirement years
% (0 iterations if N_r <= 1 -- i.e. if there is only 1 year or less of retirement)
for t = 1:N_r-1
    
    % Extract optimal k values
    % (Note that k should already bounded by the range of kgrid by the dynamic optimization solver)
    k_t = koptss(:,:,t);
    
    % Find indices of nearest values in kgrid series
    jk_lt = ones(size(k_t));
    for elem = 1:length(k_t(:))
        jk_lt(elem) = find(kgrid(1:end-1) <= k_t(elem), 1, 'last');
    end
    jk_gt = jk_lt + 1;
    
    % Calculate linear weights for lower and upper nearest values
    wk_lt = (kgrid(jk_gt) - k_t) ./ (kgrid(jk_gt) - kgrid(jk_lt));
    wk_gt = 1 - wk_lt;
    
    % Calculate distributions for next time step
    for ib = 1:nb
        for ik = 1:nk
            dist_r(jk_lt(ik,ib), ib, t+1) = dist_r(jk_lt(ik,ib), ib, t+1) + wk_lt(ik,ib) * dist_r(ik,ib,t);
            dist_r(jk_gt(ik,ib), ib, t+1) = dist_r(jk_gt(ik,ib), ib, t+1) + wk_gt(ik,ib) * dist_r(ik,ib,t);
        end
    end
    
end


% Adjust distributions based on demographics
dist_w = repmat(shiftdim(mu2_idem(T_past+1:Tr_eff), -2), [nk,nz,nb]) .* dist_w;
dist_r = repmat(shiftdim(mu2_idem(Tr_eff+1:T_eff ), -1), [nk,   nb]) .* dist_r;


% Calculate aggregates
Kalive    = [sum(reshape(kopt  (:,:,:,1:N_w) .* dist_w, [], N_w), 1), ...
             sum(reshape(koptss(:,  :,1:N_r) .* dist_r, [], N_r), 1)];

Kdead     = Kalive .* mu3_idem(T_past+1:T_eff) ./ mu2_idem(T_past+1:T_eff);

ELab      = sum(reshape(labopt(:,:,:,1:N_w) .* repmat(reshape(z(:,T_past+1:Tr_eff,idem), [1,nz,1,N_w]), [nk,1,nb,1]) .* dist_w, [], N_w), 1);

Lab       = sum(reshape(labopt(:,:,:,1:N_w) .* dist_w, [], N_w), 1);

Lfpr      = sum(reshape((labopt(:,:,:,1:N_w) >= 0.01) .* dist_w, [], N_w), 1);

Fedincome = [sum(reshape(fedincome  (:,:,:,1:N_w) .* dist_w, [], N_w), 1), ...
             sum(reshape(fedincomess(:,  :,1:N_r) .* dist_r, [], N_r), 1)];

Fedit     = [sum(reshape(fitax      (:,:,:,1:N_w) .* dist_w, [], N_w), 1), ...
             sum(reshape(fitaxss    (:,  :,1:N_r) .* dist_r, [], N_r), 1)];

SSrev     = [sum(reshape(fsstax     (:,:,:,1:N_w) .* dist_w, [], N_w), 1), ...
             sum(reshape(fsstaxss   (:,  :,1:N_r) .* dist_r, [], N_r), 1)];

Fcaptax   = [sum(reshape(fcaptax    (:,:,:,1:N_w) .* dist_w, [], N_w), 1), ...
             sum(reshape(fcaptaxss  (:,  :,1:N_r) .* dist_r, [], N_r), 1)];

SSexp     = sum(reshape(repmat(ss_benefit(:,1)', [nk,1,N_r]) .* dist_r, [], N_r), 1);


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



