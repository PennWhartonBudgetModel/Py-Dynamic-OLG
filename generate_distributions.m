%%
% Solve the dynamic optimization problem for a cohort then generate distributions and cohort aggregates using the optimal decision values.
% 
%%


function [W, dist, cohort] ...
          ...
            = generate_distributions(...
                beta, gamma, sigma, T_life, T_work, T_model_opt, T_model_dist, startyear, ...
                zs, ztrans, kv, bv, nz, nk, nb, idem, ...
                mpci, rpci, captaxshare, sstaxcredit, surv, taucap, taucapgain, sstaxs, ssincmaxs, ...
                beqs, wages, capshares, debtshares, caprates, govrates, totrates, expsubsidys, qtobin, qtobin0, Vbeq, ...
                deduction_coefs, pit_coefs, ssbenefits, ...
                dist0, mu2_idem, mu3_idem, ...
                W_static, dist_static) %#codegen


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
assert( isa(startyear,      'double') && (size(startyear,       1) == 1         ) && (size(startyear,       2) == 1             ) );

assert( isa(zs,              'double') && (size(zs,               1) <= nz_max    ) && (size(zs,               2) <= T_life_max    ) && (size(zs, 3) == 2) );
assert( isa(ztrans,           'double') && (size(ztrans,            1) <= nz_max    ) && (size(ztrans,            2) <= nz_max        ) );
assert( isa(kv,          'double') && (size(kv,           1) <= nk_max    ) && (size(kv,           2) == 1             ) );
assert( isa(bv,          'double') && (size(bv,           1) <= nb_max    ) && (size(bv,           2) == 1             ) );
assert( isa(nz,             'double') && (size(nz,              1) == 1         ) && (size(nz,              2) == 1             ) );
assert( isa(nk,             'double') && (size(nk,              1) == 1         ) && (size(nk,              2) == 1             ) );
assert( isa(nb,             'double') && (size(nb,              1) == 1         ) && (size(nb,              2) == 1             ) );
assert( isa(idem,           'double') && (size(idem,            1) == 1         ) && (size(idem,            2) == 1             ) );

assert( isa(mpci,           'double') && (size(mpci,            1) == 1         ) && (size(mpci,            2) == 1             ) );
assert( isa(rpci,           'double') && (size(rpci,            1) == 1         ) && (size(rpci,            2) == 1             ) );
assert( isa(captaxshare,  'double') && (size(captaxshare,   1) == 1         ) && (size(captaxshare,   2) == 1             ) );
assert( isa(sstaxcredit,    'double') && (size(sstaxcredit,     1) == 1         ) && (size(sstaxcredit,     2) == 1             ) );
assert( isa(surv,           'double') && (size(surv,            1) == 1         ) && (size(surv,            2) <= T_life_max    ) );
assert( isa(taucap,        'double') && (size(taucap,         1) == 1         ) && (size(taucap,         2) == 1             ) );
assert( isa(taucapgain,    'double') && (size(taucapgain,     1) == 1         ) && (size(taucapgain,     2) == 1             ) );
assert( isa(sstaxs,         'double') && (size(sstaxs,          1) == 1         ) && (size(sstaxs,          2) <= T_model_max   ) );
assert( isa(ssincmaxs,         'double') && (size(ssincmaxs,          1) == 1         ) && (size(ssincmaxs,          2) <= T_model_max   ) );

assert( isa(beqs,           'double') && (size(beqs,            1) == 1         ) && (size(beqs,            2) <= T_model_max   ) );
assert( isa(wages,          'double') && (size(wages,           1) == 1         ) && (size(wages,           2) <= T_model_max   ) );
assert( isa(capshares,     'double') && (size(capshares,      1) == 1         ) && (size(capshares,      2) <= T_model_max   ) );
assert( isa(debtshares,    'double') && (size(debtshares,     1) == 1         ) && (size(debtshares,     2) <= T_model_max   ) );
assert( isa(caprates,      'double') && (size(caprates,       1) == 1         ) && (size(caprates,       2) <= T_model_max   ) );
assert( isa(govrates,      'double') && (size(govrates,       1) == 1         ) && (size(govrates,       2) <= T_model_max   ) );
assert( isa(totrates,      'double') && (size(totrates,       1) == 1         ) && (size(totrates,       2) <= T_model_max   ) );
assert( isa(expsubsidys,   'double') && (size(expsubsidys,    1) == 1         ) && (size(expsubsidys,    2) <= T_model_max   ) );
assert( isa(qtobin,        'double') && (size(qtobin,         1) == 1         ) && (size(qtobin,         2) == 1             ) );
assert( isa(qtobin0,       'double') && (size(qtobin0,        1) == 1         ) && (size(qtobin0,        2) == 1             ) );
assert( isa(Vbeq,           'double') && (size(Vbeq,            1) <= nk_max    ) && (size(Vbeq,            2) == 1             ) );

assert( isa(deduction_coefs,          'double') && (size(deduction_coefs,           1) == 1         ) && (size(deduction_coefs,           2) == 3            ) );
assert( isa(pit_coefs,              'double') && (size(pit_coefs,               1) == 1         ) && (size(pit_coefs,               2) <= 10            ) );
assert( isa(ssbenefits,     'double') && (size(ssbenefits,      1) <= nb_max    ) && (size(ssbenefits,      2) <= T_model_max   ) );

assert( isa(dist0,          'double') && (size(dist0,           1) <= nz_max    ) && (size(dist0,           2) <= nk_max        ) && (size(dist0,         3) <= nb_max) && (size(dist0,         4) <= T_life_max ) );
assert( isa(mu2_idem,       'double') && (size(mu2_idem,        1) == 1         ) && (size(mu2_idem,        2) <= T_life_max    ) );
assert( isa(mu3_idem,       'double') && (size(mu3_idem,        1) == 1         ) && (size(mu3_idem,        2) <= T_life_max    ) );

assert( isa(W_static,  'double') && (size(W_static,   1) <= nz_max    ) && (size(W_static,   2) <= nk_max        ) && (size(W_static, 3) <= nb_max) && (size(W_static, 4) <= T_life_max ) );
assert( isa(dist_static,    'double') && (size(dist_static,     1) <= nz_max    ) && (size(dist_static,     2) <= nk_max        ) && (size(dist_static,   3) <= nb_max) && (size(dist_static,   4) <= T_model_max) );



%% Dynamic optimization

% Define dynamic aggregate generation flag
isdynamic = isempty(dist_static);

% Find number of past years, effective living years, and effective working years
T_past = max(-startyear     , 0);
S_life = max(T_life - T_past, 0);
S_work = max(T_work - T_past, 0);

% Initialize optimal decision value arrays
V           = zeros(nz,nk,nb,S_life+1);  % (1 extra time slice for initialization of backward induction)

K        = zeros(nz,nk,nb,S_life);
W      = zeros(nz,nk,nb,S_life);
B        = zeros(nz,nk,nb,S_life);

INC   = zeros(nz,nk,nb,S_life);
PIT       = zeros(nz,nk,nb,S_life);
SST      = zeros(nz,nk,nb,S_life);
CIT     = zeros(nz,nk,nb,S_life);
BEN    = zeros(nz,nk,nb,S_life);

% Specify settings for dynamic optimization subproblems
optim_options = optimset('Display', 'off', 'TolFun', 1e-4, 'TolX', 1e-4);


% Solve dynamic optimization problem through backward induction
for t = S_life:-1:1
    
    % Determine age and year, bounded by projection period
    age  = t + T_past;
    year = max(1, min(age + startyear, T_model_opt));
    
    % Extract annual parameters
    ssbenefit         = ssbenefits(:,  year);
    ssincmax    = ssincmaxs        (year);
    sstax      = sstaxs        (year);
    beq         = beqs          (year);
    wage        = wages         (year);
    caprate    = caprates     (year);
    govrate    = govrates     (year);
    capshare   = capshares    (year);
    debtshare  = debtshares   (year);
    totrate    = totrates     (year);
    expsubsidy = expsubsidys  (year);
    
    
    for ib = 1:nb
        for ik = 1:nk
            
            if (t > S_work)
                
                % Calculate expected value curve using values for next time step
                EV = (1-surv(age))*repmat(Vbeq, [1,nb]) + surv(age)*beta*reshape(V(1,:,:,t+1), [nk,nb]);
                
                % Calculate available resources and tax terms
                inc_w = ssbenefit(ib);
                [resources, inc, pit, ~, cit] ...
                    ...
                    = calculate_resources(inc_w, kv(ik), ...
                                          capshare, caprate, debtshare, govrate, captaxshare, taucap, taucapgain, expsubsidy,...
                                          deduction_coefs, pit_coefs, mpci, rpci, 0, 0, year, qtobin, qtobin0, ...
                                          totrate, beq, sstaxcredit);
                
                if isdynamic
                    
                    % Call retirement age value function to set parameters
                    value_retirement([], kv, resources, EV(:,ib), sigma, gamma);
                    
                    % Solve dynamic optimization subproblem
                    [k, v] = fminsearch(@value_retirement, kv(ik), optim_options);
                    
                    % Record values
                    V    (:,ik,ib,t) = -v     ;
                    K (:,ik,ib,t) = k      ;
                    
                end
                
                W   (:,ik,ib,t) = 0          ;
                B     (:,ik,ib,t) = bv(ib)  ;
                
                INC(:,ik,ib,t) = inc    ;
                PIT    (:,ik,ib,t) = pit       ;
                SST   (:,ik,ib,t) = 0          ;
                CIT  (:,ik,ib,t) = cit       ;
                BEN (:,ik,ib,t) = inc_w;
                
            else
                
                for iz = 1:nz
                    
                    % Calculate expected value curve using values for next time step
                    EV = (1-surv(age))*repmat(Vbeq, [1,nb]) + surv(age)*beta*reshape(sum(repmat(ztrans(iz,:)', [1,nk,nb]) .* V(:,:,:,t+1), 1), [nk,nb]);
                    
                    % Calculate effective wage
                    wage_eff = wage * zs(iz,age,idem);
                    
                    % Call resource calculation function to set parameters
                    calculate_resources([], kv(ik), ...
                                        capshare, caprate, debtshare, govrate, captaxshare, taucap, taucapgain, expsubsidy, ...
                                        deduction_coefs, pit_coefs, mpci, rpci, sstax, ssincmax, year, qtobin, qtobin0, ...
                                        totrate, beq, 0);
                    
                    if isdynamic

                        % Call working age value function and average earnings calculation function to set parameters
                        value_working([], kv, bv, wage_eff, EV, sigma, gamma)
                        calculate_b  ([], age, bv(ib), ssincmax)
                        
                        % Solve dynamic optimization subproblem
                        w0 = 0.5;
                        k0 = max(kv(ik), 0.1 * wage_eff * w0);   % (Assumes taxation will not exceed 90% of labor income)
                        
                        [x, v] = fminsearch(@value_working, [k0, w0], optim_options);
                        
                        k = x(1);
                        w = x(2);
                        
                        % Record values
                        V (iz,ik,ib,t) = -v ;
                        K (iz,ik,ib,t) = k  ;
                        
                    else
                        w = W_static(iz,ik,ib,t);
                    end
                    
                    % Calculate tax terms for optimal decision values
                    inc_w = wage_eff * w;
                    [~, inc, pit, sst, cit] = calculate_resources(inc_w);
                    
                    W    (iz,ik,ib,t) = w;
                    B      (iz,ik,ib,t) = calculate_b(inc_w);
                    
                    INC (iz,ik,ib,t) = inc;
                    PIT     (iz,ik,ib,t) = pit   ;
                    SST    (iz,ik,ib,t) = sst  ;
                    CIT   (iz,ik,ib,t) = cit   ;
                    BEN  (iz,ik,ib,t) = 0      ;
                    
                end
                
            end
            
        end
    end
end



%% Distribution generation

if isdynamic
    
    % Find number of distribution years
    T_dist = min(startyear+T_life, T_model_dist) - max(startyear, 0);
    
    % Initialize distributions
    dist = zeros(nz,nk,nb,T_dist);
    dist(:,:,:,1) = dist0(:,:,:,T_past+1);
    
    % Find distributions through forward propagation
    for t = 1:T_dist-1
        
        % Extract optimal k and b values
        k_t = K(:,:,:,t);
        b_t = B(:,:,:,t);
        
        % Find indices of nearest values in kv and bv series
        jk_lt = ones(size(k_t));
        for elem = 1:length(k_t(:))
            jk_lt(elem) = find(kv(1:end-1) <= k_t(elem), 1, 'last');
        end
        jk_gt = jk_lt + 1;
        
        jb_lt = ones(size(b_t));
        for elem = 1:length(b_t(:))
            jb_lt(elem) = find(bv(1:end-1) <= b_t(elem), 1, 'last');
        end
        jb_gt = jb_lt + 1;
        
        % Calculate linear weights for nearest values
        wk_lt = (kv(jk_gt) - k_t) ./ (kv(jk_gt) - kv(jk_lt));
        wk_gt = 1 - wk_lt;
        
        wb_lt = (bv(jb_gt) - b_t) ./ (bv(jb_gt) - bv(jb_lt));
        wb_gt = 1 - wb_lt;
        
        for jz = 1:nz
            
            % Perform productivity transformation
            dist_step = repmat(ztrans(:,jz), [1,nk,nb]) .* dist(:,:,:,t);
            
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

cohort.k_alive = sum(reshape(K  (:,:,:,1:T_dist) .* dist, [], T_dist), 1);
cohort.k_dead  = sum(reshape(K  (:,:,:,1:T_dist) .* dist, [], T_dist), 1) .* mu3_idem(T_past+(1:T_dist)) ./ mu2_idem(T_past+(1:T_dist));
cohort.w_eff   = sum(reshape(W  (:,:,:,1:T_dist) .* repmat(reshape(zs(:,T_past+(1:T_dist),idem), [nz,1,1,T_dist]), [1,nk,nb,1]) .* dist, [], T_dist), 1);
cohort.w       = sum(reshape(W  (:,:,:,1:T_dist) .* dist, [], T_dist), 1);
cohort.lfpr    = sum(reshape((W(:,:,:,1:T_dist) >= 0.01) .* dist, [], T_dist), 1);
cohort.inc     = sum(reshape(INC(:,:,:,1:T_dist) .* dist, [], T_dist), 1);
cohort.pit     = sum(reshape(PIT(:,:,:,1:T_dist) .* dist, [], T_dist), 1);
cohort.sst     = sum(reshape(SST(:,:,:,1:T_dist) .* dist, [], T_dist), 1);
cohort.cit     = sum(reshape(CIT(:,:,:,1:T_dist) .* dist, [], T_dist), 1);
cohort.ben     = sum(reshape(BEN(:,:,:,1:T_dist) .* dist, [], T_dist), 1);


end




% Resource and tax calculation function
function [resources, inc, pit, sst, cit] ...
    ...
    = calculate_resources(inc_w, kv_ik_, ...
                          capshare_, caprate_, debtshare_, govrate_, captaxshare_, taucap_, taucapgain_, expsubsidy_, ...
                          deduction_coefs_, pit_coefs_, mpci_, rpci_, sstax_, ssincmax_, year_, qtobin_, qtobin0_, ...
                          totrate_, beq_, sstaxcredit_) %#codegen

% Enforce function inlining for C code generation
coder.inline('always');

% Define parameters as persistent variables
persistent initialized
persistent kv_ik
persistent capshare
persistent caprate
persistent debtshare
persistent govrate
persistent captaxshare
persistent taucap
persistent taucapgain
persistent expsubsidy
persistent deduction_coefs
persistent pit_coefs
persistent mpci
persistent rpci
persistent sstax
persistent ssincmax
persistent year
persistent qtobin
persistent qtobin0
persistent totrate
persistent beq
persistent sstaxcredit

% Initialize parameters
if isempty(initialized)
    
    kv_ik        = 0;
    capshare       = 0;
    caprate        = 0;
    debtshare      = 0;
    govrate        = 0;
    captaxshare   = 0;
    taucap         = 0;
    taucapgain     = 0;
    expsubsidy     = 0;
    deduction_coefs = 0;
    pit_coefs       = 0;
    mpci            = 0;
    rpci            = 0;
    sstax          = 0;
    ssincmax        = 0;
    year            = 0;
    qtobin         = 0;
    qtobin0        = 0;
    totrate        = 0;
    beq             = 0;
    sstaxcredit     = 0;
    
    initialized     = true;
    
end

% Set parameters if provided
if (nargin > 1)
    
    kv_ik        = kv_ik_;
    capshare       = capshare_;
    caprate        = caprate_;
    debtshare      = debtshare_;
    govrate        = govrate_;
    captaxshare   = captaxshare_;
    taucap         = taucap_;
    taucapgain     = taucapgain_;
    expsubsidy     = expsubsidy_;
    deduction_coefs = deduction_coefs_;
    pit_coefs       = pit_coefs_;
    mpci            = mpci_;
    rpci            = rpci_;
    sstax          = sstax_;
    ssincmax        = ssincmax_;
    year            = year_;
    qtobin         = qtobin_;
    qtobin0        = qtobin0_;
    totrate        = totrate_;
    beq             = beq_;
    sstaxcredit     = sstaxcredit_;
    
    if isempty(inc_w), return, end
    
end

% Calculate taxable income
inc       = (rpci/mpci)*max(0, capshare*(caprate-1)*kv_ik*(1-captaxshare) + debtshare*(govrate-1)*kv_ik + (1-sstaxcredit)*inc_w);
deduction = max(0, deduction_coefs*inc.^[0; 1; 0.5]);
inc_eff   = max(inc - deduction, 0);
inc       = (mpci/rpci)*inc;

% Calculate personal income tax (PIT)
pit = (mpci/rpci)*pit_coefs(1)*(inc_eff - (inc_eff.^(-pit_coefs(2)) + (pit_coefs(3))).^(-1/pit_coefs(2)));

% Calculate Social Security tax
sst = sstax*min(inc_w, ssincmax);

% Calculate capital income tax (CIT)
cit = capshare*kv_ik*(taucap*((caprate-1) - expsubsidy)*captaxshare + taucapgain*(year == 1)*(qtobin - qtobin0)/qtobin0);

% Calculate available resources
resources = totrate*kv_ik + inc_w - (pit + sst + cit) + beq + (year == 1)*kv_ik*capshare*(qtobin - qtobin0)/qtobin0;

end




% Retirement age value function
function v = value_retirement(k, kv_, resources_, EV_ib_, sigma_, gamma_)

% Enforce function inlining for C code generation
coder.inline('always');

% Define parameters as persistent variables
persistent initialized
persistent kv
persistent resources
persistent EV_ib
persistent sigma
persistent gamma

% Initialize parameters
if isempty(initialized)
    
    kv           = 0;
    resources       = 0;
    EV_ib           = 0;
    sigma           = 0;
    gamma           = 0;
    
    initialized = true;
    
end

% Set parameters if provided
if (nargin > 1)
    
    kv           = kv_;
    resources       = resources_;
    EV_ib           = EV_ib_;
    sigma           = sigma_;
    gamma           = gamma_;
    
    return
    
end

% Calculate consumption
consumption = resources - k;

% Perform bound checks
if (kv(1) <= k) && (k <= kv(end)) && (0 <= consumption)
    
    % Calculate value
    v = interp1(kv, EV_ib, k, 'linear') + (consumption^(gamma*(1-sigma)))/(1-sigma);
    
    % Negate for minimization and force to scalar for C code generation
    v = -v(1);
    
else
    v = Inf;
end


end




% Working age value function
function v = value_working(x, kv_, bv_, wage_eff_, EV_, sigma_, gamma_)

% Enforce function inlining for C code generation
coder.inline('always');

% Define parameters as persistent variables
persistent initialized
persistent kv
persistent bv
persistent wage_eff
persistent EV
persistent sigma
persistent gamma

% Initialize parameters
if isempty(initialized)
    
    kv           = 0;
    bv           = 0;
    wage_eff        = 0;
    EV              = 0;
    sigma           = 0;
    gamma           = 0;
    
    initialized     = true;
    
end

% Set parameters if provided
if (nargin > 1)
    
    kv           = kv_;
    bv           = bv_;
    wage_eff        = wage_eff_;
    EV              = EV_;
    sigma           = sigma_;
    gamma           = gamma_;
    
    return
    
end

% Define decision variables and perform bound checks
k   = x(1);
w = x(2);

inc_w = wage_eff * w;
b     = calculate_b(inc_w);

if ~((0 <= w) && (w <= 1) ...
     && (kv(1) <= k) && (k <= kv(end)) ...
     && (bv(1) <= b) && (b <= bv(end)))
    
    v = Inf;
    return
    
end

% Calculate available resources
resources = calculate_resources(inc_w);

% Calculate consumption and perform bound check
consumption = resources - k;

if ~(0 <= consumption)
    v = Inf;
    return
end

% Calculate value
v = interp2(kv', bv, EV', k, b, 'linear') + (1/(1-sigma))*((consumption^gamma)*((1-w)^(1-gamma)))^(1-sigma);

% Negate for minimization and force to scalar for C code generation
v = -v(1);

end




% Average earnings calculation function
function [b] = calculate_b(inc_w, age_, bv_ib_, ssincmax_)

% Define parameters as persistent variables
persistent initialized
persistent ssincmax
persistent age
persistent bv_ib

% Initialize parameters
if isempty(initialized)
    
    ssincmax        = 0;
    age             = 0;
    bv_ib        = 0;
    
    initialized     = true;
    
end

% Set parameters if provided
if (nargin > 1)
    
    ssincmax        = ssincmax_;
    age             = age_;
    bv_ib        = bv_ib_;
    
    return
    
end

% Enforce function inlining for C code generation
coder.inline('always');

b = (1/age)*(bv_ib*(age-1) + min(inc_w, ssincmax));

end



