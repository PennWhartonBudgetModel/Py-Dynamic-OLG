%%
% Solve the dynamic optimization problem for a cohort, generating population distributions and aggregates using the optimal decision values.
% 
%%


function [LAB, DIST, Cohort] = solve_cohort(...
             startyear, T_life, T_work, T_model_opt, T_model_dist, nz, nk, nb, idem, zs, transz, ks, bs, beta, gamma, sigma, surv, V_beq, mu2_idem, mu3_idem, ...
             mpci, rpci, sstaxcredit, ssbenefits, sstaxs, ssincmaxs, deduc_coefs, pit_coefs, captaxshare, taucap, taucapgain, qtobin, qtobin0, ...
             beqs, wages, capshares, debtshares, caprates, govrates, totrates, expsubs, ...
             DIST0, LAB_static, DIST_static) %#codegen


% Define argument properties for C code generation
T_life_max  = 100;
T_model_max =  75;
nz_max      =  50;
nk_max      =  50;
nb_max      =  50;

assert( isa(startyear,      'double') && (size(startyear,       1) == 1         ) && (size(startyear,       2) == 1             ) );
assert( isa(T_life,         'double') && (size(T_life,          1) == 1         ) && (size(T_life,          2) == 1             ) );
assert( isa(T_work,         'double') && (size(T_work,          1) == 1         ) && (size(T_work,          2) == 1             ) );
assert( isa(T_model_dist,   'double') && (size(T_model_dist,    1) == 1         ) && (size(T_model_dist,    2) == 1             ) );
assert( isa(T_model_opt,    'double') && (size(T_model_opt,     1) == 1         ) && (size(T_model_opt,     2) == 1             ) );
assert( isa(nz,             'double') && (size(nz,              1) == 1         ) && (size(nz,              2) == 1             ) );
assert( isa(nk,             'double') && (size(nk,              1) == 1         ) && (size(nk,              2) == 1             ) );
assert( isa(nb,             'double') && (size(nb,              1) == 1         ) && (size(nb,              2) == 1             ) );
assert( isa(idem,           'double') && (size(idem,            1) == 1         ) && (size(idem,            2) == 1             ) );
assert( isa(zs,             'double') && (size(zs,              1) <= nz_max    ) && (size(zs,              2) <= T_life_max    ) && (size(zs, 3) == 2) );
assert( isa(transz,         'double') && (size(transz,          1) <= nz_max    ) && (size(transz,          2) <= nz_max        ) );
assert( isa(ks,             'double') && (size(ks,              1) <= nk_max    ) && (size(ks,              2) == 1             ) );
assert( isa(bs,             'double') && (size(bs,              1) <= nb_max    ) && (size(bs,              2) == 1             ) );
assert( isa(beta,           'double') && (size(beta,            1) == 1         ) && (size(beta,            2) == 1             ) );
assert( isa(gamma,          'double') && (size(gamma,           1) == 1         ) && (size(gamma,           2) == 1             ) );
assert( isa(sigma,          'double') && (size(sigma,           1) == 1         ) && (size(sigma,           2) == 1             ) );
assert( isa(surv,           'double') && (size(surv,            1) == 1         ) && (size(surv,            2) <= T_life_max    ) );
assert( isa(V_beq,          'double') && (size(V_beq,           1) <= nk_max    ) && (size(V_beq,           2) == 1             ) );
assert( isa(mu2_idem,       'double') && (size(mu2_idem,        1) == 1         ) && (size(mu2_idem,        2) <= T_life_max    ) );
assert( isa(mu3_idem,       'double') && (size(mu3_idem,        1) == 1         ) && (size(mu3_idem,        2) <= T_life_max    ) );

assert( isa(mpci,           'double') && (size(mpci,            1) == 1         ) && (size(mpci,            2) == 1             ) );
assert( isa(rpci,           'double') && (size(rpci,            1) == 1         ) && (size(rpci,            2) == 1             ) );
assert( isa(sstaxcredit,    'double') && (size(sstaxcredit,     1) == 1         ) && (size(sstaxcredit,     2) == 1             ) );
assert( isa(ssbenefits,     'double') && (size(ssbenefits,      1) <= nb_max    ) && (size(ssbenefits,      2) <= T_model_max   ) );
assert( isa(sstaxs,         'double') && (size(sstaxs,          1) == 1         ) && (size(sstaxs,          2) <= T_model_max   ) );
assert( isa(ssincmaxs,      'double') && (size(ssincmaxs,       1) == 1         ) && (size(ssincmaxs,       2) <= T_model_max   ) );
assert( isa(deduc_coefs,    'double') && (size(deduc_coefs,     1) == 1         ) && (size(deduc_coefs,     2) == 3             ) );
assert( isa(pit_coefs,      'double') && (size(pit_coefs,       1) == 1         ) && (size(pit_coefs,       2) == 3             ) );
assert( isa(captaxshare,    'double') && (size(captaxshare,     1) == 1         ) && (size(captaxshare,     2) == 1             ) );
assert( isa(taucap,         'double') && (size(taucap,          1) == 1         ) && (size(taucap,          2) == 1             ) );
assert( isa(taucapgain,     'double') && (size(taucapgain,      1) == 1         ) && (size(taucapgain,      2) == 1             ) );
assert( isa(qtobin,         'double') && (size(qtobin,          1) == 1         ) && (size(qtobin,          2) == 1             ) );
assert( isa(qtobin0,        'double') && (size(qtobin0,         1) == 1         ) && (size(qtobin0,         2) == 1             ) );

assert( isa(beqs,           'double') && (size(beqs,            1) == 1         ) && (size(beqs,            2) <= T_model_max   ) );
assert( isa(wages,          'double') && (size(wages,           1) == 1         ) && (size(wages,           2) <= T_model_max   ) );
assert( isa(capshares,      'double') && (size(capshares,       1) == 1         ) && (size(capshares,       2) <= T_model_max   ) );
assert( isa(debtshares,     'double') && (size(debtshares,      1) == 1         ) && (size(debtshares,      2) <= T_model_max   ) );
assert( isa(caprates,       'double') && (size(caprates,        1) == 1         ) && (size(caprates,        2) <= T_model_max   ) );
assert( isa(govrates,       'double') && (size(govrates,        1) == 1         ) && (size(govrates,        2) <= T_model_max   ) );
assert( isa(totrates,       'double') && (size(totrates,        1) == 1         ) && (size(totrates,        2) <= T_model_max   ) );
assert( isa(expsubs,        'double') && (size(expsubs,         1) == 1         ) && (size(expsubs,         2) <= T_model_max   ) );

assert( isa(DIST0,          'double') && (size(DIST0,           1) <= nz_max    ) && (size(DIST0,           2) <= nk_max        ) && (size(DIST0,       3) <= nb_max) && (size(DIST0,       4) <= T_life_max ) );
assert( isa(LAB_static,     'double') && (size(LAB_static,      1) <= nz_max    ) && (size(LAB_static,      2) <= nk_max        ) && (size(LAB_static,  3) <= nb_max) && (size(LAB_static,  4) <= T_life_max ) );
assert( isa(DIST_static,    'double') && (size(DIST_static,     1) <= nz_max    ) && (size(DIST_static,     2) <= nk_max        ) && (size(DIST_static, 3) <= nb_max) && (size(DIST_static, 4) <= T_model_max) );



%% Dynamic optimization

% Define dynamic aggregate generation flag
isdynamic = isempty(DIST_static);

% Find number of past years, effective living years, and effective working years
T_past = max(-startyear,      0);
S_life = max(T_life - T_past, 0);
S_work = max(T_work - T_past, 0);

% Initialize optimal decision value arrays
V   = zeros(nz,nk,nb,S_life+1); % Utility

K   = zeros(nz,nk,nb,S_life);   % Savings
LAB = zeros(nz,nk,nb,S_life);   % Labor level
B   = zeros(nz,nk,nb,S_life);   % Average earnings

INC = zeros(nz,nk,nb,S_life);   % Taxable income
PIT = zeros(nz,nk,nb,S_life);   % Personal income tax
SST = zeros(nz,nk,nb,S_life);   % Social Security tax
CIT = zeros(nz,nk,nb,S_life);   % Corporate income tax
BEN = zeros(nz,nk,nb,S_life);   % Social Security benefits

% Specify settings for dynamic optimization subproblems
optim_options = optimset('Display', 'off', 'TolFun', 1e-4, 'TolX', 1e-4);


% Solve dynamic optimization problem through backward induction
for t = S_life:-1:1
    
    % Determine age and year, bounded by projection period
    age  = t + T_past;
    year = max(1, min(age + startyear, T_model_opt));
    
    % Extract annual parameters
    ssbenefit  = ssbenefits(:, year);
    ssincmax   = ssincmaxs    (year);
    sstax      = sstaxs       (year);
    beq        = beqs         (year);
    wage       = wages        (year);
    caprate    = caprates     (year);
    govrate    = govrates     (year);
    capshare   = capshares    (year);
    debtshare  = debtshares   (year);
    totrate    = totrates     (year);
    expsub     = expsubs      (year);
    
    
    for ib = 1:nb
        for ik = 1:nk
            
            if (t > S_work)
                
                % Calculate expected value curve using values for next time step
                EV = (1-surv(age))*repmat(V_beq, [1,nb]) + surv(age)*beta*reshape(V(1,:,:,t+1), [nk,nb]);
                
                % Calculate available resources and tax terms
                labinc = ssbenefit(ib);
                [resources, inc, pit, ~, cit] = calculate_resources(...
                    labinc, ks(ik), year, ...
                    mpci, rpci, sstaxcredit, 0, 0, deduc_coefs, pit_coefs, captaxshare, taucap, taucapgain, qtobin, qtobin0, ...
                    beq, capshare, debtshare, caprate, govrate, totrate, expsub);
                
                if isdynamic
                    
                    % Call retirement age value function to set parameters
                    value_retirement([], ks, resources, EV(:,ib), sigma, gamma);
                    
                    % Solve dynamic optimization subproblem
                    [k, v] = fminsearch(@value_retirement, ks(ik), optim_options);
                    
                    % Record values
                    V(:,ik,ib,t) = -v;
                    K(:,ik,ib,t) = k ;
                    
                end
                
                LAB(:,ik,ib,t) = 0     ;
                B  (:,ik,ib,t) = bs(ib);
                
                INC(:,ik,ib,t) = inc   ;
                PIT(:,ik,ib,t) = pit   ;
                SST(:,ik,ib,t) = 0     ;
                CIT(:,ik,ib,t) = cit   ;
                BEN(:,ik,ib,t) = labinc;
                
            else
                
                for iz = 1:nz
                    
                    % Calculate expected value curve using values for next time step
                    EV = (1-surv(age))*repmat(V_beq, [1,nb]) + surv(age)*beta*reshape(sum(repmat(transz(iz,:)', [1,nk,nb]) .* V(:,:,:,t+1), 1), [nk,nb]);
                    
                    % Calculate effective wage
                    wage_eff = wage * zs(iz,age,idem);
                    
                    % Call resource calculation function to set parameters
                    calculate_resources(...
                        [], ks(ik), year, ...
                        mpci, rpci, 0, sstax, ssincmax, deduc_coefs, pit_coefs, captaxshare, taucap, taucapgain, qtobin, qtobin0, ...
                        beq, capshare, debtshare, caprate, govrate, totrate, expsub);
                    
                    if isdynamic

                        % Call working age value function and average earnings calculation function to set parameters
                        value_working([], ks, bs, wage_eff, EV, sigma, gamma)
                        calculate_b  ([], age, bs(ib), ssincmax)
                        
                        % Solve dynamic optimization subproblem
                        lab0 = 0.5;
                        k0   = max(ks(ik), 0.1 * wage_eff * lab0);   % (Assumes taxation will not exceed 90% of labor income)
                        
                        [x, v] = fminsearch(@value_working, [k0, lab0], optim_options);
                        
                        k   = x(1);
                        lab = x(2);
                        
                        % Record values
                        V(iz,ik,ib,t) = -v;
                        K(iz,ik,ib,t) = k ;
                        
                    else
                        lab = LAB_static(iz,ik,ib,t);
                    end
                    
                    % Calculate tax terms for optimal decision values
                    labinc = wage_eff * lab;
                    [~, inc, pit, sst, cit] = calculate_resources(labinc);
                    
                    LAB(iz,ik,ib,t) = lab;
                    B  (iz,ik,ib,t) = calculate_b(labinc);
                    
                    INC(iz,ik,ib,t) = inc;
                    PIT(iz,ik,ib,t) = pit;
                    SST(iz,ik,ib,t) = sst;
                    CIT(iz,ik,ib,t) = cit;
                    BEN(iz,ik,ib,t) = 0  ;
                    
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
    DIST = zeros(nz,nk,nb,T_dist);
    DIST(:,:,:,1) = DIST0(:,:,:,T_past+1);
    
    % Find distributions through forward propagation
    for t = 1:T_dist-1
        
        % Extract optimal k and b values
        k_t = K(:,:,:,t);
        b_t = B(:,:,:,t);
        
        % Find indices of nearest values in ks and bs series
        jk_lt = ones(size(k_t));
        for elem = 1:length(k_t(:))
            jk_lt(elem) = find(ks(1:end-1) <= k_t(elem), 1, 'last');
        end
        jk_gt = jk_lt + 1;
        
        jb_lt = ones(size(b_t));
        for elem = 1:length(b_t(:))
            jb_lt(elem) = find(bs(1:end-1) <= b_t(elem), 1, 'last');
        end
        jb_gt = jb_lt + 1;
        
        % Calculate linear weights for nearest values
        wk_lt = (ks(jk_gt) - k_t) ./ (ks(jk_gt) - ks(jk_lt));
        wk_gt = 1 - wk_lt;
        
        wb_lt = (bs(jb_gt) - b_t) ./ (bs(jb_gt) - bs(jb_lt));
        wb_gt = 1 - wb_lt;
        
        for jz = 1:nz
            
            % Perform productivity transformation
            DIST_step = repmat(transz(:,jz), [1,nk,nb]) .* DIST(:,:,:,t);
            
            % Calculate distributions for next time step
            for elem = 1:numel(DIST_step)
                DIST(jz, jk_lt(elem), jb_lt(elem), t+1) = DIST(jz, jk_lt(elem), jb_lt(elem), t+1) + wb_lt(elem) * wk_lt(elem) * DIST_step(elem);
                DIST(jz, jk_gt(elem), jb_lt(elem), t+1) = DIST(jz, jk_gt(elem), jb_lt(elem), t+1) + wb_lt(elem) * wk_gt(elem) * DIST_step(elem);
                DIST(jz, jk_lt(elem), jb_gt(elem), t+1) = DIST(jz, jk_lt(elem), jb_gt(elem), t+1) + wb_gt(elem) * wk_lt(elem) * DIST_step(elem);
                DIST(jz, jk_gt(elem), jb_gt(elem), t+1) = DIST(jz, jk_gt(elem), jb_gt(elem), t+1) + wb_gt(elem) * wk_gt(elem) * DIST_step(elem);
            end
            
        end
        
    end
    
    % Adjust distributions based on demographics
    DIST = repmat(shiftdim(mu2_idem(T_past+(1:T_dist)), -2), [nz,nk,nb,1]) .* DIST;
    
else
    
    DIST = DIST_static;
    T_dist = size(DIST, 4);
    
end



%% Aggregate generation

Cohort.assets  = sum(reshape(K  (:,:,:,1:T_dist) .* DIST, [], T_dist), 1) .* (1 + (mu3_idem(T_past+(1:T_dist)) ./ mu2_idem(T_past+(1:T_dist))));
Cohort.beqs    = sum(reshape(K  (:,:,:,1:T_dist) .* DIST, [], T_dist), 1) .* (0 + (mu3_idem(T_past+(1:T_dist)) ./ mu2_idem(T_past+(1:T_dist))));
Cohort.labeffs = sum(reshape(LAB(:,:,:,1:T_dist) .* repmat(reshape(zs(:,T_past+(1:T_dist),idem), [nz,1,1,T_dist]), [1,nk,nb,1]) .* DIST, [], T_dist), 1);
Cohort.labs    = sum(reshape(LAB(:,:,:,1:T_dist) .* DIST, [], T_dist), 1);
Cohort.lfprs   = sum(reshape((LAB(:,:,:,1:T_dist) >= 0.01) .* DIST, [], T_dist), 1);
Cohort.incs    = sum(reshape(INC(:,:,:,1:T_dist) .* DIST, [], T_dist), 1);
Cohort.pits    = sum(reshape(PIT(:,:,:,1:T_dist) .* DIST, [], T_dist), 1);
Cohort.ssts    = sum(reshape(SST(:,:,:,1:T_dist) .* DIST, [], T_dist), 1);
Cohort.cits    = sum(reshape(CIT(:,:,:,1:T_dist) .* DIST, [], T_dist), 1);
Cohort.bens    = sum(reshape(BEN(:,:,:,1:T_dist) .* DIST, [], T_dist), 1);


end




% Resource and tax calculation function
function [resources, inc, pit, sst, cit] = calculate_resources(...
             labinc, ks_ik_, year_, ...
             mpci_, rpci_, sstaxcredit_, sstax_, ssincmax_, deduc_coefs_, pit_coefs_, captaxshare_, taucap_, taucapgain_, qtobin_, qtobin0_, ...
             beq_, capshare_, debtshare_, caprate_, govrate_, totrate_, expsub_) %#codegen

% Enforce function inlining for C code generation
coder.inline('always');

% Define parameters as persistent variables
persistent initialized

persistent ks_ik
persistent year
persistent mpci
persistent rpci
persistent sstaxcredit
persistent sstax
persistent ssincmax
persistent deduc_coefs
persistent pit_coefs
persistent captaxshare
persistent taucap
persistent taucapgain
persistent qtobin
persistent qtobin0
persistent beq
persistent capshare
persistent debtshare
persistent caprate
persistent govrate
persistent totrate
persistent expsub

% Initialize parameters
if isempty(initialized)
    
    ks_ik           = 0;
    year            = 0;
    mpci            = 0;
    rpci            = 0;
    sstaxcredit     = 0;
    sstax           = 0;
    ssincmax        = 0;
    deduc_coefs     = 0;
    pit_coefs       = 0;
    captaxshare     = 0;
    taucap          = 0;
    taucapgain      = 0;
    qtobin          = 0;
    qtobin0         = 0;
    beq             = 0;
    capshare        = 0;
    debtshare       = 0;
    caprate         = 0;
    govrate         = 0;
    totrate         = 0;
    expsub          = 0;
    
    initialized     = true;
    
end

% Set parameters if provided
if (nargin > 1)
    
    ks_ik           = ks_ik_      ;
    year            = year_       ;
    mpci            = mpci_       ;
    rpci            = rpci_       ;
    sstaxcredit     = sstaxcredit_;
    sstax           = sstax_      ;
    ssincmax        = ssincmax_   ;
    deduc_coefs     = deduc_coefs_;
    pit_coefs       = pit_coefs_  ;
    captaxshare     = captaxshare_;
    taucap          = taucap_     ;
    taucapgain      = taucapgain_ ;
    qtobin          = qtobin_     ;
    qtobin0         = qtobin0_    ;
    beq             = beq_        ;
    capshare        = capshare_   ;
    debtshare       = debtshare_  ;
    caprate         = caprate_    ;
    govrate         = govrate_    ;
    totrate         = totrate_    ;
    expsub          = expsub_     ;
    
    if isempty(labinc), return, end
    
end

% Calculate taxable income
inc     = (rpci/mpci)*max(0, capshare*(caprate-1)*ks_ik*(1-captaxshare) + debtshare*(govrate-1)*ks_ik + (1-sstaxcredit)*labinc);
deduc   = max(0, deduc_coefs*inc.^[0; 1; 0.5]);
inc_eff = max(inc - deduc, 0);
inc     = (mpci/rpci)*inc;

% Calculate personal income tax (PIT)
pit = (mpci/rpci)*pit_coefs(1)*(inc_eff - (inc_eff.^(-pit_coefs(2)) + (pit_coefs(3))).^(-1/pit_coefs(2)));

% Calculate Social Security tax
sst = sstax*min(labinc, ssincmax);

% Calculate capital income tax (CIT)
cit = capshare*ks_ik*(taucap*((caprate-1) - expsub)*captaxshare + taucapgain*(year == 1)*(qtobin - qtobin0)/qtobin0);

% Calculate available resources
resources = totrate*ks_ik + labinc - (pit + sst + cit) + beq + (year == 1)*ks_ik*capshare*(qtobin - qtobin0)/qtobin0;

end




% Retirement age value function
function v = value_retirement(k, ks_, resources_, EV_ib_, sigma_, gamma_)

% Enforce function inlining for C code generation
coder.inline('always');

% Define parameters as persistent variables
persistent initialized

persistent ks
persistent resources
persistent EV_ib
persistent sigma
persistent gamma

% Initialize parameters
if isempty(initialized)
    
    ks          = 0;
    resources   = 0;
    EV_ib       = 0;
    sigma       = 0;
    gamma       = 0;
    
    initialized = true;
    
end

% Set parameters if provided
if (nargin > 1)
    
    ks          = ks_       ;
    resources   = resources_;
    EV_ib       = EV_ib_    ;
    sigma       = sigma_    ;
    gamma       = gamma_    ;
    
    return
    
end

% Calculate consumption
consumption = resources - k;

% Perform bound checks
if (ks(1) <= k) && (k <= ks(end)) && (0 <= consumption)
    
    % Calculate value
    v = interp1(ks, EV_ib, k, 'linear') + (consumption^(gamma*(1-sigma)))/(1-sigma);
    
    % Negate for minimization and force to scalar for C code generation
    v = -v(1);
    
else
    v = Inf;
end


end




% Working age value function
function v = value_working(x, ks_, bs_, wage_eff_, EV_, sigma_, gamma_)

% Enforce function inlining for C code generation
coder.inline('always');

% Define parameters as persistent variables
persistent initialized

persistent ks
persistent bs
persistent wage_eff
persistent EV
persistent sigma
persistent gamma

% Initialize parameters
if isempty(initialized)
    
    ks          = 0;
    bs          = 0;
    wage_eff    = 0;
    EV          = 0;
    sigma       = 0;
    gamma       = 0;
    
    initialized = true;
    
end

% Set parameters if provided
if (nargin > 1)
    
    ks          = ks_      ;
    bs          = bs_      ;
    wage_eff    = wage_eff_;
    EV          = EV_      ;
    sigma       = sigma_   ;
    gamma       = gamma_   ;
    
    return
    
end

% Define decision variables and perform bound checks
k   = x(1);
lab = x(2);

labinc = wage_eff * lab;
b      = calculate_b(labinc);

if ~((0 <= lab) && (lab <= 1) ...
     && (ks(1) <= k) && (k <= ks(end)) ...
     && (bs(1) <= b) && (b <= bs(end)))
    
    v = Inf;
    return
    
end

% Calculate available resources
resources = calculate_resources(labinc);

% Calculate consumption and perform bound check
consumption = resources - k;

if ~(0 <= consumption)
    v = Inf;
    return
end

% Calculate value
v = interp2(ks', bs, EV', k, b, 'linear') + (1/(1-sigma))*((consumption^gamma)*((1-lab)^(1-gamma)))^(1-sigma);

% Negate for minimization and force to scalar for C code generation
v = -v(1);

end




% Average earnings calculation function
function [b] = calculate_b(labinc, age_, bs_ib_, ssincmax_)

% Define parameters as persistent variables
persistent initialized

persistent ssincmax
persistent age
persistent bs_ib

% Initialize parameters
if isempty(initialized)
    
    ssincmax    = 0;
    age         = 0;
    bs_ib       = 0;
    
    initialized = true;
    
end

% Set parameters if provided
if (nargin > 1)
    
    ssincmax    = ssincmax_;
    age         = age_     ;
    bs_ib       = bs_ib_   ;
    
    return
    
end

% Enforce function inlining for C code generation
coder.inline('always');

b = (bs_ib*(age-1) + min(labinc, ssincmax))/age;

end



