%%
% Solve dynamic optimization problem for a cohort.
% 
%%


function [V, OPT] = solve_cohort(V0, LAB_static, isdynamic, ...
                        nz, nk, nb, T_past, T_shift, T_active, T_work, T_model, zs_idem, transz, kv, bv, beta, gamma, sigma, surv, V_beq, ...
                        modelunit_dollars, ...
                        sstaxcredit, ssbenefits, sstaxs, ssincmaxs, ...
                        tax_thresholds, tax_burden, tax_rates, ... 
                        captaxshare, taucap, taucapgain, qtobin, qtobin0, ...
                        beqs, wages, capshares, caprates, govrates, totrates, expsubs) %#codegen


%% Argument verification

nz_max          = 50;
nk_max          = 50;
nb_max          = 50;
T_max           = 100;
nthresholds_max = 20;

assert( isa(V0          , 'double'  ) && (size(V0           , 1) <= nz_max  ) && (size(V0           , 2) <= nk_max  ) && (size(V0           , 3) <= nb_max  ) );
assert( isa(LAB_static  , 'double'  ) && (size(LAB_static   , 1) <= nz_max  ) && (size(LAB_static   , 2) <= nk_max  ) && (size(LAB_static   , 3) <= nb_max  ) && (size(LAB_static   , 4) <= T_max   ) );
assert( isa(isdynamic   , 'logical' ) && (size(isdynamic    , 1) == 1       ) && (size(isdynamic    , 2) == 1       ) );

assert( isa(nz          , 'double'  ) && (size(nz           , 1) == 1       ) && (size(nz           , 2) == 1       ) );
assert( isa(nk          , 'double'  ) && (size(nk           , 1) == 1       ) && (size(nk           , 2) == 1       ) );
assert( isa(nb          , 'double'  ) && (size(nb           , 1) == 1       ) && (size(nb           , 2) == 1       ) );
assert( isa(T_past      , 'double'  ) && (size(T_past       , 1) == 1       ) && (size(T_past       , 2) == 1       ) );
assert( isa(T_shift     , 'double'  ) && (size(T_shift      , 1) == 1       ) && (size(T_shift      , 2) == 1       ) );
assert( isa(T_active    , 'double'  ) && (size(T_active     , 1) == 1       ) && (size(T_active     , 2) == 1       ) );
assert( isa(T_work      , 'double'  ) && (size(T_work       , 1) == 1       ) && (size(T_work       , 2) == 1       ) );
assert( isa(T_model     , 'double'  ) && (size(T_model      , 1) == 1       ) && (size(T_model      , 2) == 1       ) );
assert( isa(zs_idem     , 'double'  ) && (size(zs_idem      , 1) <= nz_max  ) && (size(zs_idem      , 2) <= T_max   ) );
assert( isa(transz      , 'double'  ) && (size(transz       , 1) <= nz_max  ) && (size(transz       , 2) <= nz_max  ) );
assert( isa(kv          , 'double'  ) && (size(kv           , 1) <= nk_max  ) && (size(kv           , 2) == 1       ) );
assert( isa(bv          , 'double'  ) && (size(bv           , 1) <= nb_max  ) && (size(bv           , 2) == 1       ) );
assert( isa(beta        , 'double'  ) && (size(beta         , 1) == 1       ) && (size(beta         , 2) == 1       ) );
assert( isa(gamma       , 'double'  ) && (size(gamma        , 1) == 1       ) && (size(gamma        , 2) == 1       ) );
assert( isa(sigma       , 'double'  ) && (size(sigma        , 1) == 1       ) && (size(sigma        , 2) == 1       ) );
assert( isa(surv        , 'double'  ) && (size(surv         , 1) == 1       ) && (size(surv         , 2) <= T_max   ) );
assert( isa(V_beq       , 'double'  ) && (size(V_beq        , 1) <= nk_max  ) && (size(V_beq        , 2) == 1       ) );

assert( isa(modelunit_dollars, 'double'  ) && (size(modelunit_dollars, 1) == 1       ) && (size(modelunit_dollars, 2) == 1       ) );

assert( isa(sstaxcredit , 'double'  ) && (size(sstaxcredit  , 1) == 1       ) && (size(sstaxcredit  , 2) == 1       ) );
assert( isa(ssbenefits  , 'double'  ) && (size(ssbenefits   , 1) <= nb_max  ) && (size(ssbenefits   , 2) <= T_max   ) );
assert( isa(sstaxs      , 'double'  ) && (size(sstaxs       , 1) == 1       ) && (size(sstaxs       , 2) <= T_max   ) );
assert( isa(ssincmaxs   , 'double'  ) && (size(ssincmaxs    , 1) == 1       ) && (size(ssincmaxs    , 2) <= T_max   ) );

assert( isa(tax_thresholds, 'double' ) && (size(tax_thresholds, 1) == 1     ) && (size(tax_thresholds, 2) <= nthresholds_max ) );
assert( isa(tax_burden    , 'double' ) && (size(tax_burden    , 1) == 1     ) && (size(tax_burden    , 2) <= nthresholds_max ) );
assert( isa(tax_rates     , 'double' ) && (size(tax_rates     , 1) == 1     ) && (size(tax_rates     , 2) <= nthresholds_max ) );

assert( isa(captaxshare , 'double'  ) && (size(captaxshare  , 1) == 1       ) && (size(captaxshare  , 2) == 1       ) );
assert( isa(taucap      , 'double'  ) && (size(taucap       , 1) == 1       ) && (size(taucap       , 2) == 1       ) );
assert( isa(taucapgain  , 'double'  ) && (size(taucapgain   , 1) == 1       ) && (size(taucapgain   , 2) == 1       ) );
assert( isa(qtobin      , 'double'  ) && (size(qtobin       , 1) == 1       ) && (size(qtobin       , 2) == 1       ) );
assert( isa(qtobin0     , 'double'  ) && (size(qtobin0      , 1) == 1       ) && (size(qtobin0      , 2) == 1       ) );

assert( isa(beqs        , 'double'  ) && (size(beqs         , 1) == 1       ) && (size(beqs         , 2) <= T_max   ) );
assert( isa(wages       , 'double'  ) && (size(wages        , 1) == 1       ) && (size(wages        , 2) <= T_max   ) );
assert( isa(capshares   , 'double'  ) && (size(capshares    , 1) == 1       ) && (size(capshares    , 2) <= T_max   ) );
assert( isa(caprates    , 'double'  ) && (size(caprates     , 1) == 1       ) && (size(caprates     , 2) <= T_max   ) );
assert( isa(govrates    , 'double'  ) && (size(govrates     , 1) == 1       ) && (size(govrates     , 2) <= T_max   ) );
assert( isa(totrates    , 'double'  ) && (size(totrates     , 1) == 1       ) && (size(totrates     , 2) <= T_max   ) );
assert( isa(expsubs     , 'double'  ) && (size(expsubs      , 1) == 1       ) && (size(expsubs      , 2) <= T_max   ) );



%% Dynamic optimization

% Initialize utility and optimal decision value arrays
V       = zeros(nz,nk,nb,T_active);   % Utility

OPT.K   = zeros(nz,nk,nb,T_active);   % Savings
OPT.LAB = zeros(nz,nk,nb,T_active);   % Labor level
OPT.B   = zeros(nz,nk,nb,T_active);   % Average earnings

OPT.INC = zeros(nz,nk,nb,T_active);   % Taxable income
OPT.PIT = zeros(nz,nk,nb,T_active);   % Personal income tax
OPT.SST = zeros(nz,nk,nb,T_active);   % Social Security tax
OPT.CIT = zeros(nz,nk,nb,T_active);   % Corporate income tax
OPT.BEN = zeros(nz,nk,nb,T_active);   % Social Security benefits

% Initialize forward-looking utility values
V_step = V0;

% Specify settings for dynamic optimization subproblems
optim_options = optimset('Display', 'off', 'TolFun', 1e-4, 'TolX', 1e-4);


% Solve dynamic optimization problem through backward induction
for t = T_active:-1:1
    
    % Determine age and year, bounded by modeling period
    age  = t + T_past;
    year = min(t + T_shift, T_model);
    
    % Extract parameters for current year
    ssbenefit  = ssbenefits(:, year);
    sstax      = sstaxs       (year);
    ssincmax   = ssincmaxs    (year);
    beq        = beqs         (year);
    wage       = wages        (year);
    caprate    = caprates     (year);
    govrate    = govrates     (year);
    capshare   = capshares    (year);
    totrate    = totrates     (year);
    expsub     = expsubs      (year);
    
    
    for ib = 1:nb
        for ik = 1:nk
            
            if (age > T_work)
                
                % Calculate available resources and tax terms
                labinc = ssbenefit(ib);
                [resources, inc, pit, ~, cit] = calculate_resources(labinc, kv(ik), year, ...
                    modelunit_dollars, ...
                    sstaxcredit, 0, 0, ...
                    tax_thresholds, tax_burden, tax_rates, ... 
                    captaxshare, taucap, taucapgain, qtobin, qtobin0, ...
                    beq, capshare, caprate, govrate, totrate, expsub);
                
                if isdynamic
                    
                    % Calculate expected value curve using forward-looking utility values
                    EV = (1-surv(age))*repmat(V_beq, [1,nb]) + surv(age)*beta*reshape(V_step(1,:,:), [nk,nb]);
                    
                    % Call retirement age value function to set parameters
                    value_retirement([], kv, resources, EV(:,ib), sigma, gamma);
                    
                    % Solve dynamic optimization subproblem
                    [k, v] = fminsearch(@value_retirement, kv(ik), optim_options);
                    
                    % Record utility and optimal decision values
                    V    (:,ik,ib,t) = -v;
                    OPT.K(:,ik,ib,t) = k ;
                    
                end
                
                OPT.LAB(:,ik,ib,t) = 0     ;
                OPT.B  (:,ik,ib,t) = bv(ib);
                
                OPT.INC(:,ik,ib,t) = inc   ;
                OPT.PIT(:,ik,ib,t) = pit   ;
                OPT.SST(:,ik,ib,t) = 0     ;
                OPT.CIT(:,ik,ib,t) = cit   ;
                OPT.BEN(:,ik,ib,t) = labinc;
                
            else
                
                for iz = 1:nz
                    
                    % Calculate effective wage
                    wage_eff = wage * zs_idem(iz,age);
                    
                    % Call resource calculation function to set parameters
                    calculate_resources([], kv(ik), year, ...
                        modelunit_dollars, ...
                        0, sstax, ssincmax, ...
                        tax_thresholds, tax_burden, tax_rates, ... 
                        captaxshare, taucap, taucapgain, qtobin, qtobin0, ...
                        beq, capshare, caprate, govrate, totrate, expsub);
                    
                    if isdynamic
                        
                        % Calculate expected value curve using forward-looking utility values
                        EV = (1-surv(age))*repmat(V_beq, [1,nb]) + surv(age)*beta*reshape(sum(repmat(transz(iz,:)', [1,nk,nb]) .* V_step, 1), [nk,nb]);
                        
                        % Call working age value function and average earnings calculation function to set parameters
                        value_working([], kv, bv, wage_eff, EV, sigma, gamma)
                        calculate_b  ([], age, bv(ib), ssincmax)
                        
                        % Solve dynamic optimization subproblem
                        lab0 = 0.5;
                        k0   = max(kv(ik), 0.1 * wage_eff * lab0);   % (Assumes taxation will not exceed 90% of labor income)
                        
                        [x, v] = fminsearch(@value_working, [k0, lab0], optim_options);
                        
                        k   = x(1);
                        lab = x(2);
                        
                        % Record utility and optimal decision values
                        V    (iz,ik,ib,t) = -v;
                        OPT.K(iz,ik,ib,t) = k ;
                        
                    else
                        lab = LAB_static(iz,ik,ib,t);
                    end
                    
                    labinc = wage_eff * lab;
                    [~, inc, pit, sst, cit] = calculate_resources(labinc);
                    
                    OPT.LAB(iz,ik,ib,t) = lab;
                    OPT.B  (iz,ik,ib,t) = calculate_b(labinc);
                    
                    OPT.INC(iz,ik,ib,t) = inc;
                    OPT.PIT(iz,ik,ib,t) = pit;
                    OPT.SST(iz,ik,ib,t) = sst;
                    OPT.CIT(iz,ik,ib,t) = cit;
                    OPT.BEN(iz,ik,ib,t) = 0  ;
                    
                end
                
            end
            
        end
    end
    
    % Update forward-looking utility values
    V_step = V(:,:,:,t);
    
end


end




% Resource and tax calculation function
function [resources, inc, pit, sst, cit] = calculate_resources(labinc, kv_ik_, year_, ...
             modelunit_dollars_, ...
             sstaxcredit_, sstax_, ssincmax_, ...
             tax_thresholds_, tax_burden_, tax_rates_, ... 
             captaxshare_, taucap_, taucapgain_, qtobin_, qtobin0_, ...
             beq_, capshare_, caprate_, govrate_, totrate_, expsub_) %#codegen

% Enforce function inlining for C code generation
coder.inline('always');

% Define parameters as persistent variables
persistent kv_ik year ...
           modelunit_dollars ...
           sstaxcredit sstax ssincmax ...
           tax_thresholds tax_burden tax_rates ... 
           captaxshare taucap taucapgain qtobin qtobin0 ...
           beq capshare caprate govrate totrate expsub ...
           capgain ...
           initialized

% Initialize parameters for C code generation
if isempty(initialized)
    kv_ik = 0; year = 0;
    modelunit_dollars = 0; 
    sstaxcredit = 0; sstax = 0; ssincmax = 0; 
    tax_thresholds = 0; tax_burden = 0; tax_rates = 0;
    captaxshare = 0; taucap = 0; taucapgain = 0; qtobin = 0; qtobin0 = 0;
    beq = 0; capshare = 0; caprate = 0; govrate = 0; totrate = 0; expsub = 0;
    capgain = 0;
    initialized = true;
end

% Set parameters if provided
if (nargin > 1)
    kv_ik = kv_ik_; year = year_;
    modelunit_dollars = modelunit_dollars_; 
    sstaxcredit = sstaxcredit_; sstax = sstax_; ssincmax = ssincmax_; 
    tax_thresholds = tax_thresholds_; tax_burden = tax_burden_; tax_rates = tax_rates_;
    captaxshare = captaxshare_; taucap = taucap_; taucapgain = taucapgain_; qtobin = qtobin_; qtobin0 = qtobin0_;
    beq = beq_; capshare = capshare_; caprate = caprate_; govrate = govrate_; totrate = totrate_; expsub = expsub_;
    
    % Pre-calculate the percent cap gain (adjusted for realization)
    capgain = 0.25*(year == 1)*(qtobin - qtobin0)/qtobin; 
    
    if isempty(labinc), return, end
end


% Calculate taxable income in dollars
%   We do not allow negative incomes
inc     = (1/modelunit_dollars)*max(0, capshare*caprate*kv_ik*(1-captaxshare) + (1-capshare)*govrate*kv_ik + (1-sstaxcredit)*labinc);

% Calculate personal income tax
%       Expect equal-size vectors with tax_thresholds(1)=0
%       tax_rates apply for income between tax_thresholds(i-1) and tax_thresholds(i)
%       tax_burden are pre-calculated total tax liability at tax_thresholds
%       pit_dollar is income tax in dollars
bracket     = find(tax_thresholds <= inc, 1, 'last');
bracket     = bracket(1);   % Force to scalar for C code generation
pit_dollar  = tax_burden(bracket) + tax_rates(bracket)*(inc - tax_thresholds(bracket));
pit         = modelunit_dollars*pit_dollar;     % Convert to model units

% Calculate Social Security tax
sst = sstax * min(labinc, ssincmax);

% Calculate corporate income tax
cit = capshare*kv_ik*(taucap*(caprate - expsub)*captaxshare + taucapgain*capgain);

% Calculate available resources
resources = (1 + totrate)*kv_ik + labinc - (pit + sst + cit) + beq + kv_ik*capshare*capgain;

end




% Retirement age value function
function v = value_retirement(k, kv_, resources_, EV_ib_, sigma_, gamma_)

% Enforce function inlining for C code generation
coder.inline('always');

% Define parameters as persistent variables
persistent kv resources EV_ib sigma gamma ...
           initialized

% Initialize parameters for C code generation
if isempty(initialized)
    kv = 0; resources = 0; EV_ib = 0; sigma = 0; gamma = 0;
    initialized = true;
end

% Set parameters if provided
if (nargin > 1)
    kv = kv_; resources = resources_; EV_ib = EV_ib_; sigma = sigma_; gamma = gamma_;
    if isempty(k), return, end
end


% Calculate consumption
consumption = resources - k;

% Perform bound checks
if (kv(1) <= k) && (k <= kv(end)) && (0 <= consumption)
    
    % Calculate utility
    v = interp1(kv, EV_ib, k, 'linear') + (consumption^(gamma*(1-sigma)))/(1-sigma);
    
    % Negate utility for minimization and force to scalar for C code generation
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
persistent kv bv wage_eff EV sigma gamma ...
           initialized

% Initialize parameters for C code generation
if isempty(initialized)
    kv = 0; bv = 0; wage_eff = 0; EV = 0; sigma = 0; gamma = 0;
    initialized = true;
end

% Set parameters if provided
if (nargin > 1)
    kv = kv_; bv = bv_; wage_eff = wage_eff_; EV = EV_; sigma = sigma_; gamma = gamma_;
    if isempty(x), return, end
end


% Define decision variables and perform bound checks
k   = x(1);
lab = x(2);

labinc = wage_eff * lab;
b      = calculate_b(labinc);

if ~((0 <= lab) && (lab <= 1) ...
     && (kv(1) <= k) && (k <= kv(end)) ...
     && (bv(1) <= b) && (b <= bv(end)))
    
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

% Calculate utility
v = interp2(kv', bv, EV', k, b, 'linear') + (1/(1-sigma))*((consumption^gamma)*((1-lab)^(1-gamma)))^(1-sigma);

% Negate utility for minimization and force to scalar for C code generation
v = -v(1);

end




% Average earnings calculation function
function [b] = calculate_b(labinc, age_, bv_ib_, ssincmax_)

% Enforce function inlining for C code generation
coder.inline('always');

% Define parameters as persistent variables
persistent age bv_ib ssincmax ...
           initialized

% Initialize parameters for C code generation
if isempty(initialized)
    age = 0; bv_ib = 0; ssincmax = 0;
    initialized = true;
end

% Set parameters if provided
if (nargin > 1)
    age = age_; bv_ib = bv_ib_; ssincmax = ssincmax_;
    if isempty(labinc), return, end
end


% Calculate average earnings
b = (bv_ib*(age-1) + min(labinc, ssincmax))/age;

end



