%%
% Solve dynamic optimization problem for a cohort.
% 
%%


function [OPT] = solve_cohort(V0, LAB_static, isdynamic, ...
                        nz, nk, nb, T_past, T_shift, T_active, T_work, T_model, ... 
                        zs_idem, transz, kv, bv, beta, gamma, sigma, surv, ...
                        bequest_phi_1, bequest_phi_2, bequest_phi_3, ...
                        modelunit_dollar, ...
                        sstaxcredit, ssbenefits, ssincmaxs, ...
                        sstax_brackets, sstax_burdens, sstax_rates, ...
                        pittax_brackets, pittax_burdens, pittax_rates, ... 
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

assert( isa(bequest_phi_1, 'double' ) && (size(bequest_phi_1, 1) == 1 ) && (size(bequest_phi_1, 2) == 1 ) );
assert( isa(bequest_phi_2, 'double' ) && (size(bequest_phi_2, 1) == 1 ) && (size(bequest_phi_2, 2) == 1 ) );
assert( isa(bequest_phi_3, 'double' ) && (size(bequest_phi_3, 1) == 1 ) && (size(bequest_phi_3, 2) == 1 ) );

assert( isa(modelunit_dollar, 'double' ) && (size(modelunit_dollar, 1) == 1 ) && (size(modelunit_dollar, 2) == 1 ) );

assert( isa(sstaxcredit , 'double'  ) && (size(sstaxcredit  , 1) == 1       ) && (size(sstaxcredit  , 2) == 1       ) );
assert( isa(ssbenefits  , 'double'  ) && (size(ssbenefits   , 1) <= nb_max  ) && (size(ssbenefits   , 2) == 1       ) );
assert( isa(ssincmaxs   , 'double'  ) && (size(ssincmaxs    , 1) == 1       ) && (size(ssincmaxs    , 2) <= T_max   ) );

assert( isa(sstax_brackets  , 'double' ) && (size(sstax_brackets  , 1) <= T_max ) && (size(sstax_brackets  , 2) <= nthresholds_max ) );
assert( isa(sstax_burdens   , 'double' ) && (size(sstax_burdens   , 1) <= T_max ) && (size(sstax_burdens   , 2) <= nthresholds_max ) );
assert( isa(sstax_rates     , 'double' ) && (size(sstax_rates     , 1) <= T_max ) && (size(sstax_rates     , 2) <= nthresholds_max ) );

assert( isa(pittax_brackets  , 'double' ) && (size(pittax_brackets  , 1) <= T_max ) && (size(pittax_brackets  , 2) <= nthresholds_max ) );
assert( isa(pittax_burdens   , 'double' ) && (size(pittax_burdens   , 1) <= T_max ) && (size(pittax_burdens   , 2) <= nthresholds_max ) );
assert( isa(pittax_rates     , 'double' ) && (size(pittax_rates     , 1) <= T_max ) && (size(pittax_rates     , 2) <= nthresholds_max ) );

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
OPT.V   = zeros(nz,nk,nb,T_active);   % Utility

OPT.K   = zeros(nz,nk,nb,T_active);   % Savings
OPT.LAB = zeros(nz,nk,nb,T_active);   % Labor level
OPT.B   = zeros(nz,nk,nb,T_active);   % Average earnings

OPT.INC = zeros(nz,nk,nb,T_active);   % Taxable income
OPT.PIT = zeros(nz,nk,nb,T_active);   % Personal income tax
OPT.SST = zeros(nz,nk,nb,T_active);   % Social Security tax
OPT.CIT = zeros(nz,nk,nb,T_active);   % Corporate income tax
OPT.BEN = zeros(nz,nk,nb,T_active);   % Social Security benefits
OPT.CON = zeros(nz,nk,nb,T_active);   % Consumption

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
    ssincmax   = ssincmaxs    (year);
    beq        = beqs         (year);
    wage       = wages        (year);
    caprate    = caprates     (year);
    govrate    = govrates     (year);
    capshare   = capshares    (year);
    totrate    = totrates     (year);
    expsub     = expsubs      (year);
    
    sst_brackets    = sstax_brackets    (year, :);
    sst_burdens     = sstax_burdens     (year, :);
    sst_rates       = sstax_rates       (year, :);
    
    pit_brackets    = pittax_brackets   (year, :);
    pit_burdens     = pittax_burdens    (year, :);
    pit_rates       = pittax_rates      (year, :);
    
    % Pre-calculate for speed and conciseness
    bequest_p_1   = beta * (1-surv(age))* bequest_phi_1;
    
    for ib = 1:nb
        for ik = 1:nk
            
            if (age > T_work)
                
                % Calculate available resources and tax terms
                ssinc = ssbenefits(ib);
                [resources, inc, pit, ~, cit] = calculate_resources(0, kv(ik), year, ...
                    modelunit_dollar, ...
                    ssinc, sstaxcredit, ...
                    sst_brackets, sst_burdens, sst_rates, ...
                    pit_brackets, pit_burdens, pit_rates, ... 
                    captaxshare, taucap, taucapgain, qtobin, qtobin0, ...
                    beq, capshare, caprate, govrate, totrate, expsub);
                
                if isdynamic
                    
                    % Calculate expected value conditional on living using forward-looking 
                    %   utility values. Pre-multiply by prob. survival and
                    %   beta to save on computation.
                    EV = surv(age)*beta*reshape(V_step(1,:,:), [nk,nb]);
                    
                    % Call retirement age value function to set parameters
                    value_retirement([], kv, resources, EV(:,ib), ... 
                                bequest_p_1, bequest_phi_2, bequest_phi_3, ...
                                sigma, gamma);
                    
                    % Solve dynamic optimization subproblem
                    [k, v] = fminsearch(@value_retirement, kv(ik), optim_options);
                    
                    % Checks -> only work in the absence of mex file!
                    assert( ~isinf(v)   , 'v is inf')
                    assert( k <= kv(end), 'k is too big!')

                    % Record utility and optimal decision values
                    OPT.V(:,ik,ib,t) = -v;
                    OPT.K(:,ik,ib,t) = k ;
                    
                else
                    
                    k = kv(ik);
                    
                end
                
                OPT.LAB(:,ik,ib,t) = 0     ;
                OPT.B  (:,ik,ib,t) = bv(ib);
                
                OPT.INC(:,ik,ib,t) = inc   ;
                OPT.PIT(:,ik,ib,t) = pit   ;
                OPT.SST(:,ik,ib,t) = 0     ;
                OPT.CIT(:,ik,ib,t) = cit   ;
                OPT.BEN(:,ik,ib,t) = ssinc;
                OPT.CON(:,ik,ib,t) = resources - k;
                
            else
                
                for iz = 1:nz
                    
                    % Calculate effective wage
                    wage_eff = wage * zs_idem(iz,age);
                    
                    % Call resource calculation function to set parameters
                    calculate_resources([], kv(ik), year, ...
                        modelunit_dollar, ...
                        0, 0, ...
                        sst_brackets, sst_burdens, sst_rates, ...
                        pit_brackets, pit_burdens, pit_rates, ... 
                        captaxshare, taucap, taucapgain, qtobin, qtobin0, ...
                        beq, capshare, caprate, govrate, totrate, expsub);
                    
                    if isdynamic
                        
                        % Calculate expected value conditional on living using forward-looking 
                        %   utility values. Pre-multiply by prob. survival and
                        %   beta to save on computation.
                        EV = surv(age)*beta*reshape(sum(repmat(transz(iz,:)', [1,nk,nb]) .* V_step, 1), [nk,nb]);
                        
                        % Call working age value function and average earnings calculation function to set parameters
                        value_working([], kv, bv, wage_eff, EV, ...
                                bequest_p_1, bequest_phi_2, bequest_phi_3, ...
                                sigma, gamma);
                        calculate_b  ([], age, bv(ib), bv(end), ssincmax);
                        
                        % Solve dynamic optimization subproblem
                        lab0 = 0.5;
                        k0   = max(kv(ik), min(kv(end), 0.1 * wage_eff * lab0));   % Assumes taxation will not exceed 90% of labor income and at the same time forces k to be in the grid
                        
                        [x, v] = fminsearch(@value_working, [k0, lab0], optim_options);

                        k   = x(1);
                        lab = x(2);       
                        
                        % Checks -> only work in the absence of mex file!
                        assert( ~isinf(v)   , 'v is inf')
                        assert( k <= kv(end), 'k is too big!')
                        
                        % Record utility and optimal decision values
                        OPT.V(iz,ik,ib,t) = -v;
                        OPT.K(iz,ik,ib,t) = k ;
                        
                    else
                        k = kv(ik);
                        lab = LAB_static(iz,ik,ib,t);
                    end
                    
                    labinc = wage_eff * lab;
                    [resources, inc, pit, sst, cit] = calculate_resources(labinc);
                    
                    OPT.LAB(iz,ik,ib,t) = lab;
                    OPT.B  (iz,ik,ib,t) = calculate_b(labinc);
                    
                    OPT.INC(iz,ik,ib,t) = inc;
                    OPT.PIT(iz,ik,ib,t) = pit;
                    OPT.SST(iz,ik,ib,t) = sst;
                    OPT.CIT(iz,ik,ib,t) = cit;
                    OPT.BEN(iz,ik,ib,t) = 0  ;
                    OPT.CON(iz,ik,ib,t) = resources - k;
                    
                end
                
            end
            
        end
    end
    
    % Update forward-looking utility values
    V_step = OPT.V(:,:,:,t);
    
end


end




% Resource and tax calculation function
function [resources, inc, pit, sst, cit] = calculate_resources(labinc, ...
             kv_ik_, year_, ...
             modelunit_dollar_, ...
             ssinc_, sstaxcredit_, ...
             sst_thresholds_, sst_burdens_, sst_rates_, ...
             pit_thresholds_, pit_burdens_, pit_rates_, ... 
             captaxshare_, taucap_, taucapgain_, qtobin_, qtobin0_, ...
             beq_, capshare_, caprate_, govrate_, totrate_, expsub_) %#codegen

% Enforce function inlining for C code generation
coder.inline('always');

% Define parameters as persistent variables
persistent kv_ik year ...
           modelunit_dollar ...
           ssinc sstaxcredit ...
           sst_thresholds sst_burdens sst_rates ...
           pit_thresholds pit_burdens pit_rates ... 
           captaxshare taucap taucapgain qtobin qtobin0 ...
           beq capshare caprate govrate totrate expsub ...
           capgain ...
           initialized

% Initialize parameters for C code generation
if isempty(initialized)
    kv_ik = 0; year = 0;
    modelunit_dollar = 0; 
    ssinc = 0; sstaxcredit = 0; 
    sst_thresholds = 0; sst_burdens = 0; sst_rates = 0;
    pit_thresholds = 0; pit_burdens = 0; pit_rates = 0;
    captaxshare = 0; taucap = 0; taucapgain = 0; qtobin = 0; qtobin0 = 0;
    beq = 0; capshare = 0; caprate = 0; govrate = 0; totrate = 0; expsub = 0;
    capgain = 0;
    initialized = true;
end

% Set parameters if provided
if (nargin > 1)
    kv_ik = kv_ik_; year = year_;
    modelunit_dollar = modelunit_dollar_; 
    ssinc = ssinc_; sstaxcredit = sstaxcredit_; 
    sst_thresholds = sst_thresholds_; sst_burdens = sst_burdens_; sst_rates = sst_rates_;
    pit_thresholds = pit_thresholds_; pit_burdens = pit_burdens_; pit_rates = pit_rates_;
    captaxshare = captaxshare_; taucap = taucap_; taucapgain = taucapgain_; qtobin = qtobin_; qtobin0 = qtobin0_;
    beq = beq_; capshare = capshare_; caprate = caprate_; govrate = govrate_; totrate = totrate_; expsub = expsub_;
    
    % Pre-calculate the percent cap gain (adjusted for realization)
    capgain = 0*(year == 1)*(qtobin - qtobin0)/qtobin; 
    
    if isempty(labinc), return, end
end


% Calculate taxable income in dollars
%   We do not allow negative incomes
inc         = max( 0,...
            capshare*caprate*kv_ik*(1-captaxshare) + (1-capshare)*govrate*kv_ik ...
            + (1-sstaxcredit)*ssinc + labinc...
            );
pit_dollar  = find_tax_liability( inc/modelunit_dollar, pit_thresholds, pit_burdens, pit_rates );
pit         = modelunit_dollar*pit_dollar;     % Convert to model units
% TBD:  INC should be in modelunit_dollars for consistency
inc         = inc/modelunit_dollar;

% Calculate Social Security tax from wage income in dollars
sst_dollar  = find_tax_liability( labinc/modelunit_dollar, sst_thresholds, sst_burdens, sst_rates );
sst         = modelunit_dollar*sst_dollar;

% Calculate corporate income tax
cit         = capshare*kv_ik*(taucap*(caprate - expsub)*captaxshare + taucapgain*capgain);

% Calculate available resources
resources   = (1 + totrate)*kv_ik + labinc + ssinc - (pit + sst + cit) + beq + kv_ik*capshare*capgain;

end




% Retirement age value function
function v = value_retirement(k, kv_, resources_, EV_ib_,... 
                    bequest_p_1_, bequest_phi_2_, bequest_phi_3_, ...
                    sigma_, gamma_)

% Enforce function inlining for C code generation
coder.inline('always');

% Define parameters as persistent variables
persistent kv resources EV_ib ...
            bequest_p_1 bequest_phi_2 bequest_phi_3 ...
            sigma gamma ...
            initialized

% Initialize parameters for C code generation
if isempty(initialized)
    kv = 0; resources = 0; EV_ib = 0; 
    bequest_p_1 = 0; bequest_phi_2 = 0; bequest_phi_3 = 0;
    sigma = 0; gamma = 0;
    initialized = true;
end

% Set parameters if provided
if (nargin > 1)
    kv = kv_; resources = resources_; EV_ib = EV_ib_; 
    bequest_p_1 = bequest_p_1_; bequest_phi_2 = bequest_phi_2_; bequest_phi_3 = bequest_phi_3_; 
    sigma = sigma_; gamma = gamma_;
    if isempty(k), return, end
end


% Calculate consumption
consumption = resources - k;

% Perform bound checks
if (kv(1) <= k) && (k <= kv(end)) && (0 <= consumption)
    
    % Residual value of bequest.
    % NOTE: (1) bequest is assets chosen for next period,
    %       (2) bequest_p_1 is beta*prob_death*bequest_phi_1
    value_bequest = bequest_p_1 * (1 + k/bequest_phi_2)^(1-bequest_phi_3);
    
    % Calculate utility
    v = (consumption^(gamma*(1-sigma)))/(1-sigma) ... % flow utility
        + interp1(kv, EV_ib, k, 'linear')         ... % continuation value of life
        + value_bequest                           ;   % value of bequest
    
    % Negate utility for minimization and force to scalar for C code generation
    v = -v(1);
    
else
    v = Inf;
end

end




% Working age value function
function v = value_working(x, kv_, bv_, wage_eff_, EV_, ...
                    bequest_p_1_, bequest_phi_2_, bequest_phi_3_, ...
                    sigma_, gamma_)

% Enforce function inlining for C code generation
coder.inline('always');

% Define parameters as persistent variables
persistent kv bv wage_eff EV ...
            bequest_p_1 bequest_phi_2 bequest_phi_3 ...
            sigma gamma ...
            initialized

% Initialize parameters for C code generation
if isempty(initialized)
    kv = 0; bv = 0; wage_eff = 0; EV = 0; 
    bequest_p_1 = 0; bequest_phi_2 = 0; bequest_phi_3 = 0;
    sigma = 0; gamma = 0;
    initialized = true;
end

% Set parameters if provided
if (nargin > 1)
    kv = kv_; bv = bv_; wage_eff = wage_eff_; EV = EV_; 
    bequest_p_1 = bequest_p_1_; bequest_phi_2 = bequest_phi_2_; bequest_phi_3 = bequest_phi_3_; 
    sigma = sigma_; gamma = gamma_;
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

% Residual value of bequest.
% NOTE: (1) bequest is assets chosen for next period,
%       (2) bequest_p_1 is beta*prob_death*bequest_phi_1
value_bequest = bequest_p_1 * (1 + k/bequest_phi_2)^(1-bequest_phi_3);
    
% Calculate utility
v = (((consumption^gamma)*((1-lab)^(1-gamma)))^(1-sigma))/(1-sigma) ... % flow utility
    + interp2(kv', bv, EV', k, b, 'linear')                         ... % continuation value of life
    + value_bequest                                                 ;   % value of bequest
 
% Negate utility for minimization and force to scalar for C code generation
v = -v(1);

end




% Average earnings calculation function
function [b] = calculate_b(labinc, age_, bv_ib_, bv_nb_, ssincmax_)

% Enforce function inlining for C code generation
coder.inline('always');

% Define parameters as persistent variables
persistent age bv_ib bv_nb ssincmax ...
           initialized

% Initialize parameters for C code generation
if isempty(initialized)
    age = 0; bv_ib = 0; bv_nb = 0; ssincmax = 0;
    initialized = true;
end

% Set parameters if provided
if (nargin > 1)
    age = age_; bv_ib = bv_ib_; bv_nb = bv_nb_; ssincmax = ssincmax_;
    if isempty(labinc), return, end
end


% Calculate average earnings and caps it to the maximum taxable earnings - bv(end)
b = min(bv_nb, (bv_ib*(age-1) + min(labinc, ssincmax))/age);

end


%%
%  Helper function to find tax liability from brackets & rates
%       NOTE:   income, burdens, & brackets are in US dollars; 
%           calculated liability is in also in US dollars.
function [tax] = find_tax_liability( income, brackets, burdens, rates )

    % Enforce function inlining for C code generation
    coder.inline('always');

    %       Expect equal-size vectors with brackets(1)=0
    %       rates apply for income between brackets(i-1) and brackets(i)
    %       burdens(i) are pre-calculated total tax liability at brackets(i)
    thebracket  = find(brackets <= income, 1, 'last');
    thebracket  = thebracket(1);   % Force to scalar for C code generation
    tax         = burdens(thebracket) + rates(thebracket)*(income - brackets(thebracket));

end % find_tax_liability



