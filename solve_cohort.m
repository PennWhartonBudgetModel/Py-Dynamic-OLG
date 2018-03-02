%%
% Solve dynamic optimization problem for a cohort.
% 
%%

function [OPT] ...
    = solve_cohort( ...
        V0, LAB_static, isdynamic, ...
        nz, ns, nb, T_past, T_shift, T_active, T_work, T_model, ... 
        zs, transz, sv, bv, beta, gamma, sigma, surv, ...
        bequest_phi_1, bequest_phi_2, bequest_phi_3, ...
        ssbenefits, ssincmins, ssincmaxs, sswageindexes, ...
        sstax_brackets, sstax_burdens, sstax_rates, ...
        pittax_sscredit, pittax_brackets, pittax_burdens, pittax_rates, ... 
        captax_shares, captax_pref_rates, capgain_taxrates, ...
        beqs, ...
        wages, ...
        fund_equityshares, ...
        equityfund_dividendrates, equityfund_prices, ...
        bondfund_dividendrates, bondfund_prices ...
    ) %#codegen


%% Argument verification

nz_max          = 50;
ns_max          = 50;
nb_max          = 50;
T_max           = 100;
nbrackets_max   = 20;

assert( isa(V0          , 'double'  ) && (size(V0           , 1) <= nz_max  ) && (size(V0           , 2) <= ns_max  ) && (size(V0           , 3) <= nb_max  ) );
assert( isa(LAB_static  , 'double'  ) && (size(LAB_static   , 1) <= nz_max  ) && (size(LAB_static   , 2) <= ns_max  ) && (size(LAB_static   , 3) <= nb_max  ) && (size(LAB_static   , 4) <= T_max   ) );
assert( isa(isdynamic   , 'logical' ) && (size(isdynamic    , 1) == 1       ) && (size(isdynamic    , 2) == 1       ) );

assert( isa(nz          , 'double'  ) && (size(nz           , 1) == 1       ) && (size(nz           , 2) == 1       ) );
assert( isa(ns          , 'double'  ) && (size(ns           , 1) == 1       ) && (size(ns           , 2) == 1       ) );
assert( isa(nb          , 'double'  ) && (size(nb           , 1) == 1       ) && (size(nb           , 2) == 1       ) );
assert( isa(T_past      , 'double'  ) && (size(T_past       , 1) == 1       ) && (size(T_past       , 2) == 1       ) );
assert( isa(T_shift     , 'double'  ) && (size(T_shift      , 1) == 1       ) && (size(T_shift      , 2) == 1       ) );
assert( isa(T_active    , 'double'  ) && (size(T_active     , 1) == 1       ) && (size(T_active     , 2) == 1       ) );
assert( isa(T_work      , 'double'  ) && (size(T_work       , 1) == 1       ) && (size(T_work       , 2) == 1       ) );
assert( isa(T_model     , 'double'  ) && (size(T_model      , 1) == 1       ) && (size(T_model      , 2) == 1       ) );
assert( isa(zs          , 'double'  ) && (size(zs           , 1) <= nz_max  ) && (size(zs           , 2) <= T_max   ) );
assert( isa(transz      , 'double'  ) && (size(transz       , 1) <= nz_max  ) && (size(transz       , 2) <= nz_max  ) && (size(transz        , 3) <= T_max  ) );
assert( isa(sv          , 'double'  ) && (size(sv           , 1) <= ns_max  ) && (size(sv           , 2) == 1       ) );
assert( isa(bv          , 'double'  ) && (size(bv           , 1) <= nb_max  ) && (size(bv           , 2) == 1       ) );
assert( isa(beta        , 'double'  ) && (size(beta         , 1) == 1       ) && (size(beta         , 2) == 1       ) );
assert( isa(gamma       , 'double'  ) && (size(gamma        , 1) == 1       ) && (size(gamma        , 2) == 1       ) );
assert( isa(sigma       , 'double'  ) && (size(sigma        , 1) == 1       ) && (size(sigma        , 2) == 1       ) );
assert( isa(surv        , 'double'  ) && (size(surv         , 1) == 1       ) && (size(surv         , 2) <= T_max   ) );

assert( isa(bequest_phi_1, 'double' ) && (size(bequest_phi_1, 1) == 1 ) && (size(bequest_phi_1, 2) == 1 ) );
assert( isa(bequest_phi_2, 'double' ) && (size(bequest_phi_2, 1) == 1 ) && (size(bequest_phi_2, 2) == 1 ) );
assert( isa(bequest_phi_3, 'double' ) && (size(bequest_phi_3, 1) == 1 ) && (size(bequest_phi_3, 2) == 1 ) );

assert( isa(ssincmins    , 'double'  ) && (size(ssincmins    , 1) <= T_max   ) && (size(ssincmins    , 2) == 1       ) );
assert( isa(ssincmaxs    , 'double'  ) && (size(ssincmaxs    , 1) <= T_max   ) && (size(ssincmaxs    , 2) == 1       ) );
assert( isa(sswageindexes, 'double'  ) && (size(sswageindexes, 1) <= T_max   ) && (size(sswageindexes, 2) == 1       ) );

assert( isa(ssbenefits      , 'double' ) && (size(ssbenefits      , 1) <= T_max ) && (size(ssbenefits      , 2) <= nb_max        ) );

assert( isa(sstax_brackets  , 'double' ) && (size(sstax_brackets  , 1) <= T_max ) && (size(sstax_brackets  , 2) <= nbrackets_max ) );
assert( isa(sstax_burdens   , 'double' ) && (size(sstax_burdens   , 1) <= T_max ) && (size(sstax_burdens   , 2) <= nbrackets_max ) );
assert( isa(sstax_rates     , 'double' ) && (size(sstax_rates     , 1) <= T_max ) && (size(sstax_rates     , 2) <= nbrackets_max ) );

assert( isa(pittax_sscredit  , 'double' ) && (size(pittax_sscredit  , 1) == 1     ) && (size(pittax_sscredit  , 2) == 1             ) );
assert( isa(pittax_brackets  , 'double' ) && (size(pittax_brackets  , 1) <= T_max ) && (size(pittax_brackets  , 2) <= nbrackets_max ) );
assert( isa(pittax_burdens   , 'double' ) && (size(pittax_burdens   , 1) <= T_max ) && (size(pittax_burdens   , 2) <= nbrackets_max ) );
assert( isa(pittax_rates     , 'double' ) && (size(pittax_rates     , 1) <= T_max ) && (size(pittax_rates     , 2) <= nbrackets_max ) );

assert( isa(captax_shares     , 'double'  ) && (size(captax_shares     , 1) <= T_max   ) && (size(captax_shares     , 2) == 1       ) );
assert( isa(captax_pref_rates , 'double'  ) && (size(captax_pref_rates , 1) <= T_max   ) && (size(captax_pref_rates , 2) == 1       ) );
assert( isa(capgain_taxrates  , 'double'  ) && (size(capgain_taxrates  , 1) <= T_max   ) && (size(capgain_taxrates  , 2) == 1       ) );

assert( isa(beqs    , 'double'  ) && (size(beqs     , 1) == 1       ) && (size(beqs     , 2) <= T_max   ) );
assert( isa(wages   , 'double'  ) && (size(wages    , 1) == 1       ) && (size(wages    , 2) <= T_max   ) );

assert( isa(fund_equityshares        , 'double'  ) && (size(fund_equityshares        , 1) == 1       ) && (size(fund_equityshares        , 2) <= T_max   ) );
assert( isa(equityfund_dividendrates , 'double'  ) && (size(equityfund_dividendrates , 1) == 1       ) && (size(equityfund_dividendrates , 2) <= T_max   ) );
assert( isa(equityfund_prices        , 'double'  ) && (size(equityfund_prices        , 1) == 1       ) && (size(equityfund_prices        , 2) <= T_max   ) );
assert( isa(bondfund_dividendrates   , 'double'  ) && (size(bondfund_dividendrates   , 1) == 1       ) && (size(bondfund_dividendrates   , 2) <= T_max   ) );
assert( isa(bondfund_prices          , 'double'  ) && (size(bondfund_prices          , 1) == 1       ) && (size(bondfund_prices          , 2) <= T_max   ) );


%% Dynamic optimization

% Initialize utility and optimal decision value arrays
OPT.V   = zeros(nz,ns,nb,T_active);   % Utility

OPT.K   = zeros(nz,ns,nb,T_active);   % Savings
OPT.LAB = zeros(nz,ns,nb,T_active);   % Labor level
OPT.B   = zeros(nz,ns,nb,T_active);   % Average earnings

OPT.INC = zeros(nz,ns,nb,T_active);   % Taxable income
OPT.PIT = zeros(nz,ns,nb,T_active);   % Personal income tax
OPT.SST = zeros(nz,ns,nb,T_active);   % Social Security tax
OPT.CIT = zeros(nz,ns,nb,T_active);   % Capital income tax
OPT.BEN = zeros(nz,ns,nb,T_active);   % Social Security benefits
OPT.CON = zeros(nz,ns,nb,T_active);   % Consumption

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
    ssincmax    = ssincmaxs    (year);
    ssincmin    = ssincmins    (year);
    sswageindex = 1; % TBD: substitute by sswageindexes(year);
    ssbenefit   = ssbenefits   (year, :);
    
    beq         = beqs                  (year);
    capshare    = fund_equityshares     (year);
    
    % Prices for current year
    govrate          = bondfund_dividendrates   (year);
    wage             = wages                    (year);
    caprate          = equityfund_dividendrates (year);
    equityfund_price = equityfund_prices        (year);
    bondfund_price   = 1; % TBD: This should be passed in
    
    % THIS IS TEMP
    equityfund_prevprice = equityfund_prices(max(year-1),1);
    capgain_share = capshare*((equityfund_price - equityfund_prevprice)/equityfund_price);
    sstaxcredit = pittax_sscredit;
    
    % Taxes
    capgain_taxrate     = capgain_taxrates      (year);
    captax_share        = captax_shares         (year);
    captax_pref_rate    = captax_pref_rates     (year);
   
    sst_brackets    = sstax_brackets    (year, :);
    sst_burdens     = sstax_burdens     (year, :);
    sst_rates       = sstax_rates       (year, :);
    
    pit_brackets    = pittax_brackets   (year, :);
    pit_burdens     = pittax_burdens    (year, :);
    pit_rates       = pittax_rates      (year, :);
        
    
    % Pre-calculate for speed and conciseness
    bequest_p_1   = beta * (1-surv(age))* bequest_phi_1;
    bequest_p_2   = bequest_phi_2;
    bequest_p_3   = bequest_phi_3;
    
    for ib = 1:nb
        for is = 1:ns
            
            if (age > T_work)
                
                % Calculate available resources and tax terms
                ssinc = ssbenefit(ib);
                
                [resources, inc, pit, ~, cit] = calculate_resources( ...
                    0, ...
                    sv(is), equityfund_price, bondfund_price, ...
                    ssinc, sstaxcredit, ...
                    sst_brackets, sst_burdens, sst_rates, ...
                    pit_brackets, pit_burdens, pit_rates, ... 
                    captax_share, captax_pref_rate, capgain_share, ...
                    beq, capshare, caprate, govrate ...
                );
                
                if isdynamic
                    
                    % Calculate expected value conditional on living using forward-looking 
                    %   utility values. Pre-multiply by prob. survival and
                    %   beta to save on computation.
                    EV = surv(age)*beta*reshape(V_step(1,:,:), [ns,nb]);
                    
                    % Solve dynamic optimization subproblem
                    [s, v] = fminsearch( ...
                        @(s) value_retirement( ...
                            s, sv, resources, EV(:,ib), ... 
                            bequest_p_1, bequest_p_2, bequest_p_3, ...
                            sigma, gamma ...
                        ), sv(is), optim_options ...
                    );
                    
                    % Checks -> only work in the absence of mex file!
                    assert( ~isinf(v)   , 'v is inf')
                    assert( s <= sv(end), 's (s_next) is too big!')

                    % Record utility and optimal decision values
                    OPT.V(:,is,ib,t) = -v;
                    OPT.K(:,is,ib,t) = s ;
                    
                else
                    
                    s = sv(is);
                    
                end
                
                OPT.LAB(:,is,ib,t) = 0     ;
                OPT.B  (:,is,ib,t) = bv(ib);
                
                OPT.INC(:,is,ib,t) = inc   ;
                OPT.PIT(:,is,ib,t) = pit   ;
                OPT.SST(:,is,ib,t) = 0     ;
                OPT.CIT(:,is,ib,t) = cit   ;
                OPT.BEN(:,is,ib,t) = ssinc ;
                OPT.CON(:,is,ib,t) = resources - s;
                
            else
                
                % Create local instance of average earnings calculation function with fixed parameters
                calculate_b_ = @(labinc) calculate_b( ...
                    labinc, age, bv(ib), ...
                    ssincmin, ssincmax, sswageindex ...
                );
                
                % Create local instance of resource calculation function with fixed parameters
                calculate_resources_ = @(labinc) calculate_resources( ...
                    labinc, ...
                    sv(is), equityfund_price, bondfund_price, ...
                    0, 0, ...
                    sst_brackets, sst_burdens, sst_rates, ...
                    pit_brackets, pit_burdens, pit_rates, ... 
                    captax_share, captax_pref_rate, capgain_share, ...
                    beq, capshare, caprate, govrate ...
                );
                
                for iz = 1:nz
                    
                    % Calculate effective wage
                    wage_eff = wage * zs(iz,age);
                    
                    if isdynamic
                        
                        % Calculate expected value conditional on living using forward-looking 
                        %   utility values. Pre-multiply by prob. survival and
                        %   beta to save on computation.
                        EV = surv(age)*beta*reshape(sum(repmat(transz(iz,:,age)', [1,ns,nb]) .* V_step, 1), [ns,nb]);
                        
                        % Solve dynamic optimization subproblem
                        lab0 = 0.5;
                        s0   = max(sv(is), min(sv(end), 0.1 * wage_eff * lab0));   % Assumes taxation will not exceed 90% of labor income and at the same time forces k to be in the grid
                        
                        [x, v] = fminsearch( ...
                            @(x) value_working( ...
                                x, sv, bv, wage_eff, EV, ...
                                bequest_p_1, bequest_p_2, bequest_p_3, ...
                                sigma, gamma, ...
                                calculate_b_, calculate_resources_ ...
                            ), [s0, lab0], optim_options ...
                        );
                        
                        s   = x(1);
                        lab = x(2);       
                        
                        % Checks -> only work in the absence of mex file!
                        assert( ~isinf(v)   , 'v is inf')
                        assert( s <= sv(end), 'k (k_next) is too big!')
                        
                        % Record utility and optimal decision values
                        OPT.V(iz,is,ib,t) = -v;
                        OPT.K(iz,is,ib,t) = s ;
                        
                    else
                        s = sv(is);
                        lab = LAB_static(iz,is,ib,t);
                    end
                    
                    labinc = wage_eff * lab;
                    [resources, inc, pit, sst, cit] = calculate_resources_(labinc);
                    
                    OPT.LAB(iz,is,ib,t) = lab;
                    OPT.B  (iz,is,ib,t) = calculate_b_(labinc);
                    
                    OPT.INC(iz,is,ib,t) = inc;
                    OPT.PIT(iz,is,ib,t) = pit;
                    OPT.SST(iz,is,ib,t) = sst;
                    OPT.CIT(iz,is,ib,t) = cit;
                    OPT.BEN(iz,is,ib,t) = 0  ;
                    OPT.CON(iz,is,ib,t) = resources - s;
                    
                end
                
            end
            
        end
    end
    
    % Update forward-looking utility values
    V_step = OPT.V(:,:,:,t);
    
end


end


% Retirement age value function
function [v] ...
    = value_retirement( ...
        k, kv, resources, EV_ib, ... 
        bequest_p_1, bequest_phi_2, bequest_phi_3, ...
        sigma, gamma ...
    )

    % Enforce function inlining for C code generation
    coder.inline('always');

    % Calculate consumption
    consumption = resources - k;

    % Perform bound checks
    if (kv(1) <= k) && (k <= kv(end)) && (0 <= consumption)

        % Residual value of bequest.
        % NOTE: (1) bequest is assets chosen for next period,
        %       (2) bequest_p_1 is beta*prob_death*bequest_phi_1
        value_bequest = bequest_p_1 * (1 + k/bequest_phi_2)^(1-bequest_phi_3);

        % Calculate utility
        v = (consumption^(gamma*(1-sigma)))*(1/(1-sigma))     ... % flow utility
            + interp1(kv, EV_ib, k, 'linear')                 ... % continuation value of life
            + value_bequest                                   ;   % value of bequest

        % Negate utility for minimization and force to scalar for C code generation
        v = -v(1);

    else
        v = Inf;
    end

end


% Working age value function
function [v] ...
    = value_working( ...
        x, kv, bv, wage_eff, EV, ...
        bequest_p_1, bequest_p_2, bequest_p_3, ...
        sigma, gamma, ...
        calculate_b_, calculate_resources_ ...
    )

    % Enforce function inlining for C code generation
    coder.inline('always');

    % Define decision variables and perform bound checks
    k   = x(1);
    lab = x(2);

    labinc = wage_eff * lab;

    if ~((0 <= lab) && (lab <= 1) ...
         && (kv(1) <= k) && (k <= kv(end)))

        v = Inf;
        return

    end

    b = calculate_b_(labinc);

    % Calculate available resources
    resources = calculate_resources_(labinc);

    % Calculate consumption and perform bound check
    consumption = resources - k;

    if ~(0 <= consumption)
        v = Inf;
        return
    end

    % Residual value of bequest.
    % NOTE: (1) bequest is assets chosen for next period,
    %       (2) bequest_p_1 is beta*prob_death*bequest_phi_1
    value_bequest = bequest_p_1 * (1 + k/bequest_p_2)^(1-bequest_p_3);

    % Calculate utility
    v = (((consumption^gamma)*((1-lab)^(1-gamma)))^(1-sigma))*(1/(1-sigma))     ... % flow utility
        + interp2(kv', bv, EV', k, b, 'linear')                                 ... % continuation value of life
        + value_bequest                                                         ;   % value of bequest

    % Negate utility for minimization and force to scalar for C code generation
    v = -v(1);

end


% Average earnings calculation function
function [b] ...
    = calculate_b( ...
        labinc, age, bv_ib, ...
        ssincmin, ssincmax, sswageindex ...
    )

    % Enforce function inlining for C code generation
    coder.inline('always');

    % Calculate average earnings, cap them to a policy ceiling, and deflate
    % them at each period using Market.priceindices.cohort_wages(:,i) = sswageindex
    % for household of cohort i. 
    if labinc > (ssincmin + 10*eps)
        b = (bv_ib*(age-1) + sswageindex*min(labinc, ssincmax)) / age - 10*eps;
    else
        b = bv_ib;
    end

end


% Resource and tax calculation function
function [resources, inc, pit, sst, cit] ...
    = calculate_resources( ...
        labinc, ...
        fund_shares, equityfund_price, bondfund_price, ...
        ssinc, sstaxcredit, ...
        sst_brackets, sst_burdens, sst_rates, ...
        pit_brackets, pit_burdens, pit_rates, ... 
        captax_share, captax_pref_rate, capgain_share, ...
        beq, capshare, caprate, govrate ...
    )

    % Enforce function inlining for C code generation
    coder.inline('always');

    % Calculate wealth in equity and g'vt bonds
    equityvalue = fund_shares * capshare * equityfund_price;
    bondvalue   = fund_shares * (1-capshare) * bondfund_price;
    
    equitydividend = caprate * equityvalue;
    bonddividend   = govrate * bondvalue;
    
    % TBD: This needs to be revised
    equitycapgain  = fund_shares * equityfund_price * capgain_share;
    
    % Calculate taxable income
    %   We do not allow negative incomes
    inc = max( 0,...
          equitydividend*(1-captax_share) + bonddividend ...
          + (1-sstaxcredit)*ssinc + labinc ...
          );
    pit = find_tax_liability( inc, pit_brackets, pit_burdens, pit_rates );

    % Calculate Social Security tax from wage income
    sst = find_tax_liability( labinc, sst_brackets, sst_burdens, sst_rates );

    % Calculate preferred rates tax 
    cit = captax_share * equitydividend * captax_pref_rate; % NOTE: capgains taxed here also

    % Calculate available resources
    resources = equityvalue + bondvalue ...
                + equitydividend + bonddividend ...
                + labinc + ssinc ...
                - (pit + sst + cit) ...
                + beq + equitycapgain;

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

