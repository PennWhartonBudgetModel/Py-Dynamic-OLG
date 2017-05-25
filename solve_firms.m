%%
% Solve dynamic optimization problem for the firms.
% 
%%


function [W, OPT] = solve_firms(W0, K_static, isdynamic, ...
                        nz, nk, T_model, transz, kv, zv,  ...
                        tfp, ...
                        modelunit_dollars, ...
                        profit_tax, investment_deduction, interest_deduction, ...
                        wages, discount_factors, govrates ) %#codegen


%% Argument verification
% TODO: FIX THIS SECTION
nz_max          = 50;
nk_max          = 50;
T_max           = 100;

assert( isa(W0          , 'double'  ) && (size(W0           , 1) <= nz_max  ) && (size(W0           , 2) <= nk_max  ) );
assert( isa(K_static    , 'double'  ) && (size(K_static     , 1) <= nz_max  ) && (size(K_static     , 2) <= nk_max  ) && (size(K_static   , 3) <= T_max   ) );
assert( isa(isdynamic   , 'logical' ) && (size(isdynamic    , 1) == 1       ) && (size(isdynamic    , 2) == 1       ) );

assert( isa(nz          , 'double'  ) && (size(nz           , 1) == 1       ) && (size(nz           , 2) == 1       ) );
assert( isa(nk          , 'double'  ) && (size(nk           , 1) == 1       ) && (size(nk           , 2) == 1       ) );
assert( isa(T_model     , 'double'  ) && (size(T_model      , 1) == 1       ) && (size(T_model      , 2) == 1       ) );
assert( isa(transz      , 'double'  ) && (size(transz       , 1) <= nz_max  ) && (size(transz       , 2) <= nz_max  ) );
assert( isa(kv          , 'double'  ) && (size(kv           , 1) <= nk_max  ) && (size(kv           , 2) == 1       ) );
assert( isa(zv          , 'double'  ) && (size(zv           , 1) <= nz_max  ) && (size(zv           , 2) == 1       ) );

assert( isa(modelunit_dollars, 'double'  ) && (size(modelunit_dollars, 1) == 1       ) && (size(modelunit_dollars, 2) == 1       ) );


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
% END TODO


%% Dynamic optimization

% Initialize valuation and optimal decision value arrays
W        = zeros(nz, nk, T_model);      % Corp value

OPT.K    = zeros(nz, nk, T_model);      % Kapital
OPT.LAB  = zeros(nz, nk, T_model);      % Labor hired

OPT.INC  = zeros(nz, nk, T_model);      % Taxable income
OPT.CIT  = zeros(nz, nk, T_model);      % Corporate income tax
OPT.BOND = zeros(nz, nk, T_model);      % Debt outstanding

% Initialize forward-looking corp values
W_step = W0;

% Specify settings for dynamic optimization subproblems
optim_options = optimset('Display', 'off', 'TolFun', 1e-4, 'TolX', 1e-4);


% Solve dynamic optimization problem through backward induction
for t = T_model:-1:1
    
    year = t; 
    
    % Extract parameters for current year
    wage            = wages             (year);
    discount_factor = discount_factors  (year);

    for ik = 1:nk

        for iz = 1:nz

                if isdynamic

                    % Calculate expected value curve using forward-looking values
                    EW = discount_factor*reshape(sum(repmat(transz(iz,:)', [1,nk]) .* V_step, 1), [nk]);

                    % Call value function to set persistent parameters
                    value_enterprise([], kv, wage_eff, EW)

                    % Solve dynamic optimization subproblem
                    k0   = kv(ik);

                    [x, v] = fminsearch(@value_enterprise, [k0], optim_options);

                    capital         = x(1);

                    % Record utility and optimal decision values
                    W    (iz,ik,t)  = -v;
                    OPT.K(iz,ik,t)  = capital ;
                end

                [labor, income, cit] = calculate_income( ...
                        capital, shock, debt, investment, new_debt, ...
                        tfp, alpha_k, alpha_n, capital_adjustment_param, ...
                        modelunit_dollars, ...
                        profit_tax, investment_deduction, interest_deduction, ...
                        wage );
                OPT.LAB(iz,ik,t) = labor;
                OPT.INC(iz,ik,t) = income;
                OPT.CIT(iz,ik,t) = cit;

        end % iz

    end % ik
    
    % Update forward-looking corp values
    W_step = W(:, :, t);
    
end % t


end % solve_firms




% Calculates firm's budget & taxes
function [labor, income, cit] = calculate_income( ...
            new_capital, new_debt, ...
            capital_, shock_, debt_, debt_coupon_,...
            tfp_, alpha_k_, alpha_n_, capital_adjustment_param_, depreciation_, ...
            profit_tax_, investment_deduction_, interest_deduction_, ...
            wage_ ) %#codegen

    % Enforce function inlining for C code generation
    coder.inline('always');

    % Define parameters as persistent variables
    persistent  capital shock debt debt_coupon ...
                tfp alpha_k alpha_n capital_adjustment_param depreciation ...
                profit_tax investment_deduction interest_deduction ...
                wage ...
                initialized;

    % Initialize parameters for C code generation
    if isempty(initialized)
        capital = 0; shock = 0; debt = 0; debt_coupon = 0;
        tfp = 0; alpha_k = 0; alpha_n = 0; capital_adjustment_param = 0; depreciation = 0;
        profit_tax = 0; investment_deduction = 0; interest_deduction = 0;
        wage = 0;
        initialized = true;
        if isempty(new_capital_), return, end
    end

    % Set parameters if provided
    if (nargin > 2)
        capital = capital_; shock = shock_; debt = debt_; debt_coupon = debt_coupon_;
        tfp = tfp_; alpha_k = alpha_k_; alpha_n = alpha_n_; capital_adjustment_param = capital_adjustment_param_; depreciation = depreciation_;
        profit_tax = profit_tax_; investment_deduction = investment_deduction_; interest_deduction = interest_deduction_;
        wage = wage_;
    end

    % Helper vars
    k_to_alpha              = capital^alpha_k;

    % Calculate optimal labor hiring choice (from FOC)
    labor                   = (wage/(alpha_n*shock*k_to_alpha))^(1/(1+alpha_n));

    net_investment          = new_capital - capital;
    investment              = net_investment + capital*depreciation;
    revenue                 = tfp*k_to_alpha*(labor^alpha_n);
    capital_adjustment_cost = (net_investment/capital)*net_investment*capital_adjustment_param;
    pretax_revenue          = revenue  ...
                                - labor*wage ...
                                - capital_adjustment_cost ...
                                - investment ...
                                - debt_coupon;

    % Calculate tax with relevant deductions
    cit                     = profit_tax*(pretax_revenue ...
                                + (1-investment_deduction)*investment ...
                                + (1-interest_deduction)*debt_coupon  );
    % After-tax income includes $ from debt issuance
    income                  = pretax_revenue - cit + (new_debt - debt);

end % calculate_income




%%
% Value function
function v = value_enterprise(x, kv_, income_, EW_)

    % Enforce function inlining for C code generation
    coder.inline('always');

    % Define parameters as persistent variables
    persistent kv income EW ...
               initialized

    % Initialize parameters for C code generation
    if isempty(initialized)
        kv = 0; income = 0; EW = 0; 
        initialized = true;
    end

    % Set parameters if provided
    if (nargin > 1)
        kv = kv_; income = income_; EW = EW_; 
        if isempty(x), return, end
    end

    % Calculate enterprise value
    v = discount_factor*interp2(kv', EW', capital, 'linear') ... 
        + income;

    % Negate value for minimization and force to scalar for C code generation
    v = -v(1);

end % value_enterprise





