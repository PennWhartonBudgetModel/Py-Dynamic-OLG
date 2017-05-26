%%
% Solve dynamic optimization problem for the firms.
% 
%%


function [W, OPT] = solve_firms(W0, K_static, DEBT_static, isdynamic, ...
                        nz, nk, nd, T_model, transz, kv, zv, dv, ...
                        tfp, alpha_k, alpha_n, capital_adjustment_param, depreciation, ...
                        profit_tax, investment_deduction, interest_deduction, ...
                        wages, discount_factors, debt_rates ) %#codegen


%% Argument verification
nz_max          = 50;
nk_max          = 50;
nd_max          = 50;
T_max           = 100;

assert( isa(W0          , 'double'  ) && (size(W0           , 1) <= nz_max  ) && (size(W0           , 2) <= nk_max  ) && (size(W0         , 3) <= nd_max  ));
assert( isa(K_static    , 'double'  ) && (size(K_static     , 1) <= nz_max  ) && (size(K_static     , 2) <= nk_max  ) && (size(K_static   , 3) <= nd_max  ) && (size(K_static   , 4) <= T_max   ) );
assert( isa(DEBT_static , 'double'  ) && (size(DEBT_static  , 1) <= nz_max  ) && (size(DEBT_static  , 2) <= nk_max  ) && (size(DEBT_static, 3) <= nd_max  ) && (size(DEBT_static, 4) <= T_max   ) );
assert( isa(isdynamic   , 'logical' ) && (size(isdynamic    , 1) == 1       ) && (size(isdynamic    , 2) == 1       ) );

assert( isa(nz          , 'double'  ) && (size(nz           , 1) == 1       ) && (size(nz           , 2) == 1       ) );
assert( isa(nk          , 'double'  ) && (size(nk           , 1) == 1       ) && (size(nk           , 2) == 1       ) );
assert( isa(nd          , 'double'  ) && (size(nd           , 1) == 1       ) && (size(nd           , 2) == 1       ) );
assert( isa(T_model     , 'double'  ) && (size(T_model      , 1) == 1       ) && (size(T_model      , 2) == 1       ) );
assert( isa(transz      , 'double'  ) && (size(transz       , 1) <= nz_max  ) && (size(transz       , 2) <= nz_max  ) );
assert( isa(kv          , 'double'  ) && (size(kv           , 1) <= nk_max  ) && (size(kv           , 2) == 1       ) );
assert( isa(zv          , 'double'  ) && (size(zv           , 1) <= nz_max  ) && (size(zv           , 2) == 1       ) );
assert( isa(dv          , 'double'  ) && (size(dv           , 1) <= nd_max  ) && (size(dv           , 2) == 1       ) );

assert( isa(tfp         , 'double'  ) && (size(tfp          , 1) == 1       ) && (size(tfp          , 2) == 1       ) );
assert( isa(alpha_k     , 'double'  ) && (size(alpha_k      , 1) == 1       ) && (size(alpha_k      , 2) == 1       ) );
assert( isa(alpha_n     , 'double'  ) && (size(alpha_n      , 1) == 1       ) && (size(alpha_n      , 2) == 1       ) );
assert( isa(capital_adjustment_cost, 'double' ) && (size(capital_adjustment_cost, 1) == 1 ) && (size(capital_adjustment_cost, 2) == 1 ) );
assert( isa(depreciation, 'double'  ) && (size(depreciation , 1) == 1       ) && (size(depreciation , 2) == 1       ) );

assert( isa(profit_tax          , 'double'  ) && (size(profit_tax          , 1) == 1 ) && (size(profit_tax          , 2) == 1 ) );
assert( isa(investment_deduction, 'double'  ) && (size(investment_deduction, 1) == 1 ) && (size(investment_deduction, 2) == 1 ) );
assert( isa(interest_deduction  , 'double'  ) && (size(interest_deduction  , 1) == 1 ) && (size(interest_deduction  , 2) == 1 ) );

assert( isa(wages           , 'double'  ) && (size(wages           , 1) == 1       ) && (size(wages           , 2) <= T_max   ) );
assert( isa(discount_factors, 'double'  ) && (size(discount_factors, 1) == 1       ) && (size(discount_factors, 2) <= T_max   ) );
% TODO: These should be firm state, not global state
assert( isa(debt_rates      , 'double'  ) && (size(debt_rates      , 1) == 1       ) && (size(debt_rates      , 2) <= T_max   ) );

%% Dynamic optimization

% Initialize valuation and optimal decision value arrays
W        = zeros(nz,nk,nd,T_model);      % Corp value

OPT.K    = zeros(nz,nk,nd,T_model);      % Kapital
OPT.DEBT = zeros(nz,nk,nd,T_model);      % Debt outstanding

OPT.LAB  = zeros(nz,nk,nd,T_model);      % Labor hired
OPT.INC  = zeros(nz,nk,nd,T_model);      % Corp distribution
OPT.CIT  = zeros(nz,nk,nd,T_model);      % Corporate income tax

% Initialize forward-looking corp values
W_step = W0;

% Specify settings for dynamic optimization subproblems
optim_options = optimset('Display', 'off', 'TolFun', 1e-4, 'TolX', 1e-4);

% Use dynamic or static method to find choice vars (capital,debt)
% We use indirection (function pointer) to avoid an "if" in the loop
if dynamic
    find_capital_debt = @(iz, ik, id, t) fminsearch(@value_enterprise, [kv(ik), dv(id)], optim_options);
else
    find_capital_debt = @(iz, ik, id, t) deal( [K_static(iz,ik,id,t) DEBT_static(iz,ik,id,t)], ...
                              value_enterprise(K_static(iz,ik,id,t), DEBT_static(iz,ik,id,t)) );
end 

% Solve dynamic optimization problem through backward induction
for t = T_model:-1:1
    
    year = t; 
    
    % Extract parameters for current year
    wage            = wages             (year);
    discount_factor = discount_factors  (year);
    debt_rate       = debt_rates        (year);

    for ik = 1:nk
        for id = 1:nd
            for iz = 1:nz
                
                % Calculate interest payment due
                %   TBD: The rate is actually a firm state, NOT a global
                %   state
                debt_coupon = debt_rate*dv(id);
                
                % Calculate expected value curve using forward-looking values
                % TODO: Fix this calculation
                EW = discount_factor*reshape(sum(repmat(transz(iz,:)', [1,nk,nd]) .* W_step, 1), [nk, nd]);

                % Call functions to set persistent parameters
                value_enterprise([], kv, kd, EW);
                calculate_income( 0, 0, ...
                    kv(ik), zv(iz), dv(id), debt_coupon,...
                    tfp, alpha_k, alpha_n, capital_adjustment_param, depreciation, ...
                    profit_tax, investment_deduction, interest_deduction, ...
                    wage );
                
                % Find choice vars (capital, debt)
                [x, v] = find_capital_debt(iz,ik,id,t);
                
% BEGIN: TEMP -- keep until indirection approach verified
%                 if isdynamic
%                     [x, v]  = fminsearch(@value_enterprise, [kv(ik), dv(id)], optim_options);
%                 else
%                     x(1)    = K_static    (iz,ik,id,t);
%                     x(2)    = DEBT_static (iz,ik,id,t);
%                     v       = value_enterprise(x(1), x(2));
%                 end % isdynamic  
% END: TEMP                 

                % Record enterprise value and optimal decision values
                W       (iz,ik,id,t) = -v;
                OPT.K   (iz,ik,id,t) = x(1);  % capital
                OPT.DEBT(iz,ik,id,t) = x(2);  % debt
                
                % Find derived vars from choice vars
                [labor, income, cit] = calculate_income( capital, debt );
                OPT.LAB(iz,ik,t) = labor;
                OPT.INC(iz,ik,t) = income;
                OPT.CIT(iz,ik,t) = cit;
            end % iz
        end % id
    end % ik
    
    % Update forward-looking enterprise values
    W_step = W(:, :, :, t);
    
end % t


end % solve_firms





%%
% Value function: State is (capital, debt)
function v = value_enterprise(x, kv_, dv_, EW_)

    % Enforce function inlining for C code generation
    coder.inline('always');

    % Define parameters as persistent variables
    persistent kv dv EW ...
               initialized

    % Initialize parameters for C code generation
    if isempty(initialized)
        kv = 0; dv = 0; EW = 0; 
        initialized = true;
    end

    % Set parameters if provided
    if (nargin > 1)
        kv = kv_; dv = dv_; EW = EW_; 
        
        % Assume if initializing, then do not continue
        return;
    end

    new_capital = x(1);
    new_debt    = x(2);
    
    [~, income, ~] = calculate_income( new_capital, new_debt );
    
    % Calculate enterprise value
    v = discount_factor*interp2(kv', dv', EW', new_capital, new_debt, 'linear') ... 
        + income;

    % Negate value for minimization and force to scalar for C code generation
    v = -v(1);

end % value_enterprise


%%
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
                labor_cost revenue ...
                initialized;

    % Initialize parameters for C code generation
    if isempty(initialized)
        capital = 0; shock = 0; debt = 0; debt_coupon = 0;
        tfp = 0; alpha_k = 0; alpha_n = 0; capital_adjustment_param = 0; depreciation = 0;
        profit_tax = 0; investment_deduction = 0; interest_deduction = 0;
        wage = 0;
        labor_cost = 0; revenue = 0;
        initialized = true;
    end

    % Set parameters if provided
    if (nargin > 2)
        capital = capital_; shock = shock_; debt = debt_; debt_coupon = debt_coupon_;
        tfp = tfp_; alpha_k = alpha_k_; alpha_n = alpha_n_; capital_adjustment_param = capital_adjustment_param_; depreciation = depreciation_;
        profit_tax = profit_tax_; investment_deduction = investment_deduction_; interest_deduction = interest_deduction_;
        wage = wage_;
        
        % Calculate derived vars: revenue and labor_cost
        %   Optimal labor hiring choice (from FOC)
        labor                   = (wage/(alpha_n*shock*(capital^alpha_k)))^(1/(1+alpha_n));
        revenue                 = tfp*(capital^alpha_k)*(labor^alpha_n);
        labor_cost              = labor*wage;

        % Assume that if initializing, we do not calculate
        return;
    end % initialize
    
    net_investment          = new_capital - capital;
    investment              = net_investment + capital*depreciation;
    capital_adjustment_cost = (net_investment/capital)*net_investment*capital_adjustment_param;
    pretax_revenue          = revenue  ...
                                - labor_cost ...
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





