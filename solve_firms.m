%%
% Solve dynamic optimization problem for the firms.
% 
%%


function OPT = solve_firms(W0, Solvent0, ...
                        K_static, DEBT_static, isdynamic, ...
                        nz, nk, nd, T_model, transz, kv, zv, dv, ...
                        tfp, alpha_k, alpha_n, k_adjustment_param, depreciation, k_default_reduction, ...
                        profit_tax, investment_deduction, interest_deduction, ...
                        wages, discount_factors ) %#codegen


%% Argument verification
nz_max          = 50;
nk_max          = 50;
nd_max          = 50;
T_max           = 100;

assert( isa(W0          , 'double'  ) && (size(W0           , 1) <= nz_max  ) && (size(W0           , 2) <= nk_max  ) && (size(W0         , 3) <= nd_max  ));
assert( isa(Solvent0    , 'double'  ) && (size(Solvent0     , 1) <= nz_max  ) && (size(Solvent0     , 2) <= nk_max  ) && (size(Solvent0   , 3) <= nd_max  ));

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
assert( isa(k_adjustment_param, 'double' ) && (size(k_adjustment_param, 1) == 1 ) && (size(k_adjustment_param, 2) == 1 ) );
assert( isa(depreciation, 'double'  ) && (size(depreciation , 1) == 1       ) && (size(depreciation , 2) == 1       ) );
assert( isa(k_default_reduction, 'double') && (size(k_default_reduction,1) == 1 ) && (size(k_default_reduction, 2) == 1 ) );

assert( isa(profit_tax          , 'double'  ) && (size(profit_tax          , 1) == 1 ) && (size(profit_tax          , 2) == 1 ) );
assert( isa(investment_deduction, 'double'  ) && (size(investment_deduction, 1) == 1 ) && (size(investment_deduction, 2) == 1 ) );
assert( isa(interest_deduction  , 'double'  ) && (size(interest_deduction  , 1) == 1 ) && (size(interest_deduction  , 2) == 1 ) );

assert( isa(wages           , 'double'  ) && (size(wages           , 1) == 1       ) && (size(wages           , 2) <= T_max   ) );
assert( isa(discount_factors, 'double'  ) && (size(discount_factors, 1) == 1       ) && (size(discount_factors, 2) <= T_max   ) );

%% Dynamic optimization

% Initialize valuation and optimal decision value arrays
OPT.W       = zeros(nz,nk,nd,T_model);      % Corp value

OPT.K       = zeros(nz,nk,nd,T_model);      % Kapital
OPT.DEBT    = zeros(nz,nk,nd,T_model);      % Debt outstanding
OPT.SOLVENT = zeros(nz,nk,nd,T_model);      % Defaulted on debt?

OPT.LAB    = zeros(nz,nk,nd,T_model);      % Labor hired
OPT.INC    = zeros(nz,nk,nd,T_model);      % Corp distribution
OPT.CIT    = zeros(nz,nk,nd,T_model);      % Corporate income tax
OPT.COUPON = zeros(nz,nk,nd,T_model);      % Debt coupon payment

% Initialize forward-looking corp values
W_step          = W0;
Solvent_step    = Solvent0;

% Specify settings for dynamic optimization subproblems
optim_options = optimset('Display', 'off', 'TolFun', 1e-4, 'TolX', 1e-4);

% Solve dynamic optimization problem through backward induction
for t = T_model:-1:1
    
    year = t; 
    
    % Extract parameters for current year
    wage                = wages             (year);
    discount_factor     = discount_factors  (year);
    interest_rate_gross = 1/discount_factor - 1;

    for ik = 1:nk
        for id = 1:nd
            for iz = 1:nz
                
                % Calculate expected value curve using forward-looking values
                EW = discount_factor*reshape(sum(repmat(transz(iz,:)', [1,nk,nd]) .* W_step, 1), [nk, nd]);
                
                % Calculate probability of solvency 
                PS = reshape(sum(repmat(transz(iz,:)', [1, nk, nd]) .* Solvent_step, 1), [nk, nd]);
                
                % Initialize functions to set persistent parameters
                value_equity([], kv, dv, EW);
                calculate_income( 0, 0, ...
                    kv(ik), zv(iz), dv(id), ...
                    tfp, alpha_k, alpha_n, k_adjustment_param, depreciation, ...
                    profit_tax, investment_deduction, interest_deduction, ...
                    PS, kv, dv, ...
                    wage );
                
                % Find choice (capital,debt)
                if isdynamic
                    k0          = kv(ik); 
                    d0          = dv(id);
                    [x, v]      = fminsearch(@value_equity, [k0, d0], optim_options);
                    capital     = x(1);
                    debt        = x(2);
                    equity      = -v;

                    is_solvent  = (equity >= 0);
                    % Are shareholders wiped out?
                    if( ~is_solvent )
                        % Are there bondholders to take-over firm?
                        if( dv(id) > 0 )
                            % Reduce amount of capital, wipe out debt,
                            % and operate firm. Bondholders are new
                            % shareholders.
                            reduced_capital = kv(ik)*k_default_reduction;
                            calculate_income( 0, 0, ...
                                reduced_capital, zv(iz), 0, ...
                                tfp, alpha_k, alpha_n, k_adjustment_param, depreciation, ...
                                profit_tax, investment_deduction, interest_deduction, ...
                                PS, kv, dv, ...
                                wage );  
                            % Re-optimize with new "state"
                            k0          = reduced_capital;
                            d0          = 0;
                            [x, v]      = fminsearch(@value_equity, [k0, d0], optim_options);
                            capital     = x(1);
                            debt        = x(2);
                            equity      = -v;
                        end % debt default
                        
                        % Is the firm (still) insolvent?
                        if( equity < 0 )
                            % Fire-sale on capital. Liquidate firm
                            % TODO: Need some way to guarantee positive
                            % equity. Ideally would also transition to
                            % different shock (from some initial
                            % distribution)
                        end % no equity
                    end % is_solvent
                else
                    capital     = K_static(iz,ik,id,t);
                    debt        = DEBT_static (iz,ik,id,t);
                    equity      = -value_equity( [capital, debt] );
                    is_solvent  = (equity >= 0);
                    % TODO: check that static is correct
                    
                end % isdynamic  
                
                % Value of firm to shareholders
                %    NOTE: This should always be > 0
                OPT.W       (iz,ik,id,t) = equity;

                % Record optimal decision values
                OPT.K      (iz,ik,id,t) = capital;
                OPT.DEBT   (iz,ik,id,t) = debt;
                OPT.SOLVENT(iz,ik,id,t) = is_solvent;
                
                % Find derived vars from choice vars
                [labor, income, cit, debt_coupon] = calculate_income( capital, debt );
                OPT.LAB   (iz,ik,id,t) = labor;
                OPT.INC   (iz,ik,id,t) = income;
                OPT.CIT   (iz,ik,id,t) = cit;
                OPT.COUPON(iz,id,id,t) = debt_coupon;
            end % iz
        end % id
    end % ik
    
    % Update forward-looking enterprise values
    W_step          = OPT.W(:, :, :, t);
    Solvent_step    = OPT.SOLVENT(:, :, :, t);
    
end % t


end % solve_firms



%%
% Value function: State is (capital, debt)
function v = value_equity(x, kv_, dv_, EW_) %#codegen

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
    %   REM: EW has already been multiplied by discount_factor
    v = interp2(kv', dv', EW', new_capital, new_debt, 'linear') ... 
        + income;

    % Negate value for minimization and force to scalar for C code generation
    v = -v(1);

end % value_enterprise


%%
% Calculates firm's budget & taxes
function [labor, income, cit, debt_coupon] = calculate_income( ...
            new_capital, new_debt, ...
            capital_, shock_, debt_, ...
            tfp_, alpha_k_, alpha_n_, k_adjustment_param_, depreciation_, ...
            profit_tax_, investment_deduction_, interest_deduction_, ...
            PS_, kv_, dv_, ...
            wage_, interest_rate_gross_ ) %#codegen

    % Enforce function inlining for C code generation
    coder.inline('always');

    % Define parameters as persistent variables
    persistent  capital shock debt ...
                tfp alpha_k alpha_n k_adjustment_param depreciation ...
                profit_tax investment_deduction interest_deduction ...
                PS kv dv ...
                wage interest_rate_gross ...
                opt_labor labor_cost revenue ...
                initialized;

    % Initialize parameters for C code generation
    if isempty(initialized)
        capital = 0; shock = 0; debt = 0; 
        tfp = 0; alpha_k = 0; alpha_n = 0; k_adjustment_param = 0; depreciation = 0;
        profit_tax = 0; investment_deduction = 0; interest_deduction = 0;
        PS = 0; kv = 0; dv = 0;
        wage = 0; interest_rate_gross = 0;
        opt_labor = 0; labor_cost = 0; revenue = 0;
        initialized = true;
    end

    % Set parameters if provided
    if (nargin > 2)
        capital = capital_; shock = shock_; debt = debt_; 
        tfp = tfp_; alpha_k = alpha_k_; alpha_n = alpha_n_; k_adjustment_param = k_adjustment_param_; depreciation = depreciation_;
        profit_tax = profit_tax_; investment_deduction = investment_deduction_; interest_deduction = interest_deduction_;
        PS = PS_; kv = kv_; dv = dv_;
        wage = wage_; interest_rate_gross = interest_rate_gross_;
        
        % Calculate derived vars: revenue and labor_cost
        %   Optimal labor hiring choice (from FOC)
        opt_labor               = (wage/(alpha_n*shock*(capital^alpha_k)))^(1/(1+alpha_n));
        revenue                 = tfp*(capital^alpha_k)*(opt_labor^alpha_n);
        labor_cost              = opt_labor*wage;

        % Assume that if initializing, we do not calculate
        return;
    end % initialize
    
    % Capital, investment, and equity
    net_investment          = new_capital - capital;
    investment              = net_investment + capital*depreciation;
    capital_adjustment_cost = (net_investment/capital)*net_investment*k_adjustment_param;
    
    % Debt and default
    probability_solvent     = interp2(kv', dv', PS', new_capital, new_debt, 'linear');
    new_debt_issue          = new_debt/interest_rate_gross * probability_solvent; 
    debt_coupon             = new_debt - new_debt_issue;
    
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
    income                  = pretax_revenue - cit + (new_debt_issue - debt);
    
    % Labor only depends on old capital
    labor                   = opt_labor;

end % calculate_income




