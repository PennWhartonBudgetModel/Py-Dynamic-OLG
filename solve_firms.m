%%
% Solve dynamic optimization problem for the firms.
% 
%%


function OPT = solve_firms(W0, Default0, ...
                        isDYNAMIC, K_static, DEBT_static, DEFAULT_static, SHUTDOWN_static, ...
                        nz, nk, nd, T_model, transz, kv, zv, dv, ...
                        tfp, alpha_k, alpha_n, k_adjustment_param, depreciation, ...
                        k_default_reduction, equity_fixed_cost, equity_linear_cost, ...
                        operation_fixed_cost, operation_linear_cost, ...
                        profit_tax, investment_deduction, interest_deduction, ...
                        wages, discount_factors ) %#codegen


%% Argument verification
nz_max          = 50;
nk_max          = 50;
nd_max          = 50;
T_max           = 100;

assert( isa(W0          , 'double'  ) && (size(W0           , 1) <= nz_max  ) && (size(W0           , 2) <= nk_max  ) && (size(W0         , 3) <= nd_max  ));
assert( isa(Default0    , 'logical' ) && (size(Default0     , 1) <= nz_max  ) && (size(Default0     , 2) <= nk_max  ) && (size(Default0   , 3) <= nd_max  ));

assert( isa(isDYNAMIC      , 'logical' ) && (size(isDYNAMIC      , 1) == 1       ) && (size(isDYNAMIC      , 2) == 1       ) );
assert( isa(K_static       , 'double'  ) && (size(K_static       , 1) <= nz_max  ) && (size(K_static       , 2) <= nk_max  ) && (size(K_static       , 3) <= nd_max  ) && (size(K_static       , 4) <= T_max   ) );
assert( isa(DEBT_static    , 'double'  ) && (size(DEBT_static    , 1) <= nz_max  ) && (size(DEBT_static    , 2) <= nk_max  ) && (size(DEBT_static    , 3) <= nd_max  ) && (size(DEBT_static    , 4) <= T_max   ) );
assert( isa(DEFAULT_static , 'logical' ) && (size(DEFAULT_static , 1) <= nz_max  ) && (size(DEFAULT_static , 2) <= nk_max  ) && (size(DEFAULT_static , 3) <= nd_max  ) && (size(DEFAULT_static , 4) <= T_max   ) );
assert( isa(SHUTDOWN_static, 'logical' ) && (size(SHUTDOWN_static, 1) <= nz_max  ) && (size(SHUTDOWN_static, 2) <= nk_max  ) && (size(SHUTDOWN_static, 3) <= nd_max  ) && (size(SHUTDOWN_static, 4) <= T_max   ) );

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

assert( isa(k_default_reduction , 'double') && (size(k_default_reduction ,1) == 1 ) && (size(k_default_reduction , 2) == 1 ) );
assert( isa(equity_fixed_cost   , 'double') && (size(equity_fixed_cost   ,1) == 1 ) && (size(equity_fixed_cost   , 2) == 1 ) );
assert( isa(equity_linear_cost  , 'double') && (size(equity_linear_cost  ,1) == 1 ) && (size(equity_linear_cost  , 2) == 1 ) );

assert( isa(operation_fixed_cost , 'double') && (size(operation_fixed_cost ,1) == 1 ) && (size(operation_fixed_cost , 2) == 1 ) );
assert( isa(operation_linear_cost, 'double') && (size(operation_linear_cost,1) == 1 ) && (size(operation_linear_cost, 2) == 1 ) );

assert( isa(profit_tax          , 'double' ) && (size(profit_tax          , 1) == 1 ) && (size(profit_tax          , 2) == 1 ) );
assert( isa(investment_deduction, 'double' ) && (size(investment_deduction, 1) == 1 ) && (size(investment_deduction, 2) == 1 ) );
assert( isa(interest_deduction  , 'double' ) && (size(interest_deduction  , 1) == 1 ) && (size(interest_deduction  , 2) == 1 ) );

assert( isa(wages           , 'double'  ) && (size(wages           , 1) == 1       ) && (size(wages           , 2) <= T_max   ) );
assert( isa(discount_factors, 'double'  ) && (size(discount_factors, 1) == 1       ) && (size(discount_factors, 2) <= T_max   ) );

%% Dynamic optimization

% Initialize equity valuation and optimal decision value arrays
OPT.W          = zeros(nz,nk,nd,T_model);      % Corp value

OPT.K          = zeros(nz,nk,nd,T_model);      % Kapital
OPT.DEBT       = zeros(nz,nk,nd,T_model);      % Debt outstanding
OPT.DEFAULT    = false(nz,nk,nd,T_model);      % Defaulted on debt?
OPT.SHUTDOWN   = false(nz,nk,nd,T_model);      % Is firm operating?

% Initialize derived values 
OPT.DIVIDEND   = zeros(nz,nk,nd,T_model);      % Equity distribution
OPT.LABOR      = zeros(nz,nk,nd,T_model);      % Labor hired
OPT.OUTPUT     = zeros(nz,nk,nd,T_model);      % Production
OPT.CIT        = zeros(nz,nk,nd,T_model);      % Corporate income tax
OPT.NEXTCOUPON = zeros(nz,nk,nd,T_model);      % Debt coupon payment for next period (if paid)
OPT.INVESTMENT = zeros(nz,nk,nd,T_model);      % Investment
OPT.COUPON     = zeros(nz,nk,nd,T_model);      % Debt coupon payment

% Initialize forward-looking corp values
W_step         = W0;
Default_step   = Default0;

% Specify settings for dynamic optimization subproblems
optim_options = optimset('Display', 'off', 'TolFun', 1e-4, 'TolX', 1e-4);

% Solve dynamic optimization problem through backward induction
for t = T_model:-1:1
    
    % Extract parameters for current year
    wage                = wages             (t);
    discount_factor     = discount_factors  (t);
    interest_rate_gross = 1/discount_factor;
    
    % Fetch t+2 values (for bond calculations)
    %    NOTE: Defaulted firms have zero bond debt (assumed to be @ index=1)
    %          The "apres-default" firms are not necessarily those
    %          who actually defaulted in the previous period. They
    %          are firms which have equity value like the firms which 
    %          defaulted.
    if (t < T_model - 1)
        discount_factor_step = discount_factors(t+1);
        W_apres_default      = OPT.W(:,:,1,t+2);
    else % last two periods
        discount_factor_step = discount_factor;
        W_apres_default      = W0(:,:,1);
    end 

    for iz = 1:nz
        % NOTE: We pre-multiply expected values by discount factors to save
        % on later multiplications. So EW = 1/1+r * E[W']

        % Calculate expected value of firms using forward-looking values
        Repay      = ~Default_step;
        EW         = reshape(sum(repmat(transz(iz,:)', [1, nk, nd]) .* (W_step .* Repay), 1), [nk, nd]) * discount_factor;

        % Calculate probability of repayment, expected value of defaulted firms
        prob_repay       = reshape(sum(repmat(transz(iz,:)', [1, nk, nd]) .* Repay      , 1), [nk, nd]);
        EW_apres_default = reshape(sum(repmat(transz(iz,:)', [1, nk]) .* W_apres_default, 1), [1 , nk]) ...
                            * discount_factor * discount_factor_step;
        
        % Initialize to set persistent parameters
        value_equity([], kv, dv, EW);
        
        for ik = 1:nk
            for id = 1:nd
                
                current_capital = kv(ik); current_debt = dv(id); current_shock = zv(iz);

                % Initialize persistent parameters with current state
                calculate_income( 0, 0, ...
                    current_capital, current_shock, current_debt, ...
                    tfp, alpha_k, alpha_n, k_adjustment_param, depreciation, ...
                    profit_tax, investment_deduction, interest_deduction, ...
                    equity_fixed_cost, equity_linear_cost, ...
                    operation_fixed_cost, operation_linear_cost, ...
                    k_default_reduction, prob_repay, EW_apres_default, kv, dv, ...
                    wage, interest_rate_gross );
                
                % Find firm's policies for capital, debt, shutdown, default
                if isDYNAMIC
                    [x, v]              = fminsearch(@value_equity, [current_capital, current_debt], optim_options);
                    new_capital         = x(1);
                    new_debt            = x(2);
                    operation_value     = -v;
                    liquidation_value   = calculate_liquidation_value( current_capital, current_debt, k_adjustment_param); 
                                        
                    is_shutdown         = false;
                    is_default          = false;
                    
                    % Set value of defaulting on debt
                    if( current_debt > 0 )
                        default_value = 0;
                    else
                        default_value = -Inf;
                    end
                    
                    %  Current shareholders pick best option
                    [equity, continuation_choice] =  ... 
                            max([liquidation_value, default_value, operation_value]);
                    
                   if( continuation_choice == 1 )
                        % Firm ceases operation and distributes.
                        % Firm has no continuation value, but shareholders
                        % get distribution of liquidation proceeds.
                        new_capital     = 0;
                        new_debt        = 0;
                        is_shutdown     = true;
                        equity          = liquidation_value;

                        OPT.DIVIDEND  (iz,ik,id,t) = liquidation_value;
                        OPT.LABOR     (iz,ik,id,t) = 0;
                        OPT.OUTPUT    (iz,ik,id,t) = 0;
                        OPT.CIT       (iz,ik,id,t) = 0;
                        OPT.NEXTCOUPON(iz,ik,id,t) = 0;  
                        OPT.INVESTMENT(iz,ik,id,t) = 0;
                    % end liquidation
                    
                    elseif (continuation_choice == 2 )
                        % Defaulted firm does not operate this period.
                        % The next period, the firm is re-established 
                        % with bondholders as new shareholders.
                        %   NOTES: 
                        %     1. We do not "repay" interest deduction.
                        %     2. Default is charged a cost based on
                        %     capital.
                        new_capital     = k_default_reduction * current_capital;
                        new_debt        = 0;
                        is_default      = true;
                        equity          = 0; % current shareholders get nothing

                        OPT.DIVIDEND  (iz,ik,id,t) = 0;
                        OPT.LABOR     (iz,ik,id,t) = 0;
                        OPT.OUTPUT    (iz,ik,id,t) = 0;
                        OPT.CIT       (iz,ik,id,t) = 0; 
                        OPT.NEXTCOUPON(iz,ik,id,t) = 0; 
                        OPT.INVESTMENT(iz,ik,id,t) = 0; 
                   end % default
                % end DYNAMIC
                
                else  % STATIC
                    new_capital         = K_static(iz,ik,id,t);
                    new_debt            = DEBT_static(iz,ik,id,t);
                    equity              = -value_equity( [current_capital, current_debt] );
                    is_default          = DEFAULT_static(iz,ik,id,t);
                    is_shutdown         = SHUTDOWN_static(iz,ik,id,t);
                end % STATIC  
                
                % Value of firm to shareholders
                %    NOTE: This should always be >= 0
                OPT.W       (iz,ik,id,t) = equity;

                % Record optimal decision values
                OPT.K       (iz,ik,id,t) = new_capital;
                OPT.DEBT    (iz,ik,id,t) = new_debt;
                OPT.DEFAULT (iz,ik,id,t) = is_default;
                OPT.SHUTDOWN(iz,ik,id,t) = is_shutdown;
                
                % Calculate derived values, unless already done by shutdown
                % or default.
                if ~(is_shutdown || is_default)
                    [dividend, labor, output, cit, next_coupon, investment] = calculate_income( new_capital, new_debt );
                    OPT.DIVIDEND  (iz,ik,id,t) = dividend;
                    OPT.LABOR     (iz,ik,id,t) = labor;
                    OPT.OUTPUT    (iz,ik,id,t) = output;
                    OPT.CIT       (iz,ik,id,t) = cit;
                    OPT.NEXTCOUPON(iz,ik,id,t) = next_coupon;
                    OPT.INVESTMENT(iz,ik,id,t) = investment;
                end 
            end % id
        end % ik
    end % iz
    
    % Update forward-looking values
    W_step          = OPT.W(:, :, :, t);
    Default_step    = OPT.DEFAULT(:, :, :, t);
    
end % t

% Solve for debt coupon
nt = T_model-1 + (T_model == 1);   % to use same loop for steady-state and transition
for t=1:nt
    for iz=1:nz
        for id=1:nd
            for ik=1:nk

                next_k_val  = OPT.K   (iz,ik,id,t);
                next_d_val  = OPT.DEBT(iz,ik,id,t);
                [~, next_k] = min(abs(kv - next_k_val)); % find closest grid point in k
                [~, next_d] = min(abs(dv - next_d_val)); % find closest grid point in debt
                next_t      = t + (T_model > 1);        % for steady-state (T_model==1) and we map within the period
                
                for next_z=1:nz
                    is_repaid = ~OPT.DEFAULT(next_z,next_k,next_d,next_t);
                    % Coupon is not paid in default
                    OPT.COUPON(next_z,next_k,next_d,next_t) ...
                           = OPT.NEXTCOUPON(iz,ik,id,t)*is_repaid;
                end % next_z
                
            end % ik
        end % id
    end % iz
end % t


end % solve_firms



%%
% Value of shareholder equity: State is (capital, debt)
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
    
    [dividend, ~, ~, ~, ~, ~] = calculate_income( new_capital, new_debt );
    
    % Calculate enterprise value
    %   REM: EW has already been multiplied by discount_factor
    v = interp2(kv', dv', EW', new_capital, new_debt, 'linear') ... 
        + dividend;

    % Negate value for minimization and force to scalar for C code generation
    v = -v(1);

end % value_enterprise


%%
% Calculates firm's budget & taxes
function [dividend, labor, output, cit, next_coupon, investment] = calculate_income( ...
            new_capital, new_debt, ...
            capital_, shock_, debt_, ...
            tfp_, alpha_k_, alpha_n_, k_adjustment_param_, depreciation_, ...
            profit_tax_, investment_deduction_, interest_deduction_, ...
            equity_fixed_cost_, equity_linear_cost_, ...
            operation_fixed_cost_, operation_linear_cost_, ...
            k_default_reduction_, prob_repay_, EW_default_, kv_, dv_, ...
            wage_, interest_rate_gross_ ) %#codegen

    % Enforce function inlining for C code generation
    coder.inline('always');

    % Define parameters as persistent variables
    persistent  capital shock debt ...
                tfp alpha_k alpha_n k_adjustment_param depreciation ...
                profit_tax investment_deduction interest_deduction ...
                equity_fixed_cost equity_linear_cost ... 
                operation_fixed_cost operation_linear_cost ...
                k_default_reduction prob_repay EW_default kv dv ...
                wage interest_rate_gross ...
                opt_labor revenue labor_cost capital_cost gross_profit ...
                debt_deduc ...
                initialized;

    % Initialize parameters for C code generation
    if isempty(initialized)
        capital = 0; shock = 0; debt = 0; 
        tfp = 0; alpha_k = 0; alpha_n = 0; k_adjustment_param = 0; depreciation = 0;
        profit_tax = 0; investment_deduction = 0; interest_deduction = 0;
        equity_fixed_cost = 0; equity_linear_cost = 0; 
        operation_fixed_cost = 0; operation_linear_cost = 0;
        k_default_reduction = 0; prob_repay = 0; EW_default = 0; kv = 0; dv = 0;
        wage = 0; interest_rate_gross = 0;
        
        opt_labor = 0; revenue = 0; labor_cost = 0; capital_cost = 0; gross_profit = 0;
        debt_deduc = 0;
        initialized = true;
    end

    % Set parameters if provided
    if (nargin > 2)
        capital = capital_; shock = shock_; debt = debt_; 
        tfp = tfp_; alpha_k = alpha_k_; alpha_n = alpha_n_; k_adjustment_param = k_adjustment_param_; depreciation = depreciation_;
        profit_tax = profit_tax_; investment_deduction = investment_deduction_; interest_deduction = interest_deduction_;
        equity_fixed_cost = equity_fixed_cost_; equity_linear_cost = equity_linear_cost_; 
        operation_fixed_cost = operation_fixed_cost_; operation_linear_cost = operation_linear_cost_; ...
        k_default_reduction = k_default_reduction_; prob_repay = prob_repay_; EW_default = EW_default_; kv = kv_; dv = dv_;
        wage = wage_; interest_rate_gross = interest_rate_gross_;
        
        % Calculate derived vars: revenue and labor_cost
        %   Optimal labor hiring choice (from FOC)
        opt_labor       = (wage/(alpha_n*shock*tfp*(capital^alpha_k)))^(1/(alpha_n-1));
        revenue         = tfp*shock*(capital^alpha_k)*(opt_labor^alpha_n);
        labor_cost      = opt_labor*wage;
        capital_cost    = operation_linear_cost*capital + operation_fixed_cost;
        gross_profit    = revenue - labor_cost - capital_cost;
        
        % Helper variable for debt calc
        debt_deduc      = interest_deduction*profit_tax/interest_rate_gross; 
        % Assume that if initializing, we do not calculate
        return;
    end % initialize
    
    % Labor and output depend only on state
    labor           = opt_labor;
    output          = revenue;

    % Capital, investment
    investment      = new_capital - (1-depreciation)*capital;
    
    % Debt and default
    P_repay         = interp2(kv', dv', prob_repay', new_capital, new_debt, 'linear');
    default_capital = new_capital * k_default_reduction;  % rem: default debt = 0
    value_default   = interp1(kv, EW_default, default_capital, 'linear'); 
    value_default   = value_default(1); % This is to fix MEX issue w/ dimensionality discovery.
    debt_principal  = (new_debt*(1 - debt_deduc)*P_repay + value_default*(1-P_repay))...
                       / ...
                      (interest_rate_gross + (debt_deduc/(1-debt_deduc))*P_repay);
    next_coupon     = new_debt - debt_principal;
    % REM: We cannot find current coupon (not in state), but the g'vt
    % "prepays" the interest deduction subsidy in the current period for
    % next period coupon. (This is lost in default).
 
    % Free cashflow before taxes and any financing costs
    pretax_profit   = gross_profit  ...
                        - capital_adjustment_cost( capital, investment, k_adjustment_param ) ...
                        - investment;

    % Calculate corp tax with relevant deductions -- investment & debt interest
    cit             = profit_tax*(pretax_profit ...
                        + (1-investment_deduction)* investment ...
                        + (1-interest_deduction  )* next_coupon );
    
    % After-tax distribution includes $ from debt issuance
    dividend        = pretax_profit - cit + (new_debt - debt);
    
    % If issuing equity (i.e. negative dividend) there are additional costs
    if( dividend < 0 )
        dividend = dividend - equity_fixed_cost + equity_linear_cost*dividend;        
    end
    
end % calculate_income


%%  
%   Capital adjustment cost function
function cost = capital_adjustment_cost( capital, investment, param)
    coder.inline('always');
    
    cost = (investment/capital)*investment*param;
end

%%  
%   Liquidation function
% It assumes full uninstall of capital and payment of
% current debt.
function v = calculate_liquidation_value( capital, debt, k_param)
    coder.inline('always');

    capital_sale_value = max(0, capital ... 
                        - capital_adjustment_cost(capital, -capital, k_param));
    v  = capital_sale_value - debt; 
end 


