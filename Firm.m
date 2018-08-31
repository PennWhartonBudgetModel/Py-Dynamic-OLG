%%
% Firms sector
%    -- currently contains only Single Firm production
%

classdef Firm < handle
    
    properties (Access = public )
        priceCapital0 = 1;
    end % public properties
    
    properties (Access = private )
        TFP                 ;           % A
        capitalShare        ;           % alpha
        depreciationRate    ;           % delta
        riskPremium         ;           % for high/low return
        interestRate        ;           % market return on capital (and debt)
        allowBusinessDebt   ;           % whether to allow borrowing
        
        ratioCorpCapital    ;           % portion of total capital which is corp (vs. pass-through)
        
        corpOtherExpensesRate   ;           % unmodelled expenses as rate on corp GDP
        corpExpensingRate       ;           % phi_exp
        corpEffectiveTaxRate    ;           % effective tax on profits
        corpStatutoryTaxRate    ;           % statutory tax rate on profits
        corpTaxbaseAdjustment   ;           % match taxbase to empirics
        corpCreditsRate         ;           % tax credits as rate on corp GDP
        corpDeductionsRate      ;           % tax deductions as rate on corp GDP
        corpInterestDeduction   ;           % phi_int
        
        otherExpensesRate   ;           % unmodelled expenses as rate on corp GDP
        expensingRate       ;           % phi_exp
        effectiveTaxRate    ;           % effective tax on profits
        statutoryTaxRate    ;           % statutory tax rate on profits
        taxbaseAdjustment   ;           % match taxbase to empirics
        creditsRate         ;           % tax credits as rate on corp GDP
        deductionsRate      ;           % tax deductions as rate on corp GDP
        interestDeduction   ;           % phi_int
        
        corpLeverageCost        ;       % nu_corp
        passLeverageCost        ;       % nu_pass
        capitalAdjustmentCost;          % eta
        
        corpInitLeverageRatio   ;       % stored value of initial leverage ratio for corps
        passInitLeverageRatio   ;       % stored value of initial leverage ratio for pass-throughs
        
        priceCapital_       ;           % cached value at time of construction
        userCostCapital0    ;
        
        capital             ;           % Capital series
        labor               ;           % Efficient labor series
        KLratio             ;           % Capital to labor ratio, ( desired not capital/labor)
        invtocaps_end       ;           % Investment/capital at time T
        
    end % properties
    
    methods
        
        %%
        %  Constructor
        %     INPUTS:   Aggregate = Dynamic or Static aggregates
        %               Market    = the prices of stuff
        %               paramsTax = ParamGenerator.tax()
        %               paramsProduction = ParamGenerator.production()
        %               interestRate = Guess at dividends rate
        function [this] = Firm( Aggregate, Market, paramsTax, paramsProduction, interestRate )
            
            this.TFP                    = paramsProduction.A;
            this.capitalShare           = paramsProduction.alpha;
            this.depreciationRate       = paramsProduction.depreciation;
            this.riskPremium            = paramsProduction.risk_premium;
            this.allowBusinessDebt      = paramsProduction.allowBusinessDebt;
            this.capitalAdjustmentCost  = paramsProduction.capitalAdjustmentCost;
            
            % TBD: read these from params
            this.ratioCorpCapital       = ones(size(paramsTax.rateCorporate));
            init_rateForPassThrough     = 0.39;  % implied tax rate for pass-through income
            
            
            this.effectiveTaxRate   = paramsTax.rateCorporate;
            this.statutoryTaxRate   = paramsTax.rateCorporateStatutory;
            this.expensingRate      = paramsTax.shareCapitalExpensing; % REVISE W/ new interface
            this.interestDeduction  = paramsTax.interestDeduction;
            this.creditsRate        = paramsTax.creditsRate;
            this.taxbaseAdjustment  = paramsTax.taxbaseAdjustment;
            this.deductionsRate     = paramsTax.deductionsRate;
            this.otherExpensesRate  = paramsProduction.otherExpensesRate;

            this.corpInitLeverageRatio  = paramsProduction.initialCorpLeverage;
            this.passInitLeverageRatio  = paramsProduction.initialPassThroughLeverage;
            
            this.setInterestRate( interestRate ); % recalculates leverage cost
            
            this.capital        = Aggregate.caps';
            this.labor          = Aggregate.labeffs';
            this.KLratio        = Market.rhos';
            this.invtocaps_end  = Market.invtocaps(end);
            
            % Calculate the price of capital (p_K, see docs)
            this.userCostCapital0 = 1 - paramsTax.init_shareCapitalExpensing .* paramsTax.init_rateCorporate ...
                                    + this.capitalAdjustmentCost .* Market.invtocaps_0;
            % TBD: Adjust userCostCapital0 to have mix of pass and corp
            this.priceCapital_    = this.priceCapital();
        end % constructor
        
        
        %%
        %  Reset the interest rate and recalculate leverage cost
        function setInterestRate( this, interestRate )
            
            if( interestRate <= 0 )
                throw(MException('LEVERAGE_COST:NOT_POSITIVE','Interest rate must be positive for leverage cost calculation.'));
            end
            
            this.interestRate   = interestRate;
            
            % Calculate or use leverage cost
            % Rem: the leverage cost is size invariant, so set capital=1
            %      also, capital from initLeverageRatio is in $, so already
            %      scaled
            this.corpLeverageCost       = this.calculateLeverageCost( this.corpInitLeverageRatio, 1);
            this.passLeverageCost       = this.calculateLeverageCost( this.passInitLeverageRatio, 1);
        end % setInterestRate
        
        %%
        % Wage level the firm is willing to pay
        %   under competitive market assumptions, this is the MP_L
        function [wages] = wageRequired( this )
            alpha = this.capitalShare;
            wages = this.TFP .* (1 - alpha) .* (this.KLratio .^ alpha);
        end % wageRequired
        
        %%
        % MPK -- Just for reporting and indexing, not a received payment
        function [r] = MPK( this )
            alpha = this.capitalShare;
            r     = this.TFP .* alpha .* ( this.KLratio .^(alpha-1) ); 
        end
        
        %%
        % Output -- Just for reporting and indexing, not a received payment
        %           Used with 'capital' input to solve steady state
        function [out] = output( this, capital, labor )
            if( nargin == 1 )
                capital = this.capital;
                labor   = this.labor;
            end
            alpha = this.capitalShare;
            out   = this.TFP .* (capital .^ alpha) .* ( labor .^(1 - alpha) ); 
        end
        
        
        %%
        %   Calculate various business payments and taxes
        function [x] = distributions( this )
            
            wage = this.wageRequired();
            
            [divs, cits, debts] = corpDividends( this, this.capital, this.KLratio, wage );
            x.corpDividends     = divs;
            x.corpTaxs          = cits;
            x.corpDebts         = debts;
            
            [divs, taxableIncome, debts] = passDividends( this, this.capital, this.KLratio, wage );
            x.passDividends     = divs;
            x.passTaxableIncome = taxableIncome;
            x.passDebts         = debts;
        end % distributions
        
        
        
        %% 
        function divs = dividends( this, capital, klRatio, wage )
            % TEMP: This should combine pass-through and corp
            %  for now, just do corp
            divs = corpDividends( this, capital, klRatio, wage );
        end % dividends
        
        
        
        %%
        function [divs, cits, debts] = corpDividends( this, capital, klRatio, wage )
            % Inputs : capital
            %          klRatio & wage 
            % Outputs: divs = dividend is the return of the corporation net of tax.
            %          cits = corporate income taxes paid by the firms
            
                
            % Labor and capital
            labor = ( capital ./ klRatio );
            
            % Total revenues
            y = this.TFP * ( (capital.^this.capitalShare) .* (labor.^(1-this.capitalShare)) );
            
            % Wage payments
            wagesPaid = wage .* labor;
            
            % Replace depreciated capital (rem: at cost to buy it)
            depreciation = this.depreciationRate .* capital .* this.priceCapital_;
            
            % Risk premium (rem: at cost to buy it)
            risk = this.riskPremium .* capital .* this.priceCapital_;
            
            % Investment
            investment = [capital(2:end) - (1 - this.depreciationRate)*capital(1:end-1); ...
                          capital(end) * this.invtocaps_end ];
            
            % Investment expensing 'subsidy'
            expensing    = max(this.expensingRate .* investment .* this.priceCapital_, 0);
            
            % Capital adjustment cost
            %    = eta/2 * (I/K)^2 * K
            adjustmentCost = (this.capitalAdjustmentCost/2) .* (investment .* investment) ./ capital;
            
            % Find optimal debt, interest tax benefit, and leverage cost
            [debts, debtCost, debtTaxBenefit]  = this.calculateDebt( capital, this.corpLeverageCost );
            
            % Combine to get net tax
            taxbase =   ( y                                           ...
                        - wagesPaid                                 ...
                        - this.interestRate .* debts                ...
                        - depreciation                              ...
                        - this.otherExpensesRate .* y               ...
                        - adjustmentCost                            ...
                        ) .* this.taxbaseAdjustment                 ...
                        - this.deductionsRate .* y                  ...
                        ;
            cits  = taxbase .* this.statutoryTaxRate            ...
                    - expensing .* this.statutoryTaxRate        ...
                    - y .* this.creditsRate                     ...
                    - debtTaxBenefit;
            % TEMP: Until we get correct inputs
            cits  = (y - wagesPaid) .* this.effectiveTaxRate    ...
                    - expensing .* this.statutoryTaxRate        ...
                    - y .* this.creditsRate                     ...
                    - debtTaxBenefit;

            % Calculate returns to owners
            divs  = y                                   ... % revenues
                  - wagesPaid                           ... % labor costs
                  - depreciation                        ... % replace depreciated capital
                  - adjustmentCost                      ... % cap adjustment cost
                  - risk                                ... % discount money lost due to risk
                  - debtCost                            ... % Cost of leverage
                  - cits;                                   % net taxes
                  
        end % corpDividends
        
        
        %%
        function [divs, taxableIncome, debts] = passDividends( this, capital, klRatio, wage )
            % Inputs : capital
            %          klRatio & wage 
            % Outputs: divs = dividend is the return of the corporation net of tax.
            %          cits = corporate income taxes paid by the firms
            
                
            % Labor and capital
            labor = ( capital ./ klRatio );
            
            % Total revenues
            y = this.TFP * ( (capital.^this.capitalShare) .* (labor.^(1-this.capitalShare)) );
            
            % Wage payments
            wagesPaid = wage .* labor;
            
            % Replace depreciated capital (rem: at cost to buy it)
            depreciation = this.depreciationRate .* capital .* this.priceCapital_;
            
            % Risk premium (rem: at cost to buy it)
            risk = this.riskPremium .* capital .* this.priceCapital_;
            
            % Investment
            investment = [capital(2:end) - (1 - this.depreciationRate)*capital(1:end-1); ...
                          capital(end) * this.invtocaps_end ];
            
            % Investment expensing 'subsidy'
            expensing    = max(this.expensingRate .* investment .* this.priceCapital_, 0);
            
            % Capital adjustment cost
            %    = eta/2 * (I/K)^2 * K
            adjustmentCost = (this.capitalAdjustmentCost/2) .* (investment .* investment) ./ capital;
            
            % Find optimal debt, interest tax benefit, and leverage cost
            [debts, debtCost, debtTaxBenefit]  = this.calculateDebt( capital, this.passLeverageCost );
            
            % Calculate returns to owners
            divs  = y                                   ... % revenues
                  - wagesPaid                           ... % labor costs
                  - depreciation                        ... % replace depreciated capital
                  - adjustmentCost                      ... % cap adjustment cost
                  - risk                                ... % discount money lost due to risk
                  - debtCost                            ... % Cost of leverage
                  ;                                   
            
            % Taxable income
            taxableIncome = divs                    ... % dividends
                          - debtTaxBenefit          ... % interest deduction
                          - expensing               ... % investment expensing
                          ;
        end % passDividends

        
        
        
        %%
        % Calculate capital gains rates from value of capital.
        function [capgains] = capitalGains( this, capital )
            if( nargin == 1 )
                priceCapital = this.priceCapital_;
            else
                priceCapital = this.priceCapital( capital );
            end
            capgains        = zeros(size(priceCapital));
            capgains(1,1)   = (priceCapital(1) - this.priceCapital0)/this.priceCapital0;
            for t = 2:length(capgains)
               capgains(t,1) = (priceCapital(t) - priceCapital(t-1))/priceCapital(t-1);
            end
        end % capitalGains
        
        
        %%
        % Calculate the price of capital
        function [price] = priceCapital( this, capital )
            if( nargin == 1 )
                if( ~isempty( this.priceCapital_ ) )
                    price = this.priceCapital_;
                    return;
                end
                capital       = this.capital;
            end
            
            % Gross Investment
            investment = [capital(2:end) - (1 - this.depreciationRate)*capital(1:end-1); ...
                          capital(end) * this.invtocaps_end ];
            invtocaps  = investment ./ capital;
            
            userCostCapitalCorp = 1 - this.expensingRate .* this.effectiveTaxRate ...
                                    + this.capitalAdjustmentCost .* invtocaps;
            userCostCapitalPass = 1;  % TBD: Fix this
            userCostCapital     = this.ratioCorpCapital .* userCostCapitalCorp      ...
                                    + (1-this.ratioCorpCapital) .* userCostCapitalPass;
                                
            price               = (this.priceCapital0/this.userCostCapital0) * userCostCapital;
        end % priceCapital

        
        
        %% 
        % Calculate K/L ratio from dividend rate
        function [KLratio, caps] = calculateKLRatio( this, fixedAfterTaxReturn, foreignTaxRates, ...
                            init_caps, labor )
            % Inputs : dividendRate = dollars received as dividend per
            %                         dollar owned of equity
            %          init_caps = capital initial guess
            %          labor = efficient units of labor (from last iteration)
            %          invtocapsT_model = last period guess of I/K
            % Outputs: KLratio that generates the dividendRate of inputs
            
            % Find dividends rate needed to fix returns to target
            dividendRate    = ( fixedAfterTaxReturn - this.capitalGains() ) ...
                              ./ (1 - foreignTaxRates);
        
            % Initialize variables
            caps    = init_caps;
            divRate = dividendRate;
            
            tolerance = 1e-12;
            err_div   = Inf;
            
            while( err_div > tolerance )
                
                % Update capital 
                %  if divRate > dividendRate 
                %    --> caps gets bigger, and divRate gets smaller
                caps(2:end) = caps(2:end) .* ((1+divRate(2:end)) ./ (1+dividendRate(2:end)) );
                
                % To prevent nonsensical results, prevent non-positive
                % capital.
                %  TBD: This may keep loop from converging under some
                %  circumstances.
                caps        = max( 1e-5, caps );
                
                
                K_by_L = caps ./ labor;
                wage   = this.TFP * (1-this.capitalShare) .* (K_by_L .^ this.capitalShare);
                
                divs = this.dividends( caps, K_by_L, wage );
                divRate = divs ./ (caps .* this.priceCapital(caps));
                
                err_div = max(abs((divRate(2:end) - dividendRate(2:end)) ./ dividendRate(2:end)));
                
                % Update dividends rate since capitalGains change with caps
                dividendRate    = ( fixedAfterTaxReturn - this.capitalGains(caps) ) ...
                              ./ (1 - foreignTaxRates);
            end % while
            
            % Calculate capital-labor ratio
            KLratio  = caps ./ labor;
        end % calculateKLRatio

        
        
        %% 
        % Calculate optimal debt 
        %   nu() = 1/nu (B/K_s)^nu , where K_s = K h p_K; h=0.1
        function [debt, debtCost, taxBenefit] = calculateDebt( this, capital, leverageCost )
           
           % If not using debt, then jump out.
           if( ~this.allowBusinessDebt )
                debt        = zeros(size(capital));
                debtCost    = zeros(size(capital));
                taxBenefit  = zeros(size(capital));
                return;
            end

            % Calculate optimal debt
            h               = 100;
            taxBenefitRate  = this.interestDeduction .* this.statutoryTaxRate .* this.interestRate;
            nu              = leverageCost;
            capital_scaled  = capital .* this.priceCapital_ .* h;
            
            d       = 1 - this.interestDeduction .* this.statutoryTaxRate;  
            x       = 1/(nu-1);
            debt    = (( taxBenefitRate ./ d ) .^ x) .* capital_scaled;
            
            debtCost    = ((1/nu) .* (debt./capital_scaled ) .^ nu) .* capital_scaled;
            taxBenefit  = debt .* taxBenefitRate;
            
        end % calculateDebt
        
        
        %% 
        % Calculate leverageCost from B/K ratio target
        %    Rem: debt and capital_value are in $
        function [nu] = calculateLeverageCost( this, debt, capital_value )
            
            h               = 100;
            taxBenefitRate  = this.interestDeduction .* this.statutoryTaxRate .* this.interestRate;
            capital_scaled  = capital_value .* h;
            
            logBK = log( debt ./ capital_scaled );
            n     = log( taxBenefitRate ) ...
                      - log( 1 - this.statutoryTaxRate .* this.interestDeduction );
            nu    = 1 + n ./ logBK;

        end % calculateLeverageCost



        %% Test
        function  [] = testMe( this )
            
            T_model = 25;
            
            % Make sample inputs
            fixedAfterTaxReturn = ones(T_model,1) .* 0.05;
            foreignTaxRates     = ones(T_model,1) .* 0.10;
            labor = ones(T_model,1);
            init_caps = ones(T_model,1) * 12;
            invtocap_Tmodel = 0.07;
            
            [klRatio tCaps]     = this.calculateKLRatio( fixedAfterTaxReturn, foreignTaxRates, ...
                            init_caps, labor )
                        
             tempCaps   = klRatio .* labor;
             tempI      = [diff(tempCaps); invtocap_Tmodel * tempCaps(T_model)];
             tempWages  = this.TFP*(1-this.capitalShare)*(klRatio.^this.capitalShare);
             
             [tempDividends, ~]  = this.dividends( tempCaps          ...
                                                     ,   tempI     ...
                                                     ,   klRatio            ...
                                                     ,   tempWages           ...
                                                     );  
              tempDivRates = (tempDividends ./ (tempCaps .* this.priceCapital));                                   

              fprintf('diff in dividends');
              tempDivRates - effectiveDividendRate

        end %testMe

    end % instance methods
    
end % classdef

