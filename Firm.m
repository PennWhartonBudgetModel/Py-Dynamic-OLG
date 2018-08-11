%%
% Firms sector
%    -- currently contains only Single Firm production
%

classdef Firm < handle
    
    properties (Constant)
        SINGLEFIRM          = 0;
        PASSTHROUGH         = 1;
    end
    
    properties 
        TFP                 ;           % A
        capitalShare        ;           % alpha
        depreciationRate    ;           % delta
        riskPremium         ;           % for high/low return
        interestRate        ;           % market return on capital (and debt)
        allowBusinessDebt   ;           % whether to allow borrowing
        
        otherExpensesRate   ;           % unmodelled expenses as rate on corp GDP
        
        expensingRate       ;           % phi_exp
        effectiveTaxRate    ;           % effective tax on profits
        statutoryTaxRate    ;           % statutory tax rate on profits
        taxbaseAdjustment   ;           % match taxbase to empirics
        creditsRate         ;           % tax credits as rate on corp GDP
        deductionsRate      ;           % tax deductions as rate on corp GDP
        interestDeduction   ;           % phi_int
        
        leverageCost        ;           % nu
        capitalAdjustmentCost;          % eta
        
        initLeverageRatio   ;           % stored value of initial leverage ratio 
        
        priceCapital        ;           % See documentation. This is p_K
        priceCapital0 = 1   ;
        
        capital             ;           % Capital series
        labor               ;           % Efficient labor series
        KLratio             ;           % Capital to labor ratio, ( desired not capital/labor)
        invtocaps           ;           % Investment/capital 
        invtocaps_0         ;           % Investment/capital at time 0 (steady state)
        firmType            ;           % from the enumeration
        
    end % properties
    
    methods
        
        %%
        %  Constructor
        %     INPUTS:   Aggregate = Dynamic or Static aggregates
        %               Market    = the prices of stuff
        %               paramsTax = ParamGenerator.tax()
        %               paramsProduction = ParamGenerator.production()
        %               interestRate = Guess at dividends rate
        function [this] = Firm( Aggregate, Market, paramsTax, paramsProduction, interestRate, firmType )
            
            this.firmType           = firmType;   
            if( ~(firmType == Firm.SINGLEFIRM || firmType == Firm.PASSTHROUGH) )
                throw(MException('Firm:firmType','firmType must be SingleFirm or PassThrough'));
            end
            
            this.TFP                    = paramsProduction.A;
            this.capitalShare           = paramsProduction.alpha;
            this.depreciationRate       = paramsProduction.depreciation;
            this.riskPremium            = paramsProduction.risk_premium;
            this.allowBusinessDebt      = paramsProduction.allowBusinessDebt;
            this.capitalAdjustmentCost  = paramsProduction.capitalAdjustmentCost;
            
            switch this.firmType
                case Firm.SINGLEFIRM
                    this.effectiveTaxRate   = paramsTax.rateCorporate;
                    this.statutoryTaxRate   = paramsTax.rateCorporateStatutory;
                    this.expensingRate      = paramsTax.shareCapitalExpensing; % REVISE W/ new interface
                    this.interestDeduction  = paramsTax.interestDeduction;
                    this.creditsRate        = paramsTax.creditsRate;
                    this.taxbaseAdjustment  = paramsTax.taxbaseAdjustment;
                    this.deductionsRate     = paramsTax.deductionsRate;
                    this.initLeverageRatio  = paramsProduction.initialCorpLeverage;
                    
                    this.otherExpensesRate  = paramsProduction.otherExpensesRate;
                    
                case Firm.PASSTHROUGH
                    this.effectiveTaxRate   = 0;  % TEMP: Should come from ParamGenerator as top marginal PIT rate
                    this.statutoryTaxRate   = this.effectiveTaxRate;
                    this.expensingRate      = paramsTax.shareCapitalExpensing; % REVISE W/ new interface
                    this.interestDeduction  = paramsTax.interestDeduction;
                    this.creditsRate        = paramsTax.creditsRate;
                    this.initLeverageRatio  = paramsProduction.initialPassThroughLeverage;
                    
                    this.otherExpensesRate  = paramsProduction.otherExpensesRate;
                    
            end
            
            this.setInterestRate( interestRate ); % recalculates leverage cost
            
            this.capital        = Aggregate.caps';
            this.labor          = Aggregate.labeffs';
            this.KLratio        = Market.rhos';
            this.invtocaps      = Market.invtocaps';
            this.invtocaps_0    = Market.invtocaps_0;
            
            % Calculate the price of capital (p_K, see docs)
            userCostCapital     = 1 - paramsTax.shareCapitalExpensing .* paramsTax.rateCorporate ...
                                    + this.capitalAdjustmentCost .* this.invtocaps;
            userCostCapital0    = 1 - paramsTax.init_shareCapitalExpensing .* paramsTax.init_rateCorporate ...
                                    + this.capitalAdjustmentCost .* this.invtocaps_0;
            this.priceCapital   = this.priceCapital0 * (userCostCapital ./ userCostCapital0);
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
            this.leverageCost       = this.calculateLeverageCost( this.initLeverageRatio, 1);
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
        function [out] = GrossOutput( this, capital, labor )
            if( nargin == 1 )
                capital = this.capital;
                labor   = this.labor;
            end
            alpha = this.capitalShare;
            out   = this.TFP .* (capital .^ alpha) .* ( labor .^(1 - alpha) ); 
        end
        
        
        %%
        function [divs, cits, debts] = dividends( this, capital, klRatio, wage )
            % Inputs : capital
            %          klRatio & wage 
            % Outputs: divs = dividend is the return of the corporation net of tax.
            %          cits = corporate income taxes paid by the firms
            
            if( nargin == 1 )
                capital = this.capital;
                klRatio = this.KLratio;
                wage = this.wageRequired();
            end
                
            % Labor and capital
            labor = ( capital ./ klRatio );
            
            % Total revenues
            y = this.TFP * ( (capital.^this.capitalShare) .* (labor.^(1-this.capitalShare)) );
            
            % Wage payments
            wagesPaid = wage .* labor;
            
            % Replace depreciated capital (rem: at cost to buy it)
            depreciation = this.depreciationRate .* capital .* this.priceCapital;
            
            % Risk premium (rem: at cost to buy it)
            risk = this.riskPremium .* capital .* this.priceCapital;
            
            % Investment
            investment = [capital(2:end) - (1 - this.depreciationRate)*capital(1:end-1); ...
                          capital(end) * this.invtocaps(end) ];
            
            % Investment expensing 'subsidy'
            expensing    = max(this.expensingRate .* investment .* this.priceCapital, 0);
            
            % Capital adjustment cost
            %    = eta/2 * (I/K)^2 * K
            adjustmentCost = (this.capitalAdjustmentCost/2) .* (investment .* investment) ./ capital;
            
            % Find optimal debt, interest tax benefit, and leverage cost
            [debts, debtCost, debtTaxBenefit]  = this.calculateDebt( capital );
            
            % Combine to get net tax
            if( this.firmType == Firm.SINGLEFIRM )
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
            else
                cits  = 0;
            end

            % Calculate returns to owners
            divs  = y                                   ... % revenues
                  - wagesPaid                           ... % labor costs
                  - depreciation                        ... % replace depreciated capital
                  - adjustmentCost                      ... % cap adjustment cost
                  - risk                                ... % discount money lost due to risk
                  - debtCost                            ... % Cost of leverage
                  - cits;                                   % net taxes
                  
        end % dividends
        
        
        %%
        % Calculate capital gains rates from value of capital.
        % Throw error if this will prevent open economy calculation
        function [capgains] = capitalGains( this )
            capgains        = zeros(size(this.priceCapital));
            capgains(1,1)   = (this.priceCapital(1) - this.priceCapital0)/this.priceCapital0;
            for t = 2:length(capgains)
               capgains(t,1) = (this.priceCapital(t) - this.priceCapital(t-1))/this.priceCapital(t-1);
            end
            
            % Set arbitrary upper limit to prevent model not solving for
            % open economy. Note: first period does not matter, since
            % foreigners cannot invest yet, so we do not fix r_world
            if( any( capgains(2:end) > (this.depreciationRate - 0.01) ) )
                error('MODEL ERROR! Capital gains are too high to allow OPEN economy to solve.');
            end
        end % capitalGains
        
        
        %% 
        % Calculate K/L ratio from dividend rate
        function [KLratio, caps] = calculateKLRatio( this, dividendRate, ...
                            init_caps, labor, invtocapsT_model )
            % Inputs : dividendRate = dollars received as dividend per
            %                         dollar owned of equity
            %          init_caps = capital initial guess
            %          labor = efficient units of labor (from last iteration)
            %          invtocapsT_model = last period guess of I/K
            % Outputs: KLratio that generates the dividendRate of inputs
            
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
                
                K_by_L = caps ./ labor;
                wage   = this.TFP * (1-this.capitalShare) .* (K_by_L .^ this.capitalShare);
                
                divs = this.dividends( caps, K_by_L, wage );
                divRate = divs ./ (caps .* this.priceCapital);
                
                err_div = max(abs((divRate(2:end) - dividendRate(2:end)) ./ dividendRate(2:end)));
                
            end % while
            
            % Calculate capital-labor ratio
            KLratio  = caps ./ labor;
        end % calculateKLRatio

        
        
        
        %% 
        % Calculate optimal debt 
        %   nu() = 1/nu (B/K_s)^nu , where K_s = K h p_K; h=0.1
        function [debt, debtCost, taxBenefit] = calculateDebt( this, capital )
           
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
            nu              = this.leverageCost;
            capital_scaled  = capital .* this.priceCapital .* h;
            
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
            effectiveDividendRate = ones(T_model,1) .* 0.05;
            labor = ones(T_model,1);
            init_caps = ones(T_model,1) * 12;
            invtocap_Tmodel = 0.07;
            
            [klRatio tCaps]     = this.calculateKLRatio( effectiveDividendRate   , ...
                                        init_caps, labor        , ...
                                        invtocap_Tmodel );
                        
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

