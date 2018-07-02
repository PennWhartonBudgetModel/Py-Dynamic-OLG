%%
% Firms sector
%    -- currently contains only Single Firm production
%

classdef Firm < handle
    
    properties (Constant)
        SINGLEFIRM          = 0;
        PASSTHROUGH         = 1;
        MULTIFIRM           = 2;
    end
    
    properties 
        TFP                 ;           % A
        capitalShare        ;           % alpha
        depreciationRate    ;           % delta
        riskPremium         ;           % for high/low return
        
        expensingRate       ;           % phi_exp
        effectiveTaxRate    ;           % effective tax on profits
        statutoryTaxRate    ;           % statutory tax rate on profits
        interestDeduction   ;           % phi_int
        leverageCost        ;           % nu
        
        initLeverageRatio   ;           % stored value of initial leverage ratio 
        debtTaxBenefitRate  ;           % pre-calculated tax value of another $1 of debt
        
        priceCapital        ;           % See documentation. This is p_K
        priceCapital0 = 1   ;
        
        firmType            ;           % from the enumeration
        
    end % properties
    
    methods
        
        %%
        %  Constructor
        %     INPUTS:   paramsTax = ParamGenerator.tax()
        %               paramsProduction = ParamGenerator.production()
        %               interestRate = Guess at dividends rate
        function [this] = Firm( paramsTax, paramsProduction, interestRate, firmType )
            
            this.firmType           = firmType;   
            if( ~(firmType == Firm.SINGLEFIRM || firmType == Firm.PASSTHROUGH) )
                throw(MException('Firm:firmType','firmType must be SingleFirm or PassThrough'));
            end
            
            this.TFP                = paramsProduction.A;
            this.capitalShare       = paramsProduction.alpha;
            this.depreciationRate   = paramsProduction.depreciation;
            this.riskPremium        = paramsProduction.risk_premium;
            
            switch this.firmType
                case Firm.SINGLEFIRM
                    this.effectiveTaxRate   = paramsTax.rateCorporate;
                    this.statutoryTaxRate   = paramsTax.rateCorporateStatutory;
                    this.expensingRate      = paramsTax.shareCapitalExpensing; % REVISE W/ new interface
                    this.interestDeduction  = paramsTax.interestDeduction;
                    this.initLeverageRatio  = paramsProduction.initialCorpLeverage;
                case Firm.PASSTHROUGH
                    this.effectiveTaxRate   = 0;  % TEMP: Should come from ParamGenerator as top marginal PIT rate
                    this.statutoryTaxRate   = this.effectiveTaxRate;
                    this.expensingRate      = paramsTax.shareCapitalExpensing; % REVISE W/ new interface
                    this.interestDeduction  = paramsTax.interestDeduction;
                    this.initLeverageRatio  = paramsProduction.initialPassThroughLeverage;
            end
            
            this.findLeverageCost( interestRate );
            
            % Calculate the price of capital (p_K, see docs)
            this.priceCapital   = this.priceCapital0 * (paramsTax.qtobin ./ paramsTax.qtobin0);
        end % constructor
        
        
        %%
        %  Reset the interest rate and recalculate leverage cost
        function findLeverageCost( this, interestRate )
            
            if( interestRate <= 0 )
                throw(MException('LEVERAGE_COST:NOT_POSITIVE','Interest rate must be positive for leverage cost calculation.'));
            end
            
            % Calculate or use leverage cost
            % Rem: the leverage cost is size invariant, so set capital=1
            %      also, capital from initLeverageRatio is in $, so already
            %      scaled
            this.debtTaxBenefitRate = (this.interestDeduction .* this.statutoryTaxRate .* interestRate);
            this.leverageCost       = this.calculateLeverageCost( this.initLeverageRatio, 1);
        end % resetInterestRate
        
        
        %%
        function [divs, cits, debts] = dividends( this, capital, invtocapsT_model, klRatio, wage )
            % Inputs : capital
            %          investment = GROSS physical investment --> K' - (1-d)K
            %          klRatio & wage (which are consistent in the current iteration)
            % Outputs: divs = dividend is the return of the corporation net of tax.
            %          cits = corporate income taxes payed by the firms
            
            % Labor
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
                          capital(end)*invtocapsT_model ];
            
            % Investment expensing 'subsidy'
            expensing    = this.expensingRate .* investment .* this.priceCapital;
            
            % Find optimal debt, interest tax benefit, and leverage cost
            [debts, debtCost, debtTaxBenefit]  = this.calculateDebt( capital );
            
            % Combine to get net tax
            if( this.firmType == Firm.SINGLEFIRM )
                cits  = (y - wagesPaid) .* this.effectiveTaxRate  ...
                        - expensing .* this.statutoryTaxRate      ...
                        - debtTaxBenefit;
            else
                cits  = 0;
            end

            % Calculate returns to owners
            divs  = y                                   ... % revenues
                  - wagesPaid                           ... % labor costs
                  - depreciation                        ... % replace depreciated capital
                  - risk                                ... % discount money lost due to risk
                  - debtCost                            ... % Cost of leverage
                  - cits;                                   % net taxes
                  
        end % dividends
        
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
                
                investment = [caps(2:end) - (1 - this.depreciationRate)*caps(1:end-1); ...
                              caps(end)*invtocapsT_model ];
                
                % Update capital 
                %  if divRate > dividendRate 
                %    --> caps gets bigger, and divRate gets smaller
                caps(2:end) = caps(2:end) .* ((1+divRate(2:end)) ./ (1+dividendRate(2:end)) );
                
                K_by_L = caps ./ labor;
                wage   = this.TFP * (1-this.capitalShare) .* (K_by_L .^ this.capitalShare);
                
                divs = this.dividends( caps, invtocapsT_model, K_by_L, wage );
                divRate = divs ./ (caps .* this.priceCapital);
                
                err_div = max(abs((divRate(2:end) - dividendRate(2:end)) ./ dividendRate(2:end)));
                
            end % while
            
            % Calculate capital-labor ratio
            KLratio  = caps ./ labor;
        end % calculateKLRatio

        
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
        
        
        %% 
        % Calculate K/L ratio from dividend rate
        function [KLratio, tempCaps] = old_calculateKLRatio( this, dividendRate, labor, invtocapsT_model )
            % Inputs : dividendRate = dollars received as dividend per
            %                         dollar owned of equity
            %          labor = efficient units of labor (from last iteration)
            %          invtocapsT_model = last period guess of I/K
            % Outputs: KLratio that generates the exact dividendRate of inputs
            
            % Initialize variables
            T_model            = length(dividendRate);
            caps               = zeros(T_model+1,1);
            invtocaps          = zeros(T_model,1);
            invtocaps(T_model) = invtocapsT_model;
            
            % Find K at T_model
            % Total expensing divided by capital (expensing subsidy rate)
            exptocaps      = this.expensingRate(T_model) * this.corpTaxRate(T_model) ...
                             * invtocaps(T_model);
            % MPK
            MPK            = ( dividendRate(T_model) + this.depreciationRate ...
                                + this.riskPremium - exptocaps     ...
                              ) * this.priceCapital(T_model) / (1 - this.corpTaxRate(T_model))            ;
            % KLratio
            KLratioT_model = ( MPK / (this.TFP * this.capitalShare) ) ^ (1/(this.capitalShare-1));
            % Capital at T_model
            caps(T_model)  = KLratioT_model * labor(T_model);
            
            % Find capital sequence by backward induction numerically solving
            % the polynomial: Psi*k(t) - Gamma*k(t)^alpha - Phi(k(t+1)) = 0 
            for t = T_model-1:-1:1
                
                % Constants
                Psi = ( dividendRate(t) + this.depreciationRate + this.riskPremium ...
                       + this.expensingRate(t) * this.corpTaxRate(t) * (1 - this.depreciationRate) ... 
                       ) * this.priceCapital(t) ...
                      / ( this.TFP * this.capitalShare * (1 - this.corpTaxRate(t)) );
                Phi = caps(t+1) * ...
                      ( ( this.expensingRate(t) * this.corpTaxRate(t) * this.priceCapital(t) ) ...
                      / (  this.TFP * this.capitalShare * (1 - this.corpTaxRate(t)) ));                
                Gamma = labor(t)^(1 - this.capitalShare);
                
                % Find polynomial root
                % caps(t+1) is used as the guess
                [x, ~, exitflag] = fzero(@(x) (x - (Gamma/Psi)*(x^this.capitalShare) - Phi/Psi), caps(t+1) );
                
                if ( exitflag ~= 1)
                    error('No root to the polynomial.')
                end
                
                % Save time t capital
                caps(t) = x;
               
            end
            
            % Calculate capital-labor ratio
            KLratio = caps(1:T_model) ./ labor;
            tempCaps = caps(1:T_model);
        end % calculateKLRatio

        
        %% 
        % Calculate optimal debt 
        %   nu() = 1/nu (B/K_s)^nu , where K_s = K h p_K; h=0.1
        function [debt, debtCost, taxBenefit] = calculateDebt( this, capital )
            
            h              = 100;
            taxBenefitRate = this.debtTaxBenefitRate;
            nu             = this.leverageCost;
            capital_scaled = capital .* this.priceCapital .* h;
            
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
            capital_scaled  = capital_value .* h;
            
            logBK = log( debt ./ capital_scaled );
            n     = log( this.debtTaxBenefitRate ) ...
                      - log( 1 - this.statutoryTaxRate .* this.interestDeduction );
            nu    = 1 + n ./ logBK;

        end % calculateLeverageCost


        
    end % instance methods
    
end % classdef

