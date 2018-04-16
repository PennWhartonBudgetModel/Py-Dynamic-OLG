%%
% Firms sector
%    -- currently contains only Single Firm production
%

classdef Firm
    
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
        corpTaxRate         ;           % tax on corp. profits
        
        priceCapital        ;           % See documentation. This is p_K
        priceCapital0 = 1   ;
        
        firmType            ;           % from the enumeration
        
    end % properties
    
    methods
        
        function [this] = Firm( scenario, firmType )
            
            this.firmType           = firmType;   
            if( ~(firmType == Firm.SINGLEFIRM || firmType == Firm.PASSTHROUGH) )
                throw(MException('Firm:firmType','firmType must be SingleFirm or PassThrough'));
            end
            
            prod                    = ParamGenerator.production( scenario );
            this.TFP                = prod.A;
            this.capitalShare       = prod.alpha;
            this.depreciationRate   = prod.depreciation;
            this.riskPremium        = prod.risk_premium;
            
            tax                 = ParamGenerator.tax( scenario );
            this.expensingRate  = tax.shareCapitalExpensing;        
        	this.corpTaxRate    = tax.rateCorporate;     
            
            % Calculate the price of capital (p_K, see docs)
            this.priceCapital   = this.priceCapital0 * (tax.qtobin ./ tax.qtobin0);
        end % constructor
        
        
        %%
        function [divs, cits] = dividends( this, capital, investment, klRatio, wage )
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
            
            % Investment expensing 'subsidy'
            expensing    = this.expensingRate .* investment .* this.priceCapital;
            
            % Combine to get net tax
            if( this.firmType == Firm.SINGLEFIRM )
                cits  = (y - wagesPaid - expensing) .* this.corpTaxRate;
            else
                cits  = 0;
            end

            % Calculate returns to owners
            divs  = y                                   ... % revenues
                  - wagesPaid                           ... % labor costs
                  - depreciation                        ... % replace depreciated capital
                  - risk                                ... % discount money lost due to risk
                  - cits;                                   % net taxes
                  
        end % dividends
        
        %% 
        % Calculate K/L ratio from dividend rate
        function [KLratio] = calculateKLRatio( this, dividendRate, labor, invtocapsT_model )
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
                       ) * this.priceCapital(T_model) ...
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
            
        end % calculateKLRatio
        
    end % instance methods
    
end % classdef

