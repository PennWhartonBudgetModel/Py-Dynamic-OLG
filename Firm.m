%%
% Firms sector
%    -- currently contains only Single Firm production
%

classdef Firm
    
    properties
        TFP                 ;           % A
        capitalShare        ;           % alpha
        depreciationRate    ;           % delta
        riskPremium         ;           % for high/low return
        
        expensingRate       ;            % phi_exp
        corpTaxRate         ;           % tax on corp. profits
        
        priceCapital        ;           % See documentation. This is p_K
        priceCapital0 = 1   ;
        
    end % properties
    
    methods
        
        function [this] = Firm( scenario )
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
            
            % Investment expensing 'subsidy'
            expensing    = this.expensingRate .* investment .* this.priceCapital;
            
            % Combine to get net tax
            cits  = (y - wagesPaid - expensing) .* this.corpTaxRate;

            % Calculate returns to owners
            divs  = y                                   ... % revenues
                  - wagesPaid                           ... % labor costs
                  - depreciation                        ... % replace depreciated capital
                  - cits;                                   % net taxes
                  
        end % dividends
        
        %% 
        % Calculate K/L ratio from dividend rate
        function [KLratio] = calculateKLRatio( this, dividendRate, labor, invtocapsT_model )
            
            % Initialize variables
            T_model            = length(dividendRate);
            caps               = zeros(T_model+1,1);
            invtocaps          = zeros(T_model,1);
            invtocaps(T_model) = invtocapsT_model;
            
            % Find K at T_model
            % Total expensing divided by capital (expensing subsidy rate)
            exptocaps      = this.expensingRate(T_model) * this.corpTaxRate(T_model) ...
                             * invtocaps(T_model) * this.priceCapital(T_model);
            % MPK
            MPK            = ( dividendRate(T_model) + this.depreciationRate ...
                                * this.priceCapital(T_model) - exptocaps     ...
                              ) / (1 - this.corpTaxRate(T_model))            ;
            % KLratio
            KLratioT_model = ( MPK / (this.TFP * this.capitalShare) ) ^ (1/(this.capitalShare-1));
            % Capital at T_model
            caps(T_model)  = KLratioT_model * labor(T_model);
            
            % Find capital sequence by backward induction numerically solving
            % the polynomial: Psi*k(t) - Gamma*k(t)^alpha - Phi(k(t+1)) = 0 
            for t = T_model-1:-1:1
                
                % Constants
                Psi = ( dividendRate(t) + this.depreciationRate * this.priceCapital(t) + ...
                        this.expensingRate(t) * this.corpTaxRate(t) * this.priceCapital(t) * ...
                        (1 - this.depreciationRate) ) ... 
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

