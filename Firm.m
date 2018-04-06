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
        function [KLratio] = calculateKLRatio( this, dividendRate, capital, investment )
            
            expensingSubsidyRate = this.expensingRate .* this.corpTaxRate ...
                                .* (investment ./ capital) ...
                                .* this.priceCapital;
            % Calculate MPK 
            r  = ( dividendRate ...
                   + this.depreciationRate .* this.priceCapital ...
                   - expensingSubsidyRate ...
                  ) ./ (1 - this.corpTaxRate);
              
            % Calculate K/L ratio from MPK
            KLratio = ( r ./ (this.TFP * this.capitalShare) ) .^ (1/(this.capitalShare-1));
        end % calculateKLRatio
        
    end % instance methods
    
end % classdef

