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
        
        expensingRate                   % phi_exp
        corpTaxRate                     % tax on corp. profits
        
        priceCapital                    % See documentation. This is p_K
        priceCapital0 = 1;
        
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
            this.priceCapital   = tax.qtobin ./ tax.qtobin0;
        end % constructor
        
        
        % dividend is the return of the corporation 
        % net of tax.
        %   NOTE: investment is GROSS physical investment --> K' - (1-d)K
        function [divs, cits] = dividends( this, capital, labor, investment, wage )
            
            % Total revenues
            y = this.TFP * ( capital.^this.capitalShare .* labor.^(1-this.capitalShare));
            
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
        
    end % instance methods
    
end % classdef

