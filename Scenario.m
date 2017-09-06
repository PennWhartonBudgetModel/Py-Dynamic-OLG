%%
%   Scenarios are runs of the model
%     
classdef Scenario
    
    properties
        ID;
        baselineScenario;
        economy;
        
        beta;
        gamma;
        sigma;
        modelunit_dollar;
        
        taxplan;
        gcut;
        legal_scale;
        prem_legal;
        amnesty;
        deportation;
    end
    
    methods
        
        % Constructors
        function this = Scenario()
            
            this.ID                 = [];
            this.baselineScenario   = [];
            this.economy            = 'steady';
            
            % these params are from 7/17/2017 
            %   with current policy, 0.5 labor elasticity, 0.5 savings
            %   elasticity
            this.beta               = 1.076897;
            this.sigma              = 7.063790;
            this.gamma              = 0.589655;
            this.modelunit_dollar   = 0.000042;
                 
            this.taxplan            = 'base';
            this.gcut               = 0.00;
            this.legal_scale        = 1.0;
            this.prem_legal         = 1.000;
            this.amnesty            = 0.00;
            this.deportation        = 0.00;
        
        end % Scenario()
        
        
        function flag = isBase(this)
            flag = isempty(this.baselineScenario);
        end %isBase
        
        
        function obj = Clone(this)
            obj = Scenario();
            obj.ID                  = this.ID;
            obj.baselineScenario    = this.baselineScenario;
            obj.economy             = this.economy;
            
            obj.beta                = this.beta;
            obj.sigma               = this.sigma;
            obj.gamma               = this.gamma;
            obj.modelunit_dollar    = this.modelunit_dollar;
                 
            obj.taxplan             = this.taxplan;
            obj.gcut                = this.gcut;
            obj.legal_scale         = this.legal_scale;
            obj.prem_legal          = this.prem_legal;
            obj.amnesty             = this.amnesty;
            obj.deportation         = this.deportation;
        end % Clone
        
        
        % Generate tags for baseline and counterfactual definitions
        %     NOTE: These tags are not guaranteed to be unique across
        %     Scenarios. These are for Development. For Production, we use
        %     the Scenario.ID
        function [basedef_tag, counterdef_tag] = generate_tags(this)

            s = [];
            s = [s, sprintf('%.3f' , this.beta)            ];
            s = [s, sprintf('_%.3f', this.gamma)           ];
            s = [s, sprintf('_%.2f', this.sigma)           ];
            s = [s, sprintf('_%e'  , this.modelunit_dollar)];
            basedef_tag = s;
            
            if( isempty(this.baselineScenario) )
                counterdef_tag = 'baseline';
            else
                s = [];
                s = [s, sprintf('%s'    , this.taxplan)        ];
                s = [s, sprintf('_%+.2f', this.gcut)           ];
                s = [s, sprintf('_%.1f' , this.legal_scale)    ];
                s = [s, sprintf('_%.3f' , this.prem_legal)     ];
                s = [s, sprintf('_%.2f' , this.amnesty)        ];
                s = [s, sprintf('_%.2f' , this.deportation)    ];
                counterdef_tag = s;
            end
        end % generate_tags

        
    end % methods
    
end % Scenario

