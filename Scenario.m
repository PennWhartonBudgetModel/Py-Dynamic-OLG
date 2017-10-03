%%
% Scenario definition for dynamic model execution.
%
%%
classdef Scenario
    
    properties (SetAccess = protected )
        
        % Economy ('steady', 'open', or 'closed')
        economy;
        
        % Preference parameters
        beta;
        gamma;
        sigma;
        modelunit_dollar;
        bequest_phi_1;
        
        % Government expenditure cut
        gcut;
        
        % Immigration parameters
        legal_scale;
        prem_legal;
        amnesty;
        deportation;
        
        % Tax parameters
        base_brackets;
        has_buffet_rule;
        has_agi_surcharge_5m;
        corporate_tax_rate;
        has_double_standard_deduction;
        has_limit_deductions;
        no_amt;
        has_expand_child_credit;
        no_aca_income_tax;
        
    end
    
    
    properties (Constant, Access = private )
    
        % Define list of required parameters
        req_params = {...
            'economy'           ;
            'beta'              ;
            'sigma'             ;
            'gamma'             ;
            'modelunit_dollar'  ;
            'bequest_phi_1'     };
        
        % Specify default values for optional parameters
        def_params = struct(...
            'gcut'                          , 0.0               , ...
            'legal_scale'                   , 1.0               , ...
            'prem_legal'                    , 1.0               , ...
            'amnesty'                       , 0.0               , ...
            'deportation'                   , 0.0               , ...
            'base_brackets'                 , 'CurrentPolicy'   , ...
            'has_buffet_rule'               , false             , ...
            'has_agi_surcharge_5m'          , false             , ...
            'corporate_tax_rate'            , 0.35              , ...
            'has_double_standard_deduction' , false             , ...
            'has_limit_deductions'          , false             , ...
            'no_amt'                        , false             , ...
            'has_expand_child_credit'       , false             , ...
            'no_aca_income_tax'             , false             );
        
    end
    
    
    methods
        
        % Constructor
        function [this] = Scenario(params)
            
            assert(isa(params, 'struct'), 'Scenario constructor expects structure of parameter values.');
            
            % Set required parameters
            for o = Scenario.req_params'
                
                assert(isfield(params, o{1}) && ~isempty(params.(o{1})), ...
                    sprintf('Scenario constructor requires nonempty <%s> parameter.', o{1}));
                
                this.(o{1}) = params.(o{1});
                
            end
            
            % Set optional parameters, using defaults where unspecified
            for o = fieldnames(Scenario.def_params)'
                
                if isfield(params, o{1})
                    this.(o{1}) = params.(o{1});
                else
                    this.(o{1}) = Scenario.def_params.(o{1});
                end
                
            end
            
            % Validate economy
            assert(any(strcmp(this.economy, {'steady', 'open', 'closed'})), ...
                '<economy> parameter must be either ''steady'', ''open'', or ''closed''.');
            
        end
        
        
        % Derived property: Get if isCurrentPolicy
        function [flag] = isCurrentPolicy(this)
             % isCurrentPolicy if no deviation from def_params  
            for f = fieldnames(Scenario.def_params)'
                this_one = (this.(f{1})); def_one = (Scenario.def_params.(f{1}));
                switch (class(this_one))
                    case 'double'
                        if (this_one ~= def_one )
                            flag = false;
                            return;
                        end
                    case 'char'
                        if( ~strcmp(this_one, def_one) )
                            flag = false;
                            return;
                        end
                end
            end
            flag = true;
        end
        
        
        % Generate corresponding current policy scenario
        function [scenario] = currentPolicy(this)
            params          = this.getParams();
            for f = fieldnames(Scenario.def_params)'
                params.(f{1}) = Scenario.def_params.(f{1});
            end
            scenario        = Scenario(params);
        end
        
        
        % Generate corresponding steady state scenario
        function [scenario] = steady(this)
            params          = this.getParams();
            params.economy  = 'steady';
            scenario        = Scenario(params);
        end
        
        
        % Generate corresponding open economy scenario
        function [scenario] = open(this)
            params          = this.getParams();
            params.economy  = 'open';
            scenario        = Scenario(params);
        end
        
        
        % Generate corresponding closed economy scenario
        function [scenario] = closed(this)
            params          = this.getParams();
            params.economy  = 'closed';
            scenario        = Scenario(params);
        end
        
        
        % Generate tags for baseline and counterfactual definitions
        function [basedef_tag, counterdef_tag] = generate_tags(this)

            basedef_tag = [...
                sprintf('_%.3f' , this.beta             ),...
                sprintf('_%.3f' , this.gamma            ), ...
                sprintf('_%.2f' , this.sigma            ), ...
                sprintf('_%e'   , this.modelunit_dollar ), ...
                sprintf('_%.3f' , this.bequest_phi_1    )];
            basedef_tag(1) = '';
            
            if this.isCurrentPolicy()
                counterdef_tag = 'baseline';
            else
                counterdef_tag = [...
                    sprintf('_%+.2f', this.gcut                             ), ...
                    sprintf('_%.1f' , this.legal_scale                      ), ...
                    sprintf('_%.2f' , this.prem_legal                       ), ...
                    sprintf('_%.2f' , this.amnesty                          ), ...
                    sprintf('_%.2f' , this.deportation                      ), ...
                    sprintf('_%s'   , this.base_brackets                    ), ...
                    sprintf('_%u'   , this.has_buffet_rule                  ), ...
                    sprintf('_%u'   , this.has_agi_surcharge_5m             ), ...
                    sprintf('_%.2f' , this.corporate_tax_rate               ), ...
                    sprintf('_%u'   , this.has_double_standard_deduction    ), ...
                    sprintf('_%u'   , this.has_limit_deductions             ), ...
                    sprintf('_%u'   , this.no_amt                           ), ...
                    sprintf('_%u'   , this.has_expand_child_credit          ), ...
                    sprintf('_%u'   , this.no_aca_income_tax                )];
                counterdef_tag(1) = '';
            end
        end
        
        
        
        
        function [params] = getParams(this)
            params = struct();
            for f = [ Scenario.req_params ; fieldnames(Scenario.def_params) ]'
                params.(f{1}) = this.(f{1});
            end
        end
        
        function [flag] = isEqual(this, scenario)
            flag = isequal(this.getParams(), scenario.getParams());
        end
        
    end
    
    
end

