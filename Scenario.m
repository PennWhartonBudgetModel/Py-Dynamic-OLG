%%
% Scenario definition for dynamic model execution.
%
%%
classdef Scenario
    
    properties (GetAccess = public, SetAccess = immutable)
        
        % Economy ('steady', 'open', or 'closed')
        economy;
        
        % Preference parameters
        beta;
        gamma;
        sigma;
        bequest_phi_1;

        % Core economy parameters
        is_low_return;
        modelunit_dollar;
        
        % Government expenditure shift
        expenditure_shift;
        
        % Immigration parameters
        legal_scale;
        prem_legal;
        amnesty;
        deportation;
        
        % Tax parameters
        base_brackets;
        has_buffet_rule;
        has_double_standard_deduction;
        has_limit_deductions;
        has_expand_child_credit;
        no_aca_income_tax;
        corporate_tax_rate;
        has_special_pass_through_rate;
        has_immediate_expensing;
        has_repeal_corporate_expensing;
        
        % Social Security parameters
        SSNRAPolicy;
        HasLimitedBenefits;
        SSTaxPolicy;
        SSBenefitsPolicy;
        SSBenefitsAccrualPolicy;
        
        % Identifier tags
        basedeftag
        counterdeftag
        
    end
    
    
    properties (Constant, Access = private)
    
        % Define list of required parameters
        req_params = {...
            'economy'           ;
            'beta'              ;
            'gamma'             ;
            'sigma'             ;
            'bequest_phi_1'     ;
            'is_low_return'     ;
            'modelunit_dollar'  ;
            };
        
        % Specify default values for optional parameters
        def_params = struct(...
            'expenditure_shift'             , 0.0               , ...
            'legal_scale'                   , 1.0               , ...
            'prem_legal'                    , 1.0               , ...
            'amnesty'                       , 0.0               , ...
            'deportation'                   , 0.0               , ...
            'base_brackets'                 , 'CurrentPolicy'   , ...
            'has_buffet_rule'               , false             , ...
            'has_double_standard_deduction' , false             , ...
            'has_limit_deductions'          , false             , ...
            'has_expand_child_credit'       , false             , ...
            'no_aca_income_tax'             , false             , ...
            'corporate_tax_rate'            , 0.35              , ...
            'has_special_pass_through_rate' , false             , ...
            'has_immediate_expensing'       , false             , ...
            'has_repeal_corporate_expensing', false             , ...
            'SSNRAPolicy'                   , 'CurrentPolicy'   , ...
            'HasLimitedBenefits'            , true              , ...
            'SSTaxPolicy'                   , 'CurrentPolicy'   , ...
            'SSBenefitsPolicy'              , 'CurrentPolicy'   , ...
            'SSBenefitsAccrualPolicy'       , 'CurrentPolicy'     ...
        );

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
            
            % Validate that there are no unused parameters
            for o = fieldnames(params)'
                
                if ( ~isfield(Scenario.def_params, o{1})     ...
                    && ~any(strcmp(o{1}, Scenario.req_params)) )
                    error( 'Field <%s> does not match any Scenario fields.', o{1} );
                end
                
            end
            
            % Validate economy
            assert(any(strcmp(this.economy, {'steady', 'open', 'closed'})), ...
                '<economy> parameter must be either ''steady'', ''open'', or ''closed''.');
            
            
            
            % Generate identifier tags for baseline and counterfactual definitions
            this.basedeftag = strjoin({...
                sprintf('%.3f'  , this.beta                 ),...
                sprintf('%.3f'  , this.gamma                ), ...
                sprintf('%.2f'  , this.sigma                ), ...
                sprintf('%.1f'  , this.modelunit_dollar*1e6 ), ...
                sprintf('%u'    , this.is_low_return        ), ...
                sprintf('%.3f'  , this.bequest_phi_1        )}, '_');
            
            if this.isCurrentPolicy()
                this.counterdeftag = 'currentpolicy';
            else
                this.counterdeftag = strjoin({...
                    sprintf('%.2f'  , abs(this.expenditure_shift)           ), ...
                    sprintf('%.1f'  , this.legal_scale                      ), ...
                    sprintf('%.1f'  , this.prem_legal                       ), ...
                    sprintf('%.1f'  , this.amnesty                          ), ...
                    sprintf('%.1f'  , this.deportation                      ), ...
                    sprintf('%s'    , this.base_brackets                    ), ...
                    sprintf('%u'    , this.has_buffet_rule                  ), ...
                    sprintf('%u'    , this.has_double_standard_deduction    ), ...
                    sprintf('%u'    , this.has_limit_deductions             ), ...
                    sprintf('%u'    , this.has_expand_child_credit          ), ...
                    sprintf('%u'    , this.no_aca_income_tax                ), ...
                    sprintf('%.2f'  , this.corporate_tax_rate               ), ...
                    sprintf('%u'    , this.has_special_pass_through_rate    ), ...
                    sprintf('%u'    , this.has_immediate_expensing          ), ...
                    sprintf('%u'    , this.has_repeal_corporate_expensing   )}, '_');
            end
            
        end
        
        
        % Identify if scenario is equivalent to another scenario
        %   Parameter representations in tags determine precision for equivalency evaluation
        function [flag] = isEquivalent(this, scenario)
            flag = strcmp(this.economy      , scenario.economy      ) ...
                && strcmp(this.basedeftag   , scenario.basedeftag   ) ...
                && strcmp(this.counterdeftag, scenario.counterdeftag);
        end
        
        
        % Identify if scenario represents current policy
        %   Current policy identified by default values for all optional parameters
        function [flag] = isCurrentPolicy(this)
            flag = true;
            for o = fieldnames(Scenario.def_params)'
                if ~isequal(this.(o{1}), Scenario.def_params.(o{1}))
                    flag = false;
                    return;
                end
            end
        end
        
        
        % Identify if scenario represents steady state
        function [flag] = isSteady(this)
            flag = strcmp(this.economy, 'steady');
        end
        
        
        % Identify if scenario represents open economy
        function [flag] = isOpen(this)
            flag = strcmp(this.economy, 'open');
        end
        
        
        % Identify if scenario represents closed economy
        function [flag] = isClosed(this)
            flag = strcmp(this.economy, 'closed');
        end
        
        
        % Get all scenario parameters
        function [params] = getParams(this)
            params = struct();
            for f = [ Scenario.req_params ; fieldnames(Scenario.def_params) ]'
                params.(f{1}) = this.(f{1});
            end
        end
        
        
        % Generate corresponding current policy scenario
        function [scenario] = currentPolicy(this)
            params = this.getParams();
            for f = fieldnames(Scenario.def_params)'
                params.(f{1}) = Scenario.def_params.(f{1});
            end
            scenario = Scenario(params);
        end
        
        
        % Generate corresponding steady state scenario
        function [scenario] = steady(this)
            params = this.getParams();
            params.economy = 'steady';
            scenario = Scenario(params);
        end
        
        
        % Generate corresponding open economy scenario
        function [scenario] = open(this)
            params = this.getParams();
            params.economy = 'open';
            scenario = Scenario(params);
        end
        
        
        % Generate corresponding closed economy scenario
        function [scenario] = closed(this)
            params = this.getParams();
            params.economy = 'closed';
            scenario = Scenario(params);
        end
        
        
    end
    
    
end

