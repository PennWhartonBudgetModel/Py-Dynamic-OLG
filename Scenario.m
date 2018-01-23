%%
% Scenario definition for dynamic model execution.
%
%%
classdef Scenario
    
    properties (GetAccess = public, SetAccess = immutable)
        
        % Identifier tags
        basedeftag
        counterdeftag
        
        %% REQUIRED parameters
        
        % Economy ('steady', 'open', or 'closed')
        economy;
        
        % Preference parameters
        beta;
        gamma;
        sigma;
        bequest_phi_1;

        % Conversion to match US dollar amounts
        modelunit_dollar;

        % Core economy parameters
        IsLowReturn;

        % Timing
        TransitionFirstYear;
        TransitionLastYear;
        

        % OPTIONAL policy parameters
        
        % Government expenditure shift
        ExpenditureShift;
        
        % Immigration parameters
        legal_scale;
        prem_legal;
        amnesty;
        deportation;
        
        % Tax parameters
        BaseBrackets;
        HasBuffetRule;
        HasDoubleStandardDeduction;
        HasLimitDeductions;
        HasExpandChildCredit;
        NoACAIncomeTax;
        CorporateTaxRate;
        HasSpecialPassThroughRate;
        HasImmediateExpensing;
        HasRepealCorporateExpensing;
        
        % Social Security parameters
        SSNRAPolicy;
        SSTaxPolicy;
        SSBenefitsPolicy;
        SSBenefitsAccrualPolicy;
        
    end
    
    
    properties (Constant, Access = private)
    
        % Define list of required parameters
        req_params = {...
            'economy'               ;
            'beta'                  ;
            'gamma'                 ;
            'sigma'                 ;
            'bequest_phi_1'         ;
            'modelunit_dollar'      ;
            'IsLowReturn'           ;
            'TransitionFirstYear'   ;
            'TransitionLastYear'    ;   
            };
        
        % Specify default values for optional parameters
        def_params = struct(...
            'ExpenditureShift'              , 0.0               , ...
            'legal_scale'                   , 1.0               , ...
            'prem_legal'                    , 1.0               , ...
            'amnesty'                       , 0.0               , ...
            'deportation'                   , 0.0               , ...
            'BaseBrackets'                  , 'CurrentPolicy'   , ...
            'HasBuffetRule'                 , false             , ...
            'HasDoubleStandardDeduction'    , false             , ...
            'HasLimitDeductions'            , false             , ...
            'HasExpandChildCredit'          , false             , ...
            'NoACAIncomeTax'                , false             , ...
            'CorporateTaxRate'              , 0.35              , ...
            'HasSpecialPassThroughRate'     , false             , ...
            'HasImmediateExpensing'         , false             , ...
            'HasRepealCorporateExpensing'   , false             , ...
            'SSNRAPolicy'                   , 'CurrentPolicy'   , ...
            'SSTaxPolicy'                   , 'CurrentPolicy'   , ...
            'SSBenefitsPolicy'              , 'CurrentPolicy'   , ...
            'SSBenefitsAccrualPolicy'       , 'CurrentPolicy'     ...
        );

    end
    
    
    methods (Access = public)
        
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
            this.basedeftag = '';
            for o = Scenario.req_params'
                if( ~strcmp( o{1}, 'economy' ) )
                    this.basedeftag = strcat( this.basedeftag, '_', num2str(this.(o{1})) );
                end
            end
            if this.isCurrentPolicy()
                this.counterdeftag = 'currentpolicy';
            else
                for o = fieldnames(Scenario.def_params)'
                    this.counterdeftag = strcat( this.counterdeftag, '_', num2str(this.(o{1})) );
                end
            end
            
        end % Scenario constructor
        
        
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

    
    
    methods (Access = private)
        
    end

    
end % Scenario


