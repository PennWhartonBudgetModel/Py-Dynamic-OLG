%%
%   Scenarios are runs of the model
%       We currently define "current policy" by defaults
%       in the constructor. In the future, we could pass in the "current
%       policy".
%
classdef Scenario
    
    properties (SetAccess = protected )
        ID;
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
        
    end % properties

    properties (Constant, Access = private )
    
        % These are _currently_ required fields
        req_fields = { 'beta', 'sigma', 'gamma', 'modelunit_dollar', 'economy' };
        
        % These are _current_ defaults for other fields
        def_fields = struct(    'ID'            , []        ...
                            ,   'taxplan'       , 'base'    ...
                            ,   'gcut'          , 0.0       ...
                            ,   'legal_scale'   , 1.0       ...
                            ,   'prem_legal'    , 1.0       ...
                            ,   'amnesty'       , 0.0       ...
                            ,   'deportation'   , 0.0       ...
                            );
    end % properties
    
    methods
        
        % Constructor
        function this = Scenario(props)
            
            if( ~isa( props, 'struct' ) )
                error( 'Scenario constructor expects struct of init variables.')
            end 
            
            % Check for required fields
            for f = Scenario.req_fields
                if( ~isfield( props, f ) )
                    error( 'Scenario constructor requires <%s> field.', string(f) );
                end
                if( isempty(props.(f{1})) )
                    error( 'Field <%s> cannot be empty.', string(f) );
                end
            end
            if( isempty(strmatch(props.economy, {'steady', 'open', 'closed'} )))
                    error( '<economy> field must be either "steady", "open", or "closed"' );
            end
            
            % Set defaults and then possibly overwrite  
            for f = fieldnames(Scenario.def_fields)'
                this.(f{1}) = Scenario.def_fields.(f{1});
            end
            
            % Set fields as passed 
            for f = fieldnames(props)'
                this.(f{1}) = props.(f{1});
            end
        
            % TODO: Warn if some fields from props not used
        end % Scenario()
        
        
        % Derived property: Get if isCurrentPolicy
        function flag = isCurrentPolicy(this)
             % isCurrentPolicy if no deviation from def_fields  
            for f = fieldnames(Scenario.def_fields)'
                this_one = (this.(f{1})); def_one = (Scenario.def_fields.(f{1}));
                if( length(this_one) ~= length(def_one) ) % Needed since Matlab2016 does not overload == for strings
                    flag = false;
                    return;
                elseif (this_one ~= def_one )
                    flag = false;
                    return;
                end
            end
            flag = true;
        end % isCurrentPolicy
        
        
        % Return cloned Scenario with override to current policy
        function obj = currentPolicy(this)
            obj = this.Clone();
            for f = fieldnames(Scenario.def_fields)'
                obj.(f{1}) = Scenario.def_fields.(f{1});
            end
        end % currentPolicy
        
        
        % Return cloned Scenario with 'steady' economy
        function obj = steady(this)
            obj = this.Clone();
            obj.economy = 'steady';
        end % steady
        
        
        % Return cloned Scenario with 'open' economy
        function obj = open(this)
            obj = this.Clone();
            obj.economy = 'open';
        end % open
        
        
        % Return cloned Scenario with 'closed' economy
        function obj = closed(this)
            obj = this.Clone();
            obj.economy = 'closed';
        end % open
        
        
        % Generate tags for baseline and counterfactual definitions
        %     NOTE: These tags are not guaranteed to be unique across
        %     Scenarios. These are for Development. For Production, we use
        %     the Scenario.ID
        function [basedef_tag, counterdef_tag] = generate_tags(this)

            str = [];
            str = [str, sprintf('%.3f' , this.beta)            ];
            str = [str, sprintf('_%.3f', this.gamma)           ];
            str = [str, sprintf('_%.2f', this.sigma)           ];
            str = [str, sprintf('_%e'  , this.modelunit_dollar)];
            basedef_tag = str;
            
            if( this.isCurrentPolicy )
                counterdef_tag = 'baseline';
            else
                str = [];
                str = [str, sprintf('%s'    , this.taxplan)        ];
                str = [str, sprintf('_%+.2f', this.gcut)           ];
                str = [str, sprintf('_%.1f' , this.legal_scale)    ];
                str = [str, sprintf('_%.3f' , this.prem_legal)     ];
                str = [str, sprintf('_%.2f' , this.amnesty)        ];
                str = [str, sprintf('_%.2f' , this.deportation)    ];
                counterdef_tag = str;
            end
        end % generate_tags

        
    end % methods
    
    methods (Access = private )
        
        function obj = Clone(this)
            % Make "empty" object to fill in
            params = struct();
            for f = Scenario.req_fields
                params.(f{1}) = 'open';  % trick to skip checks in constructor
            end
            obj = Scenario(params);
            obj.ID                  = this.ID;
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
        
    end % methods
end % Scenario

