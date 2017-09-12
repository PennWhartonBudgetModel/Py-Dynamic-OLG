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
        bequest_phi_1;
        
        taxplan;
        gcut;
        legal_scale;
        prem_legal;
        amnesty;
        deportation;
        
    end % properties

    properties (Constant, Access = private )
    
        % These are _currently_ required fields
        req_fields = { 'beta', 'sigma', 'gamma', 'modelunit_dollar', 'bequest_phi_1', 'economy' };
        
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
            %   Rem: This will fail if non-matching (extra) field
            for f = fieldnames(props)'
                this.(f{1}) = props.(f{1});
            end
             
        end % Scenario()
        
        
        % Derived property: Get if isCurrentPolicy
        function flag = isCurrentPolicy(this)
             % isCurrentPolicy if no deviation from def_fields  
            for f = fieldnames(Scenario.def_fields)'
                this_one = (this.(f{1})); def_one = (Scenario.def_fields.(f{1}));
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
            params          = this.getParams();
            params.economy  = 'steady';
            obj             = Scenario(params);
        end % steady
        
        
        % Return cloned Scenario with 'open' economy
        function obj = open(this)
            params          = this.getParams();
            params.economy  = 'open';
            obj             = Scenario(params);
        end % open
        
        
        % Return cloned Scenario with 'closed' economy
        function obj = closed(this)
            params          = this.getParams();
            params.economy  = 'closed';
            obj             = Scenario(params);
        end % open
        
        
        % Generate tags for baseline and counterfactual definitions
        %     NOTE: These tags are not guaranteed to be unique across
        %     Scenarios. These are for Development. For Production, we use
        %     the Scenario.ID
        function [basedef_tag, counterdef_tag] = generate_tags(this)

            basedef_tag = [     sprintf('%.3f' , this.beta)             ...
                            ,   sprintf('_%.3f', this.gamma)            ...     
                            ,   sprintf('_%.2f', this.sigma)            ...
                            ,   sprintf('_%e'  , this.modelunit_dollar) ...
                            ,   sprintf('_%.3f', this.bequest_phi_1)    ...
                          ];
            
            if( this.isCurrentPolicy )
                counterdef_tag = 'baseline';
            else
                counterdef_tag  = [     sprintf('%s'    , this.taxplan)        ...
                                    ,   sprintf('_%+.2f', this.gcut)           ...
                                    ,   sprintf('_%.1f' , this.legal_scale)    ...
                                    ,   sprintf('_%.3f' , this.prem_legal)     ...
                                    ,   sprintf('_%.2f' , this.amnesty)        ...
                                    ,   sprintf('_%.2f' , this.deportation)    ...
                                  ];
            end
        end % generate_tags


        function params = getParams(this)
            params = struct();
            % Create and copy all fields
            for f = Scenario.req_fields
                params.(f{1}) = this.(f{1});
            end
            for f = fieldnames(Scenario.def_fields)'
                params.(f{1}) = this.(f{1});
            end
        end % getParams
        
    end % methods
end % Scenario

