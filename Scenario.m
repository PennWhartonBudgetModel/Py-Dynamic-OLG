%%
% Scenario definition for dynamic model execution.
%
%%
classdef Scenario
    
    properties (SetAccess = protected )
        
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
        
    end
    
    
    properties (Constant, Access = private )
    
        % Define list of required fields
        req_fields = {...
            'economy'           ;
            'beta'              ;
            'sigma'             ;
            'gamma'             ;
            'modelunit_dollar'  ;
            'bequest_phi_1'     };
        
        % Specify default values for optional fields
        def_fields = struct(...
            'taxplan'       , 'base'    , ...
            'gcut'          , 0.0       , ...
            'legal_scale'   , 1.0       , ...
            'prem_legal'    , 1.0       , ...
            'amnesty'       , 0.0       , ...
            'deportation'   , 0.0       );
        
    end
    
    
    methods
        
        % Constructor
        function [this] = Scenario(props)
            
            if( ~isa( props, 'struct' ) )
                error( 'Scenario constructor expects structure of initial variables.' );
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
            if(  ~any(strcmp(props.economy, {'open', 'steady', 'closed'})) )
                error( '<economy> field must be either ''steady'', ''open'', or ''closed''' );
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
             
        end
        
        
        % Derived property: Get if isCurrentPolicy
        function [flag] = isCurrentPolicy(this)
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
        end
        
        
        % Generate analogous current policy scenario
        function [scenario] = currentPolicy(this)
            params          = this.getParams();
            for f = fieldnames(Scenario.def_fields)'
                params.(f{1}) = Scenario.def_fields.(f{1});
            end
            scenario        = Scenario(params);
        end
        
        
        % Generate analogous steady state scenario
        function [scenario] = steady(this)
            params          = this.getParams();
            params.economy  = 'steady';
            scenario        = Scenario(params);
        end
        
        
        % Generate analogous open economy scenario
        function [scenario] = open(this)
            params          = this.getParams();
            params.economy  = 'open';
            scenario        = Scenario(params);
        end
        
        
        % Generate analogous closed economy scenario
        function [scenario] = closed(this)
            params          = this.getParams();
            params.economy  = 'closed';
            scenario        = Scenario(params);
        end
        
        
        % Generate tags for baseline and counterfactual definitions
        function [basedef_tag, counterdef_tag] = generate_tags(this)

            basedef_tag = [...
                sprintf('%.3f'  , this.beta             ),...
                sprintf('_%.3f' , this.gamma            ), ...
                sprintf('_%.2f' , this.sigma            ), ...
                sprintf('_%e'   , this.modelunit_dollar ), ...
                sprintf('_%.3f' , this.bequest_phi_1    )];
            
            if( this.isCurrentPolicy )
                counterdef_tag = 'baseline';
            else
                counterdef_tag = [...
                    sprintf('%s'    , this.taxplan      ), ...
                    sprintf('_%+.2f', this.gcut         ), ...
                    sprintf('_%.1f' , this.legal_scale  ), ...
                    sprintf('_%.3f' , this.prem_legal   ), ...
                    sprintf('_%.2f' , this.amnesty      ), ...
                    sprintf('_%.2f' , this.deportation  )];
            end
        end


        function [params] = getParams(this)
            params = struct();
            for f = [ Scenario.req_fields ; fieldnames(Scenario.def_fields) ]'
                params.(f{1}) = this.(f{1});
            end
        end
        
    end
    
    
end

