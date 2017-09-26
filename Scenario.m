%%
%   Scenarios are runs of the model
%       We currently define "current policy" by defaults
%       in the constructor. In the future, we could pass in the "current
%       policy".
%
classdef Scenario
    
    properties (SetAccess = protected )
        ID;
        batchID;
        useDynamicBaseline;
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
        def_fields = struct(    'ID'                , []        ...
                            ,   'batchID'           , []        ...
                            ,   'useDynamicBaseline', []        ...
                            ,   'taxplan'           , 'base'    ...
                            ,   'gcut'              , 0.0       ...
                            ,   'legal_scale'       , 1.0       ...
                            ,   'prem_legal'        , 1.0       ...
                            ,   'amnesty'           , 0.0       ...
                            ,   'deportation'       , 0.0       ...
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
            obj = Scenario(this.getParams);
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
        
    end % instance methods
    
    methods (Static)
        
        %% 
        %  Read a batch of Scenarios from the DB
        %  
        function [scenarios] = fetch_batch( batchID )
            if( (nargin < 1) || (isempty(batchID)) )
                error( '<batchID> is required.' );
            end
            
            % Add JDBC driver to Matlab Java path
            javaaddpath(fullfile(Environment.source(), 'jar', 'sqljdbc41.jar'));
            
            % Establish connection with development database
            %   TODO: Put these params into Environment
            connection = database('second_chart', 'pwbm', 'HbXk86rabjehD2AN', ...
                                  'Vendor', 'Microsoft SQL Server', 'AuthType', 'Server', ...
                                  'Server', 'ppi-slcsql.wharton.upenn.edu', 'PortNumber', 49170);
            
            % Get batch scenarios from Scenario table
            sql         = sprintf( 'EXEC p_ScenarioBatch %u, ''D'' ', batchID );
            o           = connection.exec(sql);
            r           = o.fetch();
            dataset     = cell2struct( r.Data, o.columnnames(true), 2);
            o.close();
            connection.close();
            
            % Preload calibration matrix for inversion
            [~, f_invert] = modelCalibrator.invert();
            
            scenarios = [];
            for i = 1:size(dataset)
                
                % Invert from elasticities to beta,gamma,sigma, etc.
                savelas = dataset(i).SavingsElasticity;
                labelas = dataset(i).LaborElasticity;
                params  = f_invert( struct( 'savelas', savelas, 'labelas', labelas ) );
                
                % TODO: IMPORTANT!
                % bequest_phi_1 is missing.
                
                canAdd   = true;
                id       = dataset(i).ID;
                openecon = dataset(i).OpenEconomy;
                if( openecon == 1 )
                    economy = 'open';
                elseif( openecon == 0 ) 
                    economy = 'closed';
                else  % Cannot do anything with the convex combo Scenarios
                    canAdd = false;
                end
                            
                if( canAdd )
                    t = struct( 'ID'                , id                            ...
                            ,   'batchID'           , batchID                       ...
                            ,   'useDynamicBaseline', dataset(i).UseDynamicBaseline ...
                            ,   'economy'           , economy                       ...
                            ,   'beta'              , params.beta                   ...
                            ,   'gamma'             , params.gamma                  ...
                            ,   'sigma'             , params.sigma                  ...
                            ,   'modelunit_dollar'  , params.modelunit_dollar       ... 
                            ,   'bequest_phi_1'     , 0                             ... % TEMPORARY
                            ,   'gcut'              , -dataset(i).ExpenditureShift  ... % REM: Inconsistent definition
                            ,   'taxplan'           , dataset(i).TaxPlan            ...
                            );
                    scenarios = [scenarios, Scenario(t)];
                else
                    fprintf( 'Skipping ScenarioID=%u\n', id );
                end
            end
        end % fetch_batch
        
    end % static methods
    
end % Scenario

