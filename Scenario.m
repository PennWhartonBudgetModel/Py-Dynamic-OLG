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
            %   1. Make string of concatenated params
            %   2. Hash the string down to 120 chars
            tag = ''; 
            for o = Scenario.req_params'
                if( ~strcmp( o{1}, 'economy' ) )
                    tag = strcat( tag, '_', num2str(this.(o{1})) );
                end
            end
            this.basedeftag = Scenario.compactifyTag( tag );
            if this.isCurrentPolicy()
                this.counterdeftag = 'currentpolicy';
            else
                tag = '';
                for o = fieldnames(Scenario.def_params)'
                    tag = strcat( tag, '_', num2str(this.(o{1})) );
                end
                this.counterdeftag = Scenario.compactifyTag( tag );
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
        
        
        % Check if the Scenario has been solved and stored to files
        function [flag] = isSolved(this)
            work_dir = PathFinder.getWorkingDir(this);
            % TBD
            %  1. have paramTargets.mat hold solved info
            %  2. read it and see that it matches this scenario and is
            %  solved
            
        end % isSolved
        
        
        
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
        
        
        %%
        % Write out Scenario series for output in single CSV
        %     tag : name of file
        function [] = export(this, tag)
            
            % TBD: This should call the load() method
            
            % Identify scenario working directory
            workingdir = PathFinder.getWorkingDir(this);
            
            % Load dynamic and static variables
            Dynamic = load(fullfile(workingdir, 'dynamics.mat'));
            if this.isCurrentPolicy()
                Static = Dynamic;
            else
                Static = load(fullfile(workingdir, 'statics.mat'));
            end
            Market = load(fullfile(workingdir, 'market.mat'));
            
            % Specify data series years
            years  = (this.TransitionFirstYear : this.TransitionLastYear)';

            Dynamic.outvars = struct( ...
            'ssts'              , 'PayrollTax'              , ...
            'caprevs'           , 'CapitalTax'              , ...
            'cits_domestic'     , 'CIT_Domestic'            , ...
            'cits_foreign'      , 'CIT_Foreign'             , ...
            'caps_domestic'     , 'Capital_Domestic'        , ...
            'caps_foreign'      , 'Capital_Foreign'         , ...
            'debts_domestic'    , 'Debt_Domestic'           , ...
            'debts_foreign'     , 'Debt_Foreign'            , ...
            'outs'              , 'Output'                  , ...
            'bens'              , 'SocialSecurityBenefits'  , ...
            'caps'              , 'Capital'                 , ...
            'labeffs'           , 'EfficientLabor'          , ...
            'labs'              , 'Labor'                   , ...
            'labincs'           , 'LaborIncome'             , ...
            'capincs'           , 'CapitalIncome'           , ...
            'pops'              , 'PopulationHouseholds'      ...
            ); 
    
            Static.outvars = Dynamic.outvars;
        
            Market.outvars = struct( ...
            'caprates'          , 'CapitalReturn'           , ...
            'wages'             , 'WageLevel'                 ...
            ); 
        
            % Helper function to convert series
            function [f] = toOutputValue( outvar, outval )
                % TBD: Get these values from interface?
                %  Match population size in 2017 since pops=1 in 2017
                %  Conversion to dollar aggregates is since modelunit_dollar is
                %  currently targeted to GDP/HH
                HH_2017     = 126.22 * 1e06;
                DOLLAR      = (1/this.modelunit_dollar) * HH_2017;
                
                c = 1;
                switch outvar
                    case 'labeffs'  
                    case 'labs'     
                    case 'caprates'
                    case 'wages'
                        c = 1;
                    case 'pops'
                        c = HH_2017;
                    otherwise
                        c = DOLLAR;
                end
                f = outval * c; 
            end % toOutputValue
            
            fid = fopen(fullfile(PathFinder.getDataSeriesOutputDir(), strcat(tag, '.csv')), 'w');

            % Write header to file 
            fprintf(fid, 'Year' );
            for o = fieldnames( Dynamic.outvars )'
                fprintf( fid, ',%s', Dynamic.outvars.(o{1}) );
            end 
            for o = fieldnames( Market.outvars )'
                fprintf( fid, ',%s', Market.outvars.(o{1}) );
            end 
            for o = fieldnames( Static.outvars )'
                fprintf( fid, ',STATIC.%s', Static.outvars.(o{1}) );
            end                 
            fprintf(fid, '\n');
            
            % Write values, year is first
            for t = 1:length(years)-1
                fprintf(fid, '%u', years(t));
                for o = fieldnames( Dynamic.outvars )'
                    val = Dynamic.(o{1})(t);
                    fprintf( fid, ',%f', toOutputValue( o{1}, val ) );
                end 
                for o = fieldnames( Market.outvars )'
                    val = Market.(o{1})(t);
                    fprintf( fid, ',%f', toOutputValue( o{1}, val ) );
                end 
                for o = fieldnames( Static.outvars )'
                    val = Static.(o{1})(t);
                    fprintf( fid, ',%f', toOutputValue( o{1}, val ) );
                end 
                fprintf(fid, '\n');
            end
            
            fclose(fid);
            
        end % export
        
        
    end % instance methods, public
    
    
    
    methods (Static, Access = private)
        
        function [newtag] = compactifyTag( tag )
            
            TAG_SIZE    = 120;
            % Allow chars ASCII 65-90 and 97-122 only
            % Remap ones that fall between into numbers (0-6)
            MIN_CHAR1   = 65;   MAX_CHAR1   = 90;
            MIN_CHAR2   = 97;   MAX_CHAR2   = 122;
            CHAR0       = 48; 
            newtag      = tag( 1:min(length(tag), TAG_SIZE) );
            for i=TAG_SIZE:length(tag)
                d = mod(i, TAG_SIZE) + 1;
                c = mod(newtag(d) + tag(i), MAX_CHAR2 - MIN_CHAR1) + MIN_CHAR1;
                if ( (c > MAX_CHAR1) && (c < MIN_CHAR2) )
                    c = c - (MAX_CHAR1 - CHAR0);
                end
                newtag(d) = c;
            end
            
        end % compactifyTag
        
    end

    
end % Scenario


