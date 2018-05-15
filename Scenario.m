%%
% Scenario definition for dynamic model execution.
%
%%
classdef Scenario
    
    properties (GetAccess = public, SetAccess = immutable)
        
        % Identifier tags
        basedeftag;         % level 1 in dir structure
        counterdeftag;      % level 2 in dir structure
        economytag;         % level 3 in dir structure
        comparisontag;      % REM: This is for isEquivalent()
        
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
        
        % Immigration parameters
        legal_scale;
        prem_legal;
        amnesty;
        deportation;
        
        % Tax parameters
        TaxCode;

        % Government expenditures
        OutlaysPolicy;
        

        % Social Security parameters
        TaxRate;
        TaxMax;
        DonutHole;
        COLA;
        PIA;
        NRA;
        CreditEarningsAboveTaxMax;
        FirstYear;
        GradualChange;
        
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
            'legal_scale'               , 1.0               , ...
            'prem_legal'                , 1.0               , ...
            'amnesty'                   , 0.0               , ...
            'deportation'               , 0.0               , ...
            'TaxCode'                   , 'CurrentPolicy'   , ...
            'OutlaysPolicy'             , 'CurrentPolicy'   , ...
            'TaxRate'                   , 0                 , ...
            'TaxMax'                    , 0                 , ...
            'DonutHole'                 , 0                 , ...
            'COLA'                      , 0                 , ...
            'PIA'                       , 0                 , ...
            'NRA'                       , 0                 , ...
            'CreditEarningsAboveTaxMax' , 0                 , ...
            'FirstYear'                 , -1                , ...
            'GradualChange'             , -1                ...
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
                    warning('Scenario:constructor:noParameterMatch', 'Field <%s> does not match any Scenario fields.', o{1});
                end
                
            end
            
            % Validate economy
            assert(any(strcmp(this.economy, {'steady', 'open', 'closed'})), ...
                '<economy> parameter must be either ''steady'', ''open'', or ''closed''.');
            
            
            
            % Generate identifier tags for baseline and counterfactual definitions
            %   1. Make string of concatenated params
            %   2. Hash the string down to 120 chars
            % NOTE: comparisontag is built for isEquivalent 
            tag = '';
            for o = Scenario.req_params'
                if ~any(strcmp(o{1}, {'economy', 'TransitionLastYear'}))
                    tag = strcat(tag, '_', num2str(this.(o{1})));
                end
            end
            this.basedeftag = Scenario.compactifyTag(tag);
            
            if this.isCurrentPolicy()
                this.counterdeftag = 'currentpolicy';
            else
                tag = '';
                for o = fieldnames(Scenario.def_params)'
                    tag = strcat(tag, '_', num2str(this.(o{1})));
                end
                this.counterdeftag = Scenario.compactifyTag(tag);
            end
            
            tag = this.economy;
            if ~strcmp(this.economy, 'steady')
                tag = strcat(tag, num2str(this.TransitionLastYear));
            end
            this.economytag = tag;
            
            this.comparisontag = strcat(this.basedeftag, this.counterdeftag, this.economytag);
            
            
        end % Scenario constructor
        
        
        % Identify if scenario is equivalent to another scenario
        %   Parameter representations in tags determine precision for equivalency evaluation
        %   NOTE: For economy=steady, we ignore TransitionLastYear
        function [flag] = isEquivalent(this, scenario)
            flag = strcmp(this.comparisontag, scenario.comparisontag );
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
            flag = exist(fullfile(PathFinder.getWorkingDir(this), 'solved'), 'file');
        end % isSolved
        
        
        
        % Get all scenario parameters
        function [params] = getParams(this)
            params = struct();
            for f = [ Scenario.req_params ; fieldnames(Scenario.def_params) ]'
                params.(f{1}) = this.(f{1});
            end
        end
        
        
        % Human-readable description
        function [desc] = shortDescription(this)
            
            desc = '[';
            switch( this.economy) 
            case 'steady'
                desc = strcat(desc, 'Steady-state');
            case {'open', 'closed'}
                str1 = [upper(this.economy(1)), this.economy(2:end)];
                desc = strcat(desc, str1 );
            end
            if this.isCurrentPolicy()
                desc = strcat(desc, ' economy - baseline');
            else
                desc = strcat(desc, ' economy - counterfactual');
            end
            
            T_model = ParamGenerator.timing(this).T_model;
            desc = sprintf( '%s]', desc );
         	desc = sprintf( '%s\n \t%-25s= %u'   , desc, 'T_model'     , T_model               ); 
        	desc = sprintf( '%s\n \t%-25s= %u'   , desc, 'IsLowReturn' , this.IsLowReturn      ); 
        	desc = sprintf( '%s\n \t%-25s= %7.8f', desc, 'Beta'        , this.beta             ); 
        	desc = sprintf( '%s\n \t%-25s= %7.8f', desc, 'Gamma'       , this.gamma            ); 
        	desc = sprintf( '%s\n \t%-25s= %7.8f', desc, 'Sigma'       , this.sigma            );    
        	desc = sprintf( '%s\n \t%-25s= %e'   , desc, 'Model$'      , this.modelunit_dollar );    
 
        end %shortDescription
         
        
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
            
            % Load dynamic and static variables
            Dynamic = load(fullfile(PathFinder.getWorkingDir(this), 'dynamics.mat'));
            Market  = load(fullfile(PathFinder.getWorkingDir(this), 'market.mat'));
            if this.isCurrentPolicy()
                Static       = Dynamic;
                StaticMarket = Market;
            else
                Static       = load(fullfile(PathFinder.getWorkingDir(this), 'statics.mat'));
                StaticMarket = load(fullfile(PathFinder.getWorkingDir(this.currentPolicy()), 'market.mat' ));
            end
            
            % Specify data series years
            switch this.economy
                case 'steady'
                    years = [this.TransitionFirstYear - 1];
                otherwise
                    years  = (this.TransitionFirstYear : this.TransitionLastYear-1)';
            end

            Dynamic.outvars = struct( ...
            'ssts'              , 'PayrollTax'              , ...
            'corpTaxs'          , 'CapitalTax'              , ...
            'cits'              , 'HH_PreferredRatesTax'    , ...
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
            'investment'        , 'Investment'              , ...
            'pops'              , 'PopulationHouseholds'      ...
            ); 
            Market.outvars = struct( ...
            'MPKs'               , 'MPK'                     , ...
            'equityFundDividends', 'EquityFundDividends'     , ...
            'bondFundDividends'  , 'BondFundDividends'       , ...
            'wages'              , 'WageLevel'                 ...
            ); 
            Static.outvars       = Dynamic.outvars;
            StaticMarket.outvars = Market.outvars;
        
            % Prepare file to which to write
            if( ~exist( PathFinder.getSeriesOutputDir(), 'dir' ) )
                mkdir( PathFinder.getSeriesOutputDir() );
            end
            outputfilename  = fullfile(PathFinder.getSeriesOutputDir(), strcat(tag, '.csv'));
            fid             = fopen(outputfilename, 'w');
            if( fid == -1 )
                throw(MException('Scenario:export','Cannot open file (%s) for output.', outputfilename ));
            end 

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
            for o = fieldnames( StaticMarket.outvars )'
                fprintf( fid, ',STATIC.%s', StaticMarket.outvars.(o{1}) );
            end                 
            fprintf(fid, '\n');
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            % Helper function to convert series
            function [f] = toOutputValue( outvar, outval )
                % TBD: Get these values from interface?
                %  Match population size in 2017 since pops=1 in 2017
                %  Conversion to dollar aggregates is since modelunit_dollar is
                %  currently targeted to GDP/HH
                HH_2017     = 193.7 * 1e06; % rem: age 21+
                DOLLAR      = (1/this.modelunit_dollar) * HH_2017;
                
                c = 1;
                switch outvar
                    case 'labeffs'  
                    case 'labs'     
                    case 'MPKs'
                    case 'wages'
                    case 'equityFundDividends'
                    case 'bondFundDividends'
                        c = 1;
                    case 'pops'
                        c = HH_2017;
                    otherwise
                        c = DOLLAR;
                end
                f = outval * c; 
            end % toOutputValue
            
            %% Helper function to print values
            function [] = printValues( s )
                for o = fieldnames( s.outvars )'
                    val = s.(o{1})(t);
                    fprintf( fid, ',%f', toOutputValue( o{1}, val ) );
                end 
            end % printValues
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Write values, year is first
            for t = 1:length(years)
                fprintf(fid, '%u', years(t));
                printValues( Dynamic        );
                printValues( Market         );
                printValues( Static         );
                printValues( StaticMarket   );
                fprintf(fid, '\n');
            end
            
            fclose(fid);
            
        end % export

        %       Writes an already-solved scenario's optimal decision rules
        %       and productivity transitions to file in a long format.
        function [] = writeTransitionMatrix(scenario)
            
            working_dir = PathFinder.getWorkingDir(scenario);
            load(fullfile(working_dir, 'decisions.mat'), 'OPTs');
            
            % melt decision rules into long format
            % in steady state, no time dimension and there are 4 columns
            % on transition path, a time dimension, and 5 columns
            if scenario.isSteady()
                [productivity_index, savings_index, earnings_index, age_index] = ind2sub(size(OPTs.SAVINGS), (1:numel(OPTs.SAVINGS))');
                decision_index = [productivity_index, savings_index, earnings_index, age_index, zeros(size(age_index))];
            else
                [productivity_index, savings_index, earnings_index, age_index, time_index] = ind2sub(size(OPTs.SAVINGS), (1:numel(OPTs.SAVINGS))');
                decision_index = [productivity_index, savings_index, earnings_index, age_index, time_index];
            end

            decision_rules = [OPTs.CONSUMPTION(:), OPTs.SAVINGS(:), OPTs.LABOR(:), OPTs.AVG_EARNINGS(:), OPTs.TAXABLE_INC(:)];

            % melt productivity values
            productivity_values = ParamGenerator.grids(scenario).zs;
            [productivity_index, age_index] = ind2sub(size(productivity_values), (1:numel(productivity_values))');
            productivity_values = productivity_values(:);
            productivity_values_index = uint8([productivity_index, age_index]);
            
            % melt productivity transitions
            z_transitions = ParamGenerator.grids(scenario).transz;
            [productivity_index, productivity_prime_index, age_index] = ind2sub(size(z_transitions), (1:numel(z_transitions))');
            productivity_transitions = z_transitions(:);
            productivity_transitions_index = uint8([productivity_index, productivity_prime_index, age_index]);

            % write to file
            outputdir = PathFinder.getTransitionMatrixOutputDir();
            
            % create output folder if it does not exist
            if exist(outputdir, 'file') ~= 7
                mkdir(outputdir)
            end
            
            % check for whether scenario output subfolder exists
            % if it does, then this is a duplicate writing out
            if exist(fullfile(outputdir, scenario.basedeftag, scenario.counterdeftag), 'file') == 7
                return
            end
            
            % check if map file exists, create it if it does not
            if exist(fullfile(outputdir, 'map.csv'), 'file') ~= 2
                filehandle = fopen(fullfile(outputdir, 'map.csv'), 'w');
                fprintf(filehandle, strjoin(fieldnames(scenario), ','));
                fprintf(filehandle, '\n');
                fclose(filehandle);
            end

            % append scenario info to map file by writing out to text file
            % then loading text file back in
            values = struct2table(scenario.getParams());
            writetable(values, '.temp.txt', 'WriteVariableNames', false);
            text = fileread('.temp.txt');
            delete('.temp.txt');
            filehandle = fopen(fullfile(outputdir, 'map.csv'), 'a');
            fprintf(filehandle, [scenario.basedeftag, ',', scenario.counterdeftag, ',', text]);
            fclose(filehandle);
            
            % delete 
            % store all output in subfolder
            outputdir = fullfile(outputdir, scenario.basedeftag, scenario.counterdeftag);
            
            % load wages
            load(fullfile(working_dir, 'market.mat'), 'wages');
            
            % load distribution
            load(fullfile(working_dir, 'distribution.mat'), 'DIST');
            
            % melt distribution
            [productivity_index, savings_index, earnings_index, age_index, status_index, year_index] = ind2sub(size(DIST), (1:numel(DIST))');
            distribution_values = DIST(:);
            distribution_index = uint8([productivity_index, savings_index, earnings_index, age_index, status_index, year_index]);

            % create a folder to store output
            mkdir(outputdir)

            h5create(fullfile(outputdir, 'data.hdf5'), '/decision_rules', size(decision_rules), 'ChunkSize', size(decision_rules), 'Deflate', 9);
            h5write(fullfile(outputdir, 'data.hdf5'), '/decision_rules', decision_rules);

            h5create(fullfile(outputdir, 'data.hdf5'), '/decision_index', size(decision_index), 'ChunkSize', size(decision_index), 'Deflate', 9);
            h5write(fullfile(outputdir, 'data.hdf5'), '/decision_index', decision_index);
            
            h5create(fullfile(outputdir, 'data.hdf5'), '/productivity_values', size(productivity_values), 'ChunkSize', size(productivity_values), 'Deflate', 9);
            h5write(fullfile(outputdir, 'data.hdf5'), '/productivity_values', productivity_values);
            
            h5create(fullfile(outputdir, 'data.hdf5'), '/productivity_values_index', size(productivity_values_index), 'ChunkSize', size(productivity_values_index), 'Deflate', 9);
            h5write(fullfile(outputdir, 'data.hdf5'), '/productivity_values_index', productivity_values_index);
            
            h5create(fullfile(outputdir, 'data.hdf5'), '/earnings_values', size(ParamGenerator.grids(scenario).bv), 'ChunkSize', size(ParamGenerator.grids(scenario).bv), 'Deflate', 9);
            h5write(fullfile(outputdir, 'data.hdf5'), '/earnings_values', ParamGenerator.grids(scenario).bv);
            
            h5create(fullfile(outputdir, 'data.hdf5'), '/savings_values', size(ParamGenerator.grids(scenario).kv), 'ChunkSize', size(ParamGenerator.grids(scenario).kv), 'Deflate', 9);
            h5write(fullfile(outputdir, 'data.hdf5'), '/savings_values', ParamGenerator.grids(scenario).kv);
            
            h5create(fullfile(outputdir, 'data.hdf5'), '/productivity_transitions', size(productivity_transitions), 'ChunkSize', size(productivity_transitions), 'Deflate', 9);
            h5write(fullfile(outputdir, 'data.hdf5'), '/productivity_transitions', productivity_transitions);

            h5create(fullfile(outputdir, 'data.hdf5'), '/productivity_transitions_index', size(productivity_transitions_index), 'ChunkSize', size(productivity_transitions_index), 'Deflate', 9);
            h5write(fullfile(outputdir, 'data.hdf5'), '/productivity_transitions_index', productivity_transitions_index);

            h5create(fullfile(outputdir, 'data.hdf5'), '/distribution_values', size(distribution_values), 'ChunkSize', size(distribution_values), 'Deflate', 9);
            h5write(fullfile(outputdir, 'data.hdf5'), '/distribution_values', distribution_values);

            h5create(fullfile(outputdir, 'data.hdf5'), '/distribution_index', size(distribution_index), 'ChunkSize', size(distribution_index), 'Deflate', 9);
            h5write(fullfile(outputdir, 'data.hdf5'), '/distribution_index', distribution_index);
            
            h5create(fullfile(outputdir, 'data.hdf5'), '/wages', size(wages), 'ChunkSize', size(wages), 'Deflate', 9);
            h5write(fullfile(outputdir, 'data.hdf5'), '/wages', wages);

            
        end % writeTransitionMatrix
        
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


