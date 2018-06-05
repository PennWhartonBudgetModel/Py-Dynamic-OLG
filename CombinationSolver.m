%%
% Dynamic model combination solver and dynamic series generator.
%
%%
classdef CombinationSolver


properties (Constant)
    
    % Define scenario directory path
    scenariodir = fullfile(PathFinder.getSourceDir(), 'Scenarios');
    
    % Define current policy and counterfactual scenario file path generators
    currentpolicyfile  = @(i) fullfile(CombinationSolver.scenariodir, sprintf('currentpolicy%05d.mat' , i));
    counterfactualfile = @(i) fullfile(CombinationSolver.scenariodir, sprintf('counterfactual%05d.mat', i));
    
end



methods (Static)
    
    
    
    % Solve all combinations of inputs and outputs
    function [] = solveAll(output_parameters, input_parameters)
        
        if ~exist('output_parameters', 'var'), output_parameters = []; end
        if ~exist('input_parameters' , 'var'), input_parameters  = []; end
        
        [currentpolicys, counterfactuals] = CombinationSolver.generateScenarios(output_parameters, input_parameters);
        
        for i = 1:length(currentpolicys ), CombinationSolver.solveCurrentPolicy (i); end
        for i = 1:length(counterfactuals), CombinationSolver.solveCounterfactual(i); end
        
        CombinationSolver.checkTerminations();
        
        CombinationSolver.generateSeries();
        
    end
    
    
    
    
    % Convert an output scenario, represented as a structure, into an executable scenario
    function [scenario] = convertOutputScenario(s)
        
        % Identify economy openness, specifying closed economy for partially open economies to generate both open and closed economies
        if s.OpenEconomy == 1
            economy = 'open';
        else
            economy = 'closed';
        end
        
        % Invert elasticities
        inverse = ParamGenerator.invert(struct( ...
            'is_low_return'  , s.IsLowReturn     , ...
            'labelas'        , s.LaborElasticity ...
        ));
        
        % Initialize scenario parameter structure with required parameters
        parameters = struct( ...
            'economy'               , economy                   , ...
            'beta'                  , inverse.beta              , ...
            'gamma'                 , inverse.gamma             , ...
            'sigma'                 , inverse.sigma             , ...
            'bequest_phi_1'         , 0                         , ...
            'modelunit_dollar'      , inverse.modelunit_dollar  , ...
            'IsLowReturn'           , s.IsLowReturn             , ...
            'TransitionFirstYear'   , 2018                      , ...
            'TransitionLastYear'    , 2018 + 2                 , ...
            'ClosureYear'           , 2018 + 2                 ...
        );
        
        % Add optional parameters as available
        for field_ = fieldnames(s)', field = field_{1};
            if ~isfield(parameters, field)
                parameters.(field) = s.(field);
            end
        end
        
        % Construct scenario from parameter structure
        warningid = 'Scenario:constructor:noParameterMatch';
        warning('off', warningid), scenario = Scenario(parameters); warning('on', warningid)
        
    end
    
    
    
    
    % Generate output scenarios and executable scenarios
    % 
    %   output_parameters specify scenario combinatorics to be applied by DynamicModel
    % 
    %     output_parameters.OpenEconomy     = { 0.0, 0.5, 1.0 };
    %     output_parameters.LaborElasticity = { 0.25, 0.75 };
    % 
    %   input_parameters specify filtering on scenario combinatorics provided by input interfaces
    % 
    %     input_parameters.TaxCode = { 'CurrentPolicy' };
    %     input_parameters.TaxRate = { 0 };
    % 
    function [currentpolicys, counterfactuals] = generateScenarios(output_parameters, input_parameters)
        
        
        % Generate input sets using input interface map files, applying parameter value filters if specified
        taxcalculator_s  = generateInputSet(fullfile(PathFinder.getTaxCalculatorInputDir() , 'Map.csv'));
        oasicalculator_s = generateInputSet(fullfile(PathFinder.getOASIcalculatorInputDir(), 'map.csv'));
        
        function [input_s] = generateInputSet(mapfile)
            
            input_s  = table2struct(readtable(mapfile, 'ReadRowNames', true))';
            
            if exist('input_filters', 'var') && ~isempty(input_parameters)
                for o_ = fieldnames(input_parameters)'
                    if ischar(input_parameters.(o_{1}){1})
                        f = @strcmp; g = @cell;
                    else
                        f = @eq; g = @cell2mat;
                    end
                    if isfield(input_s, o_{1})
                        input_s = input_s(arrayfun( ...
                            @(s) any( f(s.(o_{1}), g(input_parameters.(o_{1}))) ), input_s ...
                        ));
                    end
                end
            end
            
        end
        
        
        
        % Define default output parameter value sets
        output_parameters0.OpenEconomy          = {0, 1};
        output_parameters0.UseDynamicBaseline   = {true, false};
        output_parameters0.LaborElasticity      = {0.5};
        output_parameters0.IsLowReturn          = {true};
        
        % Identify output parameter value sets, applying defaults where unspecified
        if exist('output_parameters', 'var') && ~isempty(output_parameters)
            for o = fieldnames(output_parameters0)'
                if ~isfield(output_parameters, o{1})
                    output_parameters.(o{1}) = output_parameters0.(o{1});
                end
            end
        else
            output_parameters = output_parameters0;
        end
        
        
        
        % Generate all combinations of input sets and output sets to determine full set of output scenarios
        s_sets = [ ...
            cellfun(@(o) struct(o, output_parameters.(o)), fieldnames(output_parameters)', 'UniformOutput', false), ...
            {taxcalculator_s } , ...
            {oasicalculator_s} ...
        ];
        
        output_scenarios = struct();
        for s_set_ = s_sets, s_set = s_set_{1};
            new_output_scenarios = [];
            for option = s_set
                option_scenarios = output_scenarios;
                for o = fieldnames(option)'
                    [option_scenarios.(o{1})] = deal(option.(o{1}));
                end
                new_output_scenarios = [new_output_scenarios, option_scenarios]; %#ok<AGROW>
            end
            output_scenarios = new_output_scenarios;
        end
        
        if isempty(output_scenarios), warning('No output scenarios. Scenario set generation terminated.'), return, end
        
        
        
        % Generate map file for output scenarios, removing duplicates
        map = unique(struct2table(output_scenarios), 'stable');
        output_scenarios = table2struct(map)';
        
        map.Properties.DimensionNames = {'id', 'ScenarioParameters'};
        map.Properties.RowNames = arrayfun(@(i) sprintf('%u', i), 1:height(map), 'UniformOutput', false);
        
        outputdir = PathFinder.getSeriesOutputDir();
        if exist(outputdir, 'dir'), rmdir(outputdir, 's'), end, mkdir(outputdir)
        
        writetable(map, fullfile(outputdir, 'map.csv'), 'WriteRowNames', true);
        
        
        
        % Define function to remove empty entries from a cell array
        compress = @(c) c(~cellfun(@isempty, c));
                
        % Convert output scenarios to executable scenarios
        scenarios = arrayfun(@CombinationSolver.convertOutputScenario, output_scenarios, 'UniformOutput', false);
        
        % Remove duplicate scenarios
        for i = 1:length(scenarios)
            for j = i+1:length(scenarios)
                if scenarios{j}.isEquivalent(scenarios{i})
                    scenarios{i} = [];
                    break;
                end
            end
        end
        scenarios = compress(scenarios);
        
        % Remove open economy scenarios with a corresponding closed economy scenario
        for i = 1:length(scenarios)
            if scenarios{i}.isOpen()
                scenario_closed = scenarios{i}.closed();
                for j = 1:length(scenarios)
                    if ( i ~= j && ~isempty(scenarios{j}) && scenarios{j}.isEquivalent(scenario_closed) )
                        scenarios{i} = [];
                        break;
                    end
                end
            end
        end
        scenarios = compress(scenarios);
        
        % Separate scenarios into current policy scenarios and counterfactual scenarios
        currentpolicys  = cell(size(scenarios));
        counterfactuals = cell(size(scenarios));
        
        for i = 1:length(scenarios)
            if scenarios{i}.isCurrentPolicy()
                currentpolicys{i} = scenarios{i};
            else
                currentpolicys{i} = scenarios{i}.currentPolicy();
                counterfactuals{i} = scenarios{i};
            end
            for j = 1:i-1
                if (~isempty(currentpolicys{j}) && currentpolicys{i}.isEquivalent(currentpolicys{j}))
                    currentpolicys{i} = [];
                    break;
                end
            end
        end
        
        currentpolicys  = compress(currentpolicys );
        counterfactuals = compress(counterfactuals);
        
        % Clear or create scenario directory
        if exist(CombinationSolver.scenariodir, 'dir'), rmdir(CombinationSolver.scenariodir, 's'), end, mkdir(CombinationSolver.scenariodir)
        
        % Save current policy scenarios
        for i = 1:length(currentpolicys)
            scenario = currentpolicys{i}; %#ok<NASGU>
            save(CombinationSolver.currentpolicyfile(i), 'scenario');
        end
        
        % Save counterfactual scenarios
        for i = 1:length(counterfactuals)
            scenario = counterfactuals{i}; %#ok<NASGU>
            save(CombinationSolver.counterfactualfile(i), 'scenario');
        end
        
    end
    
    
    
    
    % Solve dynamic model for a scenario stored in a scenario file
    function [] = solve(scenariofile)
        
        % Load scenario
        s = load(scenariofile);
        scenario = s.scenario;
        
        % Solve dynamic model for scenario
        workingdir = ModelSolver.solve(scenario);
        
        % Extract and save solver termination condition
        iterations = csvread(fullfile(workingdir, 'iterations.csv'));
        
        termination.iter = iterations(end, 1);
        termination.eps  = iterations(end, 2);
        
        termination = termination; %#ok<ASGSL,NASGU>
        save(scenariofile, '-append', 'termination');
        
    end
    
    function [] = solveCurrentPolicy(i)
        CombinationSolver.solve(CombinationSolver.currentpolicyfile(i));
    end
    
    function [] = solveCounterfactual(i)
        CombinationSolver.solve(CombinationSolver.counterfactualfile(i));
    end
    
    
    
    
    % Check solver termination conditions for scenarios in scenario directory
    function [] = checkTerminations()
        
        % Get all scenario files
        scenariofiles = arrayfun(@(s) fullfile(s.folder, s.name), ...
            [ dir(fullfile(CombinationSolver.scenariodir, 'currentpolicy*.mat' )) ;
              dir(fullfile(CombinationSolver.scenariodir, 'counterfactual*.mat')) ], ...
            'UniformOutput', false ...
        );
        nscenario = length(scenariofiles);
        
        % Initialize cell array of termination conditions
        terminations = cell(0,6);
        
        for i = 1:nscenario
            
            % Load scenario along with solver termination condition
            s = load(scenariofiles{i});
            
            % Determine scenario index as listed in scenario directory
            [~, index] = fileparts(scenariofiles{i});
            
            % Store termination condition, setting placeholder values if missing
            if ~isfield(s, 'termination'), s.termination = struct('iter', Inf, 'eps', Inf); end
            terminations = [terminations; ...
                { index                     , ...
                  s.scenario.basedeftag     , ...
                  s.scenario.counterdeftag  , ...
                  s.scenario.economy        , ...
                  s.termination.iter        , ...
                  s.termination.eps         } ...
            ]; %#ok<AGROW>
            
        end
        
        % Sort termination conditions by increasing iterations and error terms
        [~, sortinds] = sortrows(cell2mat(terminations(:,5:6)));
        
        % Save termination conditions to csv file
        fid = fopen(fullfile(CombinationSolver.scenariodir, 'terminations.csv'), 'w');
        fprintf(fid, 'ScenarioIndex,BaselineDefinition,CounterfactualDefinition,Economy,TerminationIteration,TerminationErrorTerm\n');
        for i = 1:nscenario, fprintf(fid, '%s,%s,%s,%s,%d,%0.4f\n', terminations{sortinds(i),:}); end
        fclose(fid);
        
    end
    
    
    
    
    % Generate dynamic series
    function [] = generateSeries()
        
        % Identify data series output directory
        outputdir = PathFinder.getSeriesOutputDir();
        
        % Generate interface dependencies file
        fid = fopen(fullfile(outputdir, 'dependencies.csv'), 'w');
        fprintf(fid, 'Component,Interface,Version\n');
        for r = PathFinder.getInputSet()
            fprintf(fid, '%s,%s,%s\n', r{1}{1}, r{1}{2}, r{1}{3});
        end
        fclose(fid);
        
        
        
        % Get output scenarios from map file
        output_scenarios = table2struct(readtable(fullfile(outputdir, 'map.csv')))';
        
        
        
        % Define mapping from series names to dynamic model variable names and static series sources
        series_names = struct( ...
            'GDP'                                   , struct('var_name', 'outs'         , 'static_source', 'projections')   ...
        ,   'GNP'                                   , struct('var_name', 'GNP'          , 'static_source', 'projections') ...
        ,   'CapitalServices'                       , struct('var_name', 'caps'         , 'static_source', 'projections')   ...
        ,   'LaborInput'                            , struct('var_name', 'labeffs'      , 'static_source', 'projections')   ...
        ,   'WagesAndSalaries'                      , struct('var_name', 'labincs'      , 'static_source', 'projections')   ...
        ,   'CompensationOfEmployees'               , struct('var_name', 'labincs'      , 'static_source', 'projections')   ...
        ,   'Employment'                            , struct('var_name', 'labs'         , 'static_source', 'projections')   ...
        ,   'GrossPrivateDomesticInvestment'        , struct('var_name', 'investment'   , 'static_source', 'projections')   ...
        ,   'ConsumptionOfFixedCapital'             , struct('var_name', 'caps'         , 'static_source', 'projections')   ...
        ... % BUDGET
        ,   'OutlaysDiscretionary'                  , struct('var_name', 'nonindexed'   , 'static_source', 'projections')    ...
        ,   'OutlaysMedicare'                       , struct('var_name', 'nonindexed'   , 'static_source', 'projections')    ...
        ,   'OutlaysMedicaid'                       , struct('var_name', 'nonindexed'   , 'static_source', 'projections')    ...
        ,   'OutlaysFederalRetirement'              , struct('var_name', 'nonindexed'   , 'static_source', 'projections')    ...
        ,   'OutlaysVeteransPrograms'               , struct('var_name', 'nonindexed'   , 'static_source', 'projections')    ...
        ,   'OutlaysOtherPrograms'                  , struct('var_name', 'nonindexed'   , 'static_source', 'projections')    ...
        ,   'OutlaysOffsettingReceipts'             , struct('var_name', 'nonindexed'   , 'static_source', 'projections')    ...
        ,   'OutlaysIncomeSecurity'                 , struct('var_name', 'nonindexed'   , 'static_source', 'taxcalculator')  ...
        ,   'OutlaysSocialSecurity'                 , struct('var_name', 'bens'         , 'static_source', 'oasicalculator') ...
        ,   'RevenuesPayrollTaxSocialSecurity'      , struct('var_name', 'ssts'         , 'static_source', 'oasicalculator') ...
        ,   'RevenuesPayrollTaxExSocialSecurity'    , struct('var_name', 'labincs'      , 'static_source', 'taxcalculator')  ...
        ,   'RevenuesIndividualIncomeTax'           , struct('var_name', 'pits'         , 'static_source', 'taxcalculator')  ...
        ,   'RevenuesCorporateIncomeTax'            , struct('var_name', 'corpTaxs'     , 'static_source', 'taxcalculator')  ...
        ,   'RevenuesEstateAndGiftTaxes'            , struct('var_name', 'caps'         , 'static_source', 'taxcalculator')  ...
        ,   'RevenuesExciseTaxes'                   , struct('var_name', 'outs'         , 'static_source', 'taxcalculator')  ...
        ,   'RevenuesCustomsDuties'                 , struct('var_name', 'GNP'          , 'static_source', 'taxcalculator')  ...
        ,   'RevenuesMiscellaneousReceipts'         , struct('var_name', 'nonindexed'   , 'static_source', 'taxcalculator')  ...
        ,   'SS_Cost'                               , struct('var_name', 'bens'         , 'static_source', 'oasicalculator') ...
        ,   'SS_NonInterestIncome'                  , struct('var_name', 'ssts'         , 'static_source', 'oasicalculator') ...
        );
        
        % Define function to construct dynamic scaling series for a specific scenario
        function scaling_series = constructScalingSeries(scenario, useDynamicBaseline)
            
            % Identify scenario working directory
            workingdir = PathFinder.getWorkingDir(scenario);
            
            % Load dynamic and static variables
            Dynamic = load(fullfile(workingdir, 'dynamics.mat'));
            if scenario.isCurrentPolicy()
                Static = Dynamic;
            else
                Static = load(fullfile(workingdir, 'statics.mat'));
            end
            
            % Load additional variables if performing dynamic baseline scaling
            if useDynamicBaseline
                Dynamic_currentPolicyOpen   = load(fullfile(PathFinder.getWorkingDir(scenario.currentPolicy().open()  ), 'dynamics.mat'));
                Dynamic_currentPolicyClosed = load(fullfile(PathFinder.getWorkingDir(scenario.currentPolicy().closed()), 'dynamics.mat'));
            else
                Dynamic_currentPolicyOpen   = {};
                Dynamic_currentPolicyClosed = {};
            end
            
            % Iterate over series names
            for o_ = fieldnames(series_names)'
                
                series_name = o_{1};
                var_name = series_names.(series_name).var_name;
                
                % Extract dynamic and static variable values and handle special series
                v_Dynamic = fetch_series( Dynamic, var_name );
                v_Static  = fetch_series( Static,  var_name );
                
                % Scale static variable values if performing dynamic baseline scaling
                if useDynamicBaseline
                    delta    = fetch_series( Dynamic_currentPolicyOpen, var_name ) ...
                                ./ fetch_series( Dynamic_currentPolicyClosed, var_name );
                    v_Static = v_Static .* delta;
                    v_Static(isnan(v_Static)) = 0;
                end
                
                % Calculate and store dynamic scaling series
                v_scale = v_Dynamic ./ v_Static;
                v_scale(isnan(v_scale)) = 1;
                scaling_series.(series_name) = v_scale';
                
            end
            
        end
        
        % Iterate over output scenarios
        n_scenarios = length(output_scenarios);
        n_failed = 0;
        for output_scenario = output_scenarios
            
            % Identify corresponding executable scenario
            scenario = CombinationSolver.convertOutputScenario(output_scenario);
            
            try
                
                % Construct dynamic scaling series, using convex combinations of open and closed economy series for partially open economy scenarios
                if (output_scenario.OpenEconomy == 0 || output_scenario.OpenEconomy == 1)
                    scaling_series = constructScalingSeries(scenario, output_scenario.UseDynamicBaseline);
                else
                    scaling_series_open   = constructScalingSeries(scenario.open()  , output_scenario.UseDynamicBaseline);
                    scaling_series_closed = constructScalingSeries(scenario.closed(), output_scenario.UseDynamicBaseline);
                    for o = fieldnames(series_names)'
                        scaling_series.(o{1}) = scaling_series_open  .(o{1})*(0 + output_scenario.OpenEconomy) ...
                                              + scaling_series_closed.(o{1})*(1 - output_scenario.OpenEconomy);
                    end
                end
                
                
                % Read static series from input interfaces
                c = fieldnames(series_names);
                first_year = scenario.TransitionFirstYear;
                last_year  = first_year + length(scaling_series.(c{1})) - 1;
                
                projections_file = fullfile(PathFinder.getProjectionsInputDir(), 'Projections.csv');
                static_series.projections = InputReader.read_series(projections_file, 'Year', first_year, last_year);
                
                taxcalculator_id = InputReader.find_input_scenario_id(fullfile(PathFinder.getTaxCalculatorInputDir(), 'Map.csv'), scenario);
                taxcalculator_file = fullfile(PathFinder.getTaxCalculatorInputDir(), strcat('Aggregates_', taxcalculator_id, '.csv'));
                static_series.taxcalculator = InputReader.read_series(taxcalculator_file, 'Year', first_year, last_year);
                
                oasicalculator_id = InputReader.find_input_scenario_id(fullfile(PathFinder.getOASIcalculatorInputDir(), 'map.csv'), scenario);
                oasicalculator_file = fullfile(PathFinder.getOASIcalculatorInputDir(), strcat('aggregates_', oasicalculator_id, '.csv'));
                static_series.oasicalculator = InputReader.read_series(oasicalculator_file, 'year', first_year, last_year);
                
                
                % Scale static series with scaling series to generate dynamic series
                for o = fieldnames(series_names)'
                    dynamic_series.(o{1}) = static_series.(series_names.(o{1}).static_source).(o{1}) .* scaling_series.(o{1});
                end
                
                
                % Write series to file
                series_table = struct2table(dynamic_series);
                
                series_table.Properties.DimensionNames = {'Year', 'Series'};
                series_table.Properties.RowNames = arrayfun(@(t) sprintf('%u', t), first_year:last_year, 'UniformOutput', false);
                
                writetable(series_table, fullfile(outputdir, sprintf('series_%u.csv', output_scenario.id)), 'WriteRowNames', true);
                
            catch e
                
                fprintf( 'Error with output: %s \n', e.message );
                n_failed = n_failed + 1;
                
            end
            
        end
        
        % Report on output generation success
        if ~n_failed
            fprintf('Successfully completed data series output generation for all %d scenarios.\n', n_scenarios);
        else
            fprintf('Failed to complete data series output generation for %d of %d scenarios.\n', n_failed, n_scenarios);
        end
        
    end
    
    
end

end



%%%  Helper function to create 'fake' series
%%   REM: Aggregate must have 'outs' series
function [series] = fetch_series( Aggregate, var_name )
    T_model = length( Aggregate.outs ); 
    switch var_name
        case 'nonindexed'
            series = [ones(1,10), Aggregate.outs(11:end)];
            series = series(1:T_model);     % truncate in case too long
        otherwise
            series = Aggregate.(var_name);
    end
end % fetch_series



