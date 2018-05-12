%%
% Dynamic model combination solver and data series output generator.
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
        
        [currentpolicys, counterfactuals] = CombinationSolver.generateScenarioSet(output_parameters, input_parameters);
        
        for i = 1:length(currentpolicys ), CombinationSolver.solveCurrentPolicy (i); end
        for i = 1:length(counterfactuals), CombinationSolver.solveCounterfactual(i); end
        
        CombinationSolver.checkTerminations();
        
        CombinationSolver.generateDataSeries();
        
    end
    
    
    
    
    % Convert an output scenario, represented as a structure, into a dynamic model scenario
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
            'TransitionLastYear'    , 2018 + 25                 ...
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
    
    
    
    
    % Generate set of executable scenarios
    % 
    %   output_parameters.OpenEconomy     = { 0.0, 0.5, 1.0 };
    %   output_parameters.LaborElasticity = { 0.25, 0.75 };
    % 
    %   input_parameters.TaxCode = { 'CurrentPolicy' };
    %   input_parameters.TaxRate = { 0 };
    % 
    function [currentpolicys, counterfactuals] = generateScenarioSet(output_parameters, input_parameters)
        
        
        % Generate input sets using input interface map files, applying parameter value filters if specified
        taxcalculator_s  = generateInputSet(fullfile(PathFinder.getTaxCalculatorInputDir() , 'Map.csv'), 'ID'               , 'taxcalculator' );
        oasicalculator_s = generateInputSet(fullfile(PathFinder.getOASIcalculatorInputDir(), 'map.csv'), 'id_OASIcalculator', 'oasicalculator');
        
        function [input_s] = generateInputSet(mapfile, idcolumn, tag)
            
            map_ = readtable(mapfile);
            map_.Properties.VariableNames{idcolumn} = ['id_', tag];
            input_s  = table2struct(map_)';

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
        
        outputdir = PathFinder.getDataSeriesOutputDir();
        if exist(outputdir, 'dir'), rmdir(outputdir, 's'), end, mkdir(outputdir)
        
        writetable(map, fullfile(outputdir, 'map.csv'), 'WriteRowNames', true);
        
        
        
        % Define function to remove empty entries from a cell array
        compress = @(c) c(~cellfun(@isempty, c));
                
        % Convert output scenarios to dynamic model scenarios
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
    
    
    
    
    % Generate data series output
    function [] = generateDataSeries()
        
        % Identify data series output directory
        outputdir = PathFinder.getDataSeriesOutputDir();
        
        % Generate interface dependencies file
        fid = fopen(fullfile(outputdir, 'dependencies.csv'), 'w');
        fprintf(fid, 'Component,Interface,Version\n');
        for r = PathFinder.getInputSet()
            fprintf(fid, '%s,%s,%s\n', r{1}{1}, r{1}{2}, r{1}{3});
        end
        fclose(fid);
        
        
        
        % Get output scenarios from map file
        output_scenarios = table2struct(readtable(fullfile(outputdir, 'map.csv')))';
        
        
        
        % Define mapping from dynamic model variable names to data series IDs
        dataseries_ids = struct( ...
            'labpits'           , [102]         , ...
            'ssts'              , [103]         , ...
            'caprevs'           , [110]         , ...
            'outs'              , [6, 199]      , ...
            'bens'              , [202]         , ...
            'caps'              , [1]           , ...
            'labeffs'           , [24]          , ...
            'lfprs'             , [25]          , ...
            'labincs'           , [28]          , ...
            'capincs'           , [29]          , ...
            'nonindexed'        , [299]         ...
        ); %#ok<NBRAK>
        
        
        
        % Specify standard data series years
        first_year  = 2017;
        last_year   = 2090;
        years       = (first_year : last_year)';
        
        % Define function to construct data series aligned to standard years for a specific scenario
        function dataseries = constructDataSeries(scenario, useDynamicBaseline)
            
            % Identify scenario working directory
            workingdir = PathFinder.getWorkingDir(scenario);
            
            % Load dynamic and static variables
            Dynamic = load(fullfile(workingdir, 'dynamics.mat'));
            if scenario.isCurrentPolicy()
                Static = Dynamic;
            else
                Static = load(fullfile(workingdir, 'statics.mat'));
            end
            
            % Determine timing parameters
            timing      = ParamGenerator.timing(scenario);
            
            % Determine number of years to pre-pad variable values
            nprepad     = timing.TransitionFirstYear - first_year;
            
            % Determine number of variable values to be trimmed or post-padded
            T_model     = timing.T_model;
            nextra      = nprepad + T_model - length(years);
            ntrim       = +max(nextra, 0);
            npostpad    = -min(nextra, 0);
            
            % Load additional variables if performing dynamic baseline scaling
            if useDynamicBaseline
                Dynamic_currentPolicyOpen   = load(fullfile(PathFinder.getWorkingDir(scenario.currentPolicy().open()  ), 'dynamics.mat'));
                Dynamic_currentPolicyClosed = load(fullfile(PathFinder.getWorkingDir(scenario.currentPolicy().closed()), 'dynamics.mat'));
            else
                Dynamic_currentPolicyOpen   = {};
                Dynamic_currentPolicyClosed = {};
            end
            
            % Construct data series
            for o_ = fieldnames(dataseries_ids)'
                
                % Extract dynamic and static variable values and handle special series
                switch o_{1}
                    case 'nonindexed'
                        v_Dynamic = [ones(1,10), Dynamic.outs(11:T_model)]; 
                        v_Static  = [ones(1,10), Dynamic.outs(11:T_model)]; 
                    otherwise
                        v_Dynamic = Dynamic.(o_{1});
                        v_Static  = Static. (o_{1});
                end
                
                % Scale static variable values if performing dynamic baseline scaling
                if useDynamicBaseline
                    v_Static = v_Static .* Dynamic_currentPolicyOpen.(o_{1}) ./ Dynamic_currentPolicyClosed.(o_{1});
                    v_Static(isnan(v_Static)) = 0;
                end
                
                % Consolidate, shift, trim, and pad dynamic and static variable values to form data series
                dataseries.(o_{1}) = [ ones(1, nprepad), v_Dynamic(1:end-ntrim), ones(1, npostpad) ;
                                       ones(1, nprepad), v_Static( 1:end-ntrim), ones(1, npostpad) ]';
                
            end
            
        end
        
        
        
        % Iterate over output scenarios
        n_scenarios = length(output_scenarios);
        n_failed = 0;
        for output_scenario = output_scenarios
            
            % Identify corresponding dynamic model scenario
            scenario = CombinationSolver.convertOutputScenario(output_scenario);
            
            try
                
                % Construct data series, using linear combinations of open and closed economy data series for partially open economy scenarios
                if (output_scenario.OpenEconomy == 0 || output_scenario.OpenEconomy == 1)
                    dataseries = constructDataSeries(scenario, output_scenario.UseDynamicBaseline);
                else
                    dataseries_open   = constructDataSeries(scenario.open()  , output_scenario.UseDynamicBaseline);
                    dataseries_closed = constructDataSeries(scenario.closed(), output_scenario.UseDynamicBaseline);
                    for o = fieldnames(dataseries_ids)'
                        dataseries.(o{1}) = dataseries_open  .(o{1})*(output_scenario.OpenEconomy    ) ...
                                          + dataseries_closed.(o{1})*(1 - output_scenario.OpenEconomy);
                    end
                end
                
                % Write data series to output files
                for o = fieldnames(dataseries_ids)'
                    for series_id = dataseries_ids.(o{1})
                        csvfile = fullfile(outputdir, sprintf('%u-%u.csv', output_scenario.id, series_id));
                        fid = fopen(csvfile, 'w'); fprintf(fid, 'Year,Dynamic,Static\n'); fclose(fid);
                        dlmwrite(csvfile, [years, dataseries.(o{1})], '-append')
                    end
                end
                
            catch
                
                % Increment failure counter
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