%%
% Dynamic model bundle solver and data series output generator.
%
%%
classdef BundleSolver


properties (Constant)
    
    % Define scenario directory path
    scenariodir = fullfile(PathFinder.getSourceDir(), 'Scenarios');
    
    % Define current policy and counterfactual scenario file path generators
    currentpolicyfile  = @(i) fullfile(BundleSolver.scenariodir, sprintf('currentpolicy%05d.mat' , i));
    counterfactualfile = @(i) fullfile(BundleSolver.scenariodir, sprintf('counterfactual%05d.mat', i));
    
end



methods (Static)
    
    
    
    % Solve a bundle
    function [] = solveBundle(bundle_id)
        
        [currentpolicys, counterfactuals] = BundleSolver.generateScenarioSet(bundle_id);
        
        for i = 1:length(currentpolicys ), BundleSolver.solveCurrentPolicy (i); end
        for i = 1:length(counterfactuals), BundleSolver.solveCounterfactual(i); end
        
        BundleSolver.generateDataSeries(bundle_id);
        
    end
    
    
    
    
    % Get bundle scenarios
    function [bundle_scenarios] = readBundle(bundle_id)
        
        % Read bundle scenarios file
        bundle_scenarios = readtable(fullfile(PathFinder.getBundleDir(bundle_id), 'scenarios.csv'));
        
        % Filter scenarios for those with dynamic model type
        bundle_scenarios = bundle_scenarios(strcmp(bundle_scenarios.ModelType, 'D'), :);
        
        % Keep only scenario parameters relevant to dynamic model
        bundle_scenarios = bundle_scenarios(:, {...
            'UseDynamicBaseline'            , ...
            'OpenEconomy'                   , ...
            'LaborElasticity'               , ...
            'IsLowReturn'                   , ...
            'ExpenditureShift'              , ...
            'BaseBrackets'                  , ...
            'HasBuffetRule'                 , ...
            'HasDoubleStandardDeduction'    , ...
            'HasLimitDeductions'            , ...
            'HasExpandChildCredit'          , ...
            'NoACAIncomeTax'                , ...
            'CorporateTaxRate'              , ...
            'HasSpecialPassThroughRate'     , ...
            'HasImmediateExpensing'         , ...
            'HasRepealCorporateExpensing'   ...
        });
        
        % Remove duplicate scenarios
        bundle_scenarios = unique(bundle_scenarios);
        
    end
    
    
    
    
    % Convert a bundle scenario, represented as a structure, into a dynamic model scenario
    function [scenario] = convertBundleScenario(s)
        
        % Identify economy openness, specifying closed economy for partially open economies to generate both open and closed economies
        if s.OpenEconomy == 1
            economy = 'open';
        else
            economy = 'closed';
        end
        
        % Invert elasticities
        inverse = ParamGenerator.invert(struct(...
            'is_low_return'  , s.IsLowReturn     , ...
            'labelas'        , s.LaborElasticity ...
        ));
        
        % Initialize scenario parameter structure with required parameters
        parameters = struct(...
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
        
        % Add optional parameters, assuming dynamic model parameter names match bundle scenario parameter names
        optional_parameters = {...
            'ExpenditureShift'              , ...
            'BaseBrackets'                  , ...
            'HasBuffetRule'                 , ...
            'HasDoubleStandardDeduction'    , ...
            'HasLimitDeductions'            , ...
            'HasExpandChildCredit'          , ...
            'NoACAIncomeTax'                , ...
            'CorporateTaxRate'              , ...
            'HasSpecialPassThroughRate'     , ...
            'HasImmediateExpensing'         , ...
            'HasRepealCorporateExpensing'   ...
        };
        for o = optional_parameters
            parameters.(o{1}) = s.(o{1});
        end
        
        % Construct scenario from parameter structure
        scenario = Scenario(parameters);
        
    end
    
    
    
    
    % Generate minimal set of executable scenarios for a bundle
    function [currentpolicys, counterfactuals] = generateScenarioSet(bundle_id)
        
        % Define function to remove empty entries from a cell array
        compress = @(c) c(~cellfun(@isempty, c));
        
        % Get bundle scenarios
        bundle_scenarios = BundleSolver.readBundle(bundle_id);
        
        % Convert bundle scenarios to dynamic model scenarios
        scenarios = arrayfun(@BundleSolver.convertBundleScenario, table2struct(bundle_scenarios), 'UniformOutput', false);
        
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
                for j = 1:i-1
                    if (~isempty(currentpolicys{j}) && currentpolicys{i}.isEquivalent(currentpolicys{j}))
                        currentpolicys{i} = [];
                        break;
                    end
                end
                counterfactuals{i} = scenarios{i};
            end
        end

        currentpolicys  = compress(currentpolicys );
        counterfactuals = compress(counterfactuals);
        
        % Clear or create scenario directory
        if exist(BundleSolver.scenariodir, 'dir'), rmdir(BundleSolver.scenariodir, 's'), end, mkdir(BundleSolver.scenariodir)
        
        % Save current policy scenarios
        for i = 1:length(currentpolicys)
            scenario = currentpolicys{i}; %#ok<NASGU>
            save(BundleSolver.currentpolicyfile(i), 'scenario');
        end
        
        % Save counterfactual scenarios
        for i = 1:length(counterfactuals)
            scenario = counterfactuals{i}; %#ok<NASGU>
            save(BundleSolver.counterfactualfile(i), 'scenario');
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
        BundleSolver.solve(BundleSolver.currentpolicyfile(i));
    end
    
    function [] = solveCounterfactual(i)
        BundleSolver.solve(BundleSolver.counterfactualfile(i));
    end
    
    
    
    
    % Check solver termination conditions for scenarios in scenario directory
    function [] = checkTerminations()
        
        % Get all scenario files
        scenariofiles = arrayfun(@(s) fullfile(s.folder, s.name), ...
            [ dir(fullfile(BundleSolver.scenariodir, 'currentpolicy*.mat' )) ;
              dir(fullfile(BundleSolver.scenariodir, 'counterfactual*.mat')) ], ...
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
        fid = fopen(fullfile(BundleSolver.scenariodir, 'terminations.csv'), 'w');
        fprintf(fid, 'ScenarioIndex,BaselineDefinition,CounterfactualDefinition,Economy,TerminationIteration,TerminationErrorTerm\n');
        for i = 1:nscenario, fprintf(fid, '%s,%s,%s,%s,%d,%0.4f\n', terminations{sortinds(i),:}); end
        fclose(fid);
        
    end
    
    
    
    
    % Generate data series output for a bundle
    function [] = generateDataSeries(bundle_id)
        
        % Get bundle scenarios
        bundle_scenarios = BundleSolver.readBundle(bundle_id);
        
        % Identify data series output directory
        outputdir = PathFinder.getDataSeriesOutputDir();
        
        % Create output directory if not already present
        if ~exist(outputdir, 'dir'), mkdir(outputdir), end
        
        % Load or initialize output mapping structure
        mapfile = fullfile(outputdir, 'map.csv');
        if exist(mapfile, 'file')
            map = readtable(mapfile, 'ReadRowNames', true);
        else
            map = cell2table(cell(0, size(bundle_scenarios, 2)), 'VariableNames', bundle_scenarios.Properties.VariableNames);
            map.Properties.DimensionNames = {'id', 'ScenarioParameters'};
        end
        
        % Determine number of scenarios already mapped
        n_mapped = height(map);
        
        % Add bundle scenarios to mapping structure, removing duplicates
        %   Original mapping structure assumed to be free of duplicates
        map = unique([map; bundle_scenarios], 'stable');
        
        % Define mapping from dynamic model variable names to data series IDs
        dataseries_ids = struct(...
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
        
        % Iterate over new bundle scenarios
        i = n_mapped;
        n_new = height(map) - n_mapped; n_failed = 0;
        
        while i < height(map)
            
            i = i+1;
            
            map.Properties.RowNames{i} = sprintf('%05d', i);
            bundle_scenario = map(i,:);
            
            % Identify corresponding dynamic model scenario
            scenario = BundleSolver.convertBundleScenario(table2struct(bundle_scenario));
            
            % Identify scenario working directory
            workingdir = PathFinder.getWorkingDir(scenario);
            
            try
                
                % Load dynamic and static variables
                Dynamic = load(fullfile(workingdir, 'dynamics.mat'));
                if scenario.isCurrentPolicy()
                    Static = Dynamic;
                else
                    Static = load(fullfile(workingdir, 'statics.mat'));
                end
                
                % Specify data series years
                first_year  = 2017;
                last_year   = 2090;
                years       = (first_year : last_year)';
                timing      = ParamGenerator.timing( scenario );
                
                % Determine number of years to pre-pad variable values
                nprepad     = timing.TransitionFirstYear - first_year;
                
                % Determine number of variable values to be trimmed or post-padded
                T_model     = timing.T_model;
                nextra      = nprepad + T_model - length(years);
                ntrim       = +max(nextra, 0);
                npostpad    = -min(nextra, 0);
                
                % Load additional variables if performing dynamic baseline scaling
                if bundle_scenario.UseDynamicBaseline
                    Dynamic_currentPolicyOpen   = load(fullfile(PathFinder.getWorkingDir(scenario.currentPolicy().open()  ), 'dynamics.mat'));
                    Dynamic_currentPolicyClosed = load(fullfile(PathFinder.getWorkingDir(scenario.currentPolicy().closed()), 'dynamics.mat'));
                else
                    Dynamic_currentPolicyOpen   = {};
                    Dynamic_currentPolicyClosed = {};
                end
                
                % Write data series output files
                for o = fieldnames(dataseries_ids)'
                    
                    % Extract variable name
                    name = o{1};
                    
                    % Extract dynamic and static variable values and handle special series
                    switch name
                        case 'nonindexed'
                            v_Dynamic = [ones(1,10), Dynamic.outs(11:T_model)]; 
                            v_Static  = [ones(1,10), Dynamic.outs(11:T_model)]; 
                        otherwise
                            v_Dynamic = Dynamic.(name);
                            v_Static  = Static. (name);
                    end
                    
                    % Scale static variable values if performing dynamic baseline scaling
                    if bundle_scenario.UseDynamicBaseline
                        v_Static = v_Static .* Dynamic_currentPolicyOpen.(name) ./ Dynamic_currentPolicyClosed.(name);
                        v_Static(isnan(v_Static)) = 0;
                    end
                    
                    % Consolidate, shift, trim, and pad dynamic and static variable values to form data series
                    dataseries = [ ones(1, nprepad), v_Dynamic(1:end-ntrim), ones(1, npostpad) ;
                                   ones(1, nprepad), v_Static( 1:end-ntrim), ones(1, npostpad) ]';
                    
                    % Write data series to csv files
                    for series_id = dataseries_ids.(name)
                        csvfile = fullfile(outputdir, sprintf('%s-%u.csv', map.Properties.RowNames{i}, series_id));
                        fid = fopen(csvfile, 'w'); fprintf(fid, 'Year,Dynamic,Static\n'); fclose(fid);
                        dlmwrite(csvfile, [years, dataseries], '-append')
                    end
                    
                end
                
            catch
                
                % Increment failure counter and delete row from mapping structure
                n_failed = n_failed+1;
                map(i,:) = [];
                i = i-1;
                
            end
            
        end
        
        % Report on output generation success
        if ~n_failed
            fprintf('Successfully completed data series output generation for all %d new scenarios.\n', n_new);
        else
            fprintf('Failed to complete data series output generation for %d of %d scenarios.\n', n_failed, n_new);
        end
        
        % Write mapping structure to output directory
        if ~isempty(map), writetable(map, fullfile(outputdir, 'map.csv'), 'WriteRowNames', true); end
        
        % Generate interface dependencies file if not already present
        dependenciesfile = fullfile(outputdir, 'dependencies.csv');
        if ~exist(dependenciesfile, 'file')
            fid = fopen(dependenciesfile, 'w');
            fprintf(fid, 'Component,Interface,Version\n');
            for r = PathFinder.getInputSet()
                fprintf(fid, '%s,%s,%s\n', r{1}{1}, r{1}{2}, r{1}{3});
            end
            fclose(fid);
        end
        
    end
    
    
end


end