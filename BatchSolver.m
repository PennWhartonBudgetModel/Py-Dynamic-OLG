%%
% Dynamic model batch scenario solver and data series output generator.
%
%%
classdef BatchSolver


properties (Constant)
    
    % Define scenario directory path
    scenariodir = fullfile(PathFinder.getSourceDir(), 'Scenarios');
    
    % Define current policy and counterfactual scenario file path generators
    currentpolicyfile  = @(i) fullfile(BatchSolver.scenariodir, sprintf('currentpolicy%05d.mat' , i));
    counterfactualfile = @(i) fullfile(BatchSolver.scenariodir, sprintf('counterfactual%05d.mat', i));
    
end



methods (Static)
    
    
    % Read batch of scenarios from database
    function [scenarios, rows] = readBatch(batch)
        
        % Add JDBC driver to Matlab Java path
        javaaddpath(fullfile(PathFinder.getSourceDir(), 'jar', 'sqljdbc41.jar'));
        
        % Establish database connection
        connection = database('second_chart', 'development', 'HbXk86rabjehD2AN', ...
                              'Vendor', 'Microsoft SQL Server', 'AuthType', 'Server', ...
                              'Server', 'ppi-slcsql.wharton.upenn.edu', 'PortNumber', 49170);
        
        % Get scenario rows from database using stored procedure
        o = connection.exec( sprintf( 'EXEC p_ScenarioBatch "%s"', batch ) );
        rows = cell2struct( o.fetch().Data, o.columnnames(true), 2 );
        o.close();
        
        % Close database connection
        connection.close();
        
        % Initialize cell array of scenarios
        %   Empty values will correspond to rows unaddressable by the dynamic model
        scenarios = cell(size(rows));
        
        for i = 1:length(rows)
            
            row = rows(i);
            
            % Identify economy openness, skipping scenarios that are neither fully open nor fully closed
            switch row.OpenEconomy
                case 1, economy = 'open'  ;
                case 0, economy = 'closed';
                otherwise, continue
            end
            
            % Invert elasticities
            inverse = ParamGenerator.invert (struct(...
                'depreciation'  , row.Depreciation    , ...
                'labelas'       , row.LaborElasticity  ));
            
            % Initialize scenario parameter structure with required parameters
            params = struct(...
                    'economy'           , economy                   ...
                ,   'beta'              , inverse.beta              ...
                ,   'gamma'             , inverse.gamma             ...
                ,   'sigma'             , inverse.sigma             ...
                ,   'bequest_phi_1'     , 0                         ... % *** To be added from elasticity inversion ***
                ,   'depreciation'      , row.Depreciation          ...
                ,   'modelunit_dollar'  , inverse.modelunit_dollar  ...
                );  
            
            
            % Populate scenario parameter structure with optional parameters as available
            colmap = struct(...
                'ExpenditureShift'              , 'expenditure_shift'               , ...
                'BaseBrackets'                  , 'base_brackets'                   , ...
                'HasBuffetRule'                 , 'has_buffet_rule'                 , ...
                'HasDoubleStandardDeduction'    , 'has_double_standard_deduction'   , ...
                'HasLimitDeductions'            , 'has_limit_deductions'            , ...
                'HasExpandChildCredit'          , 'has_expand_child_credit'         , ...
                'NoACAIncomeTax'                , 'no_aca_income_tax'               , ...
                'CorporateTaxRate'              , 'corporate_tax_rate'              , ...
                'HasSpecialPassThroughRate'     , 'has_special_pass_through_rate'   , ...
                'HasImmediateExpensing'         , 'has_immediate_expensing'         , ...
                'HasRepealCorporateExpensing'   , 'has_repeal_corporate_expensing'  );
            
            for o = fieldnames(colmap)'
                colname = o{1};
                if (isfield(row, colname) && ~isempty(row.(colname)) && ~any(isnan(row.(colname))))
                    params.(colmap.(colname)) = row.(colname);
                end
            end
            
            % Construct scenario from parameter structure and store in cell array
            scenarios{i} = Scenario(params);
            
        end
        
    end
    
    
    
    
    % Define minimal set of executable scenarios for a batch
    function [] = defineScenarios(batch)
        
        % Define function to remove empty entries from a cell array
        function c_ = compress(c)
            c_ = c(~cellfun(@isempty, c));
        end
        
        % Read batch of scenarios from database
        scenarios = BatchSolver.readBatch(batch);
        scenarios = compress(scenarios);
        
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
        if exist(BatchSolver.scenariodir, 'dir'), rmdir(BatchSolver.scenariodir, 's'), end, mkdir(BatchSolver.scenariodir)
        
        % Save current policy scenarios
        for i = 1:length(currentpolicys)
            scenario = currentpolicys{i}; %#ok<NASGU>
            save(BatchSolver.currentpolicyfile(i), 'scenario');
        end
        
        % Save counterfactual scenarios
        for i = 1:length(counterfactuals)
            scenario = counterfactuals{i}; %#ok<NASGU>
            save(BatchSolver.counterfactualfile(i), 'scenario');
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
        BatchSolver.solve(BatchSolver.currentpolicyfile(i));
    end
    
    function [] = solveCounterfactual(i)
        BatchSolver.solve(BatchSolver.counterfactualfile(i));
    end
    
    
    
    
    % Check solver termination conditions for scenarios in scenario directory
    function [] = checkTerminations()
        
        % Get all scenario files
        scenariofiles = arrayfun(@(s) fullfile(s.folder, s.name), ...
            [ dir(fullfile(BatchSolver.scenariodir, 'currentpolicy*.mat' )) ;
              dir(fullfile(BatchSolver.scenariodir, 'counterfactual*.mat')) ], ...
            'UniformOutput', false);
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
                  s.termination.eps         }]; %#ok<AGROW>
            
        end
        
        % Sort termination conditions by increasing iterations and error terms
        [~, sortinds] = sortrows(cell2mat(terminations(:,5:6)));
        
        % Save termination conditions to csv file
        fid = fopen(fullfile(BatchSolver.scenariodir, 'terminations.csv'), 'w');
        fprintf(fid, 'ScenarioIndex,BaselineDefinition,CounterfactualDefinition,Economy,TerminationIteration,TerminationErrorTerm\n');
        for i = 1:nscenario, fprintf(fid, '%s,%s,%s,%s,%d,%0.4f\n', terminations{sortinds(i),:}); end
        fclose(fid);
        
    end
    
    
    
    
    % Generate data series output for a batch
    function [] = generateDataSeries(batch)
        
        % Identify data series output directory
        outputdir = PathFinder.getDataSeriesOutputDir();
        
        % Clear or create output directory
        if exist(outputdir, 'dir'), rmdir(outputdir, 's'), end, mkdir(outputdir)
        
        % Define mapping from dynamic model variable names to data series IDs
        dataseriesmap = struct(...
            'labpits'           , [102]         , ...
            'ssts'              , [103]         , ...
            'caprevs'           , [110]         , ...
            'cits_domestic'     , [111]         , ...
            'cits_foreign'      , [112]         , ...
            'caps_domestic'     , [401]         , ...
            'caps_foreign'      , [402]         , ...
            'debts_domestic'    , [403]         , ...
            'debts_foreign'     , [404]         , ...
            'outs'              , [199, 299, 6] , ...
            'bens'              , [202]         , ...
            'caps'              , [1]           , ...
            'labeffs'           , [24]          , ...
            'lfprs'             , [25]          , ...
            'labincs'           , [28]          , ...
            'capincs'           , [29]          ); %#ok<NBRAK>
        
        
        
        % Read batch of scenarios from database
        [scenarios, rows] = BatchSolver.readBatch(batch);
        
        for i = 1:length(rows)
            
            % Extract scenario from cell array of scenarios
            scenario = scenarios{i};
            if isempty(scenario), continue, end
            
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
                first_year = 2016;
                last_year  = 2089;
                years = (first_year : last_year)';
                
                % Determine number of years to shift variable values
                nshift  = ParamGenerator.timing(scenario).first_transition_year - first_year;
                
                % Determine number of variable values to be trimmed or padded
                T_model = ParamGenerator.timing(scenario).T_model;
                nextra  = nshift + T_model - length(years);
                ntrim   = +max(nextra, 0);
                npad    = -min(nextra, 0);
                
                % Write data series output files
                writeFiles(rows(i).WithoutDynamicBaseline_Tag);
                
                if ~isnan(rows(i).WithDynamicBaseline_Tag)
                    
                    % Load additional variables for dynamic baseline scaling
                    Dynamic_currentPolicyOpen   = load(fullfile(PathFinder.getWorkingDir(scenario.currentPolicy().open()  ), 'dynamics.mat'));
                    Dynamic_currentPolicyClosed = load(fullfile(PathFinder.getWorkingDir(scenario.currentPolicy().closed()), 'dynamics.mat'));
                    
                    % Write data series output files with dynamic baseline scaling
                    writeFiles(rows(i).WithDynamicBaseline_Tag, Dynamic_currentPolicyOpen, Dynamic_currentPolicyClosed);
                    
                end
                
            catch
                
                fprintf('Failed to generate data series output for row %6u of batch %s.\n', i, num2str(batch));
                
            end
            
        end
        
        
        function [] = writeFiles(tag, Dynamic_currentPolicyOpen, Dynamic_currentPolicyClosed)
            
            for name_ = fieldnames(dataseriesmap)'
                
                % Extract variable name
                name = name_{1};
                
                % Extract dynamic and static variable values
                v_Dynamic = Dynamic.(name);
                v_Static  = Static. (name);
                
                % Scale static variable values for dynamic baseline if scaling variables provided
                if (exist('Dynamic_currentPolicyOpen', 'var') && exist('Dynamic_currentPolicyClosed', 'var'))
                    v_Static = v_Static .* Dynamic_currentPolicyOpen.(name) ./ Dynamic_currentPolicyClosed.(name);
                    v_Static(isnan(v_Static)) = 0;
                end
                
                % Consolidate, shift, trim, and pad dynamic and static variable values to form data series
                dataseries = [ ones(1, nshift), v_Dynamic(1:end-ntrim), ones(1, npad) ;
                               ones(1, nshift), v_Static( 1:end-ntrim), ones(1, npad) ]';
                
                % Write data series to csv files
                for id = dataseriesmap.(name)
                    csvfile = fullfile(outputdir, sprintf('%u-%u.csv', tag, id));
                    fid = fopen(csvfile, 'w'); fprintf(fid, 'Year,Dynamic,Static\n'); fclose(fid);
                    dlmwrite(csvfile, [years, dataseries], '-append')
                end
                
            end
            
        end
        
    end
    
    
end


end