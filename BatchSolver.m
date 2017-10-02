%%
% Dynamic model batch scenario solver and data series output generator.
%
%%
classdef BatchSolver


properties (Constant)
    
    % Define scenario directory path, file path generator, and file lister
    scenariodir   = fullfile(PathFinder.getSourceDir(), 'Scenarios');
    scenariofile  = @(iscenario) fullfile(BatchSolver.scenariodir, sprintf('scenario%04d.mat', iscenario));
    scenariofiles = @() dir(fullfile(BatchSolver.scenariodir, 'scenario*.mat'));
    
end



methods (Static)
    
    
    % Read batch of scenarios from database
    function [rows] = readBatch(batch)
        
        % Add JDBC driver to Matlab Java path
        javaaddpath(fullfile(PathFinder.getSourceDir(), 'jar', 'sqljdbc41.jar'));
        
        % Establish database connection
        connection = database('second_chart', 'pwbm', 'HbXk86rabjehD2AN', ...
                              'Vendor', 'Microsoft SQL Server', 'AuthType', 'Server', ...
                              'Server', 'ppi-slcsql.wharton.upenn.edu', 'PortNumber', 49170);
        
        % Get scenario rows from database using stored procedure
        o = connection.exec( sprintf( 'EXEC p_ScenarioBatch %u', batch ) );
        rows = cell2struct( o.fetch().Data, o.columnnames(true), 2 );
        o.close();
        
        % Close database connection
        connection.close();
        
    end
    
    
    
    % Define minimal set of executable scenarios for a batch
    function [scenarios] = defineScenarios(batch)
        
        % Read batch of scenarios from database
        rows = BatchSolver.readBatch(batch);
        
        % Preload calibration grid for parameter inversion
        [~, f_invert] = ModelCalibrator.invert();
        
        % Initialize cell array of scenarios
        scenarios = cell(size(rows));
        
        % Define function to remove empty entries from cell array of scenarios
        compress = @(scenarios_) scenarios_(~cellfun(@isempty, scenarios_));
        
        for i = 1:length(rows)
            
            % Identify economy openness, skipping scenarios that are neither fully open nor fully closed
            switch (rows(i).OpenEconomy)
                case 1, economy = 'open'  ;
                case 0, economy = 'closed';
                otherwise, continue
            end
            
            % Invert elasticities with calibration grid
            params = f_invert(struct(...
                'savelas', rows(i).SavingsElasticity, ...
                'labelas', rows(i).LaborElasticity  ));
            
            % Construct dynamic model scenario
            scenario = Scenario(struct(...
                'economy'           , economy                       , ...
                'beta'              , params.beta                   , ...
                'gamma'             , params.gamma                  , ...
                'sigma'             , params.sigma                  , ...
                'modelunit_dollar'  , params.modelunit_dollar       , ... 
                'bequest_phi_1'     , 0                             , ... % To be populated from parameter inversion
                'gcut'              , -rows(i).ExpenditureShift     ));
            
            % Store scenario in cell array
            scenarios{i} = scenario;
            
        end
        scenarios = compress(scenarios);
        
        % Remove duplicate scenarios
        for i = 1:length(scenarios)
            for j = i+1:length(scenarios)
                if scenarios{j}.isEqual(scenarios{i})
                    scenarios{i} = [];
                    break;
                end
            end
        end
        scenarios = compress(scenarios);
        
        % Remove open economy scenarios with a corresponding closed economy scenario
        for i = 1:length(scenarios)
            if strcmp(scenarios{i}.economy, 'open')
                scenario_closed = scenarios{i}.closed();
                for j = 1:length(scenarios)
                    if ( i ~= j && ~isempty(scenarios{j}) && scenarios{j}.isEqual(scenario_closed) )
                        scenarios{i} = [];
                        break;
                    end
                end
            end
        end
        scenarios = compress(scenarios);
        
        % Remove current policy scenarios with a corresponding non-current policy scenario
        for i = 1:length(scenarios)
            if scenarios{i}.isCurrentPolicy()
                for j = 1:length(scenarios)
                    if ( i ~= j && ~isempty(scenarios{j}) && scenarios{j}.currentPolicy().isEqual(scenarios{i}) )
                        scenarios{i} = [];
                        break;
                    end
                end
            end
        end
        scenarios = compress(scenarios);
        
        
        % Clear or create scenario directory
        if exist(BatchSolver.scenariodir, 'dir'), rmdir(BatchSolver.scenariodir, 's'), end, mkdir(BatchSolver.scenariodir)
        
        % Save scenarios to scenario files
        for iscenario = 1:length(scenarios)
            scenario = scenarios{iscenario}; %#ok<NASGU>
            save( BatchSolver.scenariofile(iscenario), 'scenario' );
        end
        
    end
    
    
    
    % Solve dynamic model for scenario indexed in scenario directory
    function [] = solve(iscenario)
        
        % Load scenario
        s = load(BatchSolver.scenariofile(iscenario));
        scenario = s.scenario;
        
        % Solve dynamic model for scenario
        workingdir = ModelSolver.solve(scenario);
        
        % Extract and save solver termination condition
        iterations = csvread(fullfile(workingdir, 'iterations.csv'));
        
        termination.iter = iterations(end, 1);
        termination.eps  = iterations(end, 2);
        
        termination = termination; %#ok<ASGSL,NASGU>
        save(BatchSolver.scenariofile(iscenario), '-append', 'termination');
        
    end
    
    
    
    % Check solver termination conditions for scenarios in scenario directory
    function [] = checkTerminations()
        
        % Identify number of scenario files
        nscenario = length(BatchSolver.scenariofiles());
        
        % Initialize cell array of termination conditions
        terminations = cell(0,6);
        
        for iscenario = 1:nscenario
            
            % Load scenario and termination condition
            s = load(BatchSolver.scenariofile(iscenario));
            
            % Generate baseline and counterfactual definition tags
            [basedef_tag, counterdef_tag] = s.scenario.generate_tags();
            
            % Store termination condition, setting default values if missing
            if ~isfield(s, 'termination'), s.termination = struct('iter', Inf, 'eps', Inf); end
            terminations = [terminations; {iscenario, basedef_tag, counterdef_tag, s.scenario.economy, s.termination.iter, s.termination.eps}]; %#ok<AGROW>
            
        end
        
        % Sort termination conditions by increasing iterations and error terms
        [~, sortinds] = sortrows(cell2mat(terminations(:,5:6)));
        
        % Save termination conditions to csv file
        fid = fopen(fullfile(BatchSolver.scenariodir, 'terminations.csv'), 'w');
        fprintf(fid, 'ScenarioIndex,BaselineDefinition,CounterfactualDefinition,Economy,TerminationIteration,TerminationErrorTerm\n');
        for iscenario = 1:nscenario, fprintf(fid, '%d,%s,%s,%s,%d,%0.4f\n', terminations{sortinds(iscenario),:}); end
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
        rows = BatchSolver.readBatch(batch);
        
        % Preload calibration grid for parameter inversion
        [~, f_invert] = ModelCalibrator.invert();
        
        for i = 1:length(rows)
            
            % Identify economy openness, skipping scenarios that are neither fully open nor fully closed
            switch rows(i).OpenEconomy
                case 1, economy = 'open'  ;
                case 0, economy = 'closed';
                otherwise, continue
            end
            
            % Invert elasticities with calibration grid
            params = f_invert(struct(...
                'savelas', rows(i).SavingsElasticity, ...
                'labelas', rows(i).LaborElasticity  ));
            
            % Construct dynamic model scenario
            scenario = Scenario(struct(...
                'economy'           , economy                       , ...
                'beta'              , params.beta                   , ...
                'gamma'             , params.gamma                  , ...
                'sigma'             , params.sigma                  , ...
                'modelunit_dollar'  , params.modelunit_dollar       , ... 
                'bequest_phi_1'     , 0                             , ... % To be populated from parameter inversion
                'gcut'              , -rows(i).ExpenditureShift     ));
            
            
            
            % Identify scenario working directory
            workingdir = PathFinder.getWorkingDir(scenario);
            
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
            
            
            
            for name_ = fieldnames(dataseriesmap)'
                
                % Extract variable name
                name = name_{1};
                
                % Consolidate, shift, trim, and pad dynamic and static variable values to form data series
                dataseries = [ ones(1, nshift), Dynamic.(name)(1:end-ntrim), ones(1, npad) ;
                               ones(1, nshift), Static. (name)(1:end-ntrim), ones(1, npad) ]';
                
                % Save data series to csv files
                for id = dataseriesmap.(name)
                    csvfile = fullfile(outputdir, sprintf('%u-%u.csv', rows(i).WithoutDynamicBaseline_Tag, id));
                    fid = fopen(csvfile, 'w'); fprintf(fid, 'Year,Dynamic,Static\n'); fclose(fid);
                    dlmwrite(csvfile, [years, dataseries], '-append')
                end
                
            end
            
            
            
            if ~isnan(rows(i).WithDynamicBaseline_Tag)
                
                % Load additional variables for dynamic baseline scaling
                Dynamic_currentPolicyOpen   = load(fullfile(PathFinder.getWorkingDir(scenario.currentPolicy().open()  ), 'dynamics.mat'));
                Dynamic_currentPolicyClosed = load(fullfile(PathFinder.getWorkingDir(scenario.currentPolicy().closed()), 'dynamics.mat'));
                
                for name_ = fieldnames(dataseriesmap)'
                
                    % Extract variable name
                    name = name_{1};
                    
                    % Scale static variable values for dynamic baseline
                    v = Static.(name) .* Dynamic_currentPolicyOpen.(name) ./ Dynamic_currentPolicyClosed.(name);
                    v(isnan(v)) = 0;
                    
                    % Consolidate, shift, trim, and pad dynamic and scaled static variable values to form data series
                    dataseries = [ ones(1, nshift), Dynamic.(name)(1:end-ntrim), ones(1, npad) ;
                                   ones(1, nshift), v(             1:end-ntrim), ones(1, npad) ]';
                    
                    % Save data series to csv files
                    for id = dataseriesmap.(name)
                        csvfile = fullfile(outputdir, sprintf('%u-%u.csv', rows(i).WithDynamicBaseline_Tag, id));
                        fid = fopen(csvfile, 'w'); fprintf(fid, 'Year,Dynamic,Static\n'); fclose(fid);
                        dlmwrite(csvfile, [years, dataseries], '-append')
                    end
                    
                end
                
            end
            
        end
        
    end
    
    
end


end