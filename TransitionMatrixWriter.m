%%
% Saves equilibrium transitions into a csv long format suitable for
% consumption by the microsim.
% 
% Test with:
% 
% scenario = Scenario(ModelTester.test_params).currentPolicy().closed();
% ModelSolver.solve(scenario);
% TransitionMatrixWriter.writeScenario(scenario);
% 
%%

classdef TransitionMatrixWriter
    
    methods (Static)

        %% writeScenario
        %       Writes an already-solved scenario's optimal decision rules
        %       and productivity transitions to file in a long format.
        function [] = writeScenario(scenario)
            
            working_dir = PathFinder.getWorkingDir(scenario);
            load(fullfile(working_dir, 'decisions.mat'), 'OPTs');
            
            % melt decision rules into long format
            % in steady state, no time dimension and there are 4 columns
            % on transition path, a time dimension, and 5 columns
            if size(size(OPTs.K), 2) == 4
                [productivity_index, savings_index, earnings_index, age_index] = ind2sub(size(OPTs.K), (1:numel(OPTs.K))');
                index = [productivity_index, savings_index, earnings_index, age_index, zeros(size(age_index))];
            elseif size(size(OPTs.K), 2) == 5
                [productivity_index, savings_index, earnings_index, age_index, time_index] = ind2sub(size(OPTs.K), (1:numel(OPTs.K))');
                index = [productivity_index, savings_index, earnings_index, age_index, time_index];
            else
                error("Unrecognized format for decision rules")
            end

            % stack index and rules together
            decision_rules = [index, OPTs.CON(:), OPTs.K(:), OPTs.LAB(:)];

            % melt productivity values
            productivity_values = ParamGenerator.grids(scenario).zs;
            
            [productivity_index, age_index] = ind2sub(size(productivity_values), (1:numel(productivity_values))');
            productivity_values = [productivity_index, age_index, productivity_values(:)];
            
            % save earnings and savings values directly
            earnings_values = ParamGenerator.grids(scenario).kv;
            savings_values = ParamGenerator.grids(scenario).bv;

            % melt productivity transitions into long format
            z_transitions = ParamGenerator.grids(scenario).transz;
            [productivity_index, productivity_prime_index, age_index] = ind2sub(size(z_transitions), (1:numel(z_transitions))');
            split_transitions = num2cell(ParamGenerator.grids(scenario).transz);
            productivity_transitions = [productivity_index, productivity_prime_index, age_index, vertcat(split_transitions{:})];

            % write to file
            outputdir = PathFinder.getDataSeriesOutputDir();
            
            TransitionMatrixWriter.writeToOutputFormat(             ...
                fullfile(outputdir, 'decision_rules.csv'),          ...
                decision_rules,                                     ...
                [                                                   ...
                    "productivity_index",                           ...
                    "savings_index",                                ...
                    "earnings_index",                               ...
                    "age_index",                                    ...
                    "time",                                         ...
                    "consumption",                                  ...
                    "savings",                                      ...
                    "labor"]);
            
            TransitionMatrixWriter.writeToOutputFormat(             ...
                fullfile(outputdir, 'productivity_values.csv'),     ...
                productivity_values,                                ...
                ["productivity_index", "age_index", "value"]);
            
            TransitionMatrixWriter.writeToOutputFormat(             ...
                fullfile(outputdir, 'earnings_values.csv'),         ...
                earnings_values,                                    ...
                ["earnings_index", "value"]);
        
            TransitionMatrixWriter.writeToOutputFormat(             ...
                fullfile(outputdir, 'savings_values.csv'),          ...
                savings_values,                                     ...
                ["savings_index", "value"]);
                
            TransitionMatrixWriter.writeToOutputFormat(             ...
                fullfile(                                           ...
                    outputdir,                                      ...
                    'productivity_transitions.csv'),                                                  ...
                productivity_transitions,                           ...
                [                                                   ...
                    "productivity_index",                           ...
                    "productivity_prime_index",                     ...
                    "age_index",                                    ...
                    "probability"]);
            
        end
 
        %% writeToOutputFormat
        %       Wraps the writing of matrix data structures to file.
        %       For private use only.
        function [] = writeToOutputFormat(filepath, data, headers)
           
            filehandle = fopen(filepath, 'w');
            fprintf(filehandle, strjoin(headers, ','));
            fprintf(filehandle, '\n');
            fclose(filehandle);
            
            dlmwrite(filepath, data, '-append');
            
        end

    end
    
end





