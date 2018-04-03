%%
% Saves equilibrium transitions into a csv long format suitable for
% consumption by the microsim.
% 
% Test with:
% 
% TransitionMatrixWriter.test()
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
            if scenario.isSteady()
                [productivity_index, savings_index, earnings_index, age_index] = ind2sub(size(OPTs.K), (1:numel(OPTs.K))');
                decision_index = [productivity_index, savings_index, earnings_index, age_index, zeros(size(age_index))];
            else
                [productivity_index, savings_index, earnings_index, age_index, time_index] = ind2sub(size(OPTs.K), (1:numel(OPTs.K))');
                decision_index = [productivity_index, savings_index, earnings_index, age_index, time_index];
            end

            decision_rules = [OPTs.CON(:), OPTs.K(:), OPTs.LAB(:)];

            % melt productivity values
            productivity_values = ParamGenerator.grids(scenario).zs;
            [productivity_index, age_index] = ind2sub(size(productivity_values), (1:numel(productivity_values))');
            productivity_values = [productivity_index, age_index, productivity_values(:)];
            
            % melt productivity transitions
            z_transitions = ParamGenerator.grids(scenario).transz;
            [productivity_index, productivity_prime_index, age_index] = ind2sub(size(z_transitions), (1:numel(z_transitions))');
            productivity_transitions = [productivity_index, productivity_prime_index, age_index, z_transitions(:)];

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
            filehandle = fopen(fullfile(outputdir, 'map.csv'), 'a');
            fprintf(filehandle, [scenario.basedeftag, ',', scenario.counterdeftag, ',', text]);
            fprintf(filehandle, '\n');
            fclose(filehandle);
            
            % store all output in subfolder
            outputdir = fullfile(outputdir, scenario.basedeftag, scenario.counterdeftag);
            
            % create a folder to store output
            mkdir(outputdir)

            h5create(fullfile(outputdir, 'data.hdf5'), '/decision_rules', size(decision_rules), 'ChunkSize', size(decision_rules), 'Deflate', 9);
            h5write(fullfile(outputdir, 'data.hdf5'), '/decision_rules', decision_rules);
            
            h5create(fullfile(outputdir, 'data.hdf5'), '/decision_index', size(decision_index), 'ChunkSize', size(decision_index), 'Deflate', 9);
            h5write(fullfile(outputdir, 'data.hdf5'), '/decision_index', decision_index);
            
            h5create(fullfile(outputdir, 'data.hdf5'), '/productivity_values', size(productivity_values), 'ChunkSize', size(productivity_values), 'Deflate', 9);
            h5write(fullfile(outputdir, 'data.hdf5'), '/productivity_values', productivity_values);
            
            h5create(fullfile(outputdir, 'data.hdf5'), '/earnings_values', size(ParamGenerator.grids(scenario).bv), 'ChunkSize', size(ParamGenerator.grids(scenario).bv), 'Deflate', 9);
            h5write(fullfile(outputdir, 'data.hdf5'), '/earnings_values', ParamGenerator.grids(scenario).bv);
            
            h5create(fullfile(outputdir, 'data.hdf5'), '/savings_values', size(ParamGenerator.grids(scenario).kv), 'ChunkSize', size(ParamGenerator.grids(scenario).kv), 'Deflate', 9);
            h5write(fullfile(outputdir, 'data.hdf5'), '/savings_values', ParamGenerator.grids(scenario).kv);
            
            h5create(fullfile(outputdir, 'data.hdf5'), '/productivity_transitions', size(productivity_transitions), 'ChunkSize', size(productivity_transitions), 'Deflate', 9);
            h5write(fullfile(outputdir, 'data.hdf5'), '/productivity_transitions', productivity_transitions);

        end
 
        %% test
        %       Tests writing of the interface to the transition matrix
        %       output folder.
        function [] = test()
            scenario = Scenario(ModelTester.test_params).currentPolicy().closed();
            ModelSolver.solve(scenario);
            TransitionMatrixWriter.writeScenario(scenario);
        end

    end
    
end





