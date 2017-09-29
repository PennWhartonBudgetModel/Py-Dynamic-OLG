%%
% Dynamic model batch scenario solver and results exporter.
%
%%
classdef BatchSolver


properties (Constant)
    
    % Define run directory and run file path
    run_dir  = fullfile(PathFinder.getSourceDir(), 'Runs');
    run_file = @(irun) fullfile(BatchSolver.run_dir, sprintf('run%04d.mat', irun));
    
end



methods (Static)
    
    % Define full set of production runs
    function [] = define_runs(batchID)
        if( nargin < 1 )
            error( '<batchID> is required');
        end
        
        % Get worklist of Scenarios from batchID
        % Sort them to have closed econ go first -- this is 
        % to take advantage of redundancy in solving.
        scenarios_in_batch = Scenario.fetch_batch(batchID);
        [~, ind] = sort({scenarios_in_batch.economy});
        scenarios = scenarios_in_batch(ind);
        
        % Clear or create run directory
        if exist(BatchSolver.run_dir, 'dir'), rmdir(BatchSolver.run_dir, 's'), end, mkdir(BatchSolver.run_dir)
        
        % Save individual scenarios to run
        %   NOTE: We skip UseDynamicBaseline = 1, since this is only for
        %         post-processing, so would be duplicative work.
        irun = 0;
        for i = 1:numel(scenarios)
            if( ~scenarios(i).useDynamicBaseline )
                irun = irun + 1;
                run_def = scenarios(i).getParams();
                save(BatchSolver.run_file(irun), 'run_def');
            end
        end
        
    end
    
    
    % Perform production run
    function [] = run(irun)
        
        % Fetch Scenario to run
        s           = load(BatchSolver.run_file(irun));
        scenario    = Scenario(s.run_def);
        
        % Execute dynamic model solver
        save_dir = ModelSolver.solve(scenario);
        
        % Extract and save solver termination condition
        iterations = csvread(fullfile(save_dir, 'iterations.csv'));
        
        termination.iter = iterations(end, 1);
        termination.eps  = iterations(end, 2);
        
        termination = termination; %#ok<ASGSL,NASGU>
        save(BatchSolver.run_file(irun), '-append', 'termination');
        
    end
    
    
    % Check production run termination conditions
    function [] = check_terminations()
        
        % Identify production run definition files
        run_files = dir(fullfile(BatchSolver.run_dir, 'run*.mat'));
        nrun = length(run_files);
        
        % Initialize cell array of termination conditions
        terminations = cell(0,5);
        
        for irun = 1:nrun
            
            % Load run definition and termination condition
            s = load(BatchSolver.run_file(irun));
            
            % Generate baseline and counterfactual definition tags
            [basedef_tag, counterdef_tag] = ModelSolver.generate_tags(s.basedef, s.counterdef);
            
            % Store termination condition, setting default values if missing
            if ~isfield(s, 'termination'), s.termination = struct('iter', Inf, 'eps', Inf); end
            terminations = [terminations; {irun, basedef_tag, counterdef_tag, s.termination.iter, s.termination.eps}]; %#ok<AGROW>
            
        end
        
        % Sort termination conditions by increasing iterations and error terms
        [~, sortinds] = sortrows(cell2mat(terminations(:,4:5)));
        
        % Save termination conditions to csv file
        fid = fopen(fullfile(dirFinder.saveroot(), 'terminations.csv'), 'w');
        fprintf(fid, 'Run,Baseline Definition,Counterfactual Definition,Termination Iteration,Termination Error Term\n');
        for irun = 1:nrun, fprintf(fid, '%d,%s,%s,%d,%0.4f\n', terminations{sortinds(irun),:}); end
        fclose(fid);
        
    end
    
    
    %% 
    % Package production run results into csv files for front end deployment
    function [] = package_results(batchID)
        
        % Get full work list
        scenarios = Scenario.fetch_batch(batchID);
        
        % Identify production run definition files
        run_files = dir(fullfile(BatchSolver.run_dir, 'run*.mat'));
        nrun = length(run_files);
        
        % Initialize missing aggregates flag
        missing = false;
        
        for irun = 1:nrun
            fprintf('Processing run ID %6d of %6d\n', irun, nrun)
            try
                BatchSolver.export_run(irun);
            catch
                missing = true;
                continue
            end
        end
        
        fprintf('\nResults packaged into csv files:\n');
        if missing, warning('Some aggregates not found'), end
        
    end % package_results
    
    
    %%
    % Generate exports for a particular scenario
    %    useDynamicBaseline is post-processing step
    function [] = export_results(irun)
        
        % Fetch scenario to export
        s           = load(BatchSolver.run_file(irun));
        scenario    = Scenario(s.run_def);
        
        % Identify export directory
        exportdir = PathFinder.getDataSeriesOutputDir();
        
        % Clear or create export directory
        if exist(exportdir, 'dir'), rmdir(exportdir, 's'), end, mkdir(exportdir)
    
        % Specify identifier strings and numerical codes for aggregates
        codes = { 102, 'labpits'        ;
                  103, 'ssts'           ;
                  110, 'caprevs'        ;
                  111, 'cits_domestic'  ;
                  112, 'cits_foreign'   ;
                  401, 'caps_domestic'  ;
                  402, 'caps_foreign'   ;
                  403, 'debts_domestic' ;
                  404, 'debts_foreign'  ;
                  199, 'outs'           ;
                  202, 'bens'           ;
                  299, 'outs'           ;
                    1, 'caps'           ;
                    6, 'outs'           ;
                   24, 'labeffs'        ;
                   25, 'lfprs'          ;
                   28, 'labincs'        ;
                   29, 'capincs'        };
        
        % Specify projection years for csv files
        first_year = 2016;    % first data point in file, should be before transition
        years_csv  = (first_year : 2089)';
        nyear      = length(years_csv);
        
        % Specify number of years to shift results
        nshift = ParamGenerator.timing(scenario).first_transition_year - first_year;
        
        % Get dynamic baseline flag
        %   Applicable to closed economy runs only
        usedynamicbaseline = scenario.useDynamicBaseline && strcmp( scenario.economy, 'closed' );
            
        % Identify working directories
        savedir = PathFinder.getWorkingDir(scenario);
            
        % Load aggregates
        Dynamic = load(fullfile(savedir, 'dynamics.mat'));
        isbase  = scenario.isCurrentPolicy();
        if ~isbase
            Static = load(fullfile(savedir, 'statics.mat'));
        end
            
        if usedynamicbaseline
            base_scenario       = scenario.currentPolicy();
            Dynamic_open_base   = load(fullfile(PathFinder.getWorkingDir(base_scenario.open()  ), 'dynamics.mat'));
            Dynamic_closed_base = load(fullfile(PathFinder.getWorkingDir(base_scenario.closed()), 'dynamics.mat'));
        end
            
        % Find number of entries to be trimmed or padded
        T_model = ParamGenerator.timing(scenario).T_model;
        nextra  = nshift + T_model - nyear;
        ntrim   =  max(nextra, 0);
        npad    = -min(nextra, 0);
            
        % Make the bacon
        for j = 1:size(codes,1)

            % Get aggregate name
            name = codes{j,2};

            % Extract dynamic and static aggregate series
            a.dynamic = Dynamic.(name);
            if isbase
                a.static = a.dynamic;
            else
                a.static = Static.(name);
            end

            % Adjust static aggregate series if using dynamic baseline
            if usedynamicbaseline
                a.static = a.static .* Dynamic_open_base.(name) ./ Dynamic_closed_base.(name);
                a.static(isnan(a.static)) = 0;
            end

            % Consolidate, shift, trim, and pad aggregate series
            agg_series  = [ ones(nshift, 2); 
                            [a.dynamic(1:end-ntrim)', a.static(1:end-ntrim)'];
                            ones(npad,   2) ];

            % Save aggregate series to csv file
            csvfile = fullfile(exportdir, sprintf('%u-%u.csv', scenario.ID, codes{j,1}));
            fid = fopen(csvfile, 'w'); fprintf(fid, 'Year,DynamicAggregate,StaticAggregate\n'); fclose(fid);
            dlmwrite(csvfile, [years_csv, agg_series], '-append')

        end % make a series

    end
    
    
    
    % Get batch of scenarios from database
    function [scenarios] = getBatch(batch)
        
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
        
        % Preload calibration grid for parameter inversion
        [~, f_invert] = ModelCalibrator.invert();
        
        % Initialize cell array of scenarios
        c = cell(size(rows));
        
        for i = 1:length(rows)
            
            % Identify economy openness, skipping scenarios that are neither fully open nor fully closed
            switch (rows(i).OpenEconomy)
                case 1, economy = 'open'  ;
                case 0, economy = 'closed';
                otherwise, continue
            end
            
            % Invert elasticities with calibration grid
            params = f_invert(struct(...
                'savelas', rows(i).SavingsElasticity,   ...
                'labelas', rows(i).LaborElasticity      ));
            
            c{i} = Scenario(struct(...
                'economy'           , economy                       , ...
                'beta'              , params.beta                   , ...
                'gamma'             , params.gamma                  , ...
                'sigma'             , params.sigma                  , ...
                'modelunit_dollar'  , params.modelunit_dollar       , ... 
                'bequest_phi_1'     , 0                             , ... % To be populated from parameter inversion
                'gcut'              , -rows(i).ExpenditureShift     ));
            
        end
        
        % Construct array of relevant scenarios from nonempty values of cell array of scenarios
        scenarios = [c{~cellfun(@isempty, c)}]';
        
    end

end


end