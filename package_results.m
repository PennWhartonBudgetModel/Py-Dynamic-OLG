%%
% Package all results into csv files for staging.
% 
%%


function [] = package_results()

% Identify csv save directory
csv_dir = dirFinder.csv;

% Clear or create save directory
if exist(csv_dir, 'dir')
    rmdir(csv_dir, 's')
end
mkdir(csv_dir)


% Load target elasticity sets
s = load(fullfile(dirFinder.param, 'ss_inverses.mat'));
elasticity_sets = s.targets(:,[2,3]);
deep_param_sets = s.inverses;
clear('s')


% Specify identifier strings and numerical codes for aggregates
agg_codes = { 102, 'feditlab'         ;
              103, 'ssrev'            ;
              110, 'fcaprev'          ;
              111, 'domestic_fcaptax' ;
              112, 'foreign_fcaptax'  ;
              401, 'domestic_cap'     ;
              402, 'foreign_cap'      ;
              403, 'domestic_debt'    ;
              404, 'foreign_debt'    ;
              199, 'Y'                ;
              202, 'ssexp'            ;
              299, 'Y'                ;
                1, 'cap'              ;
                6, 'Y'                ;
               24, 'elab'             ;
               25, 'lfpr'             ;
               28, 'labinc'           ;
               29, 'kinc'             };
n_aggs = size(agg_codes, 1);

% Specify projection years for csv files
years_csv = (2015 : 2089)';
n_csv = length(years_csv);

% Specify number of years to shift results
% (Currently, shifting up from 2015 to 2017)
upshift = 2;


% Initialize ID value
ID = 3000000 - 1;

for gcut = [+0.00, +0.10, +0.05, -0.05]
    for plan = {'trump', 'clinton', 'ryan'}
        
        % Skip one ID value representing the static scenario
        ID = ID + 1;
        
        for openness = [0, 0.4, 0.7, 1]
            for labor_elas = [0.25, 0.50, 0.75, 1.00]
                for savings_elas = [0.25, 0.50, 0.75, 1.00]
                    
                    % Increment ID value
                    ID = ID + 1;
                    
                    if (openness == 0) || (openness == 1)
                    
                        % Find elasticity set index, which is equivalent to the deep parameter set index
                        inddeep = all(bsxfun(@eq, [labor_elas, savings_elas], elasticity_sets), 2);

                        % Get deep parameters based on set index
                        deep_params = deep_param_sets(inddeep,:);

                        beta  = deep_params(1);
                        gamma = deep_params(2);
                        sigma = deep_params(3);

                        % Identify working directoroes
                        if (openness == 1)
                            save_dir = dirFinder.open  (beta, gamma, sigma, plan{1});
                        else
                            save_dir = dirFinder.closed(beta, gamma, sigma, plan{1}, gcut);
                        end

                        % Identify iterations log file
                        logfile = fullfile(save_dir, 'iterations.txt');
                        
                        if exist(logfile, 'file')

                            % Get last line of iterations log file
                            fid = fopen(logfile);
                            while ~feof(fid)
                                lastline = fgetl(fid);
                            end
                            fclose(fid);

                            % Check convergence error against threshold
                            lastiter = sscanf(lastline, '  %2d  --  %f', [1,2]);
                            if (lastiter(2) > 0.1), continue, end
                            
                        else
                            continue
                        end
                        
                        % Load dynamic and static aggregates
                        s_dynamic = load(fullfile(save_dir, 'aggregates.mat'));
                        s_static  = load(fullfile(save_dir, 'aggregates_static.mat'));

                        % Find number of projection years
                        Tss = length(s_dynamic.kpr_total);

                        % Find number of entries to be trimmed or padded
                        trim_or_pad = upshift + Tss - n_csv;
                        n_trim =  max(trim_or_pad, 0);
                        n_pad  = -min(trim_or_pad, 0);

                        for i = 1:n_aggs
                            
                            agg_num = agg_codes{i,1};
                            agg_str = agg_codes{i,2};

                            agg_dynamic = s_dynamic.([agg_str,'_total' ]);
                            agg_static  = s_static .([agg_str,'_static']);

                            agg_series  = [ ones(upshift, 2); ...
                                            [agg_dynamic(1:end-n_trim)', agg_static(1:end-n_trim)']; ...
                                            ones(n_pad,   2) ];

                            csvwrite(fullfile(csv_dir, sprintf('%d-%d.csv', ID, agg_num)), ...
                                [years_csv, agg_series])

                        end
                        
                    end
                    
                end
            end
        end
        
    end
end

fprintf('\nResults successfully packaged into csv files:\n\t%s\n', csv_dir)

end