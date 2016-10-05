% Pete | 2016-09-28
% 
% Check convergence for closed economy transition path runs.
% 
% 

function [lastiters] = check_closed_convergence()

lastiters = cell(0,5); %#ok<*AGROW>

for inddeep = 1:16
    for plan = {'base', 'trump', 'clinton', 'ryan'}
        for gcut = [0.10, 0.05, 0.00, -0.05]
            
            % Find save directory
            deep_params = inddeep_to_params(inddeep);
            [~, save_dir] = identify_dirs('closed', deep_params(1), deep_params(2), deep_params(3), plan{1}, gcut);
            
            % Extract subdirectories
            [up_dir, gcut_sub, dec_str] = fileparts(save_dir);
            gcut_sub = [gcut_sub, dec_str];
            
            [up_dir, plan_sub] = fileparts(up_dir);
            
            [~, param_sub, dec_str] = fileparts(up_dir);
            param_sub = [param_sub, dec_str];
            
            % Identify iterations log file
            logfile = fullfile(save_dir, 'iterations.txt');
            
            if exist(logfile, 'file')
                
                % Get last line of iterations log file
                fid = fopen(logfile);
                while ~feof(fid)
                    lastline = fgetl(fid);
                end
                fclose(fid);
                
                % Extract and store information on last iteration
                lastiters = [lastiters; {param_sub, plan_sub, gcut_sub}, num2cell(sscanf(lastline, '  %2d  --  %f', [1,2]))];
                
            else
                lastiters = [lastiters; {param_sub, plan_sub, gcut_sub, NaN, NaN}];
            end
            
        end
    end
end

% Sort by final convergence error term
[~, sortinds] = sort([lastiters{:,end}], 'descend');
lastiters = lastiters(sortinds,:);

end