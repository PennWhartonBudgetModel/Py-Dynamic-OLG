%%
% Check convergence of transition path runs.
% 
%%


function [convergence] = check_convergence()

% Generate baseline definitions
basedefs = {}; %#ok<*AGROW>
for ind = 1:16
    basedefs = [basedefs, {get_basedef(ind)}];
end

% Generate counterfactual definitions
counterdefs = {struct()}; %#ok<*AGROW>
for plan = {'trump', 'clinton', 'ryan'}
    for gcut = [+0.10, +0.05, +0.00, -0.05]
        counterdefs = [counterdefs, {struct('plan', plan{1}, 'gcut', gcut)}];
    end
end

% Initialize array of convergence information
convergence = cell(0,5); %#ok<*AGROW>

% Get convergence information from save directories
for basedef = basedefs
    for counterdef = counterdefs
        for economy = {'open', 'closed'}
            
            % Find save directory
            [save_dir, basedef_tag, counterdef_tag] = dirFinder.save(economy{1}, basedef{1}, counterdef{1});
            
            % Identify iterations log file
            log = fullfile(save_dir, 'iterations.csv');
            
            % Extract and store information on last iteration
            if exist(log, 'file')
                iterations = csvread(log);
                lastiter = num2cell(iterations(end,:));
            else
                lastiter = {NaN, NaN};
            end
            convergence = [convergence; {basedef_tag, counterdef_tag, economy{1}}, lastiter];
            
        end
    end
end

% Sort runs by number of iterations
[~, sortinds] = sort([convergence{:,end-1}], 'descend');
convergence = convergence(sortinds,:);

end