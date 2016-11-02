%%
% Directory finder.
% 
% Methods:
% 
%   source()
%       Find source code directory.
% 
%   modelroot()
%       Find model root directory.
% 
%   root()
%       Find absolute root directory.
% 
%   param()
%       Find parameters directory.
% 
%   saveroot()
%       Find save root directory.
% 
%   ss(beta, gamma, sigma)
%       Find steady state save directory.
% 
%   open(beta, gamma, sigma, plan)
%       Find open economy save directory.
% 
%   closed(beta, gamma, sigma, plan, gcut)
%       Find closed economy save directory.
% 
%   csv()
%       Find csv save directory.
% 
%%


classdef dirFinder
    
methods (Static)
    
    % Find source code directory
    function [source_dir] = source()
        source_dir = fileparts(mfilename('fullpath'));
    end
    
    
    % Find model root directory
    function [modelroot_dir] = modelroot()
        modelroot_dir = fileparts(fileparts(dirFinder.source));
    end
    
    
    % Find absolute root directory
    function [root_dir] = root()
        root_dir = fileparts(dirFinder.modelroot);
    end
    
    
    % Find parameters directory
    function [param_dir] = param()
        param_dir = fullfile(dirFinder.source, 'Parameters');
    end
    
    
    % Find save root directory
    function [saveroot_dir] = saveroot()
        if dirFinder.isproductionready
            saveroot_dir = fullfile(dirFinder.modelroot, 'Output', dirFinder.get_commit_id);
        else
            saveroot_dir = dirFinder.testout;
        end
    end
    
    
    % Find steady state save directory
    function [ss_dir] = ss(beta, gamma, sigma)
        ss_dir     = fullfile( dirFinder.saveroot, ...
                               get_deep_params_tag(beta, gamma, sigma), ...
                               'ss' );
    end
    
    
    % Find open economy save directory
    function [open_dir] = open(beta, gamma, sigma, plan)
        open_dir   = fullfile( dirFinder.saveroot, ...
                               get_deep_params_tag(beta, gamma, sigma), ...
                               get_plan_tag(plan), ...
                               'open' );
    end
    
    
    % Find closed economy save directory
    function [closed_dir] = closed(beta, gamma, sigma, plan, gcut)
        if (gcut == 0), gcut = 0; end   % Enforces positive zero in the case of a negative zero
        closed_dir = fullfile( dirFinder.saveroot, ...
                               get_deep_params_tag(beta, gamma, sigma), ...
                               get_plan_tag(plan), ...
                               sprintf('closed_gcut=%+0.2f', gcut) );
    end
    
    
    % Find csv save directory
    function [csv_dir] = csv()
        timestamp = datestr(now, 'yyyy-mm-dd-HH-MM');
        if dirFinder.isproductionready
            csvroot = fullfile(dirFinder.root, 'charts', 'version2', 'tax_dynamic_scores', dirFinder.get_commit_id);
        else
            csvroot = fullfile(dirFinder.testout, 'csv');
        end
        csv_dir = fullfile(csvroot, timestamp);
    end
    
    
    % Get identifier for active Git commit
    function [commit_id] = get_commit_id()
        
        % Get commit date (%cd), author email address (%ae), and abbreviated hash (%h)
        [~, commit_log] = system('git --no-pager log -1 --format=%cd,%ae,,%h --date=iso');
        
        % Extract commit date and reformat
        commit_date   = regexp(commit_log, '.*? .*?(?= )', 'match', 'once');
        commit_date   = datestr(commit_date, 'yyyy-mm-dd-HH-MM');
        
        % Extract author username from email address
        commit_author = regexp(commit_log, '(?<=,).*?(?=@)', 'match', 'once');
        
        % Extract abbreviated hash
        commit_hash   = regexp(commit_log, '(?<=,,).*?(?=\n)', 'match', 'once');
        
        % Construct commit identifier
        commit_id = sprintf('%s-%s-%s', commit_date, commit_author, commit_hash);
        
    end
    
end


methods (Static, Access = private)
    
    % Find testing output directory
    function [testout_dir] = testout()
        testout_dir = fullfile(dirFinder.source, 'Testing');
    end
    
    
    % Identify production run by Production stage and absence of uncommitted changes
    function [flag] = isproductionready()
        
        % Get source code stage
        [~, stage] = fileparts(fileparts(dirFinder.source));
        
        if ~strcmp(stage, 'Production')
            flag = false;
        else
            
            % Check for uncommitted changes
            % (Safeguards against any unintentional changes made in Production directory after checkout)
            [~, uncommitted] = system('git status -s');
            
            flag = isempty(uncommitted);
            
        end
        
    end
    
end

end


% Get string tag corresponding to deep parameters
function [deep_params_tag] = get_deep_params_tag(beta, gamma, sigma)
    deep_params_tag = sprintf('beta=%0.3f_gamma=%0.3f_sigma=%05.2f', beta, gamma, sigma);
end


% Get string tag corresponding to plan
function [plan_tag] = get_plan_tag(plan)
    plan_tag = sprintf('trans_plan=%s', plan);
end

