%%
% Directory finder.
% 
% Methods:
% 
%   source          Find source code directory.
% 
%   modelroot       Find model root directory.
% 
%   root            Find absolute root directory.
% 
%   param           Find parameters directory.
% 
%   saveroot        Find save root directory.
% 
%   ss              Find steady state save directory.
% 
%   open            Find open economy save directory.
% 
%   closed          Find closed economy save directory.
% 
%   csv             Find csv save directory.
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
        closed_dir = fullfile( dirFinder.saveroot, ...
                               get_deep_params_tag(beta, gamma, sigma), ...
                               get_plan_tag(plan), ...
                               sprintf('closed_gcut=%+0.2f', gcut) );
    end
    
    
    % Find csv save directory
    function [csv_dir] = csv()
        if dirFinder.isproductionready
            csv_dir = fullfile(dirFinder.root, 'charts', 'version2', 'tax_dynamic_scores', dirFinder.get_commit_id);
        else
            csv_dir = fullfile(dirFinder.testout, 'csv');
        end
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
        
        % Check for uncommitted changes
        [~, uncommitted] = system('git status -s');
        
        flag = ( strcmp(stage, 'Production') && isempty(uncommitted) );
        
    end
    
    
    % Get identifier for active Git commit
    function [commit_id] = get_commit_id()
        
        % Get date, author, and abbreviated hash for active Git commit
        [~, commit_date  ] = system('git --no-pager log -1 --format=%cd --date=iso');
        [~, commit_author] = system('git --no-pager log -1 --format=%ae'           );
        [~, commit_hash  ] = system('git --no-pager log -1 --format=%h'            );

        % Strip off trailing characters
        commit_date   = commit_date  (1:end-7);     % (Time zone stripped off; only local time needed for commit identifier)
        commit_author = commit_author(1:end-1);
        commit_hash   = commit_hash  (1:end-1);

        % Extract author username from email address
        commit_author = regexp(commit_author, '.*(?=@)', 'match');
        commit_author = commit_author{1};

        % Convert date format
        commit_date = datestr(commit_date, 'yyyy-mm-dd-HH-MM');

        % Construct identifier for active Git commit
        commit_id = sprintf('%s-%s-%s', commit_date, commit_author, commit_hash);
        
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

