%%
% Helper function to identify working directories for solvers and functions that operate on solver output.
% 
% To identify just the parameters directory, use:
%   >> param_dir = identify_dirs('');
% 
%%


function [param_dir, save_dir] = identify_dirs(solver, beta, gamma, sigma, plan, gcut)

% Identify base and parameter directories
dev_dir  = fileparts(mfilename('fullpath'));
param_dir = fullfile(dev_dir, 'Parameters');

switch solver

    % Steady state
    case 'ss',     sub_dir = 'ss';
    
    % Open economy transition path
    case 'open',   sub_dir = fullfile(sprintf('trans_plan=%s', plan), 'open');
    
    % Closed economy transition path
    case 'closed', sub_dir = fullfile(sprintf('trans_plan=%s', plan), sprintf('closed_gcut=%+0.2f', gcut));
    
    otherwise, return
    
end


% Check for uncommitted changes
[~, uncommitted] = system('git status -s');

if isempty(uncommitted)
    
    % Get date, author, and abbreviated hash for active Git commit
    [~, commit_date  ] = system('git log -1 --format=%cd --date=iso');
    [~, commit_author] = system('git log -1 --format=%ae'           );
    [~, commit_hash  ] = system('git log -1 --format=%h'            );
    
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
    
    % Designate permanent root save directory
    save_root = fullfile(dev_dir, '..', 'Output', commit_id);
    
else
    
    % Designate temporary root save directory
    save_root = fullfile(dev_dir, 'Testing');
    
end

% Identify save directory
save_dir = fullfile(save_root, sprintf('beta=%0.3f_gamma=%0.3f_sigma=%05.2f', beta, gamma, sigma), sub_dir);


end