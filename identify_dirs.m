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

    % Identify active Git commit
    %   %cd = commit date
    %   %an = author name
    [~, commit_id] = system('git log -1 --format=%cd-%an --date=format:%Y-%m-%d-%H-%M-%S');
    % [~, commit_id] = system('git log -1 --format=%cd-%an-%h --date=short');   % For older versions of Git that don't support custom date formats
    
    % Strip off trailing newline character
    commit_id(end) = '';
    
    % Designate permanent root save directory
    save_root = fullfile(dev_dir, '..', 'Output', commit_id);

else
    
    % Designate temporary root save directory
    save_root = fullfile(dev_dir, 'Testing');
    
end

% Identify save directory
save_dir = fullfile(save_root, sprintf('beta=%0.3f_gamma=%0.3f_sigma=%05.2f', beta, gamma, sigma), sub_dir);


end