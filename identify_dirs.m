% Pete | 2016-09-14
% 
% Helper function to identify working directories for steady state and transition path solvers.
% 
% To identify just the parameters directory, use:
%   >> param_dir = identify_dirs('');
% 
% 


function [param_dir, save_dir] = identify_dirs(solver, beta, gamma, sigma, plan, gcut)

% Identify base and parameter directories
main_dir  = fileparts(mfilename('fullpath'));
param_dir = fullfile(main_dir, 'Parameters');

switch solver

    % Steady state
    case 'ss',     sub_dir = 'ss';
    
    % Open economy transition path
    case 'open',   sub_dir = fullfile(sprintf('trans_plan=%s', plan), 'open');
    
    % Closed economy transition path
    case 'closed', sub_dir = fullfile(sprintf('trans_plan=%s', plan), sprintf('closed_gcut=%+0.2f', gcut));
    
    otherwise, return
    
end

% Identify save directory
save_dir = fullfile(main_dir, 'Results', sprintf('beta=%0.3f_gamma=%0.3f_sigma=%05.2f', beta, gamma, sigma), sub_dir);

end