%%
% Helper function to convert index of deep parameter set as ordered in transition_sets.mat into a 1 x 3 vector of deep parameters [beta, gamma, sigma].
% 
%%


function [deep_params] = inddeep_to_params(inddeep)

% Load deep parameter set
param_dir = identify_dirs('');
s = load(fullfile(param_dir, 'ss_inverses.mat'));
deep_params = s.inverses(inddeep,:);

end