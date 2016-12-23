%%
% Construct baseline definition structure based on index of deep parameter set in ss_inverses.mat.
% 
%%


function [basedef, deep_params] = get_basedef(ind)

% Load deep parameter set
s = load(fullfile(dirFinder.param, 'ss_inverses.mat'));
deep_params = s.inverses(ind,:);

% Construct baseline definition structure
basedef = struct('beta' , deep_params(1), ...
                 'gamma', deep_params(2), ...
                 'sigma', deep_params(3));

end