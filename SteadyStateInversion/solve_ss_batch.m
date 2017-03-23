%%
% Solve a batch of steady states as defined by generate_ss_sets.m.
% 
% Arguments:
% 
%   ibatch
%       Batch index.
% 
%%


function [] = solve_ss_batch(ibatch)

addpath('..')

% Set default batch index if none provided
if ~exist('ibatch', 'var')
    ibatch = 63;
end


% Identify save directory
sets_dir = 'Sets';


% Load batch from csv
batch = csvread(fullfile(sets_dir, sprintf('batch%05d.csv', ibatch)));

% Initialize array of results
results = zeros(size(batch));

parfor j = 1:size(batch,1)
    fprintf('Solving set %2d in batch %5d.\n', j, ibatch)
    [~, results(j,:)] = solve_ss(batch(j,:), true);
end


% Save results as csv
csvwrite(fullfile(sets_dir, sprintf('results%05d.csv', ibatch)), results)


end