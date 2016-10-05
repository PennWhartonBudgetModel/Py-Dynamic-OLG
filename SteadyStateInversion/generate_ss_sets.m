% Pete | 2016-09-22
% 
% Generate a grid of deep parameter sets over which the steady state will be solved.
% 
% Arguments:
% 
%   gridsize
%       Grid size as an integer.  The same value is used for all deep parameter dimensions.
% 
%   batchsize
%       Number of sets per batch to be solved together by solve_ss_batch.m.
% 
% 


function [] = generate_ss_sets(gridsize, batchsize)

% Set default grid size if none provided
if ~exist('gridsize', 'var')
    gridsize = 10;
end

% Set default batch size if none provided
if ~exist('batchsize', 'var')
    batchsize = 4;
end


% Specify deep parameter bounds
beta_min  = 0.970;  beta_max  = 1.110;
gamma_min = 0.150;  gamma_max = 0.900;
sigma_min = 01.50;  sigma_max = 22.00;


% Generate deep parameter sets
beta_vec  = linspace(beta_min,  beta_max,  gridsize);
gamma_vec = linspace(gamma_min, gamma_max, gridsize);
sigma_vec = linspace(sigma_min, sigma_max, gridsize);

[beta_, gamma_, sigma_] = ndgrid(beta_vec, gamma_vec, sigma_vec);
param_sets = [beta_(:), gamma_(:), sigma_(:)];


% Identify save directory and clear or create
sets_dir = 'Sets';
if exist(sets_dir, 'dir')
    rmdir(sets_dir, 's')
end
mkdir(sets_dir)

% Find number of batches
nbatches = ceil(gridsize^3/batchsize);

for i = 1:nbatches
    
    % Pull batch from sets
    if (i ~= nbatches)
        batch = param_sets( (i-1)*batchsize + (1:batchsize) , :);
    else
        batch = param_sets( (i-1)*batchsize + 1 : end, :);
    end
    
    % Save batch as csv
    csvwrite(fullfile(sets_dir, sprintf('batch%05d.csv', i)), batch)
    
end


end