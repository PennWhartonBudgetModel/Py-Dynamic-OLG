%%
% Generate an inverse interpolant for the steady state solver using a grid of solutions as produced by generate_ss_sets.m and solve_ss_batch.m.
% 
%%


function [] = generate_ss_inverter()

addpath('..')


% Identify load directory
sets_dir = 'Sets';

% Identify files containing deep parameter and elasticity sets
param_files = dir(fullfile(sets_dir, 'batch*.csv'  ));
elas_files  = dir(fullfile(sets_dir, 'results*.csv'));

nbatches = length(param_files);


% Initialize arrays for sets
param_sets = cell(nbatches, 1);
elas_sets  = cell(nbatches, 1);

% Extract deep parameter and elasticity sets from files
for i = 1:nbatches
    param_sets{i} = csvread(fullfile(sets_dir, param_files(i).name));
    elas_sets {i} = csvread(fullfile(sets_dir, elas_files (i).name));
end

param_sets = cell2mat(param_sets);
elas_sets  = cell2mat(elas_sets );


% Save deep parameter and elasticity sets
save('invert_ss.mat', 'param_sets', 'elas_sets')

% Delete load directory
rmdir(sets_dir, 's')


% Discard elasticity sets with Inf or NaN values
keep = all(isfinite(elas_sets), 2);

param_sets = param_sets(keep, :);
elas_sets  = elas_sets (keep, :);

% Extract elasticities from sets
K_vec = elas_sets(:,1);
L_vec = elas_sets(:,2);
S_vec = elas_sets(:,3);

% Create inverse interpolant
method = 'nearest';     % ('nearest', 'linear', or 'natural')

beta_interp  = scatteredInterpolant(K_vec, L_vec, S_vec, param_sets(:,1), method);
gamma_interp = scatteredInterpolant(K_vec, L_vec, S_vec, param_sets(:,2), method);
sigma_interp = scatteredInterpolant(K_vec, L_vec, S_vec, param_sets(:,3), method);

invert_ss = @(elas) [beta_interp( elas(:,1), elas(:,2), elas(:,3)), ...
                     gamma_interp(elas(:,1), elas(:,2), elas(:,3)), ...
                     sigma_interp(elas(:,1), elas(:,2), elas(:,3))]; %#ok<NASGU>


% Save inverse interpolant
save('invert_ss.mat', 'invert_ss', '-append')


end