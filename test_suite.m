%%
% Note that these tests do not check for differences in results that do not affect solver iterations.
% 

%%
deep_params = inddeep_to_params(6);
solve_ss(deep_params, false)
[isdiff, ~] = system(['fc ', 'test_suite_ss_iterations.txt', ...
                      ' ',   fullfile(dirFinder.ss(deep_params), 'iterations.txt')]);
switch isdiff, case 0, fprintf('No differences'), case 1, fprintf('Differences'), otherwise, fprintf('?'), end
fprintf(' found in steady state solver iterations.\n')


%%
deep_params = inddeep_to_params(6);
solve_open(deep_params, 'base', false)
[isdiff, ~] = system(['fc ', 'test_suite_open_base_iterations.txt', ...
                      ' ',   fullfile(dirFinder.open(deep_params, 'base'), 'iterations.txt')]);
switch isdiff, case 0, fprintf('No differences'), case 1, fprintf('Differences'), otherwise, fprintf('?'), end
fprintf(' found in open economy baseline solver iterations.\n')


%%
deep_params = inddeep_to_params(6);
solve_closed(deep_params, 'base', +0.00, false)
[isdiff, ~] = system(['fc ', 'test_suite_closed_base_iterations.txt', ...
                      ' ',   fullfile(dirFinder.closed(deep_params, 'base', +0.00), 'iterations.txt')]);
switch isdiff, case 0, fprintf('No differences'), case 1, fprintf('Differences'), otherwise, fprintf('?'), end
fprintf(' found in closed economy baseline solver iterations.\n')


%%
deep_params = inddeep_to_params(6);
solve_open(deep_params, 'clinton', false)
[isdiff, ~] = system(['fc ', 'test_suite_open_clinton_iterations.txt', ...
                      ' ',   fullfile(dirFinder.open(deep_params, 'clinton'), 'iterations.txt')]);
switch isdiff, case 0, fprintf('No differences'), case 1, fprintf('Differences'), otherwise, fprintf('?'), end
fprintf(' found in open economy counterfactual solver iterations.\n')


%%
deep_params = inddeep_to_params(6);
solve_closed(deep_params, 'clinton', +0.10, false)
[isdiff, ~] = system(['fc ', 'test_suite_closed_clinton_gcut=+0.10_iterations.txt', ...
                      ' ',   fullfile(dirFinder.closed(deep_params, 'clinton', +0.10), 'iterations.txt')]);
switch isdiff, case 0, fprintf('No differences'), case 1, fprintf('Differences'), otherwise, fprintf('?'), end
fprintf(' found in closed economy counterfactual solver iterations.\n')
