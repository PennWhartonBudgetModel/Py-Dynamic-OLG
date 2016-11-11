%%
% Dynamic model test suite.
% 
%%


classdef modelTester

properties (Constant)
    
    % Define test runs
    deep_params = inddeep_to_params(6);
    
    steady_run      = @() solve_ss    (modelTester.deep_params               );
    open_base_run   = @() solve_open  (modelTester.deep_params, 'base'       );
    open_plan_run   = @() solve_open  (modelTester.deep_params, 'ryan'       );
    closed_base_run = @() solve_closed(modelTester.deep_params, 'base', +0.00);
    closed_plan_run = @() solve_closed(modelTester.deep_params, 'ryan', +0.10);
    
end

methods (Static)
    
    % Test steady state solution
    function [] = steady(update_target)
        
        [~, save_dir] = modelTester.steady_run();
        s_solution = load(fullfile(save_dir, 'solution.mat'));
        
        s_target = load('modelTester.mat');
        target = s_target.target;
        
        if (~exist('update_target', 'var') || ~update_target)
            
            fprintf('[Test results]\n')
            isdiff = compare_values(s_solution, target.steady.Solution);
            
            if ~isdiff
                fprintf('\tNo differences from target identified.\n')
            end
            fprintf('\n')
            
        else
            
            % Update target values
            target.steady.Solution = s_solution;
            save('modelTester.mat', 'target')
            fprintf('[Test target updated]\n\n')
            
        end
    end
    
    
    % Test open economy baseline dynamic and static aggregates
    function [] = open_base(update_target)
        if ~exist('update_target', 'var'), update_target = false; end
        save_dir = modelTester.open_base_run();
        modelTester.transition('open_base', save_dir, update_target)
    end
    
    
    % Test open economy counterfactual dynamic and static aggregates
    function [] = open_plan(update_target)
        if ~exist('update_target', 'var'), update_target = false; end
        save_dir = modelTester.open_plan_run();
        modelTester.transition('open_plan', save_dir, update_target)
    end
    
    
    % Test closed economy baseline dynamic and static aggregates
    function [] = closed_base(update_target)
        if ~exist('update_target', 'var'), update_target = false; end
        save_dir = modelTester.closed_base_run();
        modelTester.transition('closed_base', save_dir, update_target)
    end
    
    
    % Test closed economy counterfactual dynamic and static aggregates
    function [] = closed_plan(update_target)
        if ~exist('update_target', 'var'), update_target = false; end
        save_dir = modelTester.closed_plan_run();
        modelTester.transition('closed_plan', save_dir, update_target)
    end
    
    
    % Perform all tests
    function [] = test_all()
        update_target = false;
        modelTester.run_all(update_target);
    end
    
    
    % Update all test targets
    function [] = update_all()
        update_target = true;
        modelTester.run_all(update_target);
    end
    
end


methods (Static, Access = private)
    
    % Test transition path dynamic and static aggregates
    function [] = transition(testname, save_dir, update_target)
        
        s_dynamic = load(fullfile(save_dir, 'aggregates.mat'       ));
        s_static  = load(fullfile(save_dir, 'aggregates_static.mat'));
        
        s_target = load('modelTester.mat');
        target = s_target.target;
        
        if ~update_target
            
            fprintf('[Test results]\n')
            isdiff_dynamic = compare_values(s_dynamic, target.(testname).Dynamic);
            isdiff_static  = compare_values(s_static , target.(testname).Static );

            if (~isdiff_dynamic && ~isdiff_static)
                fprintf('\tNo differences from target identified.\n')
            end
            fprintf('\n')
        
        else
            
            % Update target values
            target.(testname).Dynamic = s_dynamic;
            target.(testname).Static  = s_static ;
            save('modelTester.mat', 'target')
            fprintf('[Test target updated]\n\n')
            
        end
        
    end
    
    
    % Perform all tests or update all test targets
    function [] = run_all(update_target)
        modelTester.steady     (update_target);
        modelTester.open_base  (update_target);
        modelTester.open_plan  (update_target);
        modelTester.closed_base(update_target);
        modelTester.closed_plan(update_target);
    end
    
end

end


% Compare values and report differences
function [isdiff] = compare_values(output, target)
    
    isdiff = false;
    valuestrs = fieldnames(target);
    
    for i = 1:length(valuestrs)
        
        valuestr = valuestrs{i};
        
        if ~isfield(output, valuestr)
            fprintf('\t%-25sNot found.\n', valuestr)
            isdiff = true;
        else
            
            diff = output.(valuestr)(:) - target.(valuestr)(:);
            
            if any(diff)
                [maxdiff, ind] = max(diff);
                dev = maxdiff * 2 / (output.(valuestr)(ind) + target.(valuestr)(ind));
                fprintf('\t%-25s%06.2f%% deviation.\n', valuestr, abs(dev*100))
                isdiff = true;
            end
            
        end
        
    end
end

