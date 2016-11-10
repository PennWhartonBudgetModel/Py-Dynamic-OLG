%%
% Dynamic model test suite.
% 
%%


classdef modelTester

properties (Constant)
    deep_params = inddeep_to_params(6);
end

methods (Static)
    
    % Test steady state solution
    function [] = steady()
        
        [~, save_dir] = solve_ss(modelTester.deep_params);
        s_solution = load(fullfile(save_dir, 'solution.mat'));
        
        s_target = load('modelTester.mat', 'steady');
        
        fprintf('\n')
        isdiff = compare_values(s_solution, s_target.steady.Solution);
        
        if ~isdiff
            fprintf('No differences found.\n')
        end
        fprintf('\n')
        
    end
    
    
    % Test open economy baseline dynamic and static aggregates
    function [] = open_base()
        
        save_dir = solve_open(modelTester.deep_params, 'base');
        modelTester.transition(save_dir, 'open_base')
        
    end
    
    
    % Test open economy counterfactual dynamic and static aggregates
    function [] = open_plan()
        
        save_dir = solve_open(modelTester.deep_params, 'ryan');
        modelTester.transition(save_dir, 'open_plan')
        
    end
    
    
    % Test closed economy baseline dynamic and static aggregates
    function [] = closed_base()
        
        save_dir = solve_closed(modelTester.deep_params, 'base', +0.00);
        modelTester.transition(save_dir, 'closed_base')
        
    end
    
    
    % Test closed economy counterfactual dynamic and static aggregates
    function [] = closed_plan()
        
        save_dir = solve_closed(modelTester.deep_params, 'ryan', +0.10);
        modelTester.transition(save_dir, 'closed_plan')
        
    end
    
end


methods (Static, Access = private)
    
    % Test transition path dynamic and static aggregates
    function [] = transition(save_dir, testname)
        
        s_dynamic = load(fullfile(save_dir, 'aggregates.mat'       ));
        s_static  = load(fullfile(save_dir, 'aggregates_static.mat'));
        
        s_target = load('modelTester.mat', testname);
        
        fprintf('\n')
        isdiff_dynamic = compare_values(s_dynamic, s_target.(testname).Dynamic);
        isdiff_static  = compare_values(s_static , s_target.(testname).Static );
        
        if (~isdiff_dynamic && ~isdiff_static)
            fprintf('No differences found.\n')
        end
        fprintf('\n')
        
    end
    
end

end


% Compare values and report differences
function [isdiff] = compare_values(test, target)
    
    isdiff = false;
    valuestrs = fieldnames(target);
    
    for i = 1:length(valuestrs)
        
        valuestr = valuestrs{i};
        
        if ~isfield(test, valuestr)
            fprintf('Difference found: %s does not exist.\n', valuestr)
            isdiff = true;
        else
            
            diff = test.(valuestr)(:) - target.(valuestr)(:);
            
            if any(diff)
                [maxdiff, ind] = max(diff);
                dev = maxdiff * 2 / (test.(valuestr)(ind) + target.(valuestr)(ind));
                fprintf('Difference found: %s exhibits %0.2f%% deviation.\n', valuestr, abs(dev*100))
                isdiff = true;
            end
            
        end
        
    end
end

