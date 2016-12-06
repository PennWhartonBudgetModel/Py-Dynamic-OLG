%%
% Dynamic model tester.
% 
%%


classdef modelTester

properties (Constant)
    
    % Specify baseline definition
    basedef = get_basedef(6);
    
end

methods (Static)
    
    
    % Test steady state solution and elasticities
    function [] = steady()
        save_dir = dynamicSolver.steady(modelTester.basedef);
        setnames = {'solution', 'elasticities'};
        test_output(save_dir, setnames);
    end
    
    
    % Test open economy baseline solution, dynamic aggregates, and static aggregates
    function [] = open_base()
        save_dir = dynamicSolver.open(modelTester.basedef);
        setnames = {'solution', 'aggregates'};
        test_output(save_dir, setnames);
    end
    
    
    % Test open economy counterfactual solution, dynamic aggregates, and static aggregates
    function [] = open_counter()
        save_dir = dynamicSolver.open(modelTester.basedef, struct('plan', 'ryan', 'gcut', +0.10));
        setnames = {'solution', 'aggregates', 'aggregates_static'};
        test_output(save_dir, setnames);
    end
    
    
    % Test closed economy baseline solution, dynamic aggregates, and static aggregates
    function [] = closed_base()
        save_dir = dynamicSolver.closed(modelTester.basedef);
        setnames = {'solution', 'aggregates'};
        test_output(save_dir, setnames);
    end
    
    
    % Test closed economy counterfactual solution, dynamic aggregates, and static aggregates
    function [] = closed_counter()
        save_dir = dynamicSolver.closed(modelTester.basedef, struct('plan', 'ryan', 'gcut', +0.10));
        setnames = {'solution', 'aggregates', 'aggregates_static'};
        test_output(save_dir, setnames);
    end
    
    
end

end

    
% Test solver output against target values
function [] = test_output(save_dir, setnames)

% Get test name from name of calling method
callstack = dbstack();
testname = regexp(callstack(2).name, '(?<=^modelTester\.).*$', 'match', 'once');

% Load target values
s = load('modelTester.mat');
target = s.target;

% Initialize match flag
ismatch = true;

fprintf('\n[Test results]\n')
for i = 1:length(setnames)
    
    % Extract output and target values by set
    setname = setnames{i};
    output.(testname).(setname) = load(fullfile(save_dir, sprintf('%s.mat', setname)));
    
    outputset = output.(testname).(setname);
    targetset = target.(testname).(setname);
    
    % Iterate over target values
    targetvaluenames = fieldnames(targetset);
    
    for j = 1:length(targetvaluenames)
        
        valuename = targetvaluenames{j};
        
        % Identify missing value
        if ~isfield(outputset, valuename)
            fprintf('\t%-25sNot found.\n', valuename)
            ismatch = false;
        else
            
            % Identify value deviation
            delta = outputset.(valuename)(:) - targetset.(valuename)(:);
            
            if any(delta)
                dev = nanmean( delta*2 ./ (outputset.(valuename)(:) + targetset.(valuename)(:)) );
                fprintf('\t%-25s%06.2f%% deviation.\n', valuename, abs(dev*100))
                ismatch = false;
            end
            
        end
        
    end
    
    % Iterate over output values
    outputvaluenames = fieldnames(outputset);
    
    for j = 1:length(outputvaluenames)
        
        valuename = outputvaluenames{j};
        
        % Identify new value
        if ~isfield(targetset, valuename)
            fprintf('\t%-25sNew.\n', valuename)
            ismatch = false;
        end
        
    end
    
end

% Check for match
if ismatch
    fprintf('\tTarget matched.\n')
else
    
    % Query user for target update
    if strcmp(input(sprintf('\n\tUpdate test target with new values? Y/[N]: '), 's'), 'Y')
        target.(testname) = output.(testname);
        save('modelTester.mat', 'target')
        fprintf('\tTarget updated.\n')
    else
        fprintf('\tTarget retained.\n')
    end
    
end
fprintf('\n')

end

