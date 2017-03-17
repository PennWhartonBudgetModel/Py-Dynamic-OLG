%%
% Dynamic model tester.
% 
%%


classdef modelTester

properties (Constant)
    
    % Define baseline for all tests
    basedef = get_basedef(6);
    
    % Define counterfactual for counterfactual tests
    counterdef = struct('taxplan'       , 'ryan', ...
                        'gcut'          , +0.10 , ...
                        'legal_scale'   , 1.5   , ...
                        'prem_legal'    , 1.117 , ...
                        'amnesty'       , 0.05  , ...
                        'deportation'   , 0.05  );
    
end

methods (Static)
    
    
    % Test steady state solution and elasticities
    function [] = steady()
        save_dir = dynamicSolver.steady(modelTester.basedef);
        setnames = {'market', 'elasticities'};
        test_output(save_dir, setnames);
    end
    
    
    % Test open economy baseline solution, dynamic aggregates, and static aggregates
    function [] = open_base()
        save_dir = dynamicSolver.open(modelTester.basedef);
        setnames = {'market', 'dynamics'};
        test_output(save_dir, setnames);
    end
    
    
    % Test open economy counterfactual solution, dynamic aggregates, and static aggregates
    function [] = open_counter()
        save_dir = dynamicSolver.open(modelTester.basedef, modelTester.counterdef);
        setnames = {'market', 'dynamics', 'statics'};
        test_output(save_dir, setnames);
    end
    
    
    % Test closed economy baseline solution, dynamic aggregates, and static aggregates
    function [] = closed_base()
        save_dir = dynamicSolver.closed(modelTester.basedef);
        setnames = {'market', 'dynamics'};
        test_output(save_dir, setnames);
    end
    
    
    % Test closed economy counterfactual solution, dynamic aggregates, and static aggregates
    function [] = closed_counter()
        save_dir = dynamicSolver.closed(modelTester.basedef, modelTester.counterdef);
        setnames = {'market', 'dynamics', 'statics'};
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

% Define function to flag issues
function flag(str)
    fprintf('\t%-15s%-20s%s\n', setname, valuename, str)
    ismatch = false;
end

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
        
        if ~isfield(outputset, valuename)
            
            % Flag missing value
            flag('Not found');
        
        elseif any(isnan(outputset.(valuename)(:)))
            
            % Flag NaN value
            flag('NaN value');
            
        else
            
            % Identify value deviation
            delta = outputset.(valuename)(:) - targetset.(valuename)(:);
            if any(delta)
                pdev = abs(nanmean(delta*2 ./ (outputset.(valuename)(:) + targetset.(valuename)(:))))*100;
                if pdev < 0.01
                    flag(sprintf('Numerical deviation'));
                else
                    flag(sprintf('%06.2f%% deviation', pdev));
                end
            end
            
        end
        
    end
    
    % Iterate over output values
    outputvaluenames = fieldnames(outputset);
    
    for j = 1:length(outputvaluenames)
        
        valuename = outputvaluenames{j};
        
        % Identify new value
        if ~isfield(targetset, valuename)
            flag('New');
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

