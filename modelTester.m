%%
% Dynamic model tester.
% 
%%


classdef modelTester

properties (Constant)
    
    test_params = struct(...
        'beta'              , 1.000000000000000     , ...
        'gamma'             , 0.900000000000000     , ...
        'sigma'             , 3.098386676965934     , ...
        'modelunit_dollar'  , 8.308176919983091e-05 , ...
        'economy'           , 'steady'              , ...
        'bequest_phi_1'     , 0                     , ...
        'taxplan'           , 'trumpB'              );
    
end
    
methods (Static)

    % Test steady state solution and elasticities
    function [] = steady()
        scenario = Scenario(modelTester.test_params).currentPolicy.steady;
        save_dir = dynamicSolver.solve(scenario);
        setnames = {'market', 'dynamics', 'paramsTargets'};
        test_output(save_dir, setnames);
    end
    
    
    % Test open economy baseline solution, dynamic aggregates, and static aggregates
    function [] = open_base()
        scenario = Scenario(modelTester.test_params).currentPolicy.open;
        save_dir = dynamicSolver.solve(scenario);
        setnames = {'market', 'dynamics'};
        test_output(save_dir, setnames);
    end
    
    
    % Test open economy counterfactual solution, dynamic aggregates, and static aggregates
    function [] = open_counter()
        scenario = Scenario(modelTester.test_params).open;
        save_dir = dynamicSolver.solve(scenario);
        setnames = {'market', 'dynamics', 'statics'};
        test_output(save_dir, setnames);
    end
    
    
    % Test closed economy baseline solution, dynamic aggregates, and static aggregates
    function [] = closed_base()
        scenario = Scenario(modelTester.test_params).currentPolicy.closed;
        save_dir = dynamicSolver.solve(scenario);
        setnames = {'market', 'dynamics'};
        test_output(save_dir, setnames);
    end
    
    
    % Test closed economy counterfactual solution, dynamic aggregates, and static aggregates
    function [] = closed_counter()
        scenario = Scenario(modelTester.test_params).closed;
        save_dir = dynamicSolver.solve(scenario);
        setnames = {'market', 'dynamics', 'statics'};
        test_output(save_dir, setnames);
    end
    
    
    %%
    %  Testing of calibrations
    function [] = calibrate_dollar( )
        % rem: gridpoint needs to contain beta,gamma,sigma
        %  so, we can just pass test_params
        modelCalibrator.calibrate_dollar( modelTester.test_params );
    end % calibrate_dollar

    
end % methods 

end % class modelTester

    
% Test solver output against target values
function [] = test_output(save_dir, setnames)

    % Get test name from name of calling method
    callstack = dbstack();
    testname = regexp(callstack(2).name, '(?<=^modelTester\.).*$', 'match', 'once');

    % Load target values
    targetfile = fullfile(fileparts(mfilename('fullpath')), 'modelTester.mat');
    s = load(targetfile); target = s.target; clear('s')

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

            elseif any(size(outputset.(valuename)) ~= size(targetset.(valuename)))
                
                % Flag for size mismatch
                flag('Size mismatch');
                
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
            save(targetfile, 'target')
            fprintf('\tTarget updated.\n')
        else
            fprintf('\tTarget retained.\n')
        end

    end
    fprintf('\n')

end % function test_output

