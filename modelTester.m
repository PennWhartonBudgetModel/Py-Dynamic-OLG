%%
% Dynamic model tester.
% 
%%


classdef modelTester

methods (Static)
    %%
    %  Definitions of Scenarios
    function scenario = getBaselineDef()
        scenario = Scenario();
        scenario.beta               = 0.9900 ;
        scenario.gamma              = 0.6672;
        scenario.sigma              = 2.5142;
        scenario.modelunit_dollar   = 4.5408e-05;
    end % getBaseDef
    
    function scenario = getCounterDef()
        basedef = modelTester.getBaselineDef();
        scenario = basedef.Clone();
        scenario.baselineScenario = basedef;
        scenario.taxplan = 'trumpB';
    end % getBaseDef

    % Test steady state solution and elasticities
    function [] = steady()
        scenario = modelTester.getBaselineDef();
        scenario.economy = 'steady';
        save_dir = dynamicSolver.solve(scenario);
        setnames = {'market', 'dynamics', 'elasticities'};
        test_output(save_dir, setnames);
        GiniTable = momentsGenerator(scenario);
    end
    
    
    % Test open economy baseline solution, dynamic aggregates, and static aggregates
    function [] = open_base()
        scenario = modelTester.getBaselineDef();
        scenario.economy = 'open';
        save_dir = dynamicSolver.solve(scenario);
        setnames = {'market', 'dynamics'};
        test_output(save_dir, setnames);
    end
    
    
    % Test open economy counterfactual solution, dynamic aggregates, and static aggregates
    function [] = open_counter()
        scenario = modelTester.getCounterDef();
        scenario.economy = 'open';
        save_dir = dynamicSolver.solve(scenario);
        setnames = {'market', 'dynamics', 'statics'};
        test_output(save_dir, setnames);
    end
    
    
    % Test closed economy baseline solution, dynamic aggregates, and static aggregates
    function [] = closed_base()
        scenario = modelTester.getBaselineDef();
        scenario.economy = 'closed';
        save_dir = dynamicSolver.solve(scenario);setnames = {'market', 'dynamics'};
        test_output(save_dir, setnames);
    end
    
    
    % Test closed economy counterfactual solution, dynamic aggregates, and static aggregates
    function [] = closed_counter()
        scenario = modelTester.getCounterDef();
        scenario.economy = 'closed';
        save_dir = dynamicSolver.solve(scenario);setnames = {'market', 'dynamics', 'statics'};
        test_output(save_dir, setnames);
    end
    
    
    %%
    %  Testing of calibrations
    function [] = calibrate_dollar( params )
        modelCalibrator.calibrate_dollar( params )
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

