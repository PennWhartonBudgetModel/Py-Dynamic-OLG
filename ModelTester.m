%%
% Dynamic model internal output tester.
%
%%
classdef ModelTester

properties (Constant)
    
    test_params = struct(...
        'economy'               , 'steady'                  ...
    ,      'beta'               , 0.9810000000000          ...
    ,      'gamma'              , 0.7500000000000          ...
    ,      'sigma'              , 1.600000000000000         ...
    ,      'bequest_phi_1'      , 0                         ...
    ,      'modelunit_dollar'   , 3.830000000e-05          ...
    ,   'IsLowReturn'           , false                     ...
    ,   'TransitionFirstYear'   , 2018                      ...
    ,   'OpennessPath'          , 'closed'              ...
    ,   'TransitionLastYear'    , (2018+10)                 ...
    ,   'ClosureYear'           , (2018+5)                  ...
    ... % POLICY PARAMS %
    ,   'TaxCode'               , 'TopRate'           );

end

methods (Static)

    % Test steady state solution and elasticities
    function [] = steady()
        ModelTester.runTest( Scenario(ModelTester.test_params).currentPolicy().steady()     ...
                            ,   {'market', 'dynamics', 'paramsTargets'}     );
    end
    
    
    % Test open economy baseline solution, dynamic aggregates, and static aggregates
    function [] = open_base()
        ModelTester.runTest( Scenario(ModelTester.test_params).currentPolicy().open()       ...
                            ,   {'market', 'dynamics'}                      );
    end
    
    
    % Test open economy counterfactual solution, dynamic aggregates, and static aggregates
    function [] = open_counter()
        ModelTester.runTest( Scenario(ModelTester.test_params).open()                       ...
                            ,   {'market', 'dynamics', 'statics'}           );
    end
    
    
    % Test closed economy baseline solution, dynamic aggregates, and static aggregates
    function [] = closed_base()
        ModelTester.runTest( Scenario(ModelTester.test_params).currentPolicy().closed()     ...
                            ,   {'market', 'dynamics'}                      );
    end
    
    
    % Test closed economy counterfactual solution, dynamic aggregates, and static aggregates
    function [] = closed_counter()
        ModelTester.runTest( Scenario(ModelTester.test_params).closed()                      ...
                            ,   {'market', 'dynamics', 'statics'}           );
    end
    
    
end % methods


methods (Static, Access = private )
    
    function [] = runTest( scenario, testset_names )
        PathFinder.setToTestingMode();
        save_dir = ModelSolver.solve(scenario, '', true);
        test_output(save_dir, testset_names);
        PathFinder.setToDevelopmentMode();
    end
    
end

end

    
% Test solver output against target values
function [] = test_output(save_dir, setnames)

    % Get test name from name of calling method (one level up)
    callstack = dbstack();
    testname = regexp(callstack(3).name, '(?<=^ModelTester\.).*$', 'match', 'once');

    % Load target values
    targetfile = fullfile(fileparts(mfilename('fullpath')), 'ModelTester.mat');
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
            
            if (~isfield(outputset, valuename) )

                % Flag missing value
                flag('Not found');
                continue;
            end
            
            if ( isstruct( outputset.(valuename) ) )
                
                % Skip checking of structs -- it is currently just
                % priceindex which does not need to be checked
                fprintf('\tSkipping %s because it is a struct.\n', valuename);
                continue;
            end
                
            if any(isnan(outputset.(valuename)(:)))

                % Flag NaN value
                flag('NaN value');
                continue;
            end
            
            if any(size(outputset.(valuename)) ~= size(targetset.(valuename)))

                % Flag for size mismatch
                flag('Size mismatch');
                continue;
            end
            
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
            
        end % end field loop

        % Identify new values, if any
        outputvaluenames = fieldnames(outputset);

        for j = 1:length(outputvaluenames)

            valuename = outputvaluenames{j};

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

end

