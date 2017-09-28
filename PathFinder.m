%%
% Global model execution settings, most importantly file paths for reading and writing data.
%
%%
classdef PathFinder


properties (Constant, Access = private)
    
    % Specify target versions for input sets
    inputversions = struct( ...
        'cbo'           , '2017-09-14', ...
        'sim'           , '2017-09-14', ...
        'taxplan'       , '2017-09-20', ...
        'calibration'   , '2017-09-20-11-12-efraim-98d1b77' ...
    );
    
end


methods (Static, Access = private)
    
    % Singleton execution mode wrapper
    %   Serves as both getter and setter
    %   Required to work around lack of non-constant static variables in Matlab
    function mode = ExecutionMode(newmode)
        
        persistent mode_;
        
        % Set to new execution mode if provided
        if (exist('newmode', 'var'))
            if (any(strcmp(newmode, {'Development', 'Production'})))
                mode_ = newmode;
            else
                error('Execution mode must be either ''Development'' or ''Production''.');
            end
        else
            % Initialize to development mode if uninitialized
            if (isempty(mode_))
                mode_ = 'Development';
            end
        end
        
        mode = mode_;
        
    end
    
    
    
    % Get HPCC root directory
    %   Assumes that the HPCC is the only non-Windows execution location
    function [hpccrootdir] = getHpccRootDir()
        if (ispc()), d = '\\hpcc.wharton.upenn.edu'; else, d = getenv('HOME'); end
        hpccrootdir = fullfile(d, 'ppi');
    end
    
    % Get input root directory on HPCC
    function [inputrootdir] = getInputRootDir()
        inputrootdir = fullfile(PathFinder.getHpccRootDir(), 'Input');
    end
    
    % Get directory for a specific input set
    function [inputdir] = getInputDir(inputset)
        inputdir = fullfile(PathFinder.getInputRootDir(), inputset, PathFinder.inputversions.(inputset) );
    end
    
    
    
    % Get unique identifying tag for active Git commit
    function [commit_tag] = getCommitTag()
        
        % Get commit date (%cd), author email address (%ae), and abbreviated hash (%h)
        [~, commit_log] = system('git --no-pager log -1 --format=%cd,%ae,,%h --date=iso');
        
        % Extract commit date and reformat
        commit_date = regexp(commit_log, '.*? .*?(?= )', 'match', 'once');
        commit_date = datestr(commit_date, 'yyyy-mm-dd-HH-MM');
        
        % Extract author username from email address
        commit_author = regexp(commit_log, '(?<=,).*?(?=@)', 'match', 'once');
        
        % Extract abbreviated hash
        commit_hash = regexp(commit_log, '(?<=,,).*?(?=\n)', 'match', 'once');
        
        % Construct commit identifier
        commit_tag = sprintf('%s-%s-%s', commit_date, commit_author, commit_hash);
        
    end
    
    
end


methods (Static, Access = public)
    
    
    % Get execution mode
    function mode = getExecutionMode()
        mode = PathFinder.ExecutionMode();
    end
    
    % Set execution mode to development
    function [] = setToDevelopmentMode()
        PathFinder.ExecutionMode('Development');
    end
    
    % Set execution mode to production
    function [] = setToProductionMode()
        
        % Check for uncommitted changes
        [~, uncommitted] = system('git status -s');
        
        if (isempty(uncommitted))
            PathFinder.ExecutionMode('Production');
        else
            error( 'Uncommitted changes found. Cannot set to production mode.' );
        end
        
    end
    
    
    
    % Get source code directory
    function [sourcedir] = getSourceDir()
        sourcedir = fileparts(mfilename('fullpath'));
    end
    
    
    
    % Get CBO parameters directory
    function [cboparamdir] = getCboParamDir()
        cboparamdir = PathFinder.getInputDir('cbo');
    end
    
    % Get microsim parameters directory
    function [simparamdir] = getSimParamDir()
        simparamdir = PathFinder.getInputDir('sim');
    end
    
    % Get tax plan parameters directory
    function [taxplanparamdir] = getTaxPlanParamDir()
        taxplanparamdir = PathFinder.getInputDir('taxplan');
    end
    
    % Get calibration grid directory
    function [calibrationdir] = getCalibrationDir()
        calibrationdir = PathFinder.getInputDir('calibration');
    end
    
    
    
    % Get directory for newly generated calibration grid
    function [newcalibrationdir] = getNewCalibrationDir()
        newcalibrationdir = fullfile(PathFinder.getInputRoot(), 'calibration', PathFinder.getCommitTag());
    end
    
    
    
    % Get raw output directory for a specific scenario
    function [savedir, basedeftag, counterdeftag] = getSaveDir(scenario)
        switch (PathFinder.ExecutionMode())
            case 'Development'
                saverootdir = fullfile(PathFinder.getSourceDir(), 'Output');
            case 'Production'
                saverootdir = fullfile(PathFinder.getHpccRootDir(), 'DynamicModel', 'Output', PathFinder.getCommitTag());
        end
        [basedeftag, counterdeftag] = scenario.generate_tags();
        savedir = fullfile(saverootdir, basedeftag, counterdeftag, scenario.economy);
    end
    
    % Get processed output directory for a specific scenario
    function [exportdir] = getExportDir(scenario)
        switch (PathFinder.ExecutionMode())
            case 'Development'
                exportdir = PathFinder.getSaveDir(scenario);
            case 'Production'
                exportdir = fullfile(PathFinder.getHpccRootDir(), 'Output', ...
                    scenario.batchID, scenario.ID, 'DynamicModel', PathFinder.committag());
        end
    end
    
end


end