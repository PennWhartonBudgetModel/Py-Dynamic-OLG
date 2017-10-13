%%
% Path finder for dynamic model source code, inputs, and outputs.
%
%%
classdef PathFinder


properties (Constant, Access = private)
    
    % Specify input interface versions, organized by component
    inputversions = struct( ...
        'CBO'           , struct('cbo'          , '2017-09-14'                          ), ...
        'Microsim'      , struct('microsim'     , '2017-09-14'                          ), ...
        'TaxCalculator' , struct('taxplan'      , '2017-10-12'                          ), ...
        'DynamicModel'  , struct('calibration'  , '2017-10-12-14-38-danielav-594b26c'   ));
    
end


methods (Static, Access = private)
    
    % Singleton execution mode wrapper
    %   Serves as both getter and setter
    %   Required to work around lack of non-constant static variables in Matlab
    function [mode] = ExecutionMode(newmode)
        
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
    %   Local path used if executing on the HPCC or on AWS through the HPCC
    %   Server path used if executing elsewhere
    function [hpccrootdir] = getHpccRootDir()
        if ~isempty(regexp(getenv('HOSTNAME'), '^(hpcc|aws)', 'once'))
            d = fullfile(filesep, 'home', 'wcit', 'data', 'projects');
        else
            d = fullfile([filesep, filesep], 'hpcc.wharton.upenn.edu');
        end
        hpccrootdir = fullfile(d, 'ppi');
    end
    
    % Get component root directory
    %   Defaults to dynamic model if no component specified
    function [componentrootdir] = getComponentRootDir(component)
        if (~exist('component', 'var') || isempty(component)), component = 'DynamicModel'; end
        componentrootdir = fullfile(PathFinder.getHpccRootDir(), component);
    end
    
    
    
    
    % Get working root directory
    function [workingrootdir] = getWorkingRootDir()
        switch (PathFinder.ExecutionMode())
            case 'Development', workingrootdir = fullfile(PathFinder.getSourceDir(), 'Working');
            case 'Production' , workingrootdir = fullfile(PathFinder.getComponentRootDir(), 'Internal', PathFinder.getCommitTag());
        end
    end
    
    
    
    
    % Get input directory
    function [inputdir] = getInputDir(component, interface)
        inputdir = fullfile(PathFinder.getComponentRootDir(component), 'Interfaces', PathFinder.inputversions.(component).(interface), interface);
    end
    
    % Get output directory
    function [outputdir] = getOutputDir(interface)
        switch (PathFinder.ExecutionMode())
            case 'Development', outputrootdir = PathFinder.getWorkingRootDir();
            case 'Production' , outputrootdir = fullfile(PathFinder.getComponentRootDir(), 'Interfaces', PathFinder.getCommitTag());
        end
        outputdir = fullfile(outputrootdir, interface);
    end
    
    
    
    
    % Get unique identifying tag for active Git commit
    function [committag] = getCommitTag()
        
        % Get commit date (%cd), author email address (%ae), and abbreviated hash (%h)
        [~, commitlog] = system('git --no-pager log -1 --format=%cd,%ae,,%h --date=iso');
        
        % Extract commit date and reformat
        commitdate = regexp(commitlog, '.*? .*?(?= )', 'match', 'once');
        commitdate = datestr(commitdate, 'yyyy-mm-dd-HH-MM');
        
        % Extract author username from email address
        commitauthor = regexp(commitlog, '(?<=,).*?(?=@)', 'match', 'once');
        
        % Extract abbreviated hash
        commithash = regexp(commitlog, '(?<=,,).*?(?=\n)', 'match', 'once');
        
        % Construct commit identifier
        committag = sprintf('%s-%s-%s', commitdate, commitauthor, commithash);
        
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
    
    
    
    
    % Get working directory for a scenario
    function [workingdir] = getWorkingDir(scenario)
        workingdir = fullfile(PathFinder.getWorkingRootDir(), ...
            scenario.basedeftag, scenario.counterdeftag, scenario.economy);
    end
    
    
    
    
    % Get CBO input directory
    function [cboinputdir] = getCboInputDir()
        cboinputdir = PathFinder.getInputDir('CBO', 'cbo');
    end
    
    % Get microsim input directory
    function [microsiminputdir] = getMicrosimInputDir()
        microsiminputdir = PathFinder.getInputDir('Microsim', 'microsim');
    end
    
    % Get tax plan input directory
    function [taxplaninputdir] = getTaxPlanInputDir()
        taxplaninputdir = PathFinder.getInputDir('TaxCalculator', 'taxplan');
    end
    
    % Get calibration input directory
    function [calibrationinputdir] = getCalibrationInputDir()
        calibrationinputdir = PathFinder.getInputDir('DynamicModel', 'calibration');
    end
    
    
    
    
    % Get calibration output directory
    function [calibrationoutputdir] = getCalibrationOutputDir()
        calibrationoutputdir = PathFinder.getOutputDir('calibration');
    end
    
    % Get data series output directory
    function [dataseriesoutputdir] = getDataSeriesOutputDir()
        dataseriesoutputdir = PathFinder.getOutputDir('dataseries');
    end
    
    
end


end