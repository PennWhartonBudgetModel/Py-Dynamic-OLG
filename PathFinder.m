%%
% Path finder for dynamic model source code, inputs, and outputs.
%
%%
classdef PathFinder


properties (Constant, Access = private)
    
    % Specify inputsets of interface versions, organized by component
    %    Production set is for Production & Development modes
    %    Testing set is for Testing mode
    inputversions = struct( ...
        'Production'    , struct( ...
            'CBO'           , struct('series'       , '2017-12-19'                          ), ...
            'Microsim'      , struct('microsim'     , '2017-12-04-14-35-njanetos-5386dcd'   ), ...
            'TaxCalculator' , struct('taxplan'      , '2017-12-18-Conference'               ), ...
            'SocialSecurity', struct('taxpolicy'    , '2017-12-14-Test'                     , ...
                                     'nrapolicy'    , '2017-12-14-Test'                     , ...
                                     'benpolicy'    , '2017-12-14-Test'                       ...
                              ), ...
            'DynamicModel'  , struct('calibration'  , '2017-10-20-05-56-efraim-f0a4c45'     )  ...
        ), ...
        'Testing'       , struct( ...
            'CBO'           , struct('series'       , '2017-12-19'                          ), ...
            'Microsim'      , struct('microsim'     , '2017-12-04-14-35-njanetos-5386dcd'   ), ...
            'TaxCalculator' , struct('taxplan'      , '2017-12-18-Conference'               ), ...
            'SocialSecurity', struct('taxpolicy'    , '2017-12-14-Test'                     , ...
                                     'nrapolicy'    , '2017-12-14-Test'                     , ...
                                     'benpolicy'    , '2017-12-14-Test'                       ...
                              ), ...
            'DynamicModel'  , struct('calibration'  , '2017-10-20-05-56-efraim-f0a4c45'     )  ...
        )  ...
        );
    
end


methods (Static, Access = private)
    
    % Singleton execution mode wrapper
    %   Serves as both getter and setter
    %   Required to work around lack of non-constant static variables in Matlab
    function [mode] = ExecutionMode(newmode)
        
        persistent mode_;
        
        % Set to new execution mode if provided
        if (exist('newmode', 'var'))
            if (any(strcmp(newmode, {'Development', 'Production', 'Testing'})))
                mode_ = newmode;
            else
                error('Execution mode must be in (''Development'', ''Testing'', ''Production'')');
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
            case 'Testing'    , workingrootdir = fullfile(PathFinder.getSourceDir(), 'Working');
            case 'Production' , workingrootdir = fullfile(PathFinder.getComponentRootDir(), 'Internal', PathFinder.getCommitTag());
        end
    end
    
    
    
    
    % Get input directory
    function [inputdir] = getInputDir(component, interface)
        % Production & Development input sets are the same
        if( strcmp( PathFinder.ExecutionMode(), 'Testing' ) )
            inputset = 'Testing';
        else
            inputset = 'Production';
        end
        version  = PathFinder.inputversions.(inputset).(component).(interface);
        inputdir = fullfile(PathFinder.getComponentRootDir(component), 'Interfaces', version, interface);
    end
    
    % Get output directory
    function [outputdir] = getOutputDir(interface)
        switch (PathFinder.ExecutionMode())
            case 'Development', outputrootdir = PathFinder.getWorkingRootDir();
            case 'Testing'    , outputrootdir = PathFinder.getWorkingRootDir();
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
    
    % Set execution mode to Testing
    function [] = setToTestingMode()
        PathFinder.ExecutionMode('Testing');
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
    
    
    % Get list of input set used for this ExecutionMode
    function [versionset] = getInputSet()
        % Production & Development input sets are the same
        if( strcmp( PathFinder.ExecutionMode(), 'Testing' ) )
            inputset = 'Testing';
        else
            inputset = 'Production';
        end
        
        versionset = {};
        for component = fieldnames(PathFinder.inputversions.(inputset))'
            for interface = fieldnames(PathFinder.inputversions.(inputset).(component{1}))'
                version  = PathFinder.inputversions.(inputset).(component{1}).(interface{1});
                versionset{end+1} = {component{1}, interface{1}, version};
            end
        end
    end %getInputSet
    
    
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
        cboinputdir = PathFinder.getInputDir('CBO', 'series');
    end
    
    % Get microsim input directory
    function [microsiminputdir] = getMicrosimInputDir()
        microsiminputdir = PathFinder.getInputDir('Microsim', 'microsim');
    end
    
    % Get SocialSecurity tax input directory
    function [ssinputdir] = getSocialSecurityTaxInputDir()
        ssinputdir = PathFinder.getInputDir('SocialSecurity', 'taxpolicy');
    end
    
    % Get SocialSecurity NRA input directory
    function [ssinputdir] = getSocialSecurityNRAInputDir()
        ssinputdir = PathFinder.getInputDir('SocialSecurity', 'nrapolicy');
    end
    
    % Get SS benefits plan input directory
    function [ssinputdir] = getSocialSecurityBenefitsInputDir()
        ssinputdir = PathFinder.getInputDir('SocialSecurity', 'benpolicy');
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