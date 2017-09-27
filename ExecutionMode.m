%%
% Global model execution settings, most importantly file paths for reading and writing data.
% 
%%
classdef ExecutionMode


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
    function mode = theExecutionMode(newmode)
        
        persistent mode_;
        
        % Set to new execution mode if provided
        if (exist('newmode', 'var') && isa(newmode, 'ExecutionMode'))
            mode_ = newmode;
        else
            % Initialize to development mode if uninitialized
            if (isempty(mode_))
                mode_ = ExecutionMode('Development');
            end
        end
        
        mode = mode_;
        
    end
    
    
    
    % Get HPCC root directory
    %   Assumes that the HPCC is the only non-Windows execution location
    function [hpccrootdir] = hpccroot()
        if (ispc()), d = '\\hpcc.wharton.upenn.edu'; else, d = getenv('HOME'); end
        hpccrootdir = fullfile(d, 'ppi');
    end
    
    % Get input root directory on HPCC
    function [inputrootdir] = inputroot()
        inputrootdir = fullfile(ExecutionMode.hpccroot(), 'Input');
    end
    
    % Get directory for a specific input set
    function [inputdir] = input(inputset)
        inputdir = fullfile(ExecutionMode.inputroot(), inputset, ExecutionMode.inputversions.(inputset) );
    end
    
    
    
    % Get unique identifying tag for active Git commit
    function [tag] = committag()
        
        % Get commit date (%cd), author email address (%ae), and abbreviated hash (%h)
        [~, commit_log] = system('git --no-pager log -1 --format=%cd,%ae,,%h --date=iso');
        
        % Extract commit date and reformat
        commit_date   = regexp(commit_log, '.*? .*?(?= )', 'match', 'once');
        commit_date   = datestr(commit_date, 'yyyy-mm-dd-HH-MM');
        
        % Extract author username from email address
        commit_author = regexp(commit_log, '(?<=,).*?(?=@)', 'match', 'once');
        
        % Extract abbreviated hash
        commit_hash   = regexp(commit_log, '(?<=,,).*?(?=\n)', 'match', 'once');
        
        % Construct commit identifier
        tag = sprintf('%s-%s-%s', commit_date, commit_author, commit_hash);
        
    end
    
    
end


methods (Static, Access = public)
    
    
    % Get current execution mode
    function mode = getCurrent()
        mode = ExecutionMode.theExecutionMode();
    end
    
    % Set current execution mode to development
    function [] = setToDevelopment()
        e = ExecutionMode('Development');
        ExecutionMode.theExecutionMode(e);
    end
    
    % Set current execution mode to production
    function [] = setToProduction()
        
        % Check for uncommitted changes
        [~, uncommitted] = system('git status -s');
        
        if (isempty(uncommitted))
            e = ExecutionMode('Production');
            ExecutionMode.theExecutionMode(e);
        else
            error( 'Uncommitted changes found. Cannot set to Production.' );
        end
        
    end
    
    
    
    % Get source code directory
    function [sourcedir] = source()
        sourcedir = fileparts(mfilename('fullpath'));
    end
    
    
    
    % Get CBO parameters directory
    function [inputdir] = cbo_param()
        inputdir = ExecutionMode.input('cbo');
    end
    
    % Get microsim parameters directory
    function [inputdir] = sim_param()
        inputdir = ExecutionMode.input('sim');
    end
    
    % Get tax plan parameters directory
    function [inputdir] = taxplan_param()
        inputdir = ExecutionMode.input('taxplan');
    end
    
    % Get calibration grid directory
    function [inputdir] = calibration()
        inputdir = ExecutionMode.input('calibration');
    end
    
    
    
    % Get directory for newly generated calibration grid
    function [newcalibrationdir] = newcalibration()
        newcalibrationdir = fullfile(ExecutionMode.inputroot(), 'calibration', ExecutionMode.committag());
    end
    
    
end




properties (SetAccess = private)
    
    name;   % 'Development' or 'Production'
    
end


methods (Access = private)
    
    % Constructor
    function this = ExecutionMode(name)
        this.name = name;
    end
    
end


methods (Access = public)
    
    % Get raw output directory for a specific scenario
    function [savedir, basedeftag, counterdeftag] = save(this, scenario)
        switch (this.name)
            case 'Development'
                saverootdir = fullfile(ExecutionMode.source(), 'Output');
            case 'Production'
                saverootdir = fullfile(ExecutionMode.hpccroot(), 'DynamicModel', 'Output', ...
                    ExecutionMode.committag(), scenario.batchID);
        end
        [basedeftag, counterdeftag] = scenario.generate_tags();
        savedir = fullfile(saverootdir, basedeftag, counterdeftag, scenario.economy);
    end
    
    % Get processed output directory for a specific scenario
    function [exportdir] = export(this, scenario)
        switch (this.name)
            case 'Development'
                exportdir = this.save(scenario);
            case 'Production'
                assert(~isempty(scenario.batchID), 'No batch ID defined for current scenario.');
                exportdir = fullfile(ExecutionMode.hpccroot(), 'Output', ...
                    scenario.batchID, scenario.ID, 'DynamicModel', ExecutionMode.committag());
        end
    end
    
end


end