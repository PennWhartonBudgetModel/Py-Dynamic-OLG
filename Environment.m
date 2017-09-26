%%
% Environment encapsulates global execution properties, such as file locations.
% 
%%
classdef Environment
    
    
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
        
        % Singleton environment wrapper
        %   Serves as both getter and setter
        %   Required to work around lack of non-constant static variables in Matlab
        function environment = theEnvironment(newenvironment)
            
            persistent environment_;
            
            % Set to new environment if provided
            if (exist('newenvironment', 'var') && isa(newenvironment, 'Environment'))
                environment_ = newenvironment;
            else
                % Initialize to development environment if uninitialized
                if (isempty(environment_))
                    environment_ = Environment('Development');
                end
            end
            
            environment = environment_;
            
        end
        
        
        
        % Get source code directory
        function [sourcedir] = source()
            sourcedir = fileparts(mfilename('fullpath'));
        end
        
        % Get HPCC root directory
        %   Assumes that the HPCC is the only non-Windows execution location
        function [hpccrootdir] = hpccroot()
            if (ispc()), d = '\\hpcc.wharton.upenn.edu'; else, d = getenv('HOME'); end
            hpccrootdir = fullfile(d, 'ppi');
        end
        
        % Get input root directory on HPCC
        function [inputrootdir] = inputroot()
            inputrootdir = fullfile(Environment.hpccroot(), 'Input');
        end
        
        % Get directory for a specific input set
        function [inputdir] = input(inputset)
            inputdir = fullfile(Environment.inputroot(), inputset, Environment.inputversions.(inputset) );
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
        
        
        % Get current environment
        function environment = getCurrent()
            environment = Environment.theEnvironment();
        end
        
        % Set current environment to development
        function [] = setToDevelopment()
            e = Environment('Development');
            Environment.theEnvironment(e);
        end
        
        % Set current environment to production
        function [] = setToProduction()
            
            % Check for uncommitted changes
            [~, uncommitted] = system('git status -s');
            
            if (isempty(uncommitted))
                e = Environment('Production');
                Environment.theEnvironment(e);
            else
                error( 'Uncommitted changes found. Cannot set to Production.' );
            end
            
        end
        
        
        
        % Get CBO parameters directory
        function [inputdir] = cbo_param()
            inputdir = Environment.input('cbo');
        end
        
        % Get microsim parameters directory
        function [inputdir] = sim_param()
            inputdir = Environment.input('sim');
        end
        
        % Get tax plan parameters directory
        function [inputdir] = taxplan_param()
            inputdir = Environment.input('taxplan');
        end
        
        % Get calibration grid directory
        function [inputdir] = calibration()
            inputdir = Environment.input('calibration');
        end
        
        
        
        % Get directory for newly generated calibration grid
        function [newcalibrationdir] = newcalibration()
            newcalibrationdir = fullfile(Environment.inputroot(), 'calibration', Environment.committag());
        end
        
        
    end
    
    
    
    
    properties (SetAccess = private)
        
        name;   % 'Development' or 'Production'
        
    end
    
    
    methods (Access = private)
        
        % Constructor
        function this = Environment(name)
            this.name = name;
        end
        
    end
    
    
    methods (Access = public)
        
        % Get raw output directory for a specific scenario
        function [savedir, basedeftag, counterdeftag] = save(this, scenario)
            switch (this.name)
                case 'Development'
                    saverootdir = fullfile(Environment.source(), 'Output');
                case 'Production'
                    saverootdir = fullfile(Environment.hpccroot(), 'DynamicModel', 'Output', Environment.committag(), scenario.batchID);
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
                    exportdir = fullfile(Environment.hpccroot(), 'Output', scenario.batchID, 'DynamicModel', Environment.committag());
            end
        end
        
    end
    

end