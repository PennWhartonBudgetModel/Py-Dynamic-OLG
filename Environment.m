%%
% Environment determines things like: 
%    file locations, production readiness, etc.
% 
%%
classdef Environment

    
    
    properties ( Constant, Access = private )
        %  NOTE: We think of a 'commit' moving from Dev->Test->Prod
        %  So, there appears to be a need for only one set of environment
        %  param dirs as it is linked to a 'commit'.
        param_dirs = ... 
            struct(     'cbo'           , '2017-09-14'      ...
                    ,   'sim'           , '2017-09-14'      ...
                    ,   'taxplan'       , '2017-09-14'      ...
                    ,   'calibration'   , '2017-08-15-12-16-danielav-193b73c' ...
                  );
    end % private static properties
    
    
    
    properties ( SetAccess = private )
        name;       % Testing, Development, Production
    end % instance properties
    
    
        
    methods (Static)

        
        % Get singleton 
        function environment = getCurrent()
            % if not set, set to default -- Development
            environment = Environment.theEnvironment([]);
            if( isempty(environment) )
                e = Environment('Development');
                environment = Environment.theEnvironment(e);
            end
        end % getCurrent
        
        
        function [] = setTesting()
            e = Environment('Testing');
            Environment.theEnvironment(e);
        end
        
        
        % Production environment 
        function [] = setProduction()
            if( Environment.isproductionready() )
                e = Environment('Production');
                Environment.theEnvironment(e);
            else
                error( 'Uncommitted changes found. Cannot set to Production.' );
            end
        end 
        
    end % public static methods
    
    
    
    methods ( Static, Access = private )
        
        % Wrapper on singleton variable
        %   This is a work-around to Matlab not being able to have 
        %   static variables that are not constants.
        %   This function serves as Getter and Setter.
        function environment = theEnvironment( env )
            persistent theSingleton;
            
            if( ~isempty( env ) )
                theSingleton = env;
            end
            environment = theSingleton;
            
        end % singleton wrapper
        
        
        % Identify production run by absence of uncommitted changes
        function [flag] = isproductionready()

            % Check for uncommitted changes
            % (Safeguards against any unintentional changes made in Production directory after checkout)
            [~, uncommitted] = system('git status -s');

            flag = isempty(uncommitted);

        end %isproductionready
        
        
        % Get identifier for active Git commit
        function [commit_tag] = get_commit_tag()
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
            commit_tag = sprintf('%s-%s-%s', commit_date, commit_author, commit_hash);
            
        end % get_git_commit_id

        
        % Fetch from Environment properties
        function [val] = prop_val(topic)
            val = Environment.param_dirs.(topic);
        end % prop_val
        
        

    end % private static methods
    
    
    
    methods (Access = private ) % Instance
        
        % Constructor
        function this = Environment( name )
            this.name = name;
        end % constructor
        
        
        % Make input dir name from static pointers struct
        function [inputdir] = input_dir(this, topic)
            inputdir    = fullfile(this.input_root(), topic, this.prop_val(topic) );
        end % input_dir
        
        
    end % private instance methods
    
    
    
    methods  % Instance 
        
        
        % Source code directory
        function [source_dir] = source(this)
            source_dir = fileparts(mfilename('fullpath'));
        end


        % Model root directory
        function [modelroot_dir] = modelroot(this)
            modelroot_dir = fileparts(fileparts(fileparts(this.source)));
        end


        % Absolute root directory
        function [root_dir] = root(this)
            root_dir = fileparts(this.modelroot);
        end


        % CBO parameters directory
        function [param_dir] = cbo_param(this)
            param_dir = this.input_dir('cbo');
        end
        
        
        % SIM parameters directory
        function [param_dir] = sim_param(this)
            param_dir = this.input_dir('sim');
        end
        
        
        % Taxplan parameters directory
        function [param_dir] = taxplan_param(this)
            param_dir = this.input_dir('taxplan');
        end

        
        % Calibration grid directory
        function [param_dir] = calibration(this)
            param_dir = this.input_dir('calibration');
        end
        
        function [newdir] = new_calibration_dir(this)
            newdir  = fullfile( this.input_root(), 'calibration', this.get_commit_tag());
        end % new calibration dir

        % Save directory for Scenario run
        %   Use Scenario IDs for Production runs
        function [save_dir, basedef_tag, counterdef_tag] = save(this, scenario)
            if( nargin < 2 )
                error( '<scenario> required.');
            end
            switch( this.name )
                case 'Production'
                    save_dir = fullfile(    this.modelroot      ...   
                                        ,   'Output'            ...
                                        ,   scenario.batchID    ...
                                        ,   this.get_commit_tag ...
                                        ,   'MAT'               ...
                                        );
                case {'Testing', 'Development'}
                    save_dir = fullfile(    this.source         ...
                                        ,   'Testing'           ...
                                        );
            end
            
            [basedef_tag, counterdef_tag] = scenario.generate_tags();
            save_dir = fullfile(    save_dir            ...
                                ,   basedef_tag         ...
                                ,   counterdef_tag      ...
                                ,   scenario.economy    ...
                                );
        end


        % CSV save directory
        function [csv_dir] = csv(this, scenario)
            if( nargin < 2 )
                error( '<scenario> required.');
            end
            switch( this.name )
                case 'Production'
                    csv_dir = fullfile(     this.modelroot      ...   
                                        ,   'Output'            ...
                                        ,   scenario.batchID    ...
                                        ,   this.get_commit_tag ...
                                        ,   'CSV'               ...
                                    );               
                case {'Testing', 'Development'}
                    csv_dir = this.save(scenario);
                otherwise
                    csv_dir = [];
            end
        end


        % Get the location of the input files, e.g. taxplans
        %  REM: Unix system has no UNC pathing
        function [input_root] = input_root(this)
            if( ispc )
                input_root = fullfile('\\hpcc.wharton.upenn.edu', 'ppi', 'Input');
            else
                input_root = fullfile( this.root(), 'Input' );
            end
        end % input

        
    end % methods
    

end % Environment



