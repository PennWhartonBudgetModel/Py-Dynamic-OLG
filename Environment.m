%%
% Environment determines things like: 
%    file locations, production readiness, etc.
% 
%%
classdef Environment

    
    
    properties ( Constant, Access = private )
        % HPCC share in UNC format -- note not used for UNIX
        HPCC_share = fullfile('\\hpcc.wharton.upenn.edu', 'ppi');
        
        %  NOTE: We think of a 'commit' moving from Dev->Test->Prod
        %  So, there appears to be a need for only one set of environment
        %  param dirs as it is linked to a 'commit'.
        param_dirs = ... 
            struct(     'cbo'           , '2017-09-14'      ...
                    ,   'sim'           , '2017-09-14'      ...
                    ,   'taxplan'       , '2017-09-20'      ...
                    ,   'calibration'   , '2017-09-20-11-12-efraim-98d1b77' ...
                  );
              
    end % private static properties
    
    
    
    properties ( SetAccess = private )
        name;       % Development, Production
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
        
        
        function [] = setDevelopment()
            e = Environment('Development');
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
            inputdir    = fullfile(this.inputroot(), topic, this.prop_val(topic) );
        end % input_dir
        
        
        % Input root on HPCC
        function [inputdir] = inputroot(this)
            inputdir    = fullfile(this.root(), 'Input');
        end % inputroot
        
    end % private instance methods
    
    
    
    methods  % Instance 
        
        
        % Source code directory
        function [source_dir] = source(this)
            source_dir = fileparts(mfilename('fullpath'));
        end

        
        % Root dir is the HPCC Share
        function rootdir = root(this)
            if( ispc )
                rootdir = fullfile(Environment.HPCC_share);
            else
                rootdir = getenv('HOME'); % on UNIX
            end
        end % root
        
        
        % Model root directory
        function [modelroot_dir] = modelroot(this)
            switch this.name
                case 'Production'
                    modelroot_dir = this.root();
                case 'Development'
                    modelroot_dir = this.source();
            end
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
            newdir  = fullfile( this.inputroot(), 'calibration', this.get_commit_tag());
        end % new calibration dir

        
        % Save directory for Scenario run
        %   Use Scenario IDs for Production runs
        function [save_dir, basedef_tag, counterdef_tag] = save(this, scenario)
            if( nargin < 2 )
                error( '<scenario> required.');
            end
            switch( this.name )
                case 'Production'
                    save_dir = fullfile(    this.modelroot()    ...   
                                        ,   'Output'            ...
                                        ,   scenario.batchID    ...
                                        ,   this.get_commit_tag ...
                                        ,   'MAT'               ...
                                        );
                case 'Development'
                    save_dir = fullfile(    this.modelroot()         ...
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
                case 'Development'
                    csv_dir = this.save(scenario);
                otherwise
                    csv_dir = [];
            end
        end



        
    end % methods
    

end % Environment



