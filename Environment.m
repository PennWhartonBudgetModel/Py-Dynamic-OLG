%%
% Environment determines things like: 
%    file locations, production readiness, etc.
% 
%%
classdef Environment

    
    
    properties ( Constant, Access = private )
        param_dirs = ... 
            struct( 'Development', ...
                    struct(     'cbo'           , '2017-09-14'      ...
                            ,   'sim'           , '2017-09-14'      ...
                            ,   'taxplan'       , '2017-09-14'      ...
                            ,   'calibration'   , '2017-08-15-12-16-danielav-193b73c' ...
                            ), ...
                    'Testing', ...
                    struct(     'cbo'           , '2017-09-14'      ...
                            ,   'sim'           , '2017-09-14'      ...
                            ,   'taxplan'       , '2017-09-14'      ...
                            ,   'calibration'   , '2017-08-15-12-16-danielav-193b73c' ...
                            ), ...
                    'Production', ...
                    struct(     'cbo'           , '2017-09-14'      ...
                            ,   'sim'           , '2017-09-14'      ...
                            ,   'taxplan'       , '2017-09-14'      ...
                            ,   'calibration'   , '2017-08-15-12-16-danielav-193b73c' ...
                            ) ...
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
                e = Environment('Production', batchID );
                Environment.theEnvironment(e);
            else
                error( 'Cannot set to Production.' );
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
        function [commit_id] = get_git_commit_id()

            try
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
                commit_id = sprintf('%s-%s-%s', commit_date, commit_author, commit_hash);
            catch ex
                commit_id = '<no-git>';
            end
            
        end % get_git_commit_id

        
    end % private static methods
    
    
    
    methods (Access = private ) % Instance
        
        % Constructor
        function this = Environment( name )
            this.name = name;
        end % constructor
        
        
        % Fetch from Environment properties
        function [val] = prop_val(this, topic)
            val = Environment.param_dirs.(this.name).(topic);
        end % prop_val
        
        
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
        

        % Save directory for Scenario run
        %   Use Scenario IDs for Production runs
        function [save_dir, basedef_tag, counterdef_tag] = save(this, scenario)
            if( strcmp(this.name, 'Production') ) 
                % NOTE: currently saving all batches, including
                % 'calibration' batch to the Output dir.
                % We can write the calibration batch to the calibration
                % Input dir.
                save_dir = fullfile(    this.modelroot      ...   
                                    ,   'Output'            ...
                                    ,   scenario.batchID    ...
                                    ,   this.get_commit_tag ...
                                    );
            else
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
            if( strcmp(this.name, 'Production') ) 
                csv_dir = fullfile(     this.modelroot      ...   
                                    ,   'Output'            ...
                                    ,   scenario.batchID    ...
                                    ,   this.get_commit_tag ...
                                    ,   'CSV'               ...
                                    );               
            else
                csv_dir = this.save(scenario);
            end
        end


        % Get the location of the input files, e.g. taxplans
        function [input_root] = input_root(this)
            input_root = '\\hpcc.wharton.upenn.edu\ppi\Input';
        end % input

        
        % Full tag of commit (code + input params)
        function [commit_tag] = get_commit_tag(this)
            commit_tag = [      Environment.get_git_commit_id() ...
                            ,   this.prop_val('cbo')            ...
                            ,   this.prop_val('sim')            ...
                            ,   this.prop_val('taxplan')        ...
                            ,   this.prop_val('calibration')    ...
                         ];
        end % get_commit_tag
        

    end % methods
    

end % Environment



