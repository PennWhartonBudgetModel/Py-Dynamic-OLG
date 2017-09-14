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
        batchID;    % BatchID only for Production
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
        
        
        function [] = setProduction( batchID )
            if( nargin < 1 )
                error( '<batchID> is required for Production.' );
            end
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
        
        % Identify production run by Production stage and absence of uncommitted changes
        function [flag] = isproductionready()

            % Get source code stage
            [~, stage] = fileparts(mfilename('fullpath'));

            if ~strcmp(stage, 'Production')
                flag = false;
            else

                % Check for uncommitted changes
                % (Safeguards against any unintentional changes made in Production directory after checkout)
                [~, uncommitted] = system('git status -s');

                flag = isempty(uncommitted);

            end
        end %isproductionready
        
    end % private static methods
    
    
    methods (Access = private )
        
        % Constructor
        function this = Environment(name, batchID )
            this.name = name;
            if( nargin == 2 && ~isempty(batchID) )
                this.batchID = batchID;
            end;
        end % constructor
        
        
        % Testing/Development output directory
        function [testout_dir] = testout(this)
            testout_dir = fullfile(this.source, 'Testing');
        end
        
        % Make input dir name from static pointers struct
        function [inputdir] = fetch_dir(this, topic)
            p_dir       = Environment.param_dirs.(this.name).(topic);
            inputdir    = fullfile(this.input(), topic, p_dir );
        end
    end % private instance methods
    
    
    methods 
        
        %% FILE LOCATIONS
        % Source code directory
        function [source_dir] = source(this)
            source_dir = fileparts(mfilename('fullpath'));
        end


        % Model root directory
        function [modelroot_dir] = modelroot(this)
            modelroot_dir = fileparts(fileparts(fileparts(Environment.source)));
        end


        % Absolute root directory
        function [root_dir] = root(this)
            root_dir = fileparts(this.modelroot);
        end


        % CBO parameters directory
        function [param_dir] = cbo_param(this)
            param_dir = this.fetch_dir('cbo');
        end
        
        % SIM parameters directory
        function [param_dir] = sim_param(this)
            param_dir = this.fetch_dir('sim');
        end
        
        % Taxplan parameters directory
        function [param_dir] = taxplan_param(this)
            param_dir = this.fetch_dir('taxplan');
        end

        % Calibration grid directory
        function [param_dir] = calibration(this)
            param_dir = this.fetch_dir('calibration');
        end
        
        % Save root directory
        function [saveroot_dir] = saveroot(this)
            if( strcmp(this.name, 'Production' ) )
                saveroot_dir = fullfile(    this.modelroot      ...
                                        ,   'Output'            ...
                                        ,   this.batchID        ...
                                        ,   this.get_commit_id  ...
                                        );
                %  has subfolders MAT and CSV
            else
                saveroot_dir = this.testout;
            end
        end


        % Save directory for Scenario run
        %   Use Scenario IDs for Production runs
        function [save_dir, basedef_tag, counterdef_tag] = save(this, scenario)
            save_dir = fullfile(this.saveroot, 'MAT' );
            if( strcmp(this.name, 'Production') ) 
                save_dir = fullfile(save_dir, scenario.ID, scenario.economy);
            else
                [basedef_tag, counterdef_tag] = scenario.generate_tags();
                save_dir = fullfile(save_dir, basedef_tag, counterdef_tag, scenario.economy);
            end
        end


        % CSV save directory
        function [csv_dir] = csv(this)
            save_dir = fullfile(this.saveroot, 'CSV' );
        end


        % Get the location of the microSIM input files, e.g. taxplans
        function [input_dir] = input(this)
            input_dir = '\\hpcc.wharton.upenn.edu\ppi\Input';
        end % input

        % Get identifier for active Git commit
        function [commit_id] = get_commit_id()

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

        end

    end % methods
    

end % Environment



