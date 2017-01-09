function [base_dir, save_dir] = identify_dirs(experiment)

% Identify the base directory
base_dir = fileparts(mfilename('fullpath'));  % 'fileparts' gives the path name, but not the file name or extension.
                                              % mfilename gives the location of the file that's being run (i.e., identify_dirs).
sub_dir = sprintf('experiment_number=%d', experiment);
                                              
save_dir = fullfile(base_dir, 'Results', sub_dir);


end

                                              
                         
