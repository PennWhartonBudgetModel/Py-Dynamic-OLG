% Pete | 2016-09-14
% 
% Load file with multiple attempts and random pauses between attempts to ameliorate parallel read-write conflicts.
% 
% 


function [s] = hardyload(filename)

% Set maximum number of attempts and maximum pause time in seconds
maxtries = 50;
maxpause = 1.0;

for itry = 1:maxtries
    try
        s = load(filename);
        break
    catch
        if (itry == maxtries)
            error('Failed to load %s after %d attempts.', filename, maxtries)
        end
        pause(rand * maxpause)
    end
end

end