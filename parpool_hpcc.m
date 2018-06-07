% Initialize parallel pool on HPCC, supporting multiple simultaneous pools
function parpool_hpcc()
    
    % Get local parallel pool cluster configuration
    config = parcluster('local');
    
    % Define, create, and set parallel pool management directory
    %   High precision current time used to avoid conflicts between multiple near-simultaneous pools
    parpool_dir = fullfile(getenv('TMP'), sprintf('%u', fix(now*1e13)));
    mkdir(parpool_dir)
    config.JobStorageLocation = parpool_dir;
    
    % Initialize parallel pool
    parpool(config)
    
end
