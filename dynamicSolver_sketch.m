%%
% Dynamic model solver.
% 
% Methods:
% 
%   steady(deep_params)
%       Solve steady state.
% 
%   open(deep_params, plan, gcut)
%       Solve open economy transition path.
% 
%   closed(deep_params, plan, gcut)
%       Solve closed economy transition path.
% 
%%


classdef dynamicSolver_sketch

properties (Constant)
    
    % Define default parameters
    
    param_dir = dirFinder.param;
    A = 10;
    
    
    % (Can use nested structures as well to group parameters)
    default = struct('param_dir', param_dir, ...
                     'A'        , A        );
    
end


methods (Static)
    
    function [] = steady(definition) %#ok<INUSD>
        
    end
    
    
    function [] = open(definition, partag)
        
        if ~exist('partag', 'var'), partag = ''; end
        dynamicSolver.transition(true, definition, partag);
        
    end
    
    
    function [] = closed(definition, partag)
        
        if ~exist('partag', 'var'), partag = ''; end
        dynamicSolver.transition(false, definition, partag);
        
    end
    
end


methods (Static, Access = private)
    
    function [] = transition(isopen, definition, partag)
        
        
        default = dynamicSolver.default;
                
        % Extract parameters from default, overwriting with parameters from definition if available
        
        % (dirFinder should generate a name based on parameters supplied in definition)
        
        isbase = strcmp(plan, 'base');
        % Enforce baseline conditions (e.g. gcut = +0.00) regardless of parameters in definition
        
        
    end
    
end

end