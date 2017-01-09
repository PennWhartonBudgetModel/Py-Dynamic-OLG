%%
% mex function builder.
% 
%%


classdef mexBuilder

methods (Static)
    
    % Build all mex functions
    function [] = build()
        mexBuilder.build_function('solve_cohort')
    end
    
end

methods (Static, Access = private)
    
    % Build mex function according to source function name
    function [] = build_function(fname)
        
        fprintf('\nBuilding mex function for %s.', fname)
        
        % Specify code generation directory
        codegen_dir = sprintf('%s_codegen', fname);
        
        % Configure code generation
        mex_cfg = coder.config('mex');
        
        mex_cfg.ExtrinsicCalls            = false;
        mex_cfg.IntegrityChecks           = false;
        mex_cfg.SaturateOnIntegerOverflow = false;
        
        % Generate mex function
        codegen('-d', codegen_dir, '-config', 'mex_cfg', '-o', fname, fname)
        
        % Clean up code generation directory
        rmdir(codegen_dir, 's')
        
        fprintf('\nBuild complete.\n\n')
        
    end
    
end

end