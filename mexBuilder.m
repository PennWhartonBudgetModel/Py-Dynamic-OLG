%%
% mex function builder.
% 
%%


classdef mexBuilder

methods (Static)
    
    % Build all mex functions
    function [] = all()
        mexBuilder.solve_cohort();
        mexBuilder.generate_distribution();
        mexBuilder.solve_firms();
    end
    
    % Build solve_cohort
    function [] = solve_cohort()
        mexBuilder.build('solve_cohort');
    end
    
    % Build generate_distribution
    function [] = generate_distribution()
        mexBuilder.build('generate_distribution');
    end

    % Build solve_firms
    function [] = solve_firms()
        mexBuilder.build('solve_firms');
    end

end

methods (Static, Access = private)
    
    % Build mex function according to source function name
    function [] = build(fname)
        
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
        
        fprintf('\nBuild complete.\n')
        
    end
    
end

end