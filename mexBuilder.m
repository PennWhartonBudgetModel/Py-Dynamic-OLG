% Pete | 2016-10-07
% 
% mex function builder.
% 
% 


classdef mexBuilder

methods (Static)
    
    % Build all mex functions
    function [] = build_all()
        mexBuilder.build_solve_dynamic_optimization
        mexBuilder.build_calculate_static_taxes
    end
    
    % Build dynamic optimization solver mex function
    function [] = build_solve_dynamic_optimization()
        build_mex('solve_dynamic_optimization')
    end
    
    % Build static aggregate tax calculator mex function
    function [] = build_calculate_static_taxes()
        build_mex('calculate_static_taxes')
    end
    
end

end


function [] = build_mex(funstr)

% Specify code generation directory
codegen_dir = sprintf('%s_codegen', funstr);

% Configure code generation
mex_cfg = coder.config('mex');

mex_cfg.ExtrinsicCalls            = false;
mex_cfg.IntegrityChecks           = false;
mex_cfg.SaturateOnIntegerOverflow = false;

% Generate mex function
codegen('-d', codegen_dir, '-config', 'mex_cfg', funstr)

% Clean up code generation directory
rmdir(codegen_dir, 's')

end