%%
% Tests for firm sector.
% 
%%


classdef firmTester

properties (Constant)
    
    % Define baseline for all tests
    basedef = struct('beta', 1.1145, ...
                     'gamma', 0.5638, ...
                     'sigma', 9.0198, ...
                     'modelunit_dollars', 3.98e-05 );
                 
    extremedef = struct('beta', 1.2, ...
                     'gamma', 0.9, ...
                     'sigma', 30.0, ...
                     'modelunit_dollars', 3.98e-05 );
    % Define counterfactual for counterfactual tests
    counterdef = struct('taxplan'       , 'ryan', ...
                        'gcut'          , +0.10 , ...
                        'legal_scale'   , 1.5   , ...
                        'prem_legal'    , 1.117 , ...
                        'amnesty'       , 0.05  , ...
                        'deportation'   , 0.05  );
    
end

methods (Static)
    
    function oneRun()
        T_model = 2; nz = 3; nk = 6; nd = 6;
        zv = [0.1, 1, 2]';
        kv = [0, 1, 2, 3, 4, 100.1]';
        dv = [0, 1, 2, 3, 4, 100.1]';
        transz = [0.1, 0.5, 0.4; 0.1, 0.5, 0.4; 0.1, 0.5, 0.4;];
        isDYNAMIC = true; 
        K_static        = zeros(nz,nk,nd,T_model); 
        DEBT_static     = zeros(nz,nk,nd,T_model); 
        DEFAULT_static  = false(nz,nk,nd,T_model); 
        SHUTDOWN_static = false(nz,nk,nd,T_model); 
        tfp = 1; alpha_k = 0.3; alpha_n = 0.6; k_adjustment_param = 8; depreciation = 0.08;
        k_default_reduction = 0.9; equity_fixed_cost = 0; equity_linear_cost = 0.02;
        operation_fixed_cost = 0; operation_linear_cost = 0;
        profit_tax = 0.35; investment_deduction = 0.6; interest_deduction = 1;
        wages            = [1.0, 1]; 
        discount_factors = [0.9, 0.9];
        
        W0       = ones (nz,nk,nd); 
        Default0 = false(nz,nk,nd);
        tic;
        OPT = solve_firms(W0, Default0, ...
                        isDYNAMIC, K_static, DEBT_static, DEFAULT_static, SHUTDOWN_static, ...
                        nz, nk, nd, T_model, transz, kv, zv, dv, ...
                        tfp, alpha_k, alpha_n, k_adjustment_param, depreciation, ...
                        k_default_reduction, equity_fixed_cost, equity_linear_cost, ...
                        operation_fixed_cost, operation_linear_cost, ...
                        profit_tax, investment_deduction, interest_deduction, ...
                        wages, discount_factors );
        
        fprintf( 'Elapsed time %f seconds\n', toc );
        t           = 1;
        k1          = sum(reshape(OPT.K         (:,:,:,t), 1, []));
        debt1       = sum(reshape(OPT.DEBT      (:,:,:,t), 1, []));
        default1    = sum(reshape(OPT.DEFAULT   (:,:,:,t), 1, [])) / (nz*nk*nd);
        shutdown1   = sum(reshape(OPT.SHUTDOWN  (:,:,:,t), 1, [])) / (nz*nk*nd);

        out1        = sum(reshape(OPT.OUTPUT    (:,:,:,t), 1, []));
        coupon1     = sum(reshape(OPT.COUPON    (:,:,:,t), 1, []));
        nextcoupon1 = sum(reshape(OPT.NEXTCOUPON(:,:,:,t), 1, []));
        dividend1   = sum(reshape(OPT.DIVIDEND  (:,:,:,t), 1, []));
        
        fprintf( 'Time: t = %u \n'     , t );
        fprintf( 'K            : %f \n', k1);
        fprintf( 'Debt         : %f \n', debt1);
        fprintf( 'Default rate : %f \n', default1);
        fprintf( 'Shutdown rate: %f \n', shutdown1);
        fprintf( '\n' );
        fprintf( 'Output       : %f \n', out1 );
        fprintf( 'Coupon       : %f \n', coupon1);
        fprintf( 'Next coupon  : %f \n', nextcoupon1);
        fprintf( 'Dividend     : %f \n', dividend1);
                    
    end % testSolve
    
end % methods
end % firmTester
   

