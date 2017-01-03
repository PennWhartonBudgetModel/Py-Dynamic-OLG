
% Generates parameters.  
% 
% Methods:
% 
%     firm()
%     
%     hh(tauchen, heathcote)
%     
%     tax(no_tax)
    
    
classdef paramGenerator
    
methods (Static)
    
    function [firm_params] = firm()


        % Determine productivity shock process
        n_prodshocks = 5;
        ar1_coef     = .9;
        var_innov    = .118;
        [prod_shocks, prod_transprob, ~] = cooper_mex(n_prodshocks,0,ar1_coef,sqrt(var_innov));
        prod_shocks   = exp(prod_shocks);


        % Other parameters
        depreciation                = .085; %#ok<*NASGU>
        adj_cost_param              = .05;
        returns_to_scale_adjustment = .85;    % Want: f(k,n) = [(k^a)*(n^(1-a))]^b = (k^(a*b))*(n^((1-a)*b)).
        capital_share               = .33;
        capital_share               = capital_share*returns_to_scale_adjustment;
        labor_share                 = (1-capital_share)*returns_to_scale_adjustment;
        discount_factor             = .95;
        nk                          = 15;
        % kgrid                       = linspace(.1,2000,nk);
        kgrid                       = logspace(log(.1)/log(10),log(2000)/log(10),nk);


        % Create structure
        firm_params.n_prodshocks    = n_prodshocks;
        firm_params.prod_transprob  = prod_transprob;
        firm_params.prod_shocks     = prod_shocks;
        firm_params.depreciation    = depreciation;
        firm_params.adj_cost_param  = adj_cost_param;
        firm_params.capital_share   = capital_share; 
        firm_params.labor_share     = labor_share; 
        firm_params.nk              = nk; 
        firm_params.kgrid           = kgrid; 
        firm_params.discount_factor = discount_factor; %#ok<*STRNU>

    end
    
    
    
    
    function [hh_params] = hh(tauchen, heathcote)


        % Determine productivity shock process

        if tauchen
            n_prodshocks = 2;
            ar1_coef     = .75;
            var_innov    = .20;
            [prod_shocks, prod_transprob, ~] = cooper_mex(n_prodshocks,0,ar1_coef,sqrt(var_innov));
            prod_shocks   = exp(prod_shocks)';
        end
        if heathcote
            n_prodshocks = 3;
            prod_shocks = [0.167, 0.839, 5.087];
            prod_transprob = zeros(n_prodshocks);
            diagonal_probabilities = [.9, .99, .9];
            prod_transprob(1,1) = diagonal_probabilities(1);
            prod_transprob(1,2) = 1-diagonal_probabilities(1);
            prod_transprob(2,1) = (1-diagonal_probabilities(2))/2;
            prod_transprob(2,2) = diagonal_probabilities(2);
            prod_transprob(2,3) = (1-diagonal_probabilities(2))/2;
            prod_transprob(3,2) = 1-diagonal_probabilities(3);
            prod_transprob(3,3) = diagonal_probabilities(3);
        end



        % Other parameters
        crra = 2;
        discount_factor = .96;
        ns = 50;
        shares_grid = [0,logspace(-1,log(10000)/log(10),ns-1)];
        % shares_grid = linspace(0,10000,ns);


        % Create structure
        hh_params.n_prodshocks     = n_prodshocks;
        hh_params.prod_transprob   = prod_transprob;
        hh_params.prod_shocks      = prod_shocks;
        hh_params.crra             = crra;
        hh_params.discount_factor  = discount_factor; 
        hh_params.ns               = ns; 
        hh_params.shares_grid      = shares_grid; %#ok<*STRNU>

    end
    
    
    
    function [taxes] = tax(no_tax) %#ok<INUSD>


        if ~exist('no_tax','var')
            % Define firm tax rates
            taxes.firm.income    = .35;
            taxes.firm.exp_share = .58;


            % Define household tax rates
            taxes.hh.dividend      = .15;
            taxes.hh.inc_prog      = .181;
            taxes.hh.inc_scale     = 1.02; 
            taxes.hh.div_inc_share = 0;  

        else
            % Define firm tax rates
            taxes.firm.income    = 0;
            taxes.firm.exp_share = 1;


            % Define household tax rates
            taxes.hh.dividend      = 0;
            taxes.hh.inc_prog      = 0;
            taxes.hh.inc_scale     = 1; 
            taxes.hh.div_inc_share = 0;
        end



    end
    
end
    
end