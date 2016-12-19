function [] = generate_tax_parameters(no_tax) %#ok<INUSD>


if ~exist('no_tax','var')
    % Define firm tax rates
    taxes.firm.income    = .35;
    taxes.firm.exp_share = .58;


    % Define household tax rates
    taxes.hh.dividend      = .15;
    taxes.hh.inc_prog      = .181;
    taxes.hh.inc_scale     = 1.02; 
    taxes.hh.div_inc_share = 0;  %#ok<STRNU>
    
else
    % Define firm tax rates
    taxes.firm.income    = 0;
    taxes.firm.exp_share = 0;


    % Define household tax rates
    taxes.hh.dividend      = 0;
    taxes.hh.inc_prog      = 0;
    taxes.hh.inc_scale     = 1; 
    taxes.hh.div_inc_share = 0;  %#ok<STRNU>
end


save('tax_parameters.mat','taxes')




end