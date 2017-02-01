% Solves the firm optimization problem, distribution, and calculates aggregates.  The code will automatically determine if it is solving a transition path by finding the maximum length of either the prices
% or policy instruments.  IMPORTANT NOTE: This procedure requires that the first element of the each aggregate state space vector is the initial steady-state.


function [capital_total, labor_total, eq_total, inv_total, adjcost_total, V_total, output_total, tax_total, dist] =...
        solve_firm_optimization( prices, taxes, firm_params, tolerance) %#ok<*INUSD>
%#codegen


%% Determine if transition path and accommodate if transition.

% Find max length in prices and taxes.
T_prices = transition_length_finder(prices);
T_taxes  = transition_length_finder(taxes);
T = max( T_prices, T_taxes );

% Give warning if price length is less than policy length
if T_prices<T_taxes
    warning('Policy length exceeds price length.  Extend price length to ensure cool down period.')
end



% Assign transition if T>1
if T>1, istransition = true; else istransition = false; end

% For transition, determine initial steady-state, terminal steady-state, and extend vectors to length T.
if istransition
    % Since "firm" is a field of taxes, must ensure that all other fields of "taxes" are included for code generating purposes.
    taxes_initial    = taxes;
    taxes_transition = taxes;
    taxes_terminal   = taxes;
    
    [prices_initial, prices_transition, prices_terminal] = transition_accommodator(prices,      T);
    [taxes_initial.firm,  taxes_transition.firm,  taxes_terminal.firm] = transition_accommodator(taxes.firm,  T);
    

end


if ~istransition
    % Solve firm optimization problem for the case of steady-state.
    Vpr = zeros( firm_params.nk, firm_params.n_prodshocks); 
    [Vopt, kopt, labopt, invopt, eqopt, adjcostopt, taxopt] = solve_firm_policies(prices, taxes, firm_params, Vpr, tolerance, T); 
    
elseif istransition
    
    % First, solve terminal steady-state.
    Vpr = zeros( firm_params.nk, firm_params.n_prodshocks); 
    [Vopt, ~, ~, ~, ~, ~, ~                                                    ] = solve_firm_policies( prices_terminal  ,...
                                                                   taxes_terminal,   firm_params, Vpr,  tolerance, 1);
    
    % Second, solve backwards from terminal steady-state
    [Vopt,    kopt,    labopt,    invopt,    eqopt,    adjcostopt,    taxopt   ] = solve_firm_policies( prices_transition,...
                                                                   taxes_transition, firm_params, Vopt, tolerance, T);
    
    % Finally, solve initial steady-state.  That gives the initial distribution.
    Vpr = zeros( firm_params.nk, firm_params.n_prodshocks); 
    [Vopt_ss, kopt_ss, labopt_ss, invopt_ss, eqopt_ss, adjcostopt_ss, taxopt_ss] = solve_firm_policies( prices_initial   ,...
                                                                   taxes_initial,    firm_params, Vpr,  tolerance, 1);
end


display('Pause here')

% Solve stationary distribution
[capital_total, labor_total, eq_total, inv_total, adjcost_total, V_total, output_total, tax_total, dist] = solve_firm_distribution(kopt, labopt,...
                                                        eqopt, invopt, adjcostopt, taxopt, Vopt, revenues, outputs, firm_params);

end






function [Vopt, kopt, labopt, invopt, eqopt, adjcostopt, taxopt] = solve_firm_policies(prices, taxes, firm_params, Vpr, tolerance, T)

    if T>1, issteadystate = false; else issteadystate = true; end

    %% Unload firm parameter structure
    n_prodshocks    = firm_params.n_prodshocks;
    prod_transprob  = firm_params.prod_transprob;
    prod_shocks     = firm_params.prod_shocks;
    depreciation    = firm_params.depreciation;
    adj_cost_param  = firm_params.adj_cost_param;
    capital_share   = firm_params.capital_share;
    labor_share     = firm_params.labor_share;
    nk              = firm_params.nk; 
    kgrid           = firm_params.kgrid; %#ok<*STRNU>


    %% Initiate arrays
    time_dimension = max(1,T-1);
    
    Vopt       = zeros( nk, n_prodshocks, time_dimension);
    kopt       = zeros( nk, n_prodshocks, time_dimension);
    labopt     = zeros( nk, n_prodshocks, time_dimension);
    invopt     = zeros( nk, n_prodshocks, time_dimension);
    eqopt      = zeros( nk, n_prodshocks, time_dimension);
    adjcostopt = zeros( nk, n_prodshocks, time_dimension);
    taxopt     = zeros( nk, n_prodshocks, time_dimension);


    % Specify settings for dynamic optimization subproblems
    optim_options = optimset('Display', 'off', 'TolFun', 1e-4, 'TolX', 1e-4);
    
    %% Solve firm optimization
    max_iter = 20000;
    iter=0;
    breakcondition = false;
    while true
        iter = iter+1;
        if issteadystate
            it = 1;
        elseif ~issteadystate
            it   = T - iter;
        end
        for ik = 1:nk
            for ip = 1:n_prodshocks
                % Presolve quantities
                capital_production  = prod_shocks(ip)*(kgrid(ik)^capital_share);
                labopt   (ik,ip,it) = prices.wage(it)./(labor_share*capital_production).^(1/(labor_share-1));
                revenue             = capital_production.*labopt(ik,ip,it)^labor_share - prices.wage(it)*labopt(ik,ip,it);
                output              = prod_shocks(ip)*(kgrid(ik)^capital_share)*(labopt(ik,ip,it)^labor_share);

                EVpr = prod_transprob(ip,:)*Vpr';
                firm_params.discount_factor = prices.rate(it);
                value_function([], revenue, kgrid(ik), EVpr, it, taxes, firm_params);
                k_opt = 0;  
                V_opt = 0;  
                [k_opt,V_opt] = fminsearch(@value_function,kgrid(ik),optim_options);
                Vopt      (ik,ip,it) = -V_opt;
                kopt      (ik,ip,it) = k_opt;
                invopt    (ik,ip,it) = k_opt - (1 - depreciation)*kgrid(ik);
                adjcostopt(ik,ip,it) = (adj_cost_param/2)*((invopt(ik,ip,it)/kgrid(ik))^2);
                taxopt    (ik,ip,it) = taxes.firm.income(it)*(revenue - adjcostopt(ik,ip,it) - taxes.firm.exp_share(it)*invopt(ik,ip,it));
                eqopt     (ik,ip,it) = revenue - invopt(ik,ip,it) - adjcostopt(ik,ip,it) - taxopt(ik,ip,it);

            end
        end
        error = max(max(max(abs(Vopt(:,:,it)-Vpr))))
        if issteadystate
            if ((error<tolerance)||(iter>max_iter)), breakcondition = true; end
            Vpr = Vopt;
        elseif ~issteadystate
            if it == 1, breakcondition = true; end
            Vpr = Vopt(:,:,it);
        end
        if breakcondition, break, end
        Vpr = squeeze(Vpr);
    end
    
end





function vf = value_function( k_prime, revenue_, k_ik_, EVpr_, it_, taxes, firm_params)

% Enforce function inlining for C code generation
coder.inline('always');

% Define parameters as persistent variables
persistent initialized   
persistent revenue
persistent k_ik
persistent EVpr
persistent it
persistent kgrid
persistent depreciation
persistent adj_cost_param 
persistent discount_factor 
persistent corp_inc_tax
persistent exp_share


% Initialize parameters
if isempty(initialized)
    revenue           = 0;
	k_ik              = 0;
	EVpr              = 0;
    it                = 0;
    kgrid             = 0;
	depreciation      = 0;
	adj_cost_param    = 0;
	discount_factor   = 0;
    corp_inc_tax      = 0;
    exp_share         = 0;
    
    initialized = true;
    
end

% Set parameters if provided
if (nargin > 1)
    revenue           = revenue_;
	k_ik              = k_ik_;
	EVpr              = EVpr_;
    it                = it_;
    kgrid             = firm_params.kgrid;
	depreciation      = firm_params.depreciation;
	adj_cost_param    = firm_params.adj_cost_param;
	discount_factor   = firm_params.discount_factor;
    corp_inc_tax      = taxes.firm.income(it);
    exp_share         = taxes.firm.exp_share(it);
    
    vf = [];
    return
end

% Check boundary condition.
if (k_prime<kgrid(1))||(k_prime>kgrid(end))
    vf = Inf;
    return
end


investment = k_prime - (1-depreciation)*k_ik;
adjustment_costs = (adj_cost_param/2)*((investment/k_ik)^2);

taxbill = corp_inc_tax*(revenue - exp_share*investment - adjustment_costs);
vf = revenue - investment - adjustment_costs - taxbill + discount_factor*interp1(kgrid,EVpr,k_prime,'linear');
vf = -vf;


end








function [capital_total, labor_total, eq_total, inv_total, adjcost_total, V_total, output_total, tax_total, dist] =...
         solve_firm_distribution(kopt, labopt, eqopt, invopt, adjcostopt, taxopt, Vopt, revenues, outputs, firm_params) 
kgrid          = firm_params.kgrid;
nk             = firm_params.nk;
n_prodshocks   = firm_params.n_prodshocks;
prod_transprob = firm_params.prod_transprob;


% Solve the distribution of households
% Find indices of nearest values in kgrid and cgrid series
jk_lt = ones(size(kopt));
for elem = 1:length(kopt(:))
    jk_lt(elem) = find(kgrid(1:end-1) <= kopt(elem), 1, 'last');
end
jk_gt = jk_lt + 1;

% Calculate linear weights for lower and upper nearest values
wk_lt = (kgrid(jk_gt) - kopt) ./ (kgrid(jk_gt) - kgrid(jk_lt));
wk_gt = 1 - wk_lt;

dist = (1/(nk*n_prodshocks)).*ones(nk, n_prodshocks);
tolerance = 10e-8; 
max_iter = 10000;
iter = 0;


while true
    iter = iter + 1;
    distpr = zeros(nk,n_prodshocks);
    for ip = 1:n_prodshocks
        % Perform productivity transformation
        dist_step = permute(repmat(prod_transprob(:,ip),...
                        [1,nk]),[2,1]) .* dist;
        % Calculate distributions for next time step
        for elem = 1:numel(dist_step)
            distpr(jk_lt(elem), ip) = distpr(jk_lt(elem), ip) +  wk_lt(elem) * dist_step(elem);
            distpr(jk_gt(elem), ip) = distpr(jk_gt(elem), ip) +  wk_gt(elem) * dist_step(elem);
        end
    end
    error = max(max(max(abs(distpr-dist))));
    if ((error<tolerance)||(iter>max_iter)), break, end
    dist = distpr;
end

capital_total = sum(dist(:).*kopt(:)      );
labor_total   = sum(dist(:).*labopt(:)    );
eq_total      = sum(dist(:).*eqopt(:)     );
inv_total     = sum(dist(:).*invopt(:)    );
adjcost_total = sum(dist(:).*adjcostopt(:));
V_total       = sum(dist(:).*Vopt(:)      );
output_total  = sum(dist(:).*outputs(:)   );
tax_total     = sum(dist(:).*taxopt(:)    );

end






function [max_fieldlength] = transition_length_finder(variable)
% This function solves the maximum length of the elements in var1
variable_fields = fieldnames(variable);
n_fields    = length(variable_fields);

max_fieldlength = 0;
for iv = 1:n_fields
    fieldlength = length( variable.( variable_fields{iv} ) );
    max_fieldlength = max( fieldlength, max_fieldlength );
end

end






function [variable_initial, variable_transition, variable_terminal] = transition_accommodator(variable, T)
% This function extends the vectors to length T+1, then divides it into 3 parts: initial, transition, and terminal.
variable_fields = fieldnames(variable);
n_fields    = length(variable_fields);

for iv = 1:n_fields
    fieldlength = length( variable.( variable_fields{iv} ) );
    
    if fieldlength>T, error('Field length exceeds transition path.'), end
    
    if fieldlength<T
        
        % Extend vector by repeating last value until the end of the steady-state.
        variable.( variable_fields{iv} ) = [variable.( variable_fields{iv} ), variable.( variable_fields{iv} )(fieldlength).*ones( 1, T - fieldlength)];
        
    end
    
    variable_initial.   ( variable_fields{iv} ) =  variable.( variable_fields{iv} )(1)  ;
    variable_transition.( variable_fields{iv} ) =  variable.( variable_fields{iv} )(2:T);
    variable_terminal.  ( variable_fields{iv} ) =  variable.( variable_fields{iv} )(T)  ;

end

end







