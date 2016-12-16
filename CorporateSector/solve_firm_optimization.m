function [capital_total, labor_total, eq_total, inv_total, adjcost_total, V_total, output_total, dist] =...
        solve_firm_optimization( prices, firm_params) %#ok<*INUSD>
%#codegen

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
Vopt        = zeros( nk, n_prodshocks);
kopt        = zeros( nk, n_prodshocks);
Vpr         = zeros( nk, n_prodshocks);
invopt      = zeros( nk, n_prodshocks);
eqopt       = zeros( nk, n_prodshocks);
adjcostopt  = zeros( nk, n_prodshocks);


%% Presolve quantities
capital_production =  (prod_shocks*(kgrid.^capital_share))';
labopt             = (prices.wage./(labor_share*capital_production)).^(1/(labor_share-1));
revenues           = capital_production.*labopt.^labor_share - prices.wage.*labopt;
outputs            = ((prod_shocks*((kgrid.^capital_share)).*(labopt'.^labor_share)))';


% Specify settings for dynamic optimization subproblems
optim_options = optimset('Display', 'off', 'TolFun', 1e-4, 'TolX', 1e-4);


%% Solve firm optimization problem
tolerance = .001;
max_iter = 2000;
iter=0;
while true
    iter = iter+1;
    for ik = 1:nk
        for ip = 1:n_prodshocks
%             EVpr = reshape( sum(bsxfun(@times, prod_transprob(ip,:), Vpr' ), 2), [nk,1] );
            EVpr = prod_transprob(ip,:)*Vpr';
            revenue = revenues(ik,ip);
            value_function([], revenue, kgrid(ik), EVpr, firm_params);
            k_opt = 0; %#ok<NASGU>
            V_opt = 0;  %#ok<NASGU>
            [k_opt,V_opt] = fminsearch(@value_function,kgrid(ik),optim_options);
            Vopt      (ik,ip) = -V_opt;
            kopt      (ik,ip) = k_opt;
            invopt    (ik,ip) = k_opt - (1-depreciation)*kgrid(ik);
            eqopt     (ik,ip) = revenue - invopt(ik,ip) - (adj_cost_param/2)*(invopt(ik,ip)^2);
            adjcostopt(ik,ip) = (adj_cost_param/2)*(invopt(ik,ip)^2);
            
        end
    end
    error = max(max(max(abs(Vopt-Vpr))));
    if ((error<tolerance)||(iter>max_iter)), break, end
    Vpr = Vopt;

end



% Solve stationary distribution
[capital_total, labor_total, eq_total, inv_total, adjcost_total, V_total, output_total, dist] = solve_firm_distribution(kopt, labopt,...
                                                        eqopt, invopt, adjcostopt, Vopt, revenues, outputs, firm_params);

end


function vf = value_function( k_prime, revenue_, k_ik_, EVpr_, firm_params)

% Enforce function inlining for C code generation
coder.inline('always');

% Define parameters as persistent variables
% Define parameters as persistent variables
persistent initialized   
persistent revenue
persistent k_ik
persistent EVpr
persistent kgrid
persistent depreciation
persistent adj_cost_param 
persistent discount_factor 


% Initialize parameters
if isempty(initialized)
    revenue           = 0;
	k_ik              = 0;
	EVpr              = 0;
    kgrid             = 0;
	depreciation      = 0;
	adj_cost_param    = 0;
	discount_factor   = 0;
    
    initialized = true;
    
end

% Set parameters if provided
if (nargin > 1)
    revenue           = revenue_;
	k_ik              = k_ik_;
	EVpr              = EVpr_;
    kgrid             = firm_params.kgrid;
	depreciation      = firm_params.depreciation;
	adj_cost_param    = firm_params.adj_cost_param;
	discount_factor   = firm_params.discount_factor;
    
    vf = [];
    return
end

% Check boundary condition.
if (k_prime<kgrid(1))||(k_prime>kgrid(end))
    vf = Inf;
    return
end


investment = k_prime - (1-depreciation)*k_ik;
adjustment_costs = (adj_cost_param/2)*(investment^2);

vf = revenue - investment - adjustment_costs + discount_factor*interp1(kgrid,EVpr,k_prime,'linear');
vf = -vf;


end





function [capital_total, labor_total, eq_total, inv_total, adjcost_total, V_total, output_total, dist] =...
         solve_firm_distribution(kopt, labopt, eqopt, invopt, adjcostopt, Vopt, revenues, outputs, firm_params) %#ok<INUSL>
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

end

