function [] = solve_firm_optimization( prices, firm_params) %#ok<*INUSD>


% Initiate arrays
Vopt   = zeros( firm_params.nk, firm_params.n_prodshocks);
kopt   = zeros( firm_params.nk, firm_params.n_prodshocks);
Vpr    = zeros( firm_params.nk, firm_params.n_prodshocks);
invopt = zeros( firm_params.nk, firm_params.n_prodshocks);
eqopt  = zeros( firm_params.nk, firm_params.n_prodshocks);

% Presolve quantities
revenues = prices.consumption*(firm_params.prod_shocks*(firm_params.kgrid.^firm_params.prod_func_param));


% Specify settings for dynamic optimization subproblems
optim_options = optimset('Display', 'off', 'TolFun', 1e-4, 'TolX', 1e-4);


% Solve firm optimization problem
tolerance = .01;
max_iter = 200;
iter=0;
while true
    iter = iter+1;
    for ik = 1:firm_params.nk
        for ip = 1:firm_params.n_prodshocks
            EVpr = reshape( sum(bsxfun(@times, firm_params.prod_transprob(ip,:), Vpr' ), 2), [nk,1] );
            revenue = revenues(ip,ik);
            value_function([],revenue,params.kgrid(ik),EVpr,firm_params);
            x_opt = 0;
            V_opt = 0; %#ok<*NASGU>
            [x_opt,V_opt] = fminsearch(@value_function,kgrid(ik),optim_options);
            Vopt  (ik,ip) = -V_opt;
            kopt  (ik,ip) = x_opt;
            invopt(ik,ip) = x_opt - (1-firm_params.depreciation)*kgrid(ik);
            eqopt (ik,ip) = revenue - invopt(ik,ip);
        end
    end
    error = max(max(max(abs(Vopt-Vpr))));
    if ((error<tolerance)||(iter>max_iter)), break, end
    Vpr = Vopt;

end


end


function vf = value_function( x_prime, revenue_, k_ik_, EVpr_, params )

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
persistent cap_adj 
persistent rate 


% Initialize parameters
if isempty(initialized)
    revenue    = 0;
	k_ik           = 0;
	EVpr           = 0;
    kgrid          = 0;
	depreciation   = 0;
	cap_adj        = 0;
	rate           = 0;
    
    initialized = true;
    
end

% Set parameters if provided
if (nargin > 1)
    revenue        = revenue_;
	k_ik           = k_ik_;
	EVpr           = EVpr_;
    kgrid          = params.kgrid;
	depreciation   = params.depreciation;
	cap_adj        = params.cap_adj;
	rate           = params.rate;
    
    vf = [];
    return
end

investment = k_prime - (1-depreciation)*k_ik;







end

