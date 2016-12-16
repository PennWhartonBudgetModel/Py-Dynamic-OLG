function [V, sopt, shares_total, consumption_total, dist] = solve_hh_optimization(hh_params, prices, dividend, tolerance)
%#codegen
% This function solves the household optimization problem over the entire state space.

% Unload structure
n_prodshocks    = hh_params.n_prodshocks;
prod_transprob  = hh_params.prod_transprob;
prod_shocks     = hh_params.prod_shocks;
crra            = hh_params.crra;  %#ok<NASGU>
discount_factor = hh_params.discount_factor;  %#ok<NASGU>
ns              = hh_params.ns; 
shares_grid     = hh_params.shares_grid; 

% Specify settings for dynamic optimization subproblems
optim_options = optimset('Display', 'off', 'TolFun', 1e-8, 'TolX', 1e-6);


% Initialize arrays
V    = zeros(ns, n_prodshocks);
sopt = zeros(ns, n_prodshocks);
copt = zeros(ns, n_prodshocks);
Vpr  = zeros(ns, n_prodshocks);

iter = 0;
maxiter = 100000;
while true
    for ia = 1:ns
        for ip = 1:n_prodshocks
%             EVpr = sum(bsxfun(@times, prod_transprob(ip,:), Vpr),2);
            EVpr = prod_transprob(ip,:)*Vpr';
            financial_resources = shares_grid(ia)*(prices.fund + dividend) + prices.wage*prod_shocks(ip);
            
            % Call valfunc to intiate parameters
            value_function([], hh_params, EVpr, financial_resources, prices);
            a0 = shares_grid(ia);
            [x_opt, V_opt, exitflag] = fminsearch(@value_function, a0, optim_options);
            if exitflag~=1, break, end
            sopt(ia,ip) = x_opt;
            copt(ia,ip) = financial_resources - prices.fund*x_opt;
            V   (ia,ip) = -V_opt;
        end
    end
    error = max(abs(V(:)-Vpr(:)));
    if (error<tolerance)||(iter>maxiter), break, end
    Vpr = V;
end

% Solve HH distribution
if n_prodshocks>1
    [shares_total, consumption_total, dist] = solve_hh_distribution(sopt, copt, hh_params);
else
    shares_total = NaN;
    consumption_total = NaN;
    dist = NaN;
end



end



function vf = value_function( s_prime, hh_params, EVpr_, financial_resources_, prices)

% Enforce function inlining for C code generation
coder.inline('always');

% Define parameters as persistent variables
persistent initialized
persistent EVpr
persistent financial_resources
persistent crra
persistent discount_factor
persistent shares_grid
persistent price_fund


% Initialize parameters
if isempty(initialized)
    EVpr                = 0;
    financial_resources = 0;
    crra                = 0;
    discount_factor     = 0;
    shares_grid         = 0;
    price_fund          = 0;
    
    initialized = true;
    
end

% Set parameters if provided
if (nargin > 1)
    EVpr                = EVpr_;
    financial_resources = financial_resources_;
    crra                = hh_params.crra;
    discount_factor     = hh_params.discount_factor;
    shares_grid         = hh_params.shares_grid;
    price_fund          = prices.fund;
    
    vf = [];
    return
end

% Check boundary conditions
if s_prime<shares_grid(1) || s_prime>shares_grid(end)
    vf = Inf;
    return
end
    

consumption = financial_resources - price_fund*s_prime;
if consumption<0
    vf = Inf;
    return
end

vf = (1/(1-crra))*(consumption^(1-crra)) + discount_factor*interp1(shares_grid,EVpr,s_prime,'linear');

vf = -vf;

end





function [shares_total, consumption_total, dist] = solve_hh_distribution(sopt, copt, hh_params)
shares_grid     = hh_params.shares_grid;
ns             = hh_params.ns;
n_prodshocks   = hh_params.n_prodshocks;
prod_transprob = hh_params.prod_transprob;


% Solve the distribution of households
% Find indices of nearest values in kgrid and cgrid series
js_lt = ones(size(sopt));
for elem = 1:length(sopt(:))
    js_lt(elem) = find(shares_grid(1:end-1) <= sopt(elem), 1, 'last');
end
js_gt = js_lt + 1;

% Calculate linear weights for lower and upper nearest values
ws_lt = (shares_grid(js_gt) - sopt) ./ (shares_grid(js_gt) - shares_grid(js_lt));
ws_gt = 1 - ws_lt;

dist = (1/(ns*n_prodshocks)).*ones(ns, n_prodshocks);
tolerance = 10e-8; 
max_iter = 10000;
iter = 0;
while true
    iter = iter + 1;
    distpr = zeros(ns,n_prodshocks);
    for ip = 1:n_prodshocks
        % Perform productivity transformation
        dist_step = permute(repmat(prod_transprob(:,ip),...
                        [1,ns]),[2,1]) .* dist;
        % Calculate distributions for next time step
        for elem = 1:numel(dist_step)
            distpr(js_lt(elem), ip) = distpr(js_lt(elem), ip) +  ws_lt(elem) * dist_step(elem);
            distpr(js_gt(elem), ip) = distpr(js_gt(elem), ip) +  ws_gt(elem) * dist_step(elem);
        end
    end
    error = max(max(max(abs(distpr-dist))));
    if ((error<tolerance)||(iter>max_iter)), break, end
    dist = distpr;
end

shares_total      = sum(dist(:).*sopt(:));
consumption_total = sum(dist(:).*copt(:));

end



            
            
            
            
            
            



