function [V, aopt, assets_total, consumption_total, dist] = solve_hh_optimization(hh_params,rate)
%#codegen
% This function solves the household optimization problem over the entire state space.

% Unload structure
n_prodshocks    = hh_params.n_prodshocks;
prod_transprob  = hh_params.prod_transprob;
prod_shocks     = hh_params.prod_shocks;
crra            = hh_params.crra;  %#ok<NASGU>
discount_factor = hh_params.discount_factor;  %#ok<NASGU>
na              = hh_params.na; 
asset_grid      = hh_params.asset_grid; 

% Specify settings for dynamic optimization subproblems
optim_options = optimset('Display', 'off', 'TolFun', 1e-8, 'TolX', 1e-6);


% Initialize arrays
V    = zeros(na, n_prodshocks);
aopt = zeros(na, n_prodshocks);
copt = zeros(na, n_prodshocks);
Vpr  = zeros(na, n_prodshocks);

tolerance = .001;
iter = 0;
maxiter = 10000;
while true
    iter = iter+1;
    for ia = 1:na
        for ip = 1:n_prodshocks
%             EVpr = sum(bsxfun(@times, prod_transprob(ip,:), Vpr),2);
            EVpr = prod_transprob(ip,:)*Vpr';
            financial_resources = rate*asset_grid(ia) + prod_shocks(ip);
            
            % Call valfunc to intiate parameters
            value_function([],hh_params,EVpr,financial_resources);
            a0 = asset_grid(ia);
            [x_opt, V_opt, exitflag] = fminsearch(@value_function, a0, optim_options);
            if exitflag~=1, break, end
            aopt(ia,ip) = x_opt;
            copt(ia,ip) = financial_resources - x_opt;
            V   (ia,ip) = -V_opt;
        end
    end
    error = max(abs(V(:)-Vpr(:)));
    if (error<tolerance)||(iter>maxiter), break, end
    Vpr = V;
end

% Solve HH distribution
if n_prodshocks>1
    [assets_total, consumption_total, dist] = solve_hh_distribution(aopt, copt, hh_params);
else
    assets_total = NaN;
    consumption_total = NaN;
    dist = NaN;
end



end



function vf = value_function( a_prime, hh_params, EVpr_, financial_resources_)

% Enforce function inlining for C code generation
coder.inline('always');

% Define parameters as persistent variables
persistent initialized
persistent EVpr
persistent financial_resources
persistent crra
persistent discount_factor
persistent asset_grid


% Initialize parameters
if isempty(initialized)
    EVpr                = 0;
    financial_resources = 0;
    crra                = 0;
    discount_factor     = 0;
    asset_grid          = 0;

    
    initialized = true;
    
end

% Set parameters if provided
if (nargin > 1)
    EVpr                = EVpr_;
    financial_resources = financial_resources_;
    crra                = hh_params.crra;
    discount_factor     = hh_params.discount_factor;
    asset_grid          = hh_params.asset_grid;
    
    vf = [];
    return
end

% Check boundary conditions
if a_prime<asset_grid(1) || a_prime>asset_grid(end)
    vf = Inf;
    return
end
    

consumption = financial_resources - a_prime;
if consumption<0
    vf = Inf;
    return
end

vf = (1/(1-crra))*(consumption^(1-crra)) + discount_factor*interp1(asset_grid,EVpr,a_prime,'spline');
vf = -vf;

end





function [assets_total, consumption_total, dist] = solve_hh_distribution(aopt, copt, hh_params)
asset_grid     = hh_params.asset_grid;
na             = hh_params.na;
n_prodshocks   = hh_params.n_prodshocks;
prod_transprob = hh_params.prod_transprob;


% Solve the distribution of households
% Find indices of nearest values in kgrid and cgrid series
ja_lt = ones(size(aopt));
for elem = 1:length(aopt(:))
    ja_lt(elem) = find(asset_grid(1:end-1) <= aopt(elem), 1, 'last');
end
ja_gt = ja_lt + 1;

% Calculate linear weights for lower and upper nearest values
wa_lt = (asset_grid(ja_gt) - aopt) ./ (asset_grid(ja_gt) - asset_grid(ja_lt));
wa_gt = 1 - wa_lt;

dist = (1/(na*n_prodshocks)).*ones(na, n_prodshocks);
tolerance = 10e-8; 
max_iter = 10000;
iter = 0;
while true
    iter = iter + 1;
    distpr = zeros(na,n_prodshocks);
    for ip = 1:n_prodshocks
        % Perform productivity transformation
        dist_step = permute(repmat(prod_transprob(:,ip),...
                        [1,na]),[2,1]) .* dist;
        % Calculate distributions for next time step
        for elem = 1:numel(dist_step)
            distpr(ja_lt(elem), ip) = distpr(ja_lt(elem), ip) +  wa_lt(elem) * dist_step(elem);
            distpr(ja_gt(elem), ip) = distpr(ja_gt(elem), ip) +  wa_gt(elem) * dist_step(elem);
        end
    end
    error = max(max(max(abs(distpr-dist))));
    if ((error<tolerance)||(iter>max_iter)), break, end
    dist = distpr;
end

assets_total      = sum(dist(:).*aopt(:));
consumption_total = sum(dist(:).*copt(:));

end



            
            
            
            
            
            



