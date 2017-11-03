%%
% Generate moments distribution for transition path
% 
%%

%% Scenarios and directories

% Params for K/Y=3 with d=0.056
% params.beta =  0.98745000;
% params.gamma = 0.75000000;
% params.sigma = 1.24000000;
% params.modelunit_dollar = 4.359874681178362e-05;
% params.depreciation = 0.056;

% Params for K/Y=3 with d=0.08
params.beta =  1.003510000000000;
params.gamma = 0.680000000000000;
params.sigma = 1.500000000000000;
params.modelunit_dollar = 4.135682750000000e-05;
params.depreciation = 0.08;

% Solve for baseline steady state
scenario   = Scenario(struct('economy', 'closed', 'beta', params.beta, 'gamma', params.gamma, ...
                             'sigma', params.sigma, 'modelunit_dollar', params.modelunit_dollar, ...
                             'depreciation', params.depreciation, 'bequest_phi_1', 0, 'base_brackets', 'C'));

sc_steady  = scenario.currentPolicy.steady;
steady_dir = PathFinder.getWorkingDir(sc_steady);
sc_dir     = PathFinder.getWorkingDir(scenario);

if scenario.economy == 'closed'
    sc_base = scenario.currentPolicy.closed;
else
    sc_base = scenario.currentPolicy.open;
end
base_dir = PathFinder.getWorkingDir(sc_base);


%% PARAMETERS

% Define time constants
s = ParamGenerator.timing( scenario );
T_life  = s.T_life;    % Total life years
T_work  = s.T_work;    % Retirement age
T_model = s.T_model;   % Transition path model years

% Define grids
s = ParamGenerator.grids( scenario );
ndem   = s.ndem;       % demographic types
ng     = s.ng;         % num groups
nz     = s.nz;         % num labor productivity shocks
zs     = s.zs;         % shocks grid (by demographic type and age)
nk     = s.nk;         % num asset points
nb     = s.nb;         % num avg. earnings points
kv     = s.kv;         % asset grid


%% Households distribution

[dist dist_static] = append_dist(nz, nk, nb, T_life, ng, T_model, ndem, steady_dir, sc_dir);

%% Import taxable income and define quintiles

[inc inc_static] = append_decisions('INC', nz, nk, nb, T_life, ng, T_model, ndem, kv, steady_dir, sc_dir);

inc_groups = zeros(T_model+1,5);
inc_index  = zeros(T_model+1,5);
sort_inc   = zeros(nz*nk*nb*T_life*ng*ndem, T_model+1);
inc_groups_static = zeros(T_model+1,5);
inc_index_static  = zeros(T_model+1,5);
sort_inc_static   = zeros(nz*nk*nb*T_life*ng*ndem, T_model+1);

[inc_groups sort_inc inc_index] = get_quintiles(inc, dist, sort_inc, inc_index, T_model);
inc_groups = params.modelunit_dollar * inc_groups;

[inc_groups_static sort_inc_static inc_index_static] = get_quintiles(inc_static, dist_static, sort_inc, inc_index, T_model);
inc_groups_static = params.modelunit_dollar * inc_groups_static;


%% Generate quintiles for other variables based on taxable income distribution

var_list = {'CIT', 'K', 'PIT'};
for var_name = var_list
        
    [var var_static] = append_decisions(var_name{1}, nz, nk, nb, T_life, ng, T_model, ndem, kv, steady_dir, sc_dir);
    quintiles.(var_name{1}) = generate_moments(var, dist, sort_inc, inc_index, T_model);
    quintiles.(strcat(var_name{1},'_static')) = generate_moments(var_static, dist_static, sort_inc_static, inc_index_static, T_model);
    
end

save_dir = PathFinder.getWorkingDir(scenario);
save(fullfile(save_dir, 'quintiles.mat' ), '-struct', 'quintiles')


%% FUNCTIONS

% Pre-appends steady state distribution
function [x x_static] = append_dist(nz, nk, nb, T_life, ng, T_model, ndem, dir_ss, dir_sc)

    x = zeros(nz,nk,nb,T_life,ng,T_model+1,ndem);
    x_static = zeros(nz,nk,nb,T_life,ng,T_model+1,ndem);

    s = load( fullfile(dir_ss, 'distribution.mat' ) );
    x       (:,:,:,:,:,1,:) = s.DIST;
    x_static(:,:,:,:,:,1,:) = s.DIST;

    s = load( fullfile(dir_sc, 'distribution.mat' ) );
    x(:,:,:,:,:,2:end,:) = s.DIST;
    s = load( fullfile(dir_sc, 'Static_distribution.mat' ) );
    x_static(:,:,:,:,:,2:end,:) = s.Static_DIST;
        
end

% Pre-appends steady state variables from all_decisions mat file
function [x x_static] = append_decisions(x_name, nz, nk, nb, T_life, ng, T_model, ndem, kv, dir_ss, dir_sc)

	if x_name == 'K'
        
        x = repmat(reshape(kv, [1,nk,1,1,1,1,1]),[nz,1,nb,T_life,ng,T_model+1,ndem]);
        x_static = repmat(reshape(kv, [1,nk,1,1,1,1,1]),[nz,1,nb,T_life,ng,T_model+1,ndem]);
        
    else
        
        x = zeros(nz,nk,nb,T_life,ng,T_model+1,ndem);
        x_static = zeros(nz,nk,nb,T_life,ng,T_model+1,ndem);
        f = @(X) repmat(reshape(X, [nz,nk,nb,T_life,1,1,ndem]), [1,1,1,1,ng,1,1]);
        s = load( fullfile(dir_ss, 'all_decisions.mat' ) );
        x(:,:,:,:,:,1,:) = f(s.(x_name));
        x_static(:,:,:,:,:,1,:) = f(s.(x_name));
        f = @(X) repmat(reshape(X, [nz,nk,nb,T_life,1,T_model,ndem]), [1,1,1,1,ng,1,1]);
        s = load( fullfile(dir_sc, 'all_decisions.mat' ) );
        x(:,:,:,:,:,2:end,:) = f(s.(x_name));
        s = load( fullfile(dir_sc, 'Static_all_decisions.mat' ) );
        x_static(:,:,:,:,:,2:end,:) = f(s.(x_name));
        
    end

end


% Get quintiles function

function [quintiles sort_x index_x] = get_quintiles(x, dist, sort_x, index_x, T_model)

    for t = 1:T_model+1

    % Vectorize variables
    x_t    = x(:,:,:,:,:,t,:);
    x_t    = x_t(:);
    dist_t = dist(:,:,:,:,:,t,:);
    dist_t = dist_t(:);

    % Sort variables
    [x_t, sort_x(:,t)] = sort(x_t);
    dist_t = dist_t(sort_x(:,t));

    % Find quintiles
    x_cum = zeros(1,5);
    q = 1;
    for quintile = 0.2:0.2:0.8

        % Taxable income distribution
        i = find(cumsum(dist_t) >= quintile,1);
        x_cum(1,q)   = sum(x_t(1:i).*dist_t(1:i));
        index_x(t,q) = i;

        % Counter
        q = q + 1;

    end

    % Top quintile
    x_cum(1,q)   = sum(x_t.*dist_t);
	index_x(t,q) = size(x_t,1);
    % Groups
    quintiles(t,:) = diff([0 x_cum(1,:)]);

end

end

% generate_moments function
function [x_groups] = generate_moments(x, dist, sortx, q_index, T_model)

    x_groups = zeros(T_model+1,5);

    for t = 1:T_model+1

        % Vectorize variables
        x_t    = x(:,:,:,:,:,t,:);
        dist_t = dist(:,:,:,:,:,t,:);
        x_t    = x_t(:);
        dist_t = dist_t(:);

        % Sort variables
        x_t    = x_t(sortx(:,t));
        dist_t = dist_t(sortx(:,t));

        % Find quintiles
        x_cum = zeros(1,5);
        q = 1;
        for quintile = 0.2:0.2:0.8

            % Asset distributed according to taxable income
            i = q_index(t,q);
            x_cum(1,q) = sum(x_t(1:i).*dist_t(1:i));
            % Counter
            q = q + 1;

        end

        % Top quintile
        x_cum(1,q)    = sum(x_t(1:end).*dist_t(1:end));
        x_groups(t,:) = diff([0 x_cum(1,:)]);

    end

end
