%%
% Generate moments distribution for the retired population in the transition path
% 
%%

% Scenarios
params.beta =  1.003341000000000;
params.gamma = 0.680000000000000;
params.sigma = 1.500000000000000;
params.modelunit_dollar = 4.135682750000000e-05;
params.depreciation = 0.08;

scenario   = Scenario(struct('economy', 'closed', 'beta', params.beta, 'gamma', params.gamma, ...
                             'sigma', params.sigma, 'modelunit_dollar', params.modelunit_dollar, ...
                             'depreciation', params.depreciation, 'bequest_phi_1', 0, 'base_brackets', 'SenCMA'));
sc_steady  = scenario.currentPolicy.steady;
if isequal(scenario, sc_steady); error('TransitionMoments scenario economy cannot be steady state.'); end
if strcmp(scenario.economy, 'open')
    sc_base = scenario.currentPolicy.open();
else
    sc_base = scenario.currentPolicy.closed();
end
if isequal(scenario, sc_base)  ; error('Scenario must use a counterfactual policy plan.'); end

% Directories
steady_dir = PathFinder.getWorkingDir(sc_steady);
sc_dir     = PathFinder.getWorkingDir(scenario);
base_dir   = PathFinder.getWorkingDir(sc_base);


%% PARAMETERS

% Define time constants
s = ParamGenerator.timing( scenario );
T_life  = s.T_life;    % Total life years
T_work  = s.T_work;    % Last working age
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
bv     = s.bv;         % SS transfers grid

%% Import total income without Social Security transfers in the baseline economy and define percentiles

[totinc_base, ~] = append_decisions('TOTINCbase', nz, nk, nb, T_life, ng, T_model, ndem, kv, zs, steady_dir, base_dir, {});
[dist_base, ~] = append_dist(nz, nk, nb, T_life, ng, T_model, ndem, steady_dir, base_dir, 1);
[~ , sort_base, index_base] = get_percentiles(totinc_base, dist_base, T_model);

totinc_base_ss = totinc_base(:,:,:,:,:,1,:);
totinc_base_ss = totinc_base_ss(sort_base(:,1));
totinc_threshold = totinc_base_ss(index_base(1,:));

%% Households distribution

[dist, dist_static] = append_dist(nz, nk, nb, T_life, ng, T_model, ndem, steady_dir, sc_dir, 0);

%% Variables

% Capital income distribution
[AINC, AINC_static] = append_decisions('AINC', nz, nk, nb, T_life, ng, T_model, ndem, kv, zs, steady_dir, sc_dir, base_dir);
groups.AINC = generate_groups(AINC, dist, sort_base, index_base, T_model);

% Total income distribution
[TOTINC, TOTINC_static] = append_decisions('TOTINC', nz, nk, nb, T_life, ng, T_model, ndem, kv, zs, steady_dir, sc_dir, base_dir);
groups.TOTINC = generate_groups(TOTINC, dist, sort_base, index_base, T_model);

%% All households - I will focus on steady state for simplicity

dist_ss = dist(:,:,:,:,:,1,:);
ainc_ss = AINC(:,:,:,:,:,1,:);
totinc_ss = TOTINC(:,:,:,:,:,1,:);

% Eliminate observation if mass is too small
for iz = 1:nz; for ik = 1:nk; for ib = 1:nb; for age = 1:T_life; for ig = 1:ng; for idem = 1:ndem
                        if dist_ss(iz,ik,ib,age,ig,1,idem) <= 1e-11
                            ainc_ss  (iz,ik,ib,age,ig,1,idem) = 0;
                            totinc_ss(iz,ik,ib,age,ig,1,idem) = 0;
                        end
end; end; end; end; end; end

[totinc_max totinc_ind] = max(totinc_ss(:))

%% Retired households - steady state analysis only

dist_ss = dist(:,:,:,T_work+1:end,:,1,:);
ainc_ss = AINC(:,:,:,T_work+1:end,:,1,:);
totinc_ss = TOTINC(:,:,:,T_work+1:end,:,1,:);

% Eliminate observation if mass is too small
for iz = 1:nz; for ik = 1:nk; for ib = 1:nb; for age = 1:33; for ig = 1:ng; for idem = 1:ndem
                        if dist_ss(iz,ik,ib,age,ig,1,idem) <= 1e-11
                            ainc_ss  (iz,ik,ib,age,ig,1,idem) = 0;
                            totinc_ss(iz,ik,ib,age,ig,1,idem) = 0;
                        end
end; end; end; end; end; end

[totinc_max totinc_ind] = max(totinc_ss(:))

mass_totinc_max = dist_ss(:);
mass_totinc_max = mass_totinc_max(totinc_ind)

karray = repmat(reshape(kv, [1,nk,1,1,1,1,1]),[nz,1,nb,T_life-T_work,ng,1,ndem]);
kvec = karray(:);
kvec(totinc_ind)

s = load( fullfile(steady_dir, 'market.mat' ) );
totrates  = s.totrates;
check = totinc_max - totrates*kvec(totinc_ind)
check = totinc_max - ainc_ss(totinc_ind)

barray = repmat(reshape(bv, [1,1,nb,1,1,1,1]), [nz,nk,1,T_life-T_work,ng,1,ndem]);
bvec = barray(:);
pension_of_max = bvec(totinc_ind)

agev = [T_work+1:T_life];
agearray = repmat(reshape(agev, [1,1,1,T_life-T_work,1,1,1]), [nz,nk,nb,1,ng,1,ndem]);
agevec = agearray(:);
agevec(totinc_ind)

%% FUNCTIONS

% Pre-appends steady state distribution
function [x, x_static] = append_dist(nz, nk, nb, T_life, ng, T_model, ndem, dir_ss, dir_sc, base)

% Inputs:  (nz, nk, nb, T_life, ng, T_model, ndem) = dimensions of DIST (number of states for each state variable)
%          dir_ss = steady state directory
%          dir_sc = directory of the scenario of interest
% Outputs: x        = dynamic distribution array with pre-appended steady state distribution
%          x_static = static distribution array with pre-appended steady state distribution

    x        = zeros(nz,nk,nb,T_life,ng,T_model+1,ndem);

    s = load( fullfile(dir_ss, 'distribution.mat' ) );
    x       (:,:,:,:,:,1,:) = s.DIST;

    if ~base
        x_static = zeros(nz,nk,nb,T_life,ng,T_model+1,ndem);
        x_static(:,:,:,:,:,1,:) = s.DIST;
        s = load( fullfile(dir_sc, 'Static_distribution.mat' ) );
        x_static(:,:,:,:,:,2:end,:) = s.Static_DIST;
    else
        x_static = {};
    end

    s = load( fullfile(dir_sc, 'distribution.mat' ) );
    x(:,:,:,:,:,2:end,:) = s.DIST;

end

% Pre-appends steady state variables from all_decisions mat file
function [x, x_static] = append_decisions(x_name, nz, nk, nb, T_life, ng, T_model, ndem, kv, zs, dir_ss, dir_sc, dir_bs)

% Inputs:  (nz, nk, nb, T_life, ng, T_model, ndem) = dimensions of DIST (number of states for each state variable)
%          x_name = name of variable of interest
%          kv     = capital vector
%          zs     = productivity shocks matrix
%          dir_ss = steady state directory
%          dir_sc = directory of the scenario of interest
%          dir_bs = directory of the baseline case of the scenario of interest
% Outputs: x        = dynamic variable array with pre-appended steady state value
%          x_static = static variable array with pre-appended steady state value

    % Assets case
    if strcmp(x_name, 'K')

        x        = repmat(reshape(kv, [1,nk,1,1,1,1,1]),[nz,1,nb,T_life,ng,T_model+1,ndem]);
        x_static = repmat(reshape(kv, [1,nk,1,1,1,1,1]),[nz,1,nb,T_life,ng,T_model+1,ndem]);

    % Labor income case
    elseif strcmp(x_name, 'LABINC')

        x        = zeros(nz,nk,nb,T_life,ng,T_model+1,ndem);
        x_static = zeros(nz,nk,nb,T_life,ng,T_model+1,ndem);

        s        = load( fullfile(dir_ss, 'market.mat' ) );
        wages_ss = s.wages;
        f = @(X) repmat(reshape(X, [nz,nk,nb,T_life,1,1,ndem]), [1,1,1,1,ng,1,1]);
        s = load( fullfile(dir_ss, 'all_decisions.mat' ) );
        x(:,:,:,:,:,1,:)        = wages_ss*f(s.LAB .* repmat(reshape(zs, [nz,1,1,T_life,1,ndem]), [1,nk,nb,1,1,1]));
        x_static(:,:,:,:,:,1,:) = wages_ss*f(s.LAB .* repmat(reshape(zs, [nz,1,1,T_life,1,ndem]), [1,nk,nb,1,1,1]));

        s     = load( fullfile(dir_sc, 'market.mat' ) );
        wages = s.wages;
        f = @(X) repmat(reshape(X, [nz,nk,nb,T_life,1,T_model,ndem]), [1,1,1,1,ng,1,1]);
        s = load( fullfile(dir_sc, 'all_decisions.mat' ) );
        x(:,:,:,:,:,2:end,:) = f(repmat(reshape(wages, [1,1,1,1,T_model,1]), [nz,nk,nb,T_life,1,ndem]) .* s.LAB .* repmat(reshape(zs, [nz,1,1,T_life,1,ndem]), [1,nk,nb,1,T_model,1]));
        s     = load( fullfile(dir_bs, 'market.mat' ) );
        wages_static = s.wages;
        s = load( fullfile(dir_sc, 'Static_all_decisions.mat' ) );
        x_static(:,:,:,:,:,2:end,:) = f(repmat(reshape(wages_static, [1,1,1,1,T_model,1]), [nz,nk,nb,T_life,1,ndem]) .*s.LAB .* repmat(reshape(zs, [nz,1,1,T_life,1,ndem]), [1,nk,nb,1,T_model,1]));

    % Assets income case
    elseif strcmp(x_name, 'AINC')

        x        = zeros(nz,nk,nb,T_life,ng,T_model+1,ndem);
        x_static = zeros(nz,nk,nb,T_life,ng,T_model+1,ndem);

        s        = load( fullfile(dir_ss, 'market.mat' ) );
        totrates_ss  = s.totrates;
        x(:,:,:,:,:,1,:)        = totrates_ss*repmat(reshape(kv, [1,nk,1,1,1,1,1]),[nz,1,nb,T_life,ng,1,ndem]);
        x_static(:,:,:,:,:,1,:) = totrates_ss*repmat(reshape(kv, [1,nk,1,1,1,1,1]),[nz,1,nb,T_life,ng,1,ndem]);

        f = @(X) repmat(reshape(X, [1,1,1,1,1,T_model,1]), [nz,nk,nb,T_life,ng,1,ndem]);
        g = @(X) repmat(reshape(X, [1,nk,1,1,1,1,1]), [nz,1,nb,T_life,ng,T_model,ndem]);
        s     = load( fullfile(dir_sc, 'market.mat' ) );
        totrates  = s.totrates;
        x(:,:,:,:,:,2:end,:) = f(totrates).* g(kv);
        s     = load( fullfile(dir_bs, 'market.mat' ) );
        totrates_static  = s.totrates;
        x_static(:,:,:,:,:,2:end,:) = f(totrates_static).* g(kv);

    % Total taxes paid case
    elseif strcmp(x_name, 'TAX')

        x        = zeros(nz,nk,nb,T_life,ng,T_model+1,ndem);
        x_static = zeros(nz,nk,nb,T_life,ng,T_model+1,ndem);

        f = @(X) repmat(reshape(X, [nz,nk,nb,T_life,1,1,ndem]), [1,1,1,1,ng,1,1]);
        s = load( fullfile(dir_ss, 'all_decisions.mat' ) );
        x(:,:,:,:,:,1,:)        = f(s.CIT + s.PIT + s.SST);
        x_static(:,:,:,:,:,1,:) = f(s.CIT + s.PIT + s.SST);

        f = @(X) repmat(reshape(X, [nz,nk,nb,T_life,1,T_model,ndem]), [1,1,1,1,ng,1,1]);
        s = load( fullfile(dir_sc, 'all_decisions.mat' ) );
        x(:,:,:,:,:,2:end,:) = f(s.CIT + s.PIT + s.SST);
        s = load( fullfile(dir_sc, 'Static_all_decisions.mat' ) );
        x_static(:,:,:,:,:,2:end,:) = f(s.CIT + s.PIT + s.SST);

    % Total income case
    elseif strcmp(x_name, 'TOTINC')

        x        = zeros(nz,nk,nb,T_life,ng,T_model+1,ndem);
        x_static = zeros(nz,nk,nb,T_life,ng,T_model+1,ndem);

        s        = load( fullfile(dir_ss, 'market.mat' ) );
        wages_ss = s.wages;
        totrates_ss = s.totrates;
        f = @(X) repmat(reshape(X, [nz,nk,nb,T_life,1,1,ndem]), [1,1,1,1,ng,1,1]);
        s = load( fullfile(dir_ss, 'all_decisions.mat' ) );
        x(:,:,:,:,:,1,:)        = f(wages_ss*s.LAB .* repmat(reshape(zs, [nz,1,1,T_life,1,ndem]), [1,nk,nb,1,1,1]) + ...
                                  totrates_ss*repmat(reshape(kv, [1,nk,1,1,1,1]),[nz,1,nb,T_life,1,ndem]));
        x_static(:,:,:,:,:,1,:) = f(wages_ss*s.LAB .* repmat(reshape(zs, [nz,1,1,T_life,1,ndem]), [1,nk,nb,1,1,1]) + ...
                                  totrates_ss*repmat(reshape(kv, [1,nk,1,1,1,1]),[nz,1,nb,T_life,1,ndem]));

        s     = load( fullfile(dir_sc, 'market.mat' ) );
        wages = s.wages;
        totrates = s.totrates;
        f = @(X) repmat(reshape(X, [nz,nk,nb,T_life,1,T_model,ndem]), [1,1,1,1,ng,1,1]);
        s = load( fullfile(dir_sc, 'all_decisions.mat' ) );
        x(:,:,:,:,:,2:end,:) = f(repmat(reshape(wages, [1,1,1,1,T_model,1]), [nz,nk,nb,T_life,1,ndem]) .* s.LAB .* repmat(reshape(zs, [nz,1,1,T_life,1,ndem]), [1,nk,nb,1,T_model,1]) + ...
                                 repmat(reshape(totrates, [1,1,1,1,T_model,1]), [nz,nk,nb,T_life,1,ndem]) .* repmat(reshape(kv, [1,nk,1,1,1,1]),[nz,1,nb,T_life,T_model,ndem]));
        s     = load( fullfile(dir_bs, 'market.mat' ) );
        wages_static = s.wages;
        totrates_static = s.totrates;
        s = load( fullfile(dir_sc, 'Static_all_decisions.mat' ) );
        x_static(:,:,:,:,:,2:end,:) = f(repmat(reshape(wages_static, [1,1,1,1,T_model,1]), [nz,nk,nb,T_life,1,ndem]) .*s.LAB .* repmat(reshape(zs, [nz,1,1,T_life,1,ndem]), [1,nk,nb,1,T_model,1]) + ...
                                        repmat(reshape(totrates_static, [1,1,1,1,T_model,1]), [nz,nk,nb,T_life,1,ndem]) .* repmat(reshape(kv, [1,nk,1,1,1,1]),[nz,1,nb,T_life,T_model,ndem]));

    % Total income with Social Security transfers case
    elseif strcmp(x_name, 'TOTINCwSS')

        x        = zeros(nz,nk,nb,T_life,ng,T_model+1,ndem);
        x_static = zeros(nz,nk,nb,T_life,ng,T_model+1,ndem);

        s        = load( fullfile(dir_ss, 'market.mat' ) );
        wages_ss = s.wages;
        totrates_ss = s.totrates;
        f = @(X) repmat(reshape(X, [nz,nk,nb,T_life,1,1,ndem]), [1,1,1,1,ng,1,1]);
        s = load( fullfile(dir_ss, 'all_decisions.mat' ) );
        x(:,:,:,:,:,1,:)        = f(wages_ss*s.LAB .* repmat(reshape(zs, [nz,1,1,T_life,1,ndem]), [1,nk,nb,1,1,1]) + ...
                                  totrates_ss*repmat(reshape(kv, [1,nk,1,1,1,1]),[nz,1,nb,T_life,1,ndem]) + s.BEN);
        x_static(:,:,:,:,:,1,:) = f(wages_ss*s.LAB .* repmat(reshape(zs, [nz,1,1,T_life,1,ndem]), [1,nk,nb,1,1,1]) + ...
                                  totrates_ss*repmat(reshape(kv, [1,nk,1,1,1,1]),[nz,1,nb,T_life,1,ndem]) + s.BEN);

        s     = load( fullfile(dir_sc, 'market.mat' ) );
        wages = s.wages;
        totrates = s.totrates;
        f = @(X) repmat(reshape(X, [nz,nk,nb,T_life,1,T_model,ndem]), [1,1,1,1,ng,1,1]);
        s = load( fullfile(dir_sc, 'all_decisions.mat' ) );
        x(:,:,:,:,:,2:end,:) = f(repmat(reshape(wages, [1,1,1,1,T_model,1]), [nz,nk,nb,T_life,1,ndem]) .* s.LAB .* repmat(reshape(zs, [nz,1,1,T_life,1,ndem]), [1,nk,nb,1,T_model,1]) + ...
                                 repmat(reshape(totrates, [1,1,1,1,T_model,1]), [nz,nk,nb,T_life,1,ndem]) .* repmat(reshape(kv, [1,nk,1,1,1,1]),[nz,1,nb,T_life,T_model,ndem]) + s.BEN);
        s     = load( fullfile(dir_bs, 'market.mat' ) );
        wages_static = s.wages;
        totrates_static = s.totrates;
        s = load( fullfile(dir_sc, 'Static_all_decisions.mat' ) );
        x_static(:,:,:,:,:,2:end,:) = f(repmat(reshape(wages_static, [1,1,1,1,T_model,1]), [nz,nk,nb,T_life,1,ndem]) .*s.LAB .* repmat(reshape(zs, [nz,1,1,T_life,1,ndem]), [1,nk,nb,1,T_model,1]) + ...
                                        repmat(reshape(totrates_static, [1,1,1,1,T_model,1]), [nz,nk,nb,T_life,1,ndem]) .* repmat(reshape(kv, [1,nk,1,1,1,1]),[nz,1,nb,T_life,T_model,ndem]) + s.BEN);

    % Total income without Social Security transfers case for baseline economy
    elseif strcmp(x_name, 'TOTINCbase')

        x        = zeros(nz,nk,nb,T_life,ng,T_model+1,ndem);
        x_static = {};

        s        = load( fullfile(dir_ss, 'market.mat' ) );
        wages_ss = s.wages;
        totrates_ss = s.totrates;
        f = @(X) repmat(reshape(X, [nz,nk,nb,T_life,1,1,ndem]), [1,1,1,1,ng,1,1]);
        s = load( fullfile(dir_ss, 'all_decisions.mat' ) );
        x(:,:,:,:,:,1,:)        = f(wages_ss*s.LAB .* repmat(reshape(zs, [nz,1,1,T_life,1,ndem]), [1,nk,nb,1,1,1]) + ...
                                  totrates_ss*repmat(reshape(kv, [1,nk,1,1,1,1]),[nz,1,nb,T_life,1,ndem]));

        s     = load( fullfile(dir_sc, 'market.mat' ) );
        wages = s.wages;
        totrates = s.totrates;
        f = @(X) repmat(reshape(X, [nz,nk,nb,T_life,1,T_model,ndem]), [1,1,1,1,ng,1,1]);
        s = load( fullfile(dir_sc, 'all_decisions.mat' ) );
        x(:,:,:,:,:,2:end,:) = f(repmat(reshape(wages, [1,1,1,1,T_model,1]), [nz,nk,nb,T_life,1,ndem]) .* s.LAB .* repmat(reshape(zs, [nz,1,1,T_life,1,ndem]), [1,nk,nb,1,T_model,1]) + ...
                                 repmat(reshape(totrates, [1,1,1,1,T_model,1]), [nz,nk,nb,T_life,1,ndem]) .* repmat(reshape(kv, [1,nk,1,1,1,1]),[nz,1,nb,T_life,T_model,ndem]));

    % All other variables (already stored in all_decisions.mat file)    
    else

        x        = zeros(nz,nk,nb,T_life,ng,T_model+1,ndem);
        x_static = zeros(nz,nk,nb,T_life,ng,T_model+1,ndem);

        f = @(X) repmat(reshape(X, [nz,nk,nb,T_life,1,1,ndem]), [1,1,1,1,ng,1,1]);
        s = load( fullfile(dir_ss, 'all_decisions.mat' ) );
        x(:,:,:,:,:,1,:)        = f(s.(x_name));
        x_static(:,:,:,:,:,1,:) = f(s.(x_name));

        f = @(X) repmat(reshape(X, [nz,nk,nb,T_life,1,T_model,ndem]), [1,1,1,1,ng,1,1]);
        s = load( fullfile(dir_sc, 'all_decisions.mat' ) );
        x(:,:,:,:,:,2:end,:) = f(s.(x_name));
        s = load( fullfile(dir_sc, 'Static_all_decisions.mat' ) );
        x_static(:,:,:,:,:,2:end,:) = f(s.(x_name));

    end

end


% Get percentiles function
function [percentiles, sort_x, index_x] = get_percentiles(x, dist, T_model)

% Inputs:  x       = array with the variable of interest
%          dist    = measure of households at each state of x
%          T_model = number of transition periods
% Outputs: percentiles = total amount of x held by each percentiles
%          sort_x    = array on how to sort x in ascending order at each period
%          index_x   = array with the cutoffs delimiting each quintile

    percentiles = zeros(T_model+1, 7);
    index_x   = zeros(T_model+1, 7);

    for t = 1:T_model+1

        % Vectorize variables
        x_t    = x(:,:,:,:,:,t,:);
        x_t    = x_t(:);
        dist_t = dist(:,:,:,:,:,t,:);
        dist_t = dist_t(:);
        dist_t = dist_t/sum(dist_t(:));

        % Sort variables
        [x_t, sort_x(:,t)] = sort(x_t);
        dist_t = dist_t(sort_x(:,t));

        % Find percentiles
        x_cum = zeros(1,7);
        q = 1;
        for percentile = [0.2, 0.4, 0.6, 0.8, 0.9, 0.95]

            % Total income with SS distribution
            i = find(cumsum(dist_t) >= percentile,1);
            x_cum(1,q)   = sum(x_t(1:i).*dist_t(1:i));
            index_x(t,q) = i;

            % Counter
            q = q + 1;

        end

        % Top percentile
        x_cum(1,q)   = sum(x_t.*dist_t);
        index_x(t,q) = size(x_t,1);

        % Groups
        percentiles(t,:) = diff([0 x_cum(1,:)]);

    end

end % get_percentiles

% Generate percentile-like groups according to how households are distributed wrt total income with SS percentiles in the baseline economy
function [x_groups] = generate_groups(x, dist, sortx, q_index, T_model)

% Inputs:  x        = array with the variable of interest
%          dist     = measure of households at each state of x
%          sort_x   = array on how to sort x in total income without SS ascending order at each period
%          index_x  = array with the cutoffs delimiting each total income without SS percentile
%          T_model  = number of transition periods
% Outputs: x_groups = total amount of x held by each 'percentile-like group'

    x_groups = zeros(T_model+1,7);

    for t = 1:T_model+1

        % Vectorize variables
        x_t    = x(:,:,:,:,:,t,:);
        dist_t = dist(:,:,:,:,:,t,:);
        x_t    = x_t(:);
        dist_t = dist_t(:);
        dist_t = dist_t/sum(dist_t(:));

        % Sort variables
        x_t    = x_t(sortx(:,2));
        dist_t = dist_t(sortx(:,2));

        % Find percentiles
        x_cum = zeros(1,7);
        q = 1;
        for percentile = [0.2, 0.4, 0.6, 0.8, 0.9, 0.95]

            % Asset distributed according to total income without SS
            i = q_index(2,q);
            x_cum(1,q) = sum(x_t(1:i).*dist_t(1:i));
            % Counter
            q = q + 1;

        end

        % Top quintile
        x_cum(1,q)    = sum(x_t(1:end).*dist_t(1:end));

        % Groups
        x_groups(t,:) = diff([0 x_cum(1,:)]);

    end

end % generate_groups


