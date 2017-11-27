%%
% Generate moments distribution for transition path
% 
%%

function [] = TransitionMoments(scenario)

    if strcmp(scenario.economy, 'steady')
        error('TransitionMoments scenario cannot be steady state.');        
    end
                         
    %% Scenarios and directories
    sc_steady  = scenario.currentPolicy.steady;
    if strcmp(scenario.economy, 'open')
        sc_base    = scenario.currentPolicy.open();
    else
        sc_base    = scenario.currentPolicy.closed();
    end
    steady_dir = PathFinder.getWorkingDir(sc_steady);
    sc_dir     = PathFinder.getWorkingDir(scenario);
    base_dir   = PathFinder.getWorkingDir(sc_base);


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

    % Define captaxshare
    s = ParamGenerator.tax( sc_steady );
    captaxshare_ss = s.captaxshare;
    s = ParamGenerator.tax( scenario );
    captaxshare = s.captaxshare;

    %% Households distribution

    [dist dist_static] = append_dist(nz, nk, nb, T_life, ng, T_model, ndem, steady_dir, sc_dir);

    %% Import taxable income and define quintiles

    [inc inc_static] = append_decisions('INC', nz, nk, nb, T_life, ng, T_model, ndem, kv, zs, captaxshare, captaxshare_ss, steady_dir, sc_dir, base_dir);

    inc_groups = zeros(T_model+1,5);
    inc_index  = zeros(T_model+1,5);
    sort_inc   = zeros(nz*nk*nb*T_life*ng*ndem, T_model+1);
    inc_groups_static = zeros(T_model+1,5);
    inc_index_static  = zeros(T_model+1,5);
    sort_inc_static   = zeros(nz*nk*nb*T_life*ng*ndem, T_model+1);

    [inc_groups sort_inc inc_index] = get_quintiles(inc, dist, sort_inc, inc_index, T_model);
    inc_groups = scenario.modelunit_dollar * inc_groups;

    [inc_groups_static sort_inc_static inc_index_static] = get_quintiles(inc_static, dist_static, sort_inc, inc_index, T_model);
    inc_groups_static = scenario.modelunit_dollar * inc_groups_static;


    %% Generate groups for other variables based on taxable income distribution

    for var_name = {'CIT', 'K', 'PIT', 'SST', 'BEN', 'LABINC', 'AINC', 'TAX', 'TOTINC', 'TOTINCwSS'}

        [var var_static] = append_decisions(var_name{1}, nz, nk, nb, T_life, ng, T_model, ndem, kv, zs, captaxshare, captaxshare_ss, steady_dir, sc_dir, base_dir);
        quintiles.(var_name{1}) = generate_moments(var, dist, sort_inc_static, inc_index_static, T_model);
        quintiles.(strcat(var_name{1},'_static')) = generate_moments(var_static, dist_static, sort_inc_static, inc_index_static, T_model);
        quintiles.(strcat(var_name{1},'_delta'))  = quintiles.(var_name{1})./ quintiles.(strcat(var_name{1},'_static'));

    end

    save(fullfile(sc_dir, 'quintiles.mat' ), '-struct', 'quintiles');
    

    header = {'year', 'INC_delta_q1', 'INC_delta_q2', 'INC_delta_q3', 'INC_delta_q4', 'INC_delta_q5', ...
                   'CIT_delta_q1', 'CIT_delta_q2', 'CIT_delta_q3', 'CIT_delta_q4', 'CIT_delta_q5', ...
                   'asset_delta_q1', 'asset_delta_q2', 'asset_delta_q3', 'asset_delta_q4', 'asset_delta_q5', ...
                   'PIT_delta_q1', 'PIT_delta_q2', 'PIT_delta_q3', 'PIT_delta_q4', 'PIT_delta_q5', ...
                   'LABINC_delta_q1', 'LABINC_delta_q2', 'LABINC_delta_q3', 'LABINC_delta_q4', 'LABINC_delta_q5', ...
                   'AINC_delta_q1', 'AINC_delta_q2', 'AINC_delta_q3', 'AINC_delta_q4', 'AINC_delta_q5', ...
                   'TAX_delta_q1', 'TAX_delta_q2', 'TAX_delta_q3', 'TAX_delta_q4', 'TAX_delta_q5', ...
                   'TOTINC_delta_q1', 'TOTINC_delta_q2', 'TOTINC_delta_q3', 'TOTINC_delta_q4', 'TOTINC_delta_q5', ...
                   'TOTINCwSS_delta_q1', 'TOTINCwSS_delta_q2', 'TOTINCwSS_delta_q3', 'TOTINCwSS_delta_q4', 'TOTINCwSS_delta_q5'};

    quint_table = table([2017:1:(2017 + T_model)]', inc_groups(:,1)./inc_groups_static(:,1), ...
                        inc_groups(:,2)./inc_groups_static(:,2), inc_groups(:,3)./inc_groups_static(:,3), ...
                        inc_groups(:,4)./inc_groups_static(:,4), inc_groups(:,5)./inc_groups_static(:,5), ...
                        quintiles.CIT_delta(:,1), quintiles.CIT_delta(:,2), ...
                        quintiles.CIT_delta(:,3), quintiles.CIT_delta(:,4), quintiles.CIT_delta(:,5), ...
                        quintiles.K_delta(:,1), quintiles.K_delta(:,2), quintiles.K_delta(:,3), ...
                        quintiles.K_delta(:,4), quintiles.K_delta(:,5), quintiles.PIT_delta(:,1),  ...
                        quintiles.PIT_delta(:,2), quintiles.PIT_delta(:,3), quintiles.PIT_delta(:,4), ...
                        quintiles.PIT_delta(:,5), quintiles.LABINC_delta(:,1),  ...
                        quintiles.LABINC_delta(:,2), quintiles.LABINC_delta(:,3), quintiles.LABINC_delta(:,4), ...
                        quintiles.LABINC_delta(:,5), quintiles.AINC_delta(:,1),  ...
                        quintiles.AINC_delta(:,2), quintiles.AINC_delta(:,3), quintiles.AINC_delta(:,4), ...
                        quintiles.AINC_delta(:,5), quintiles.TAX_delta(:,1),  ...
                        quintiles.TAX_delta(:,2), quintiles.TAX_delta(:,3), quintiles.TAX_delta(:,4), ...
                        quintiles.TAX_delta(:,5), quintiles.TOTINC_delta(:,1), quintiles.TOTINC_delta(:,2), ...
                        quintiles.TOTINC_delta(:,3), quintiles.TOTINC_delta(:,4), quintiles.TOTINC_delta(:,5), ...
                        quintiles.TOTINCwSS_delta(:,1), quintiles.TOTINCwSS_delta(:,2), quintiles.TOTINCwSS_delta(:,3), ...
                        quintiles.TOTINCwSS_delta(:,4), quintiles.TOTINCwSS_delta(:,5), ...
                        'VariableNames', header);

    writetable(quint_table, fullfile(sc_dir, 'quint_table.csv'));


    %% FUNCTIONS

    % Pre-appends steady state distribution
    function [x x_static] = append_dist(nz, nk, nb, T_life, ng, T_model, ndem, dir_ss, dir_sc)

        x        = zeros(nz,nk,nb,T_life,ng,T_model+1,ndem);
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
    function [x x_static] = append_decisions(x_name, nz, nk, nb, T_life, ng, T_model, ndem, kv, zs, captaxshare, captaxshare_ss, dir_ss, dir_sc, dir_bs)

        if strcmp(x_name, 'K')

            x        = repmat(reshape(kv, [1,nk,1,1,1,1,1]),[nz,1,nb,T_life,ng,T_model+1,ndem]);
            x_static = repmat(reshape(kv, [1,nk,1,1,1,1,1]),[nz,1,nb,T_life,ng,T_model+1,ndem]);

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


    % Get quintiles function

    function [quintiles sort_x index_x] = get_quintiles(x, dist, sort_x, index_x, T_model)

        quintiles = zeros(T_model+1, 5);
        
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

    end % get_quintiles

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
            x_t    = x_t(sortx(:,2));
            dist_t = dist_t(sortx(:,2));

            % Find quintiles
            x_cum = zeros(1,5);
            q = 1;
            for quintile = 0.2:0.2:0.8

                % Asset distributed according to taxable income
                i = q_index(2,q);
                x_cum(1,q) = sum(x_t(1:i).*dist_t(1:i));
                % Counter
                q = q + 1;

            end

            % Top quintile
            x_cum(1,q)    = sum(x_t(1:end).*dist_t(1:end));
            x_groups(t,:) = diff([0 x_cum(1,:)]);

        end

    end
    
end
