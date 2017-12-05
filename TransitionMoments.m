%%
% Generate moments distribution for the transition path of a counterfactual economy whose working folder is complete, that is, after solving for it
% 
%%
classdef TransitionMoments

    methods (Static)

    function [] = showDistribution(scenario)

        % Scenarios
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
        
        % Define tax
        s = ParamGenerator.tax( scenario );
        shareCapitalCorporate = s.shareCapitalCorporate;

        %% Import total income without Social Security transfers in the baseline economy and define percentiles

        [totinc_base, ~] = append_decisions('TOTINCbase', nz, nk, nb, T_life, ng, T_model, ndem, kv, zs, shareCapitalCorporate, steady_dir, base_dir, {});
        [dist_base, ~] = append_dist(nz, nk, nb, T_life, ng, T_model, ndem, steady_dir, base_dir, 1);
        [~ , sort_base, index_base] = get_percentiles(totinc_base, dist_base, T_model);

        %% Households distribution

        [dist, dist_static] = append_dist(nz, nk, nb, T_life, ng, T_model, ndem, steady_dir, sc_dir, 0);

        %% Generate groups for other variables based on the distribution of total income without SS in baseline economy

        for var_name = {'CIT', 'K', 'PIT', 'SST', 'BEN', 'LABINC', 'AINC', 'TAX', 'TOTINC', 'TOTINCwSS'}

            % append steady state values
            [var, var_static] = append_decisions(var_name{1}, nz, nk, nb, T_life, ng, T_model, ndem, kv, zs, shareCapitalCorporate, steady_dir, sc_dir, base_dir);

            % generate percentile-like groups (dynamic and static)
            groups.(var_name{1}) = generate_groups(var, dist, sort_base, index_base, T_model);
            groups.(strcat(var_name{1},'_static')) = generate_groups(var_static, dist_static, sort_base, index_base, T_model);

            % generate deltas
            groups.(strcat(var_name{1},'_delta'))  = groups.(var_name{1})./ groups.(strcat(var_name{1},'_static'));

        end

        save(fullfile(sc_dir, 'groups.mat' ), '-struct', 'groups');   

        % Save deltas in a spreadsheet
        header = {'year'};
        data_table = [2017:1:(2017 + T_model)]';
        for v = fieldnames( groups )'
            if contains( v{1}, '_delta' )
                for i = 1:7
                    header_name = sprintf( '%s_%i', v{1}, i );
                    header{end+1} = header_name;
                    data_table(:,end+1) = groups.(v{1})(:,i);
                end
            end
        end

        T = array2table(data_table, 'VariableNames', header);
        writetable(T, fullfile(sc_dir, 'groups_table.csv'));


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
        function [x, x_static] = append_decisions(x_name, nz, nk, nb, T_life, ng, T_model, ndem, kv, zs, shareCapitalCorporate, dir_ss, dir_sc, dir_bs)

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
                                          totrates_ss*(1 - shareCapitalCorporate)*repmat(reshape(kv, [1,nk,1,1,1,1]),[nz,1,nb,T_life,1,ndem]));
                x_static(:,:,:,:,:,1,:) = f(wages_ss*s.LAB .* repmat(reshape(zs, [nz,1,1,T_life,1,ndem]), [1,nk,nb,1,1,1]) + ...
                                          totrates_ss*(1 - shareCapitalCorporate)*repmat(reshape(kv, [1,nk,1,1,1,1]),[nz,1,nb,T_life,1,ndem]));

                s     = load( fullfile(dir_sc, 'market.mat' ) );
                wages = s.wages;
                totrates = s.totrates;
                f = @(X) repmat(reshape(X, [nz,nk,nb,T_life,1,T_model,ndem]), [1,1,1,1,ng,1,1]);
                s = load( fullfile(dir_sc, 'all_decisions.mat' ) );
                x(:,:,:,:,:,2:end,:) = f(repmat(reshape(wages, [1,1,1,1,T_model,1]), [nz,nk,nb,T_life,1,ndem]) .* s.LAB .* repmat(reshape(zs, [nz,1,1,T_life,1,ndem]), [1,nk,nb,1,T_model,1]) + ...
                                         repmat(reshape(totrates, [1,1,1,1,T_model,1]), [nz,nk,nb,T_life,1,ndem]) .* repmat(reshape(kv, [1,nk,1,1,1,1]),[nz,1,nb,T_life,T_model,ndem]) * ...
                                         (1 - shareCapitalCorporate));
                s     = load( fullfile(dir_bs, 'market.mat' ) );
                wages_static = s.wages;
                totrates_static = s.totrates;
                s = load( fullfile(dir_sc, 'Static_all_decisions.mat' ) );
                x_static(:,:,:,:,:,2:end,:) = f(repmat(reshape(wages_static, [1,1,1,1,T_model,1]), [nz,nk,nb,T_life,1,ndem]) .*s.LAB .* repmat(reshape(zs, [nz,1,1,T_life,1,ndem]), [1,nk,nb,1,T_model,1]) + ...
                                                repmat(reshape(totrates_static, [1,1,1,1,T_model,1]), [nz,nk,nb,T_life,1,ndem]) .* repmat(reshape(kv, [1,nk,1,1,1,1]),[nz,1,nb,T_life,T_model,ndem]) * ...
                                                (1 - shareCapitalCorporate));

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
                                          totrates_ss*(1 - shareCapitalCorporate)*repmat(reshape(kv, [1,nk,1,1,1,1]),[nz,1,nb,T_life,1,ndem]) + s.BEN);
                x_static(:,:,:,:,:,1,:) = f(wages_ss*s.LAB .* repmat(reshape(zs, [nz,1,1,T_life,1,ndem]), [1,nk,nb,1,1,1]) + ...
                                          totrates_ss*(1 - shareCapitalCorporate)*repmat(reshape(kv, [1,nk,1,1,1,1]),[nz,1,nb,T_life,1,ndem]) + s.BEN);

                s     = load( fullfile(dir_sc, 'market.mat' ) );
                wages = s.wages;
                totrates = s.totrates;
                f = @(X) repmat(reshape(X, [nz,nk,nb,T_life,1,T_model,ndem]), [1,1,1,1,ng,1,1]);
                s = load( fullfile(dir_sc, 'all_decisions.mat' ) );
                x(:,:,:,:,:,2:end,:) = f(repmat(reshape(wages, [1,1,1,1,T_model,1]), [nz,nk,nb,T_life,1,ndem]) .* s.LAB .* repmat(reshape(zs, [nz,1,1,T_life,1,ndem]), [1,nk,nb,1,T_model,1]) + ...
                                         repmat(reshape(totrates, [1,1,1,1,T_model,1]), [nz,nk,nb,T_life,1,ndem]) .* repmat(reshape(kv, [1,nk,1,1,1,1]),[nz,1,nb,T_life,T_model,ndem]) * ...
                                         (1 - shareCapitalCorporate) + s.BEN);
                s     = load( fullfile(dir_bs, 'market.mat' ) );
                wages_static = s.wages;
                totrates_static = s.totrates;
                s = load( fullfile(dir_sc, 'Static_all_decisions.mat' ) );
                x_static(:,:,:,:,:,2:end,:) = f(repmat(reshape(wages_static, [1,1,1,1,T_model,1]), [nz,nk,nb,T_life,1,ndem]) .*s.LAB .* repmat(reshape(zs, [nz,1,1,T_life,1,ndem]), [1,nk,nb,1,T_model,1]) + ...
                                                repmat(reshape(totrates_static, [1,1,1,1,T_model,1]), [nz,nk,nb,T_life,1,ndem]) .* repmat(reshape(kv, [1,nk,1,1,1,1]),[nz,1,nb,T_life,T_model,ndem]) * ...
                                                (1 - shareCapitalCorporate) + s.BEN);

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
                                          totrates_ss*(1 - shareCapitalCorporate)*repmat(reshape(kv, [1,nk,1,1,1,1]),[nz,1,nb,T_life,1,ndem]));

                s     = load( fullfile(dir_sc, 'market.mat' ) );
                wages = s.wages;
                totrates = s.totrates;
                f = @(X) repmat(reshape(X, [nz,nk,nb,T_life,1,T_model,ndem]), [1,1,1,1,ng,1,1]);
                s = load( fullfile(dir_sc, 'all_decisions.mat' ) );
                x(:,:,:,:,:,2:end,:) = f(repmat(reshape(wages, [1,1,1,1,T_model,1]), [nz,nk,nb,T_life,1,ndem]) .* s.LAB .* repmat(reshape(zs, [nz,1,1,T_life,1,ndem]), [1,nk,nb,1,T_model,1]) + ...
                                         repmat(reshape(totrates, [1,1,1,1,T_model,1]), [nz,nk,nb,T_life,1,ndem]) .* repmat(reshape(kv, [1,nk,1,1,1,1]),[nz,1,nb,T_life,T_model,ndem]) * ...
                                         (1 - shareCapitalCorporate));

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

    end % showDistribution

    end % methods

end %classdef
