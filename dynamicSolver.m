%%
% Dynamic model solver.
% 
%%


classdef dynamicSolver

methods (Static)
    
    
    % Solve steady state
    function [save_dir] = steady(basedef, callertag)
        
        if ~exist('callertag' , 'var'), callertag  = ''; end
        
        save_dir = dynamicSolver.solve('steady', basedef, [], callertag);
        
    end
    
    
    % Solve open economy transition path
    function [save_dir] = open(basedef, counterdef, callertag)
        
        if ~exist('counterdef', 'var'), counterdef = []; end
        if ~exist('callertag' , 'var'), callertag  = ''; end
        
        save_dir = dynamicSolver.solve('open', basedef, counterdef, callertag);
        
    end
    
    
    % Solve closed economy transition path
    function [save_dir] = closed(basedef, counterdef, callertag)
        
        if ~exist('counterdef', 'var'), counterdef = []; end
        if ~exist('callertag' , 'var'), callertag  = ''; end
        
        save_dir = dynamicSolver.solve('closed', basedef, counterdef, callertag);
        
    end
    
    
    % Generate tags for baseline and counterfactual definitions
    function [basedef_tag, counterdef_tag] = generate_tags(basedef, counterdef)
        
        % Define baseline and counterfactual parameter formats
        basedef_format    = struct( 'beta'      , '%0.3f'   , ...
                                    'gamma'     , '%0.3f'   , ...
                                    'sigma'     , '%05.2f'  );
        
        counterdef_format = struct( 'taxplan'   , '%s'      , ...
                                    'gcut'      , '%+0.2f'  );
        
        % Define function to construct tag from definition and format specifications
        function tag = construct_tag(def, format)
            strs = {};
            for field = fields(def)'
                strs = [strs, {sprintf([field{1}, '=', format.(field{1})], def.(field{1}))}]; %#ok<AGROW>
            end
            tag = strjoin(strs, '_');
        end
        
        % Generate baseline tag
        basedef_tag = construct_tag(basedef, basedef_format);
        
        % Generate counterfactual tag
        if (isempty(counterdef) || isempty(fields(counterdef)))
            counterdef_tag = 'baseline';
        else
            counterdef_tag = construct_tag(dynamicSolver.fill_default(counterdef), counterdef_format);
        end
        
    end
    
end


methods (Static, Access = private)
    
    % Generate counterfactual definition filled with default parameter values where necessary
    function [counterdef_filled] = fill_default(counterdef)
        
        % Define default counterfactual parameter values
        % (These are the values used for a baseline run)
        counterdef_filled = struct( 'taxplan' , 'base'    , ...
                                    'gcut'    , +0.00     );
        
        % Override default parameter values with values from counterfactual definition
        for field = fields(counterdef)'
            counterdef_filled.(field{1}) = counterdef.(field{1});
        end
        
    end
    
    
    % Solve dynamic model
    function [save_dir] = solve(economy, basedef, counterdef, callertag)
        
        
        %% Initialization
        
        % Start parallel pool if JVM enabled and pool does not already exist
        if usejava('jvm'), gcp; end
        
        % Unpack parameters from baseline definition
        beta  = basedef.beta ;
        gamma = basedef.gamma;
        sigma = basedef.sigma;
        
        % Identify baseline run by empty counterfactual definition
        if isempty(counterdef), counterdef = struct(); end
        isbase = isempty(fields(counterdef));
        
        % Generate counterfactual definition filled with default parameter values where necessary
        counterdef_filled = dynamicSolver.fill_default(counterdef);
        
        % Unpack parameters from filled counterfactual definition
        taxplan = counterdef_filled.taxplan;
        gcut    = counterdef_filled.gcut   ;
        
        
        % Identify working directories
        param_dir = dirFinder.param;
        [save_dir, ~, counterdef_tag] = dirFinder.save(economy, basedef, counterdef);
        
        % Append caller tag to save directory name and generate calling tag
        % (Obviates conflicts between parallel runs)
        save_dir   = [save_dir                                  , callertag];
        callingtag = [sprintf('^%s_%s', counterdef_tag, economy), callertag];
        
        % Clear or create save directory
        if exist(save_dir, 'dir'), rmdir(save_dir, 's'), end, mkdir(save_dir)
        
        
        
        %% Parameter loading
        
        % Load global parameters
        s = load(fullfile(param_dir, 'param_global.mat'));
        
        T_life = s.T_life;
        switch economy
            case 'steady'
                T_model    = 1;
                startyears = 0;
            case {'open', 'closed'}
                T_model    = s.T_model;
                startyears = (-s.T_life+1):(s.T_model-1);
        end
        nstartyears = length(startyears);
        
        nz   = s.nz;
        nk   = s.nk;
        nb   = s.nb;
        ndem = s.ndem;
        
        zs     = s.z;
        transz = s.tr_z;
        DISTz  = s.proddist(:,1)';
        ks     = s.kgrid;
        bs     = [0; s.bgrid(2:end)];
        
        A   = s.A;
        alp = s.alp;
        d   = s.d;
        
        surv = [s.surv(1:T_life-1), 0];
        V_beq = s.phi1.*((1+ks./s.phi2).^(1-s.phi3));
        
        mu2 = s.demdist_2015 * (s.Mu2/sum(s.Mu2));
        mu3 = repmat(1-surv, [ndem,1]) .* mu2;
        
        mpci = s.mpci;
        rpci = s.rpci;
        
        deducscale  = s.deduc_scale;
        sstaxcredit = s.ss_tax_cred;
        
        
        % Load social security parameters
        s = load(fullfile(param_dir, 'param_socsec.mat'));
        
        T_work     = s.NRA(1);
        ssbenefits = s.ss_benefit(:,1:T_model);
        sstaxs     = s.ss_tax(1:T_model);
        ssincmaxs  = s.taxmax(1:T_model);
        
        
        % Load CBO parameters
        s = load(fullfile(param_dir, 'param_cbo.mat'));
        
        debtout     = s.FederalDebtHeldbythePublic(1)/100;
        fedgovtnis  = s.fedgovtnis(1:T_model);
        cborates    = s.r_cbo(1:T_model);
        cbomeanrate = mean(s.r_cbo);
        
        
        % Load income tax parameters
        s = load(fullfile(param_dir, sprintf('param_inctax_%s.mat', taxplan)));
        
        deduc_coefs = [deducscale*s.avg_deduc, s.coefs(1:2)];
        pit_coefs   = [s.limit, s.X];
        
        
        % Load business tax parameters
        s = load(fullfile(param_dir, sprintf('param_bustax_%s.mat', taxplan)));
        
        captaxshare = s.cap_tax_share;
        expshare    = s.exp_share;
        taucap      = s.tau_cap;
        taucapgain  = s.tau_capgain;
        
        qtobin = 1 - taucap*expshare;
        
        s_base = load(fullfile(param_dir, 'param_bustax_base.mat'));
        taucap_base = s_base.tau_cap;
        qtobin0 = 1 - taucap_base*s_base.exp_share;
        
        
        
        %% Aggregate generation function
        
        function [Aggregate, LABs_, DISTs] = generate_aggregates(Market, DISTs_steady, LABs_static, DISTs_static)
            
            % Define dynamic aggregate generation flag
            isdynamic = isempty(LABs_static) && isempty(DISTs_static);
            
            % Set empty values for static optimal decision values and distributions if not provided
            if isdynamic
                LABs_static  = cell(nstartyears, ndem);
                DISTs_static = cell(nstartyears, ndem);
            end
            
            % Initialize optimal labor, distribution, and cohort aggregate arrays
            LABs_   = cell(nstartyears, ndem);
            DISTs   = cell(nstartyears, ndem);
            Cohorts = cell(nstartyears, ndem);
            
            % Initialize aggregates
            series = {'assets', 'beqs', 'labeffs', 'labs', 'lfprs', 'incs', 'pits', 'ssts', 'cits', 'bens'};
            for o = series, Aggregate.(o{1}) = zeros(1,T_model); end
            for a = series, Aggregate_.(a{1}) = []; end
            
            
            for idem = 1:ndem
                
                Ks   = zeros(nz,nk,nb,T_life,T_model); Ks_par   = repmat({Ks  }, [1,nstartyears]);
                LABs = zeros(nz,nk,nb,T_life,T_model); LABs_par = repmat({LABs}, [1,nstartyears]);
                Bs   = zeros(nz,nk,nb,T_life,T_model); Bs_par   = repmat({Bs  }, [1,nstartyears]);
                INCs = zeros(nz,nk,nb,T_life,T_model); INCs_par = repmat({INCs}, [1,nstartyears]);
                PITs = zeros(nz,nk,nb,T_life,T_model); PITs_par = repmat({PITs}, [1,nstartyears]);
                SSTs = zeros(nz,nk,nb,T_life,T_model); SSTs_par = repmat({SSTs}, [1,nstartyears]);
                CITs = zeros(nz,nk,nb,T_life,T_model); CITs_par = repmat({CITs}, [1,nstartyears]);
                BENs = zeros(nz,nk,nb,T_life,T_model); BENs_par = repmat({BENs}, [1,nstartyears]);
                
                
                % Package fixed dynamic optimization parameters into anonymous function
                solve_cohort_ = @(T_past, T_shift, T_active, V0, DIST0, LAB_static, DIST_static) solve_cohort(...
                    T_past, T_shift, T_active, T_work, T_model, nz, nk, nb, zs(:,:,idem), transz, ks, bs, beta, gamma, sigma, surv, V_beq, mu2(idem,:), mu3(idem,:), ...
                    mpci, rpci, sstaxcredit, ssbenefits, sstaxs, ssincmaxs, deduc_coefs, pit_coefs, captaxshare, taucap, taucapgain, qtobin, qtobin0, ...
                    Market.beqs, Market.wages, Market.capshares, Market.debtshares, Market.caprates, Market.govrates, Market.totrates, Market.expsubs, ...
                    V0, DIST0, LAB_static, DIST_static, isdynamic);
                
                
                % Define initial distributions
                switch economy
                    
                    case 'steady'
                        DIST0s = zeros(nz,nk,nb);
                        DIST0s(:,1,1) = DISTz';
                        
                    case {'open', 'closed'}
                        if isdynamic
                            DIST0s = DISTs_steady{1,idem} ./ repmat(reshape(mu2(idem, 1:T_life), [1,1,1,T_life]), [nz,nk,nb,1]);
                        else
                            DIST0s = double.empty(0,0,0,T_life);
                        end
                        
                end
                
                
                % Initialize series of terminal utility values
                V0s = zeros(nz,nk,nb,T_life+1);
                
                
                % Solve steady state / post-transition path cohort
                if isdynamic
                    
                    % Extract starting year and derive time constants
                    startyear = startyears(end);
                    T_past  = max(-startyear, 0);
                    T_shift = max(+startyear, 0);
                    
                    % Define active time as full lifetime
                    T_active = T_life;
                    
                    % Define terminal utility values
                    V0 = zeros(nz,nk,nb);
                    
                    % Extract initial distribution
                    DIST0 = DIST0s(:,:,:,T_past+1);
                    
                    % Solve dynamic optimization
                    [DISTs{end,idem}, Cohorts{end,idem}, V, ~, LABs_{end,idem}] = solve_cohort_(T_past, T_shift, T_active, V0, DIST0, [], []);
                    
                    % Define series of terminal utility values
                    V0s(:,:,:,1:T_life) = V;
                    
                end
                
                
                switch economy
                    
                    case {'steady'}
                        
                        % Add cohort aggregates to total aggregates
                        for o = series, Aggregate.(o{1}) = Aggregate.(o{1}) + sum(Cohorts{1,idem}.(o{1})); end
                    
                    case {'open', 'closed'}
                        
                        % Solve transition path cohorts
                        parfor i = 1:nstartyears
                            
                            % Extract starting year and derive time constants
                            startyear = startyears(i);
                            T_past  = max(-startyear, 0);
                            T_shift = max(+startyear, 0);
                            
                            % Define active time as life years within modeling period
                            T_active = min(startyear+T_life, T_model) - T_shift;
                            
                            % Extract terminal utility values
                            V0 = V0s(:,:,:,min(T_model-startyear, T_life)+1); %#ok<PFBNS>
                            
                            % Extract initial distribution
                            DIST0 = DIST0s(:,:,:,T_past+1); %#ok<PFBNS>
                            
                            % Solve dynamic optimization
                            [DISTs{i,idem}, Cohort, ~, K, LAB, B, INC, PIT, SST, CIT, BEN] = solve_cohort_(T_past, T_shift, T_active, V0, DIST0, LABs_static{i,idem}, DISTs_static{i,idem});
                            
                            LABs_{i,idem} = LAB;
                            
                            for t = 1:T_active
                                
                                age  = t + T_past ;
                                year = t + T_shift;
                                
                                Ks_par  {i}(:,:,:,age,year) = K  (:,:,:,t);
                                LABs_par{i}(:,:,:,age,year) = LAB(:,:,:,t);
                                Bs_par  {i}(:,:,:,age,year) = B  (:,:,:,t);
                                INCs_par{i}(:,:,:,age,year) = INC(:,:,:,t);
                                PITs_par{i}(:,:,:,age,year) = PIT(:,:,:,t);
                                SSTs_par{i}(:,:,:,age,year) = SST(:,:,:,t);
                                CITs_par{i}(:,:,:,age,year) = CIT(:,:,:,t);
                                BENs_par{i}(:,:,:,age,year) = BEN(:,:,:,t);
                                
                            end
                            
                            
                            % Align cohort aggregates with model years
                            for o = series
                                Cohorts{i,idem}.(o{1}) = zeros(1,T_model);
                                Cohorts{i,idem}.(o{1})(T_shift+(1:T_active)) = Cohorts{i,idem}.(o{1})(T_shift+(1:T_active)) + Cohort.(o{1});
                            end
                            
                        end
                        
                        
                        
                        for i = 1:nstartyears
                            Ks   = Ks   + Ks_par  {i};
                            LABs = LABs + LABs_par{i};
                            Bs   = Bs   + Bs_par  {i};
                            INCs = INCs + INCs_par{i};
                            PITs = PITs + PITs_par{i};
                            SSTs = SSTs + SSTs_par{i};
                            CITs = CITs + CITs_par{i};
                            BENs = BENs + BENs_par{i};
                        end
                        
                        
                        if isdynamic
                            
                            % DIST__ = DISTs_steady{1,idem};
                            DIST__ = DIST0s;
                            
                            year = 1;
                            lastyear = T_model;
                            
                            while (true)
                                
                                for a = series, if (length(Aggregate_.(a{1})) < year), Aggregate_.(a{1})(year) = 0; end, end
                                
                                DIST_mu2 = DIST__ .* repmat(reshape(mu2(idem,:), [1,1,1,T_life]), [nz,nk,nb,1]);
                                
                                A_.assets  = DIST_mu2 .* repmat(reshape(ks, [1,nk,1,1]), [nz,1,nb,T_life]);
                                A_.beqs    = DIST_mu2 .* Ks  (:,:,:,:,min(year, T_model)).*repmat(reshape(1-surv, [1,1,1,T_life]), [nz,nk,nb,1]);
                                A_.labeffs = DIST_mu2 .* LABs(:,:,:,:,min(year, T_model)).*repmat(reshape(zs(:,:,idem), [nz,1,1,T_life]), [1,nk,nb,1]);
                                A_.labs    = DIST_mu2 .* LABs(:,:,:,:,min(year, T_model));
                                A_.lfprs   = DIST_mu2 .* (LABs(:,:,:,:,min(year, T_model)) > 0);
                                A_.incs    = DIST_mu2 .* INCs(:,:,:,:,min(year, T_model));
                                A_.pits    = DIST_mu2 .* PITs(:,:,:,:,min(year, T_model));
                                A_.ssts    = DIST_mu2 .* SSTs(:,:,:,:,min(year, T_model));
                                A_.cits    = DIST_mu2 .* CITs(:,:,:,:,min(year, T_model));
                                A_.bens    = DIST_mu2 .* BENs(:,:,:,:,min(year, T_model));
                                
%                                 A_.assets  = DIST__ .* repmat(reshape(ks, [1,nk,1,1]), [nz,1,nb,T_life]);
%                                 A_.beqs    = DIST__ .* Ks  (:,:,:,:,min(year, T_model)).*repmat(reshape(1-surv, [1,1,1,T_life]), [nz,nk,nb,1]);
%                                 A_.labeffs = DIST__ .* LABs(:,:,:,:,min(year, T_model)).*repmat(reshape(zs(:,:,idem), [nz,1,1,T_life]), [1,nk,nb,1]);
%                                 A_.labs    = DIST__ .* LABs(:,:,:,:,min(year, T_model));
%                                 A_.lfprs   = DIST__ .* (LABs(:,:,:,:,min(year, T_model)) > 0);
%                                 A_.incs    = DIST__ .* INCs(:,:,:,:,min(year, T_model));
%                                 A_.pits    = DIST__ .* PITs(:,:,:,:,min(year, T_model));
%                                 A_.ssts    = DIST__ .* SSTs(:,:,:,:,min(year, T_model));
%                                 A_.cits    = DIST__ .* CITs(:,:,:,:,min(year, T_model));
%                                 A_.bens    = DIST__ .* BENs(:,:,:,:,min(year, T_model));
%                                 
                                for a = series, Aggregate_.(a{1})(year) = Aggregate_.(a{1})(year) + sum(A_.(a{1})(:)); end
                                
                                
                                if (year < lastyear), year = year + 1; else, break, end
                                
                                
                                DIST_next = zeros(nz,nk,nb,T_life);
                                
                                % DIST_next(:,1,1,1) = (sum(DIST__(:))/T_life) * reshape(DISTz, [nz,1,1,1]);
                                DIST_next(:,1,1,1) = reshape(DISTz, [nz,1,1,1]);
                                
                                for age = 2:T_life
                                    
                                    % Extract optimal k and b decision values
                                    k_t = Ks(:,:,:,age-1,year-1);
                                    b_t = Bs(:,:,:,age-1,year-1);
                                    
                                    % Find indices of nearest values in ks and bs series
                                    jk_lt = ones(size(k_t));
                                    for elem = 1:length(k_t(:))
                                        jk_lt(elem) = find(ks(1:end-1) <= k_t(elem), 1, 'last');
                                    end
                                    jk_gt = jk_lt + 1;
                                    
                                    jb_lt = ones(size(b_t));
                                    for elem = 1:length(b_t(:))
                                        jb_lt(elem) = find(bs(1:end-1) <= b_t(elem), 1, 'last');
                                    end
                                    jb_gt = jb_lt + 1;
                                    
                                    % Calculate linear weights for nearest values
                                    wk_lt = (ks(jk_gt) - k_t) ./ (ks(jk_gt) - ks(jk_lt));
                                    wk_gt = 1 - wk_lt;
                                    
                                    wb_lt = (bs(jb_gt) - b_t) ./ (bs(jb_gt) - bs(jb_lt));
                                    wb_gt = 1 - wb_lt;
                                    
                                    for jz = 1:nz
                                        
                                        % Apply survival and productivity transformations to cohort distribution from current year
                                        % DIST_transz = DIST__(:,:,:,age-1) * surv(age-1) .* repmat(reshape(transz(:,jz), [nz,1,1]), [1,nk,nb]);
                                        DIST_transz = DIST__(:,:,:,age-1) .* repmat(reshape(transz(:,jz), [nz,1,1]), [1,nk,nb]);
                                        
                                        % Redistribute cohort for next year according to target indices and weights
                                        for elem = 1:numel(DIST_transz)
                                            DIST_next(jz, jk_lt(elem), jb_lt(elem), age) = DIST_next(jz, jk_lt(elem), jb_lt(elem), age) + wk_lt(elem)*wb_lt(elem)*DIST_transz(elem);
                                            DIST_next(jz, jk_gt(elem), jb_lt(elem), age) = DIST_next(jz, jk_gt(elem), jb_lt(elem), age) + wk_gt(elem)*wb_lt(elem)*DIST_transz(elem);
                                            DIST_next(jz, jk_lt(elem), jb_gt(elem), age) = DIST_next(jz, jk_lt(elem), jb_gt(elem), age) + wk_lt(elem)*wb_gt(elem)*DIST_transz(elem);
                                            DIST_next(jz, jk_gt(elem), jb_gt(elem), age) = DIST_next(jz, jk_gt(elem), jb_gt(elem), age) + wk_gt(elem)*wb_gt(elem)*DIST_transz(elem);
                                        end
                                        
                                    end
                                    
                                end
                                
                                DIST__ = DIST_next;
                                
                            end
                            
                        end
                        
                        
                        
                        % Add cohort aggregates to total aggregates
                        % (Separate loop necessary due to restrictions with use of structures within parfor loops)
                        for i = 1:nstartyears
                            for o = series, Aggregate.(o{1}) = Aggregate.(o{1}) + Cohorts{i,idem}.(o{1}); end
                        end
                        
                end
                
            end
            
        end
        
        
        
        %% Static aggregate generation
        
        if ~isbase
            
            % Identify baseline generator and save directory
            base_generator = @() dynamicSolver.solve(economy, basedef, [], callingtag);
            base_dir = dirFinder.save(economy, basedef);
            
            % Load baseline market conditions, optimal decision values, and distributions
            Market = hardyload('market.mat'       , base_generator, base_dir);
            
            s      = hardyload('decisions.mat'    , base_generator, base_dir);
            LABs_static  = s.LABs ;
            
            s      = hardyload('distributions.mat', base_generator, base_dir);
            DISTs_static = s.DISTs;
            
            
            % Generate static aggregates
            % (Intermediary structure used to filter out extraneous fields)
            [Static_] = generate_aggregates(Market, {}, LABs_static, DISTs_static);
            
            for oo = {'incs', 'pits', 'ssts', 'cits', 'bens'}
                Static.(oo{1}) = Static_.(oo{1});
            end
            
            
            % Copy additional static aggregates from baseline aggregates
            Dynamic_base = hardyload('dynamics.mat', base_generator, base_dir);
            
            for oo = {'labeffs', 'caps', 'lfprs', 'labincs', 'capincs', 'outs', 'caps_domestic', 'caps_foreign', 'debts_domestic', 'debts_foreign'}
                Static.(oo{1}) = Dynamic_base.(oo{1});
            end
            
            
            % Calculate additional static aggregates
            Static.labpits       = Static.pits .* Static.labincs ./ Static.incs;
            
            Static.cits_domestic = Static.cits + Static.pits - Static.labpits;
            Static.cits_foreign  = zeros(1,T_model);
            
            Static.caprevs       = Static.cits_domestic;
            
            
            % Save static aggregates
            save(fullfile(save_dir, 'statics.mat'), '-struct', 'Static')
            
        end
        
        
        
        %% Dynamic aggregate generation
        
        switch economy
            
            case 'steady'
                
                % Load neutral market conditions as initial conditions
                Market0 = load(fullfile(param_dir, 'market0.mat'));
                
                DISTs_steady = {};
                
            case {'open', 'closed'}
                
                % Identify steady state generator and save directory
                steady_generator = @() dynamicSolver.steady(basedef, callingtag);
                steady_dir = dirFinder.save('steady', basedef);
                
                % Load steady state market conditions as initial conditions
                Market0 = hardyload('market.mat'       , steady_generator, steady_dir);
                
                % Load steady state distributions
                s       = hardyload('distributions.mat', steady_generator, steady_dir);
                
                DISTs_steady = s.DISTs;
                
        end
        
        
        switch economy
            
            case 'open'
                
                % Load government expenditure adjustment parameters
                s = load(fullfile(param_dir, 'param_gtilde.mat'));
                
                indtaxplan = cellfun(@(str) strncmp(taxplan, str, length(str)), {'base'; 'trump'; 'clinton'; 'ryan'});
                revperc = s.revenue_percent(indtaxplan,:);
                revperc = [revperc, revperc(end)*ones(1, max(T_model-length(revperc), 0))];
                
                if ~isbase
                    
                    % Calculate government expenditure adjustments
                    Dynamic_base = hardyload('dynamics.mat', base_generator, base_dir);
                    
                    Gtilde = Dynamic_base.Gtilde - gcut*s.GEXP_percent(1:T_model).*Dynamic_base.outs;
                    Ttilde = revperc.*Dynamic_base.outs - Static.pits - Static.ssts - Static.cits;
                    
                end
                
            case 'closed'
                
                % Identify open economy generator and save directory
                open_generator = @() dynamicSolver.open(basedef, counterdef, callingtag);
                open_dir = dirFinder.save('open', basedef, counterdef);
                
                % Load government expenditure adjustments
                Dynamic_open = hardyload('dynamics.mat', open_generator, open_dir);
                
                Gtilde = Dynamic_open.Gtilde;
                Ttilde = Dynamic_open.Ttilde;
                
        end
        
        
        % Define convergence tolerance and initialize error term
        tol = 1e-3;
        eps = Inf;
        
        % Initialize iteration count and set maximum number of iterations
        iter    =  0;
        itermax = 25;
        
        % Create file for logging iterations
        % (Note that this should be done after any parallel pool is started to avoid file access issues)
        iterlog = fopen(fullfile(save_dir, 'iterations.csv'), 'w');
        
        % Display header
        fprintf('\n[')
        switch economy
    
            case 'steady'
                fprintf('Steady state')
            
            case {'open', 'closed'}
                str1 = [upper(economy(1)), economy(2:end)];
                if isbase, str2 = 'baseline'; else, str2 = 'counterfactual'; end
                fprintf('%s economy %s', str1, str2)
                
        end
        fprintf(']\n')
        
        
        while (eps > tol && iter < itermax)
            
            % Increment iteration count
            iter = iter + 1;
            fprintf('\tIteration %2d  ...  ', iter)
            
            
            % Define prices
            switch economy
                
                case {'steady', 'closed'}
                    
                    if (iter == 1)
                        Market.rhos   = Market0.rhos  *ones(1,T_model);
                        Market.beqs   = Market0.beqs  *ones(1,T_model);
                        Market.assets = Market0.assets*ones(1,T_model);
                        Market.debts  = Market0.debts *ones(1,T_model);
                        Market.caps   = (Market.assets - Market.debts)/qtobin;
                    else
                        switch economy
                            case 'steady'
                                Market.rhos  = 0.5*Market.rhos + 0.5*Dynamic.rhos;
                                Market.debts = debtout*Dynamic.outs;
                            case 'closed'
                                Market.rhos  = 0.3*Market.rhos + 0.7*Dynamic.rhos;
                                Market.debts = Dynamic.debts;
                        end
                        Market.beqs   = Dynamic.beqs  ;
                        Market.assets = Dynamic.assets;
                        Market.caps   = Dynamic.caps  ;
                    end
                    
                    Market.capshares  = (Market.assets - Market.debts) ./ Market.assets;
                    Market.debtshares = 1 - Market.capshares;
                    Market.caprates   = (A*alp*(Market.rhos.^(alp-1)) - d)/qtobin;
                    switch economy
                        case 'steady', Market.govrates = cbomeanrate;
                        case 'closed', Market.govrates = cborates   ;
                    end
                    Market.totrates   = Market.capshares.*Market.caprates + Market.debtshares.*Market.govrates;
                    
                    
                case 'open'
                    
                    if (iter == 1)
                        
                        Market.assets = Market0.assets*ones(1,T_model);
                        Market.debts  = Market0.debts *ones(1,T_model);
                        Market.caps   = (Market.assets - Market.debts)/qtobin0;
                        
                        Market.capshares  = Market0.capshares *ones(1,T_model);
                        Market.debtshares = Market0.debtshares*ones(1,T_model);
                        Market.caprates   = Market0.caprates  *ones(1,T_model)*(1-taucap_base)/(1-taucap);
                        Market.govrates   = Market0.govrates  *ones(1,T_model);
                        Market.totrates   = Market0.totrates  *ones(1,T_model);
                        
                        Market.rhos   = ((qtobin*Market.caprates + d)/alp).^(1/(alp-1));
                        Market.beqs   = Market0.beqs  *ones(1,T_model);
                        
                    else
                        Market.beqs   = Dynamic.beqs;
                        Market.caps   = Dynamic.caps;
                    end
                    
            end
            
            Market.wages   = A*(1-alp)*(Market.rhos.^alp);
            Market.expsubs = [expshare * max(diff(Market.caps), 0), 0] ./ Market.caps;
            
            
            % Generate dynamic aggregates
            [Dynamic, LABs, DISTs] = generate_aggregates(Market, DISTs_steady, {}, {});
            
            
            % Calculate additional dynamic aggregates
            % (Note that open economy requires capital calculation before debt calculation while closed economy requires the reverse)
            switch economy
                
                case 'steady'
                    
                    % Calculate debt
                    Dynamic.debts = Market.debts;
                    
                    % Calculate capital and output
                    Dynamic.caps = (Dynamic.assets - Dynamic.debts)/qtobin;
                    Dynamic.outs = A*(max(Dynamic.caps, 0).^alp).*(Dynamic.labeffs.^(1-alp));
                    
                    % Calculate market clearing series
                    Dynamic.rhos = (max(Dynamic.assets - Dynamic.debts, 0)/qtobin) ./ Dynamic.labeffs;
                    clearing = Market.rhos - Dynamic.rhos;
                    
                case 'open'
                    
                    % Calculate capital and output
                    Dynamic.caps = Market.rhos .* Dynamic.labeffs;
                    Dynamic.outs = A*(max(Dynamic.caps, 0).^alp).*(Dynamic.labeffs.^(1-alp));
                    
                    Dynamic.caps_domestic = Market.capshares .* [Market0.assets, Dynamic.assets(1:end-1)];
                    Dynamic.caps_foreign  = qtobin*Dynamic.caps - Dynamic.caps_domestic;
                    
                    % Calculate debt
                    Dynamic.cits_domestic = Dynamic.cits;
                    Dynamic.cits_foreign  = taucap * Market.caprates * captaxshare .* Dynamic.caps_foreign;
                    Dynamic.cits          = Dynamic.cits_domestic + Dynamic.cits_foreign;
                    
                    if isbase
                        Gtilde = (revperc - fedgovtnis).*Dynamic.outs - Dynamic.bens;
                        Ttilde = revperc.*Dynamic.outs - Dynamic.pits - Dynamic.ssts - Dynamic.cits;
                    end
                    
                    Dynamic.revs  = Dynamic.pits + Dynamic.ssts + Dynamic.cits - Dynamic.bens;
                    Dynamic.debts = [Market0.debts, zeros(1,T_model-1)];
                    for t = 1:T_model-1
                        Dynamic.debts(t+1) = Gtilde(t) - Ttilde(t) - Dynamic.revs(t) + Dynamic.debts(t)*(1 + cborates(t));
                    end
                    Dynamic.Gtilde = Gtilde;
                    Dynamic.Ttilde = Ttilde;
                    
                    Dynamic.debts_domestic = (1 - Market.capshares) .* Dynamic.assets;
                    Dynamic.debts_foreign  = Dynamic.debts - Dynamic.debts_domestic;
                    
                    % Calculate income
                    Dynamic.labincs = Dynamic.labeffs .* Market.wages;
                    Dynamic.capincs = qtobin * Market.caprates .* Dynamic.caps;
                    
                    Dynamic.labpits = Dynamic.pits .* Dynamic.labincs ./ Dynamic.incs;
                    Dynamic.caprevs = Dynamic.cits + Dynamic.pits - Dynamic.labpits;
                    
                    % Calculate market clearing series
                    Dynamic.rhos = Market.rhos;
                    clearing = Market.beqs - Dynamic.beqs;
                    
                case 'closed'
                    
                    % Calculate debt
                    Dynamic.cits_domestic = Dynamic.cits;
                    Dynamic.cits_foreign  = zeros(1,T_model);
                    Dynamic.cits          = Dynamic.cits_domestic + Dynamic.cits_foreign;
                    
                    Dynamic.revs  = Dynamic.pits + Dynamic.ssts + Dynamic.cits - Dynamic.bens;
                    Dynamic.debts = [Market0.debts, zeros(1,T_model-1)];
                    for t = 1:T_model-1
                        Dynamic.debts(t+1) = Gtilde(t) - Ttilde(t) - Dynamic.revs(t) + Dynamic.debts(t)*(1 + cborates(t));
                    end
                    Dynamic.Gtilde = Gtilde;
                    Dynamic.Ttilde = Ttilde;
                    
                    Dynamic.debts_domestic = Dynamic.debts;
                    Dynamic.debts_foreign  = zeros(1,T_model);
                    
                    % Calculate capital and output
                    Dynamic.caps = ([(Market0.assets - Market0.debts)/qtobin0, (Dynamic.assets(1:end-1) - Dynamic.debts(2:end))/qtobin]);
                    Dynamic.outs = A*(max(Dynamic.caps, 0).^alp).*(Dynamic.labeffs.^(1-alp));
                    
                    Dynamic.caps_domestic = [qtobin0 * Dynamic.caps(1), qtobin * Dynamic.caps(2:end)];
                    Dynamic.caps_foreign  = zeros(1,T_model);
                    
                    % Calculate income
                    Dynamic.labincs = Dynamic.labeffs .* Market.wages;
                    Dynamic.capincs = qtobin * Market.caprates .* Dynamic.caps;
                    
                    Dynamic.labpits = Dynamic.pits .* Dynamic.labincs ./ Dynamic.incs;
                    Dynamic.caprevs = Dynamic.cits + Dynamic.pits - Dynamic.labpits;
                    
                    % Calculate market clearing series
                    Dynamic.rhos = (max([Market0.assets, Dynamic.assets(1:end-1)] - Dynamic.debts, 0)/qtobin) ./ Dynamic.labeffs;
                    clearing = Market.rhos - Dynamic.rhos;
                    
            end
            
            
            % Check market clearing series
            eps = max(abs(clearing));
            
            fprintf('Error term = %7.4f\n', eps)
            fprintf(iterlog, '%u,%0.4f\n', iter, eps);
            
            
        end
        fprintf('\n')
        fclose(iterlog);
        
        % Issue warning if maximum iterations reached
        if (iter == itermax), warning('Maximum iterations reached.'), end
        
        
        % Save optimal decision values and distributions for baseline
        if isbase
            save(fullfile(save_dir, 'decisions.mat'    ), 'LABs' )
            save(fullfile(save_dir, 'distributions.mat'), 'DISTs')
        end
        
        % Save final market conditions
        save(fullfile(save_dir, 'market.mat'), '-struct', 'Market')
        
        % Save dynamic aggregates
        switch economy
            case {'open', 'closed'}, save(fullfile(save_dir, 'dynamics.mat'), '-struct', 'Dynamic')
        end
        
        
        
        %% Elasticity calculation
        
        switch economy
            case 'steady'
                
                % Calculate capital to output ratio
                capout = (Dynamic.assets - Dynamic.debts) / Dynamic.outs;
                
                
                % Calculate labor elasticity
                workmass = 0;
                frisch   = 0;
                
                for jdem = 1:ndem
                    
                    LAB  = LABs {1,jdem};
                    DIST = DISTs{1,jdem};
                    
                    workind = (LAB > 0.01);
                    
                    workmass = workmass + sum(DIST(workind));
                    frisch   = frisch   + sum(DIST(workind).*(1-LAB(workind))./LAB(workind))*(1-gamma*(1-sigma))/sigma;
                    
                end
                
                labelas = frisch / workmass;
                
                
                % Calculate savings elasticity
                ratedev = 0.01;
                Market_dev = Market;
                
                Market_dev.caprates = Market.caprates * (1 + ratedev);
                Market_dev.govrates = Market.govrates * (1 + ratedev);
                Market_dev.totrates = Market.totrates * (1 + ratedev);
                
                [Dynamic_dev] = generate_aggregates(Market_dev, {}, {}, {});
                
                savelas = (Dynamic_dev.assets - Dynamic.assets) / (Dynamic.assets * ratedev);
                
                
                % Save and display elasticities
                save(fullfile(save_dir, 'elasticities.mat'), 'capout', 'labelas', 'savelas')
                
                elasticities = { 'Capital to output ratio' , capout  ;
                                 'Labor elasticity'        , labelas ;
                                 'Savings elasticity'      , savelas };
                for j = 1:size(elasticities, 1)
                    fprintf('\t%-25s= % 7.4f\n', elasticities{j,:})
                end
                fprintf('\n')
                
        end
        
        
    end
    
end

end




%%
% Check if file exists and generate before loading if necessary, handling parallel write and read conflicts.
% 
%%


function [s] = hardyload(filename, generator, save_dir)

% Check if file exists and generate if necessary
filepath = fullfile(save_dir, filename);
if ~exist(filepath, 'file')
    
    tagged_dir = generator();
    
    % Check for file again, possibly generated by a parallel run
    if ~exist(filepath, 'file')
        % Attempt to copy to save directory and continue if error encountered
        % (Errors are typically due to multiple simultaneous copy attempts, which should result in at least one successful copy)
        try copyfile(tagged_dir, save_dir); catch, end
    end
    
    % Clean up temporary directory
    rmdir(tagged_dir, 's')
    
end


% Turn on pause and store current pause state
pause0 = pause('on');

% Set maximum number of attempts and maximum pause time in seconds
maxtries = 200;
maxpause = 2.0;

for itry = 1:maxtries
    try
        % Attempt load
        s = load(filepath);
        break
    catch e
        if (itry == maxtries)
            error('Failed to load ''%s'' after %d attempts.\n%s', filepath, maxtries, e.message)
        end
        % Take a breather
        pause(rand * maxpause)
    end
end

% Reset pause state
pause(pause0)

end

