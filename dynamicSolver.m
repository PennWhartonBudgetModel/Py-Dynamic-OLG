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
        basedef_format    = struct( 'beta'          , '%.3f'    , ...
                                    'gamma'         , '%.3f'    , ...
                                    'sigma'         , '%.2f'    );
        
        counterdef_format = struct( 'taxplan'       , '%s'      , ...
                                    'gcut'          , '%+.2f'   , ...
                                    'legal_scale'   , '%.1f'    , ...
                                    'prem_legal'    , '%.3f'    , ...
                                    'amnesty'       , '%.2f'    , ...
                                    'deportation'   , '%.2f'    );
        
        % Define function to construct tag from definition and format specifications
        function tag = construct_tag(def, format)
            strs = {};
            for field = fields(def)'
                strs = [strs, {sprintf(format.(field{1}), def.(field{1}))}]; %#ok<AGROW>
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
        counterdef_filled = struct( 'taxplan'       , 'base', ...
                                    'gcut'          , +0.00 , ...
                                    'legal_scale'   , 1.0   , ...
                                    'prem_legal'    , 1.000 , ...
                                    'amnesty'       , 0.00  , ...
                                    'deportation'   , 0.00  );
        
        % Override default parameter values with values from counterfactual definition
        for field = fields(counterdef)'
            counterdef_filled.(field{1}) = counterdef.(field{1});
        end
        
    end
    
    
    % Solve dynamic model
    function [save_dir] = solve(economy, basedef, counterdef, callertag) %#ok<*FXUP>
        
        
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
        taxplan     = counterdef_filled.taxplan    ;
        gcut        = counterdef_filled.gcut       ;
        legal_scale = counterdef_filled.legal_scale;
        prem_legal  = counterdef_filled.prem_legal ;
        amnesty     = counterdef_filled.amnesty    ;
        deportation = counterdef_filled.deportation;
        
        
        % Identify working directories
        param_dir = dirFinder.param();
        [save_dir, ~, counterdef_tag] = dirFinder.save(economy, basedef, counterdef);
        
        % Append caller tag to save directory name and generate calling tag
        % (Obviates conflicts between parallel runs)
        save_dir   = [save_dir                                  , callertag];
        callingtag = [sprintf('^%s_%s', counterdef_tag, economy), callertag];
        
        % Clear or create save directory
        if exist(save_dir, 'dir'), rmdir(save_dir, 's'), end, mkdir(save_dir)
        
        
        
        %% Parameter definition
        
        % Define time constants
        T_life = 80;    % Total life years
        T_work = 47;    % Total working years
        
        switch economy
            case 'steady'
                T_model    = 1;                         % Steady state total modeling years
                startyears = 0;                         % Steady state cohort start year
            case {'open', 'closed'}
                T_model    = 24;                        % Transition path total modeling years
                startyears = (-T_life+1):(T_model-1);   % Transition path cohort start years
        end
        nstartyears = length(startyears);
        
        T_pasts   = max(-startyears, 0);                            % Life years before first model year
        T_shifts  = max(+startyears, 0);                            % Model years before first life year
        T_actives = min(startyears+T_life, T_model) - T_shifts;     % Life years within modeling period
        T_ends    = min(T_model-startyears, T_life);                % Maximum age within modeling period
        
        % Calculate productivities
        [nz, zs, transz, ndem] = calculate_productivities(T_life);
        
        % Define savings and average earnings discretization vectors
        % (Upper bound of average earnings defined as maximum possible Social Security benefit)
        f = @(lb, ub, n) lb + (ub-lb)*((0:n-1)/(n-1))'.^2;
        nk = 10; kv = f(1e-3, 120           , nk);
        nb =  5; bv = f(0   , 1.5*max(zs(:)), nb);
        
        % Define utility of bequests
        % (Currently defined to be zero for all savings levels)
        phi1 = 0; phi2 = 11.6; phi3 = 1.5;
        V_beq = phi1 * (1 + kv/phi2).^(1 - phi3);
        
        % Define output parameters
        A     = 1;      % Total factor productivity
        alpha = 0.45;   % Capital share of output
        d     = 0.085;  % Depreciation rate
        
        % Define constants relating model units to real units
        mpci = 3.25;    % Model per capita income            (To be updated; should be specific to each baseline)
        rpci = 1.19e5;  % Real per capita income in dollars  (To be updated; currently chosen to match percentage of individuals above Social Security maximum taxable earning level)
        
        % Define population growth parameters
        birth_rate   = 0.0200;                  % Annual birth rate
        legal_rate   = 0.0016 * legal_scale;    % Annual legal immigration rate
        illegal_rate = 0.0024;                  % Annual illegal immigration rate
        
        % Load age-dependent parameters
        % (To be updated; original sources outdated)
        s = load(fullfile(param_dir, 'age.mat'));
        surv         = s.surv;                  % Survival rates by age
        imm_age      = s.imm_age;               % New immigrant population distribution over age
        DISTz_age    = s.DISTz_age;             % Productivity distributions by age
        
        % Define population group index mapping
        groups = {'citizen', 'legal', 'illegal'};
        ng = length(groups);
        for ig = 1:ng, g.(groups{ig}) = ig; end
        
        % Shift legal immigrant productivity distribution according to policy parameter
        for age = 1:T_life
            v = mean(zs(:,age,:), 3);
            ztarget = sum(v .* DISTz_age(:,age,g.legal)) * prem_legal;
            p = (v(nz) - ztarget) / (v(nz)*(nz-1) - sum(v(1:nz-1)));
            DISTz_age(:,age,g.legal) = [p*ones(nz-1,1); 1 - p*(nz-1)];
        end
        
        % Define Social Security parameters
        ssthresholds = [856, 5157]*12*(mpci/rpci);  % Thresholds for earnings brackets
        ssrates      = [0.9, 0.32, 0.15];           % Marginal benefit rates for earnings brackets
        ss_scale     = 1.6;                         % Benefit scaling factor used to match total outlays as a percentage of GDP
        
        ssbenefit = [ max(min(bv, ssthresholds(1)) - 0              , 0) , ...
                      max(min(bv, ssthresholds(2)) - ssthresholds(1), 0) , ...
                      max(min(bv, Inf            ) - ssthresholds(2), 0) ] * ssrates' * ss_scale;
        
        ssbenefits  = repmat(ssbenefit          , [1,T_model]);  % Benefits
        sstaxs      = repmat(0.124              , [1,T_model]);  % Tax rates
        ssincmaxs   = repmat(1.185e5*(mpci/rpci), [1,T_model]);  % Maximum taxable earnings
        
        sstaxcredit = 0.15;     % Benefit tax credit percentage
        
        % Load CBO parameters
        % (Values originally extracted from PWBM_GrowthAdjustedBudget.xlsx)
        s = load(fullfile(param_dir, 'PWBM_GrowthAdjustedBudget.mat'));
        
        debttoout   = s.FederalDebtHeldbythePublic(1)/100;  % Debt-to-output ratio
        fedgovtnis  = s.fedgovtnis(1:T_model);              % Federal government noninterest surplus
        cborates    = s.r_cbo(1:T_model);                   % Growth-adjusted average interest rates
        cbomeanrate = mean(s.r_cbo);                        % Mean growth-adjusted interest rate over modeling period
        
        % Load tax parameters
        % (Generated by generate_tax() using data from TPC_Inputs.xlsx)
        s = load(fullfile(param_dir, 'tax.mat'));
        
        pit_coefs   = s.(taxplan).pit_coefs;
        deduc_coefs = [s.(taxplan).avg_deduc, s.(taxplan).coefs(1:2)];
        
        captaxshare = s.(taxplan).cap_tax_share;
        expshare    = s.(taxplan).exp_share;
        taucap_base = s.('base' ).tau_cap;
        taucap      = s.(taxplan).tau_cap;
        taucapgain  = s.(taxplan).tau_capgain;
        
        qtobin0 = 1 - s.('base' ).tau_cap*s.('base' ).exp_share;
        qtobin  = 1 - s.(taxplan).tau_cap*s.(taxplan).exp_share;
        
        
        
        %% Aggregate generation function
        
        function [Aggregate, LABs, DIST, pgr] = generate_aggregates(Market, DIST_steady, LABs_static, DIST_static)
            
            % Define dynamic aggregate generation flag
            isdynamic = isempty(LABs_static) || isempty(DIST_static);
            
            % Set static optimal decision values to empty values for dynamic aggregate generation
            if isdynamic, LABs_static = cell(nstartyears, ndem); end
            
            % Initialize optimal decision value arrays
            os = {'K', 'LAB', 'B', 'INC', 'PIT', 'SST', 'CIT', 'BEN'};
            for o = os, OPTs.(o{1}) = zeros(nz,nk,nb,T_life,T_model,ndem); end
            
            % Initialize array of cohort optimal labor values
            LABs = cell(nstartyears, ndem);
            
            % Initialize population distribution array
            DIST = zeros(nz,nk,nb,T_life,ng,T_model,ndem);
            
            
            for idem = 1:ndem
                
                % Package fixed dynamic optimization arguments into anonymous function
                solve_cohort_ = @(V0, LAB_static, T_past, T_shift, T_active) solve_cohort(V0, LAB_static, isdynamic, ...
                    nz, nk, nb, T_past, T_shift, T_active, T_work, T_model, zs(:,:,idem), transz, kv, bv, beta, gamma, sigma, surv, V_beq, ...
                    mpci, rpci, sstaxcredit, ssbenefits, sstaxs, ssincmaxs, deduc_coefs, pit_coefs, captaxshare, taucap, taucapgain, qtobin, qtobin0, ...
                    Market.beqs, Market.wages, Market.capshares, Market.caprates, Market.govrates, Market.totrates, Market.expsubs);
                
                
                % Initialize series of terminal utility values
                V0s = zeros(nz,nk,nb,T_life);
                
                % Solve steady state / post-transition path cohort
                if isdynamic
                    
                    % Solve dynamic optimization
                    % (Note that active time is set to full lifetime)
                    [V, OPT] = solve_cohort_(V0s(:,:,:,T_life), [], T_pasts(end), T_shifts(end), T_life);
                    
                    % Define series of terminal utility values
                    V0s(:,:,:,1:T_life-1) = V(:,:,:,2:T_life);
                    
                end
                
                
                switch economy
                    
                    case 'steady'
                        
                        % Store optimal decision values
                        for o = os, OPTs.(o{1})(:,:,:,:,1,idem) = OPT.(o{1}); end
                        LABs{1,idem} = OPT.LAB;
                        
                    case {'open', 'closed'}
                        
                        % Solve transition path cohorts
                        OPTs_cohort = cell(1,nstartyears);
                        
                        parfor i = 1:nstartyears
                            
                            % Extract terminal utility values
                            V0 = V0s(:,:,:,T_ends(i)); %#ok<PFBNS>
                            
                            % Solve dynamic optimization
                            [~, OPTs_cohort{i}] = solve_cohort_(V0, LABs_static{i,idem}, T_pasts(i), T_shifts(i), T_actives(i));
                            
                            LABs{i,idem} = OPTs_cohort{i}.LAB;
                            
                        end
                        
                        % Construct optimal decision value arrays
                        for i = 1:nstartyears
                            for t = 1:T_actives(i)
                                age  = t + T_pasts (i);
                                year = t + T_shifts(i);
                                for o = os, OPTs.(o{1})(:,:,:,age,year,idem) = OPTs_cohort{i}.(o{1})(:,:,:,t); end
                            end
                        end
                        
                end
                
                
                if isdynamic
                    
                    % Define initial population distribution and distribution generation termination conditions
                    switch economy
                        
                        case 'steady'
                            DIST_next = ones(nz,nk,nb,T_life,ng) / (nz*nk*nb*T_life*ng*ndem);
                            lastyear = Inf;
                            disttol = 1e-6;
                            
                        case {'open', 'closed'}
                            DIST_next = DIST_steady(:,:,:,:,:,1,idem);
                            lastyear = T_model;
                            disttol = -Inf;
                            
                    end
                    
                    year = 1;
                    disteps = Inf;
                    
                    while (disteps > disttol && year <= lastyear)
                        
                        % Store population distribution for current year
                        DIST_year = DIST_next;
                        DIST(:,:,:,:,:,min(year, T_model),idem) = DIST_year;
                        
                        % Extract optimal decision values for current year
                        K = OPTs.K(:,:,:,:,min(year, T_model),idem);
                        B = OPTs.B(:,:,:,:,min(year, T_model),idem);
                        
                        % Define population growth distribution
                        DIST_grow = zeros(nz,nk,nb,T_life,ng);
                        P = sum(DIST_year(:));
                        
                        DIST_grow(:,1,1,1,g.citizen) = reshape(DISTz_age(:,1,g.citizen), [nz,1,1,1     ,1]) * P * birth_rate   ;
                        DIST_grow(:,1,1,:,g.legal  ) = reshape(DISTz_age(:,:,g.legal  ), [nz,1,1,T_life,1]) * P * legal_rate   .* repmat(reshape(imm_age, [1,1,1,T_life,1]), [nz,1,1,1,1]);
                        DIST_grow(:,1,1,:,g.illegal) = reshape(DISTz_age(:,:,g.illegal), [nz,1,1,T_life,1]) * P * illegal_rate .* repmat(reshape(imm_age, [1,1,1,T_life,1]), [nz,1,1,1,1]);
                        
                        % Generate population distribution for next year
                        DIST_next = generate_distribution(DIST_year, DIST_grow, K, B, nz, nk, nb, T_life, ng, transz, kv, bv, surv);
                        
                        % Increase legal immigrant population for amnesty, maintaining distributions over productivity
                        DISTz_legal = DIST_next(:,:,:,:,g.legal) ./ repmat(sum(DIST_next(:,:,:,:,g.legal), 1), [nz,1,1,1,1]);
                        DISTz_legal(isnan(DISTz_legal)) = 1/nz;
                        
                        DIST_next(:,:,:,:,g.legal) = DIST_next(:,:,:,:,g.legal) + repmat(sum(amnesty*DIST_next(:,:,:,:,g.illegal), 1), [nz,1,1,1,1]).*DISTz_legal;
                        
                        % Reduce illegal immigrant population for amnesty and deportation
                        DIST_next(:,:,:,:,g.illegal) = (1-amnesty-deportation)*DIST_next(:,:,:,:,g.illegal);
                        
                        % Calculate age distribution convergence error
                        f = @(D) sum(sum(reshape(D, [], T_life, ng), 1), 3) / sum(D(:));
                        disteps = max(abs(f(DIST_next) - f(DIST_year)));
                        
                        % Calculate net population growth rate
                        % (Note that rate is assumed to be independent of demographic type)
                        pgr = sum(DIST_next(:)) / sum(DIST_year(:));
                        
                        year = year + 1;
                        
                    end
                    
                else
                    DIST = DIST_static;
                end
                
            end
            
            
            % Normalize steady state population distribution
            switch economy, case 'steady', DIST = DIST / sum(DIST(:)); end
            
            
            % Generate aggregates
            DIST_gs = reshape(sum(DIST, 5), [nz,nk,nb,T_life,T_model,ndem]);
            f = @(F) sum(sum(reshape(DIST_gs .* F, [], T_model, ndem), 1), 3);
            
            Aggregate.pops     = f(1);                                                                                   % Population
            Aggregate.assets   = f(repmat(reshape(kv, [1,nk,1,1,1,1]), [nz,1,nb,T_life,T_model,ndem]));                  % Assets
            Aggregate.bequests = f(OPTs.K .* repmat(reshape(1-surv, [1,1,1,T_life,1,1]), [nz,nk,nb,1,T_model,ndem]));    % Bequests
            Aggregate.labs     = f(OPTs.LAB);                                                                            % Labor
            Aggregate.labeffs  = f(OPTs.LAB .* repmat(reshape(zs, [nz,1,1,T_life,1,ndem]), [1,nk,nb,1,T_model,1]));      % Effective labor
            Aggregate.lfprs    = f(OPTs.LAB > 0.01) ./ f(1);                                                             % Labor force participation rate
            Aggregate.incs     = f(OPTs.INC);                                                                            % Income
            Aggregate.pits     = f(OPTs.PIT);                                                                            % Personal income tax
            Aggregate.ssts     = f(OPTs.SST);                                                                            % Social Security tax
            Aggregate.cits     = f(OPTs.CIT);                                                                            % Capital income tax
            Aggregate.bens     = f(OPTs.BEN);                                                                            % Social Security benefits
            
            
        end
        
        
        
        %% Static aggregate generation
        
        if ~isbase
            
            % Identify baseline generator and save directory
            base_generator = @() dynamicSolver.solve(economy, basedef, [], callingtag);
            base_dir = dirFinder.save(economy, basedef);
            
            % Load baseline market conditions, optimal labor values, and population distribution
            Market = hardyload('market.mat'      , base_generator, base_dir);
            
            s      = hardyload('decisions.mat'   , base_generator, base_dir);
            LABs_static = s.LABs;
            
            s      = hardyload('distribution.mat', base_generator, base_dir);
            DIST_static = s.DIST;
            
            
            % Generate static aggregates
            % (Intermediary structure used to filter out extraneous fields)
            [Static_] = generate_aggregates(Market, {}, LABs_static, DIST_static);
            
            for series = {'incs', 'pits', 'ssts', 'cits', 'bens'}
                Static.(series{1}) = Static_.(series{1});
            end
            
            
            % Copy additional static aggregates from baseline aggregates
            Dynamic_base = hardyload('dynamics.mat', base_generator, base_dir);
            
            for series = {'labeffs', 'caps', 'lfprs', 'labincs', 'capincs', 'outs', 'caps_domestic', 'caps_foreign', 'debts_domestic', 'debts_foreign'}
                Static.(series{1}) = Dynamic_base.(series{1});
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
                
                % Load initial conditions
                Market0 = load(fullfile(param_dir, 'market0.mat'));
                
                DIST_steady = {};
                
            case {'open', 'closed'}
                
                % Identify steady state generator and save directory
                steady_generator = @() dynamicSolver.steady(basedef, callingtag);
                steady_dir = dirFinder.save('steady', basedef);
                
                % Load steady state market conditions and dynamic aggregates
                Market0  = hardyload('market.mat'      , steady_generator, steady_dir);
                Dynamic0 = hardyload('dynamics.mat'    , steady_generator, steady_dir);
                
                % Load steady state population distribution
                s        = hardyload('distribution.mat', steady_generator, steady_dir);
                
                DIST_steady = s.DIST;
                
        end
        
        
        switch economy
            
            case 'open'
                
                % Load government expenditure adjustment parameters
                % (revenue_percent originally extracted from RevenuesPctGDP.xlsx)
                % (GEXP_percent originally derived from data in LTBO2015.xlsx -- Total Non-Interest Spending - Social Security - Medicare)
                s = load(fullfile(param_dir, 'govexp.mat'));
                
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
        
        
        % Define marketing clearing tolerance and initialize error term
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
            isinitial = iter == 1;
            
            
            % Define market conditions
            if isinitial
                Market.beqs      = Market0.beqs*ones(1,T_model);
                Market.capshares = Market0.capshares*ones(1,T_model);
                Market.expsubs   = zeros(1,T_model);
            else
                Market.beqs      = beqs;
                Market.expsubs   = [expshare * max(diff(Dynamic.caps), 0), 0] ./ Dynamic.caps;
            end
            
            switch economy
                
                case {'steady', 'closed'}
                    
                    if isinitial
                        Market.rhos      = Market0.rhos*ones(1,T_model);
                        switch economy
                            case 'steady', Market.govrates = cbomeanrate;
                            case 'closed', Market.govrates = cborates;
                        end
                    else
                        rhostep = 0.5;
                        Market.rhos      = rhostep*rhos + (1-rhostep)*Market.rhos;
                        Market.capshares = (Dynamic.assets - Dynamic.debts) ./ Dynamic.assets;
                    end
                    
                    Market.caprates = (A*alpha*(Market.rhos.^(alpha-1)) - d)/qtobin;
                    Market.totrates = Market.capshares.*Market.caprates + (1-Market.capshares).*Market.govrates;
                    
                case 'open'
                    
                    if isinitial
                        Market.caprates  = Market0.caprates*ones(1,T_model)*(1-taucap_base)/(1-taucap);
                        Market.govrates  = Market0.govrates*ones(1,T_model);
                        Market.totrates  = Market0.totrates*ones(1,T_model);
                        Market.rhos      = ((qtobin*Market.caprates + d)/alpha).^(1/(alpha-1));
                    end
                    
            end
            
            Market.wages = A*(1-alpha)*(Market.rhos.^alpha);
            
            
            % Generate dynamic aggregates
            [Dynamic, LABs, DIST, pgr] = generate_aggregates(Market, DIST_steady, {}, {});
            
            
            % Calculate additional dynamic aggregates
            % (Note that open economy requires capital calculation before debt calculation while closed economy requires the reverse)
            switch economy
                
                case 'steady'
                    
                    % Calculate debt, capital, and output
                    % (Numerical solver used due to absence of closed form solution)
                    f_debts = @(outs ) debttoout*outs;
                    f_caps  = @(debts) (Dynamic.assets - debts)/qtobin;
                    f_outs  = @(caps ) A*(max(caps, 0).^alpha).*(Dynamic.labeffs.^(1-alpha));
                    x_ = fsolve(@(x) x - [f_debts(x(3)); f_caps(x(1)); f_outs(x(2))], zeros(3,1), optimoptions('fsolve', 'Display', 'none'));
                    Dynamic.debts = x_(1);
                    Dynamic.caps  = x_(2);
                    Dynamic.outs  = x_(3);
                    
                    % Calculate market clearing series
                    rhos = max(Dynamic.caps, 0) / Dynamic.labeffs;
                    beqs = Dynamic.bequests / pgr;
                    clearing = Market.rhos - rhos;
                    
                case 'open'
                    
                    % Calculate capital and output
                    Dynamic.caps = Market.rhos .* Dynamic.labeffs;
                    Dynamic.outs = A*(max(Dynamic.caps, 0).^alpha).*(Dynamic.labeffs.^(1-alpha));
                    
                    Dynamic.caps_domestic = Market.capshares .* [Dynamic0.assets, Dynamic.assets(1:T_model-1)];
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
                    Dynamic.debts = [Dynamic0.debts, zeros(1,T_model-1)];
                    for year = 1:T_model-1
                        Dynamic.debts(year+1) = Gtilde(year) - Ttilde(year) - Dynamic.revs(year) + Dynamic.debts(year)*(1 + cborates(year));
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
                    beqs = [Market0.beqs, Dynamic.bequests(1:T_model-1) ./ Dynamic.pops(2:T_model)];
                    clearing = Market.beqs - beqs;
                    
                case 'closed'
                    
                    % Calculate debt
                    Dynamic.cits_domestic = Dynamic.cits;
                    Dynamic.cits_foreign  = zeros(1,T_model);
                    Dynamic.cits          = Dynamic.cits_domestic + Dynamic.cits_foreign;
                    
                    Dynamic.revs  = Dynamic.pits + Dynamic.ssts + Dynamic.cits - Dynamic.bens;
                    Dynamic.debts = [Dynamic0.debts, zeros(1,T_model-1)];
                    for year = 1:T_model-1
                        Dynamic.debts(year+1) = Gtilde(year) - Ttilde(year) - Dynamic.revs(year) + Dynamic.debts(year)*(1 + cborates(year));
                    end
                    Dynamic.Gtilde = Gtilde;
                    Dynamic.Ttilde = Ttilde;
                    
                    Dynamic.debts_domestic = Dynamic.debts;
                    Dynamic.debts_foreign  = zeros(1,T_model);
                    
                    % Calculate capital and output
                    Dynamic.caps = ([(Dynamic0.assets - Dynamic0.debts)/qtobin0, (Dynamic.assets(1:T_model-1) - Dynamic.debts(2:T_model))/qtobin]);
                    Dynamic.outs = A*(max(Dynamic.caps, 0).^alpha).*(Dynamic.labeffs.^(1-alpha));
                    
                    Dynamic.caps_domestic = [qtobin0 * Dynamic.caps(1), qtobin * Dynamic.caps(2:T_model)];
                    Dynamic.caps_foreign  = zeros(1,T_model);
                    
                    % Calculate income
                    Dynamic.labincs = Dynamic.labeffs .* Market.wages;
                    Dynamic.capincs = qtobin * Market.caprates .* Dynamic.caps;
                    
                    Dynamic.labpits = Dynamic.pits .* Dynamic.labincs ./ Dynamic.incs;
                    Dynamic.caprevs = Dynamic.cits + Dynamic.pits - Dynamic.labpits;
                    
                    % Calculate market clearing series
                    rhos = (max([Dynamic0.assets, Dynamic.assets(1:end-1)] - Dynamic.debts, 0)/qtobin) ./ Dynamic.labeffs;
                    beqs = [Market0.beqs, Dynamic.bequests(1:T_model-1) ./ Dynamic.pops(2:T_model)];
                    clearing = Market.rhos - rhos;
                    
            end
            
            
            % Calculate maximum error in market clearing series
            eps = max(abs(clearing));
            
            fprintf('Error term = %7.4f\n', eps)
            fprintf(iterlog, '%u,%0.4f\n', iter, eps);
            
            
        end
        fprintf('\n')
        fclose(iterlog);
        
        % Issue warning if maximum iterations reached
        if (iter == itermax), warning('Maximum iterations reached'), end
        
        
        % Save baseline optimal labor values and population distribution
        if isbase
            save(fullfile(save_dir, 'decisions.mat'   ), 'LABs')
            save(fullfile(save_dir, 'distribution.mat'), 'DIST')
        end
        
        % Save market conditions and dynamic aggregates
        save(fullfile(save_dir, 'market.mat'  ), '-struct', 'Market' )
        save(fullfile(save_dir, 'dynamics.mat'), '-struct', 'Dynamic')
        
        
        
        %% Elasticity calculation
        
        switch economy
            case 'steady'
                
                % Calculate capital to output ratio
                captoout = (Dynamic.assets - Dynamic.debts) / Dynamic.outs;
                
                
                % Calculate labor elasticity
                workmass = 0;
                frisch   = 0;
                
                for idem = 1:ndem
                    
                    LAB_idem  = LABs{1,idem};
                    DIST_idem = sum(DIST(:,:,:,:,:,1,idem), 5);
                    
                    workind = (LAB_idem > 0.01);
                    
                    workmass = workmass + sum(DIST_idem(workind));
                    frisch   = frisch   + sum(DIST_idem(workind) .* (1 - LAB_idem(workind)) ./ LAB_idem(workind)) * (1 - gamma*(1-sigma))/sigma;
                    
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
                save(fullfile(save_dir, 'elasticities.mat'), 'captoout', 'labelas', 'savelas')
                
                for label = { {'Capital to output ratio' , captoout } , ...
                              {'Labor elasticity'        , labelas  } , ...
                              {'Savings elasticity'      , savelas  } }
                    fprintf('\t%-25s= % 7.4f\n', label{1}{:})
                end
                fprintf('\n')
                
        end
        
        
    end
    
end

end




%%
% Calculate productivities.
% 
%%


function [nz, zs, transz, ndem] = calculate_productivities(T_life)
    
    function [y] = generate_normdist_shocks(nt,sigt)
        
        % This function gives the cutoff points for nt slices of a normal
        % distribution with equal weight.
        
        epsgrid = zeros(nt+1,1);    % creates grid for endpoints
        
        for i = 1:nt+1
            epsgrid(i) = sigt*norminv((i-1)/nt,0,1);
        end
        
        y = zeros(nt,1);    % grid of conditional means
        for i = 1:nt
            e1 = (epsgrid(i))/sigt;
            e2 = (epsgrid(i+1))/sigt;
            y(i) = nt*sigt*(normpdf(e1,0,1) - normpdf(e2,0,1));
        end
        
    end
    
    function [y, yprob, epsgrid] = generate_persistent_shocks(ny,const,lambda,sigep)
        
        % This function creates a Markov matrix that approximates an AR(1) process.  Written for C-code generation.
        % ny = gridsize
        % const = constant on ar1 process
        % lambda = coefficient on ar1 process
        % sigep = standard deviation of error terms
        % y(t) = mu*(1-lambda) + lambda(t-1) + eps
        % recovering mu
        
        mu = const/(1-lambda);
        sigy = sigep/(sqrt(1-lambda^2));
        epsgrid = zeros(ny+1,1);    % creates grid for endpoints
        
        for i = 1:ny+1    
            epsgrid(i) = sigy*norminv((i-1)/ny,0,1) + mu;
        end
        
        y = zeros(ny,1);    % grid of conditional means
        for i = 1:ny    
            e1 = (epsgrid(i)-mu)/sigy;
            e2 = (epsgrid(i+1)-mu)/sigy;
            y(i) = ny*sigy*(normpdf(e1,0,1) - normpdf(e2,0,1)) + mu;
        end
        
        yprob = zeros(ny,ny);
        for i = 1:ny
            for i2 = 1:ny
                ei1 = epsgrid(i);
                ei2 = epsgrid(i+1);
                ej1 = epsgrid(i2);
                ej2 = epsgrid(i2+1);
                yprob(i,i2) = (ny/(sqrt(2*pi*(sigy^2))))*quadgk(@(x) integrand(x,ej1,ej2,sigy,sigep,mu,lambda), ei1, ei2);
            end
        end
        
    end
    
    function [c2] = integrand(x,ej1,ej2,sigy,sigep,mu,lambda)
        
        nx = length(x);
        c2 = zeros(1,nx);
        
        for i1 = 1:nx
            term1 = exp(-((x(i1)-mu)^2)/(2*sigy^2));
            term2 = normcdf((ej2-mu*(1-lambda)-lambda*x(i1))/sigep,0,1);
            term3 = normcdf((ej1-mu*(1-lambda)-lambda*x(i1))/sigep,0,1);
            c2(i1) = term1*(term2-term3);
        end

    end
    
    
    % (Storesletten, Telmer, and Yaron, 2004)
    
    % Permanent shock
    np = 2;  % Number of permanent shocks
    sig_p = 0.2105^0.5;    % Standard deviation
    z_perm = generate_normdist_shocks(np,sig_p);
    
    % Persistent shock
    nz = 2;
    sigep_1 = 0.124^0.5;
    const = 0;
    pep = 0.973;    % coefficient on lagged productivity 
    sigep = 0.018^0.5;    % conditional stdev of productivity 
    [z1, tr_z1, epsgrid] = generate_persistent_shocks(nz,const,pep,sigep);    % returns the transition matrix
    proddist = zeros(1,nz);
    for i = 1:nz
        proddist(i) = normcdf(epsgrid(i+1),0,sigep_1) - normcdf(epsgrid(i),0,sigep_1);
    end
    
    % Transitory shock 
    nt = 2;  % Number of transitory shocks
    sig_t = (.0630)^(.5);
    z_trans = generate_normdist_shocks(nt,sig_t); 
    tr_t = (1/nt).*ones(nt);  % Initial distribution over transitory shocks.
    
    transz = kron(tr_z1,tr_t);    % Provides final productivity Markov transition matrix.
    
    
    zs = zeros(nz*nt, T_life, np);
    % Deterministic component of lifecycle productivity.  Estimates are from Barro and Barnes Medicaid working paper.
    for t1 = 1:T_life
        for d1 = 1:np
            shock_node = 0;
            for idio_shock = 1:nz
                for trans_shock = 1:nt
                    shock_node = shock_node+1;
                    zs(shock_node,t1,d1) = (t1+19)* 0.2891379 -0.0081467*(t1+19)^2 +  0.000105*(t1+19)^3 -5.25E-07*(t1+19)^4 -1.203521 + z1(idio_shock) + z_trans(trans_shock) + z_perm(d1);
                end
            end
        end
    end
    
    zs = max(zs, 0);
    nz = nz*nt;  % reassigning nz so that we don't have to make any other adjustments in the code.
    
    ndem = np;    % Redefining "np" as "ndem" to reflect patch for the creation of demographic types.
    
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

