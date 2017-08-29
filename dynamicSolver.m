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
        basedef_format    = struct( 'beta'              , '%.3f'    , ...
                                    'gamma'             , '%.3f'    , ...
                                    'sigma'             , '%.2f'    , ...
                                    'modelunit_dollars' , '%e'       );
        
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
        beta                = basedef.beta ;
        gamma               = basedef.gamma;
        sigma               = basedef.sigma;
        modelunit_dollars   = basedef.modelunit_dollars;
        
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
        
        % map model inputs (or outputs) to actual years
        first_transition_year  = 2018;
        
        % Identify working directories
        param_dir = dirFinder.param();
        [save_dir, ~, counterdef_tag] = dirFinder.save(economy, basedef, counterdef);
        
        % Append caller tag to save directory name and generate calling tag
        % (Obviates conflicts between parallel runs)
        save_dir   = [save_dir                                  , callertag];
        callingtag = [sprintf('^%s_%s', counterdef_tag, economy), callertag];
        
        % Clear or create save directory
        if exist(save_dir, 'dir'), rmdir(save_dir, 's'), end, mkdir(save_dir)
        
        
        
        %% PARAMETERS
        
        % Define time constants
        s       = paramGenerator.timing();
        T_life  = s.T_life;    % Total life years
        T_work  = s.T_work;    % Total working years
        T_model = s.T_model;   % Transition path model years
        switch economy
            case 'steady'
                T_model    = 1;                         % Steady state total modeling years
                startyears = 0;                         % Steady state cohort start year
            case {'open', 'closed'}
                T_model    = T_model;                   % Transition path total modeling years
                startyears = (-T_life+1):(T_model-1);   % Transition path cohort start years
        end
        nstartyears = length(startyears);
        
        T_pasts   = max(-startyears, 0);                            % Life years before first model year
        T_shifts  = max(+startyears, 0);                            % Model years before first life year
        T_actives = min(startyears+T_life, T_model) - T_shifts;     % Life years within modeling period
        T_ends    = min(T_model-startyears, T_life);                % Maximum age within modeling period
        
        
        % Discretized grids, including shock process
        %   ndem is number of permanent types (hi, low)
        %   g    is population subgroup set: 'citizen',
        %           'legal', 'illegal'
        %   NOTE: DISTz and zs are indexed by g
        s      = paramGenerator.grids( T_life, prem_legal );
        ndem   = s.ndem;      % demographic types
        g      = s.g;         % groups: citizen, legal, illegal 
        ng     = s.ng;        % num groups
        nz     = s.nz;        % num labor productivity shocks
        transz = s.transz;    % transition matrix of z's
        DISTz  = s.DISTz;     % steady-state distribution of z's
        zs     = s.zs;        % shocks grid (by demographic type and age)
        nk     = s.nk;        % num asset points
        kv     = s.kv;        % assets grid
        nb     = s.nb;        % num avg. earnings points
        bv     = s.bv;        % avg. earnings grid
        
        % Load production parameters
        s       = paramGenerator.production();
        A       = s.A;        % Total factor productivity
        alpha   = s.alpha;    % Capital share of output
        d       = s.d;        % Depreciation rate
        
        % Load population growth parameters
        % Load age-dependent parameters
        % (To be updated; original sources outdated)        s       = paramGenerator.demographics();
        s = paramGenerator.demographics();
        birth_rate      = s.birth_rate;                 % Annual birth rate
        legal_rate      = s.legal_rate * legal_scale;   % Annual legal immigration rate
        illegal_rate    = s.illegal_rate;               % Annual illegal immigration rate
        surv            = s.surv;                       % Survival probabilities by age
        imm_age         = s.imm_age;                    % Immigrants' age distribution
        
        % Load Social Security parameters
        s           = paramGenerator.social_security( modelunit_dollars, bv, T_model );
        ssbenefits  = s.ssbenefits ;    % Benefits
        sstaxs      = s.sstaxs     ;    % Tax rates
        ssincmaxs   = s.ssincmaxs  ;    % Maximum taxable earnings
        sstaxcredit = s.sstaxcredit;    % Benefit tax credit percentage
        
        
        %%  CBO interest rates, expenditures, and debt
        s               = paramGenerator.budget( first_transition_year, T_model );
        GEXP_by_GDP     = s.GEXP_by_GDP;    % Gvt expenditures as pct GDP
        debt            = s.debt;           % Gvt debt as pct gdp
        debttoout       = s.debttoout;      % Initial debt/gdp (for steady state)
        fedgovtnis      = s.fedgovtnis;     % Gvt net interest surplus (deficit)
        cborates        = s.cborates;       % Interest rates on gvt debt (from CBO)
        cbomeanrate     = s.cbomeanrate;    % Avg of cborates (for steady state)
        
        %% Tax parameters
        s                   = paramGenerator.tax( taxplan );
        tax_thresholds      = s.tax_thresholds; % Tax func is linearized, these are income thresholds 
        tax_burden          = s.tax_burden;     % Tax burden (cumulative tax) at thresholds
        tax_rates           = s.tax_rates;      % Effective marginal tax rate between thresholds
        
        captaxshare         = s.captaxshare;    % Portion of capital income taxed at preferred rates
        expshare            = s.expshare;       % Portion of investment which can be expensed
        taucap              = s.taucap;         % Capital tax rate
        taucapgain          = s.taucapgain;     % Capital gains tax rate
        
        s_base = paramGenerator.tax( 'base' );
        expshare_base = s_base.expshare;        % expshare for baseline
        taucap_base   = s_base.taucap;          % taucap for baseline
        
        qtobin0 = 1 - expshare_base*taucap_base;
        qtobin  = 1 - expshare     *taucap     ;

        % Define utility of bequests
        % (Currently defined to be zero for all savings levels)
        phi1 = 0; phi2 = 11.6; phi3 = 1.5;
        V_beq = phi1 * (1 + kv/phi2).^(1 - phi3);

        
        
        %% Aggregate generation function
        
        function [Aggregate, LABs, DIST, pgr, policy_fncs] = generate_aggregates(Market, DIST_steady, LABs_static, DIST_static)
            
            % Define dynamic aggregate generation flag
            isdynamic = isempty(LABs_static) || isempty(DIST_static);
            
            % Set static optimal decision values to empty values for dynamic aggregate generation
            if isdynamic, LABs_static = cell(nstartyears, ndem); end
            
            % Initialize optimal decision value arrays
            os = {'K', 'LAB', 'B', 'INC', 'PIT', 'SST', 'CIT', 'BEN', 'CON'};
            for o = os, OPTs.(o{1}) = zeros(nz,nk,nb,T_life,T_model,ndem); end
            
            % Initialize array of cohort optimal labor values
            LABs = cell(nstartyears, ndem);
            
            % Initialize population distribution array
            DIST = zeros(nz,nk,nb,T_life,ng,T_model,ndem);
            
            
            for idem = 1:ndem
                
                % Package fixed dynamic optimization arguments into anonymous function
                solve_cohort_ = @(V0, LAB_static, T_past, T_shift, T_active) solve_cohort(V0, LAB_static, isdynamic, ...
                    nz, nk, nb, T_past, T_shift, T_active, T_work, T_model, zs(:,:,idem), transz, kv, bv, beta, gamma, sigma, surv, V_beq, ...
                    modelunit_dollars, ...
                    sstaxcredit, ssbenefits, sstaxs, ssincmaxs, ...
                    tax_thresholds, tax_burden, tax_rates, ... 
                    captaxshare, taucap, taucapgain, qtobin, qtobin0, ...
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
                        
                        DIST_grow(:,1,1,1,g.citizen) = reshape(DISTz(:,1,g.citizen), [nz,1,1,1     ,1]) * P * birth_rate   ;
                        DIST_grow(:,1,1,:,g.legal  ) = reshape(DISTz(:,:,g.legal  ), [nz,1,1,T_life,1]) * P * legal_rate   .* repmat(reshape(imm_age, [1,1,1,T_life,1]), [nz,1,1,1,1]);
                        DIST_grow(:,1,1,:,g.illegal) = reshape(DISTz(:,:,g.illegal), [nz,1,1,T_life,1]) * P * illegal_rate .* repmat(reshape(imm_age, [1,1,1,T_life,1]), [nz,1,1,1,1]);
                        
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
            Aggregate.cons     = f(OPTs.CON);                                                                            % Consumption
            
            policy_fncs = OPTs;
            
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
                % Tax revenues as fraction of GDP are loaded from
                % single-series CSV files which contain data from TPC by
                % tax plan (base, trumpA, trumpB)
                % Input: TPCRevenues_<taxplan>.csv -- TPC estimated of tax
                %           revenues as percent GDP
                %           Format is (Year), (PctRevenues) w/ header row.
                filename            = strcat('TPCRevenues_', taxplan, '.csv');
                tax_revenue_by_GDP  = read_series(filename, first_transition_year, dirFinder.param);
                tax_revenue_by_GDP  = tax_revenue_by_GDP'; 
                if( T_model - length(tax_revenue_by_GDP) < 0 )
                    tax_revenue_by_GDP = tax_revenue_by_GDP(1:T_model);
                else
                    tax_revenue_by_GDP  = [tax_revenue_by_GDP, ...
                        tax_revenue_by_GDP(end)*ones(1, T_model-length(tax_revenue_by_GDP))];
                end
                
                if ~isbase
                    
                    % Calculate government expenditure adjustments
                    Dynamic_base = hardyload('dynamics.mat', base_generator, base_dir);
                    
                    Gtilde = Dynamic_base.Gtilde - gcut*GEXP_by_GDP(1:T_model).*Dynamic_base.outs;
                    Ttilde = tax_revenue_by_GDP.*Dynamic_base.outs - Static.pits - Static.ssts - Static.cits;
                    
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
        fprintf( 'Started at: %s \n', datetime );
        for label = { {'Beta'          , beta              } , ...
                      {'Gamma'         , gamma             } , ...
                      {'Sigma'         , sigma             } , ...
                      {'Model$'        , modelunit_dollars } }
            fprintf('\t%-25s= % 7.8f\n', label{1}{:})
        end
        
        
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
            [Dynamic, LABs, DIST, pgr, policy_fncs] = generate_aggregates(Market, DIST_steady, {}, {});
            
            
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
                        Gtilde = (tax_revenue_by_GDP - fedgovtnis).*Dynamic.outs - Dynamic.bens;
                        Ttilde = tax_revenue_by_GDP.*Dynamic.outs - Dynamic.pits - Dynamic.ssts - Dynamic.cits;
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
        
        Dynamic.is_converged = (eps <= tol);
        % Issue warning if did not converge
        if (~Dynamic.is_converged)
            warning('Model did not converge.')
        end
        
        
        % Save baseline optimal labor values and population distribution
        if isbase
            save(fullfile(save_dir, 'decisions.mat'   ), 'LABs')
            save(fullfile(save_dir, 'all_decisions.mat'   ), '-struct', 'policy_fncs')
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
                
                % Calculate $GDP/HH
                outperHH = (Dynamic.outs./Dynamic.pops)./modelunit_dollars;
                
                % Save and display elasticities
                save(fullfile(save_dir, 'elasticities.mat') ...
                    , 'captoout', 'labelas', 'savelas', 'outperHH' ...
                    , 'beta', 'gamma', 'sigma', 'modelunit_dollars' );
                
                fprintf( '\n' );
                fprintf( 'Finished at: %s\n', datetime );
                for label = { {'Beta'          , beta              } , ...
                              {'Gamma'         , gamma             } , ...
                              {'Sigma'         , sigma             } , ...
                              {'Model$'        , modelunit_dollars } }
                    fprintf('\t%-25s= % 7.8f\n', label{1}{:})
                end
                fprintf( '--------------\n' );
                for label = { {'Capital/Output'          , captoout } , ...
                              {'Labor elasticity'        , labelas  } , ...
                              {'Savings elasticity'      , savelas  } , ...
                              {'Output/HH'               , outperHH } }
                    fprintf('\t%-25s= % 7.4f\n', label{1}{:})
                end
                fprintf('\n');
                
        end
        
        
    end
    
end

end



%%
% Check if file exists and generate before loading if necessary, handling parallel write and read conflicts.
% 
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

%%
% Read a CSV file in format (Index), (Value)
%    For time series, (Index) is (Year), 
%       return a vector with first element being the value at
%       first_index (that is, first_year_transition - 1). The next value is the value at the start of
%       the transition path.
%   For other series (e.g. age_survival_probability, (Index) is (Age)
%       return the series from first_index. 
%   If first_index is empty, then return whole vector.
function [series] = read_series(filename, first_index, param_dir )

    warning( 'off', 'MATLAB:table:ModifiedVarnames' );

    % Check if file exists 
    filepath    = fullfile(param_dir, filename);
    if ~exist(filepath, 'file')
        err_msg = strcat('Cannot find file = ', strrep(filepath, '\', '\\'));
        throw(MException('read_series:FILENAME', err_msg ));
    end;
        
    T           = readtable(filepath, 'Format', '%u%f');
    indices     = table2array(T(:,1));
    vals        = table2array(T(:,2));
    
    if( isempty(first_index) )
        idx_start   = 1;
    else
        idx_start   = find( indices == first_index, 1);
    end;

    if( isempty(idx_start) )
        throw(MException('read_series:FIRSTINDEX','Cannot find first index in file.'));
    end;
    
    series      = vals(idx_start:end );

end