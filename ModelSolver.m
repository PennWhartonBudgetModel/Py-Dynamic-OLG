%%
% Dynamic model solver.
%
%%
classdef ModelSolver

methods (Static)

    % Solve dynamic model
    function [save_dir] = solve(scenario, callertag) %#ok<*FXUP>
        
        if ~exist('callertag' , 'var'), callertag  = ''; end
        economy = scenario.economy;
        
        %% Initialization
        
        % Start parallel pool if JVM enabled and pool does not already exist
        if usejava('jvm'), gcp; end
        
        % Unpack parameters from baseline definition
        beta                = scenario.beta ;
        gamma               = scenario.gamma;
        sigma               = scenario.sigma;
        modelunit_dollar    = scenario.modelunit_dollar;
        
        % Identify baseline run 
        isbase = scenario.isCurrentPolicy();
        
        % Unpack parameters from filled counterfactual definition
        expenditure_shift   = scenario.expenditure_shift;
        legal_scale         = scenario.legal_scale      ;
        prem_legal          = scenario.prem_legal       ;
        amnesty             = scenario.amnesty          ;
        deportation         = scenario.deportation      ;
        
        % Identify working directory
        save_dir = PathFinder.getWorkingDir(scenario);
        
        % Append caller tag to save directory name and generate calling tag
        %   Obviates conflicts between parallel solver calls
        save_dir   = [save_dir                                          , callertag];
        callingtag = [sprintf('^%s_%s', scenario.counterdeftag, economy), callertag];
        
        % Clear or create save directory
        if exist(save_dir, 'dir'), rmdir(save_dir, 's'), end, mkdir(save_dir)
        
        
        
        %% PARAMETERS
        
        % Define time constants
        s = ParamGenerator.timing(scenario);
        first_transition_year   = s.first_transition_year;  % map model inputs (or outputs) to actual years
        T_life                  = s.T_life;                 % Total life years
        T_work                  = s.T_work;                 % Total working years
        T_model                 = s.T_model;                % Transition path model years
        startyears              = s.startyears;             % Cohort start years as offsets to year 1
        nstartyears             = length(startyears);
        
        T_pasts   = max(-startyears, 0);                            % Life years before first model year
        T_shifts  = max(+startyears, 0);                            % Model years before first life year
        T_actives = min(startyears+T_life, T_model) - T_shifts;     % Life years within modeling period
        T_ends    = min(T_model-startyears, T_life);                % Maximum age within modeling period
        
        
        % Discretized grids, including shock process
        %   ndem is number of permanent types (hi, low)
        %   g    is population subgroup set: 'citizen',
        %           'legal', 'illegal'
        %   NOTE: DISTz and zs are indexed by g
        s = ParamGenerator.grids( scenario );
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
        s = ParamGenerator.production( scenario );
        A       = s.A;        % Total factor productivity
        alpha   = s.alpha;    % Capital share of output
        d       = s.d;        % Depreciation rate
        
        % Load population growth parameters
        % Load age-dependent parameters
        s = ParamGenerator.demographics();
        birth_rate      = s.birth_rate;                 % Annual birth rate
        legal_rate      = s.legal_rate * legal_scale;   % Annual legal immigration rate
        illegal_rate    = s.illegal_rate;               % Annual illegal immigration rate
        surv            = s.surv;                       % Survival probabilities by age
        imm_age         = s.imm_age;                    % Immigrants' age distribution
        
        % Load Social Security parameters
        s = ParamGenerator.social_security( scenario );
        sstax_brackets  = s.brackets    ;    % Payroll tax brackets (currentlaw is 0 to taxmax)
        sstax_rates     = s.rates       ;    % Payroll tax rates (currentlaw is 12.4%)
        sstax_burdens   = s.burdens     ;    % Cumulative tax liability at each bracket
        ssbenefits      = s.ssbenefits  ;    % Benefits
        sstaxcredit     = s.taxcredit   ;    % Benefit tax credit percentage
        ssincmaxs       = s.ssincmaxs   ;    % Maximum wage allowed for benefit calculation
        
        %%  Budget: CBO interest rates, expenditures, and debt
        s = ParamGenerator.budget( scenario );
        GEXP_by_GDP      = s.GEXP_by_GDP;           % Gvt expenditures as pct GDP
        debt             = s.debt;                  % Gvt debt as pct gdp
        debttoout        = s.debttoout;             % Initial debt/gdp (for steady state)
        debttoout_trans1 = s.debttoout_trans1;      % First transition debt/gdp
        fedgovtnis       = s.fedgovtnis;            % Gvt net interest surplus (deficit)
        cborates         = s.cborates;              % Interest rates on gvt debt (from CBO)
        cbomeanrate      = s.cbomeanrate;           % Avg of cborates (for steady state)
        % Tax revenue targets (for Ttilde), depend on tax plan
        tax_revenue_by_GDP = s.tax_revenue_by_GDP;
        
        %% Tax parameters
        s = ParamGenerator.tax( scenario );
        pittax_brackets     = s.tax_thresholds; % Tax func is linearized, these are income thresholds 
        pittax_burdens      = s.tax_burdens;    % Tax burden (cumulative tax) at thresholds
        pittax_rates        = s.tax_rates;      % Effective marginal tax rate between thresholds
        
        captaxshare         = s.captaxshare;            % Portion of capital income taxed at preferred rates
        expshare            = s.shareCapitalExpensing;  % Portion of investment which can be expensed
        taucap              = s.taucap;                 % Capital tax rate
        taucapgain          = s.taucapgain;             % Capital gains tax rate
        
        s_base = ParamGenerator.tax( scenario.currentPolicy );
        expshare_base = s_base.shareCapitalExpensing;   % expshare for baseline
        taucap_base   = s_base.taucap;                  % taucap for baseline
        
        qtobin0 = 1 - expshare_base*taucap_base;
        qtobin  = 1 - expshare     *taucap     ;

        % Define parameters on residual value of bequest function.
        s = ParamGenerator.bequest_motive( scenario );
        bequest_phi_1 = s.phi1;                 % phi1 reflects parent's concern about leaving bequests to her children (THIS IS THE ONE WE WANT TO CALIBRATE FOR LATER!)
        bequest_phi_2 = s.phi2;                 % phi2 measures the extent to which bequests are a luxury good
        bequest_phi_3 = s.phi3;                 % phi3 is the relative risk aversion coefficient

        
        %% Aggregate generation function
        
        function [Aggregate, LABs, DIST, OPTs, DIST_trans] = generate_aggregates(Market, DIST_steady, LABs_static, DIST_static)
            
            % Define dynamic aggregate generation flag
            isdynamic = isempty(LABs_static) || isempty(DIST_static);
            
            % Set static optimal decision values to empty values for dynamic aggregate generation
            if isdynamic, LABs_static = cell(nstartyears, ndem); end
            
            % Initialize optimal decision value arrays
            os = {'K', 'LAB', 'B', 'INC', 'PIT', 'SST', 'CIT', 'BEN', 'CON', 'V'};
            for o = os, OPTs.(o{1}) = zeros(nz,nk,nb,T_life,T_model,ndem); end
            
            % Initialize array of cohort optimal labor values
            LABs = cell(nstartyears, ndem);
            
            % Initialize population distribution array
            DIST = zeros(nz,nk,nb,T_life,ng,T_model,ndem);
            if (strcmp(scenario.economy, 'steady'))
                DIST_trans = zeros(nz,nk,nb,T_life,ng,T_model,ndem);
            else
                DIST_trans = {};
            end
            
            
            for idem = 1:ndem
                
                % Package fixed dynamic optimization arguments into anonymous function
                solve_cohort_ = @(V0, LAB_static, T_past, T_shift, T_active, T_work) solve_cohort(V0, LAB_static, isdynamic, ...
                    nz, nk, nb, T_past, T_shift, T_active, T_work, T_model, zs(:,:,idem), transz, Market.kpricescale*kv, bv, beta, gamma, sigma, surv, ...
                    bequest_phi_1, bequest_phi_2, bequest_phi_3, ...
                    modelunit_dollar, ...
                    sstaxcredit, ssbenefits, ssincmaxs, ...
                    sstax_brackets, sstax_burdens, sstax_rates, ...
                    pittax_brackets, pittax_burdens, pittax_rates, ... 
                    captaxshare, taucap, taucapgain, qtobin, qtobin0, ...
                    Market.beqs, Market.wages, Market.capshares, Market.caprates, Market.govrates, Market.totrates, Market.expsubs);
                
                
                % Initialize series of terminal utility values
                V0s = zeros(nz,nk,nb,T_life);
                
                % Solve steady state / post-transition path cohort
                if isdynamic
                    
                    % Solve dynamic optimization
                    % (Note that active time is set to full lifetime)
                    OPT = solve_cohort_(V0s(:,:,:,T_life), [], T_pasts(end), T_shifts(end), T_life, T_work(end));
                    
                    % Define series of terminal utility values
                    V0s(:,:,:,1:T_life-1) = OPT.V(:,:,:,2:T_life);
                    
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
                            OPTs_cohort{i} = solve_cohort_(V0, LABs_static{i,idem}, T_pasts(i), T_shifts(i), T_actives(i), T_work(i));
                            
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
                        DIST_next = generate_distribution(DIST_year, DIST_grow, K, B, nz, nk, nb, T_life, ng, transz, Market.kpricescale*kv, bv, surv);
                        assert(all(DIST_next(:)>=0), 'Negative mass of people at DIST_next.')
                        
                        % Increase legal immigrant population for amnesty, maintaining distributions over productivity
                        DISTz_legal = DIST_next(:,:,:,:,g.legal) ./ repmat(sum(DIST_next(:,:,:,:,g.legal), 1), [nz,1,1,1,1]);
                        DISTz_legal(isnan(DISTz_legal)) = 1/nz;
                        
                        DIST_next(:,:,:,:,g.legal) = DIST_next(:,:,:,:,g.legal) + repmat(sum(amnesty*DIST_next(:,:,:,:,g.illegal), 1), [nz,1,1,1,1]).*DISTz_legal;
                        
                        % Reduce illegal immigrant population for amnesty and deportation
                        DIST_next(:,:,:,:,g.illegal) = (1-amnesty-deportation)*DIST_next(:,:,:,:,g.illegal);
                        
                        % Calculate age distribution convergence error
                        f = @(D) sum(sum(reshape(D, [], T_life, ng), 1), 3) / sum(D(:));
                        disteps = max(abs(f(DIST_next) - f(DIST_year)));
                        
                        year = year + 1;
                        
                        switch economy, case 'steady'
                            DIST_trans(:,:,:,:,:,1,idem) = DIST_next;
                        end
                        
                    end
                    
                else
                    DIST = DIST_static;
                end
                
            end
            
            
            % Normalize steady state population distribution
            switch economy, case 'steady'
                DIST_trans = DIST_trans / sum(DIST(:));
                DIST = DIST / sum(DIST(:));
            end
            
            % Generate aggregates
            assert(all(DIST(:)>=0),'WARNING! Negative mass of people at DIST.')
            DIST_gs = reshape(sum(DIST, 5), [nz,nk,nb,T_life,T_model,ndem]);
            f = @(F) sum(sum(reshape(DIST_gs .* F, [], T_model, ndem), 1), 3);
            
            Aggregate.pops     = f(1);                                                                                   % Population
            Aggregate.assets   = f(repmat(reshape(Market.kpricescale*kv, [1,nk,1,1,1,1]), [nz,1,nb,T_life,T_model,ndem])); % Assets
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
            
        end
        
        
        
        %% Static aggregate generation
        
        if ~isbase
            
            baselineScenario = scenario.currentPolicy();
            base_generator = @() ModelSolver.solve(baselineScenario, callingtag);
            base_dir = PathFinder.getWorkingDir(baselineScenario);
            
            % Load baseline market conditions, optimal labor values, and population distribution
            Market = hardyload('market.mat'      , base_generator, base_dir);
            
            s      = hardyload('decisions.mat'   , base_generator, base_dir);
            LABs_static = s.LABs;
            
            s      = hardyload('distribution.mat', base_generator, base_dir);
            DIST_static = s.DIST;
            
            
            % Generate static aggregates
            % (Intermediary structure used to filter out extraneous fields)
            [Static, ~, Static_DIST, Static_OPTs, ~] = ...
                generate_aggregates(Market, {}, LABs_static, DIST_static);
            
            % Copy additional static aggregates from baseline aggregates
            Dynamic_base = hardyload('dynamics.mat', base_generator, base_dir);
            
            for series = {'caps', 'caps_domestic', 'caps_foreign', 'capincs', 'labincs', 'outs', 'debts_domestic', 'debts_foreign', 'Gtilde', 'Ttilde'}
                Static.(series{1}) = Dynamic_base.(series{1});
            end

            % Calculate static budgetary aggregate variables
            Static.cits_domestic = Static.cits;
            Static.cits_foreign  = taucap * Market.caprates * captaxshare .* (qtobin*Static.caps_foreign);
            Static.cits          = Static.cits_domestic + Static.cits_foreign;
            Static.revs          = Static.pits + Static.ssts + Static.cits - Static.bens;            
            Static.labpits       = Static.pits .* Static.labincs ./ Static.incs;
            Static.caprevs       = Static.cits + Static.pits - Static.labpits;

            Static.debts = [Dynamic_base.debts(1), zeros(1,T_model-1)];
            for year = 1:T_model-1
                Static.debts(year+1) = Static.Gtilde(year) - Static.Ttilde(year) - Static.revs(year) + Static.debts(year)*(1 + cborates(year));
            end
                        
            % Save static aggregates
            save(fullfile(save_dir, 'statics.mat'), '-struct', 'Static')
            save(fullfile(save_dir, 'Static_all_decisions.mat'   ), '-struct', 'Static_OPTs')
            save(fullfile(save_dir, 'Static_distribution.mat'), 'Static_DIST')        
            
        end
        
        
        
        %% Dynamic aggregate generation
        
        switch economy
            
            case 'steady'
                
                % Load initial conditions
                Market0 = struct('beqs',0.0927,'capshares',3/(3+debttoout),'rhos',6.2885);
                    % Initial guesses set as follows:
                    % capshare = (K/Y / (K/Y + D/Y)), where K/Y = captoout = 3 and D/Y = debttoout.
                    % beqs and rhos are guesses from previous code.
                
                DIST_steady = {};
                
            case {'open', 'closed'}
                
                % Make Scenario for current policy, steady state. 
                steadyBaseScenario = scenario.currentPolicy().steady();
                steady_generator = @() ModelSolver.solve(steadyBaseScenario, callingtag);
                steady_dir = PathFinder.getWorkingDir(steadyBaseScenario);
                
                % Load steady state market conditions and dynamic aggregates
                Market0  = hardyload('market.mat'      , steady_generator, steady_dir);
                Dynamic0 = hardyload('dynamics.mat'    , steady_generator, steady_dir);
                
                % Load steady state population distribution
                s        = hardyload('distribution.mat', steady_generator, steady_dir);
                
                DIST_steady = s.DIST_trans;
                
        end
        
        
        switch economy
            
            case 'open'
                
                if ~isbase
                    
                    % Calculate government expenditure adjustments
                    Dynamic_base = hardyload('dynamics.mat', base_generator, base_dir);
                    
                    Gtilde = Dynamic_base.Gtilde + expenditure_shift*GEXP_by_GDP(1:T_model).*Dynamic_base.outs;
                    Ttilde = tax_revenue_by_GDP.*Dynamic_base.outs - Static.pits - Static.ssts - Static.cits;
                    
                end
                
            case 'closed'
                
                % Make Scenario for the open economy. 
                openScenario = scenario.open();
                open_generator = @() ModelSolver.solve(openScenario, callingtag);
                open_dir = PathFinder.getWorkingDir(openScenario);
                
                % Load government expenditure adjustments
                Dynamic_open = hardyload('dynamics.mat', open_generator, open_dir);
                
                Gtilde = Dynamic_open.Gtilde;
                Ttilde = Dynamic_open.Ttilde;
                
        end
        
        
        % Define marketing clearing tolerance and initialize error term
        tol = 1e-5;
        eps = Inf;
        
        % Initialize iteration count and set maximum number of iterations
        iter    =  0;
        itermax = 30;
        
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
                      {'phi_1'         , bequest_phi_1     } , ...
                      {'depreciation'  , d                 } , ...
                      {'Model$'        , modelunit_dollar  }   ...
                    };
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
                    
                    Market.caprates = max((A*alpha*(Market.rhos.^(alpha-1)) - d), 0);
                    Market.totrates = Market.capshares.*Market.caprates + (1-Market.capshares).*Market.govrates;
                    
                case 'open'
                    
                    if isinitial
                        Market.caprates  = Market0.caprates*ones(1,T_model)*(1-taucap_base)/(1-taucap);
                        Market.govrates  = Market0.govrates*ones(1,T_model);
                        Market.totrates  = Market0.totrates*ones(1,T_model);
                    end
                    
            end
            
            Market.rhos  = ((Market.caprates + d)/(A*alpha)).^(1/(alpha-1));
            Market.wages = A*(1-alpha)*(Market.rhos.^alpha);
            
            Market.kpricescale = 1 + Market.capshares(1)*(qtobin - qtobin0)/qtobin;
            Market.qtobin0     = qtobin0;
            Market.qtobin      = qtobin;
            
            % Generate dynamic aggregates
            [Dynamic, LABs, DIST, OPTs, DIST_trans] = generate_aggregates(Market, DIST_steady, {}, {});
            
            
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
                    beqs = Dynamic.bequests / (sum(DIST_trans(:))/sum(DIST(:)));
                    clearing = Market.rhos - rhos;
                    
                    % Calculate income - THIS SHOULD BE OUTSIDE THE CASE
                    % ECONOMY LOOP!
                    Dynamic.labincs = Dynamic.labeffs .* Market.wages;
                    Dynamic.capincs = qtobin * Market.caprates .* Dynamic.caps;
                    
                    Dynamic.labpits = Dynamic.pits .* Dynamic.labincs ./ Dynamic.incs;
                    Dynamic.caprevs = Dynamic.cits + Dynamic.pits - Dynamic.labpits;

                case 'open'
                    
                    % Calculate capital and output
                    Dynamic.caps = Market.rhos .* Dynamic.labeffs;
                    Dynamic.outs = A*(max(Dynamic.caps, 0).^alpha).*(Dynamic.labeffs.^(1-alpha));
                    
                    Dynamic.caps_domestic = (Market.capshares .* [Dynamic0.assets, Dynamic.assets(1:T_model-1)])/qtobin;
                    Dynamic.caps_foreign  = Dynamic.caps - Dynamic.caps_domestic;
                    
                    % Calculate debt
                    Dynamic.cits_domestic = Dynamic.cits;
                    Dynamic.cits_foreign  = taucap * Market.caprates * captaxshare .* (qtobin*Dynamic.caps_foreign);
                    Dynamic.cits          = Dynamic.cits_domestic + Dynamic.cits_foreign;
                    
                    if isbase
                        Gtilde = (tax_revenue_by_GDP - fedgovtnis).*Dynamic.outs - Dynamic.bens;
                        Ttilde = tax_revenue_by_GDP.*Dynamic.outs - Dynamic.pits - Dynamic.ssts - Dynamic.cits;
                    end
                    
                    Dynamic.revs  = Dynamic.pits + Dynamic.ssts + Dynamic.cits - Dynamic.bens;
                    Dynamic.debts = [debttoout_trans1*Dynamic0.outs, zeros(1,T_model-1)];
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
                    Dynamic.debts = [debttoout_trans1*Dynamic0.outs, zeros(1,T_model-1)];
                    for year = 1:T_model-1
                        Dynamic.debts(year+1) = Gtilde(year) - Ttilde(year) - Dynamic.revs(year) + Dynamic.debts(year)*(1 + cborates(year));
                    end
                    Dynamic.Gtilde = Gtilde;
                    Dynamic.Ttilde = Ttilde;
                    
                    Dynamic.debts_domestic = Dynamic.debts;
                    Dynamic.debts_foreign  = zeros(1,T_model);
                    
                    % Calculate capital and output
                    Dynamic.caps = ([Dynamic0.assets, Dynamic.assets(1:end-1)] - Dynamic.debts)/qtobin;
                    Dynamic.outs = A*(max(Dynamic.caps, 0).^alpha).*(Dynamic.labeffs.^(1-alpha));
                    
                    Dynamic.caps_domestic = Dynamic.caps;
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
            
            Dynamic.revs  = Dynamic.pits + Dynamic.ssts + Dynamic.cits - Dynamic.bens;
            
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
        end
        
        % Save market conditions and dynamic aggregates
        save(fullfile(save_dir, 'market.mat'  ), '-struct', 'Market' )
        save(fullfile(save_dir, 'dynamics.mat'), '-struct', 'Dynamic')
        save(fullfile(save_dir, 'all_decisions.mat'   ), '-struct', 'OPTs')
        switch economy
            case 'steady'
                DIST = struct('DIST', DIST, 'DIST_trans', DIST_trans);
                save(fullfile(save_dir, 'distribution.mat'), '-struct', 'DIST')
            case {'open', 'closed'}
                save(fullfile(save_dir, 'distribution.mat'), 'DIST')
        end                
        
        
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
                    DIST_idem = sum(DIST.DIST(:,:,:,:,:,1,idem), 5);
                    
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
                outperHH = (Dynamic.outs./Dynamic.pops)./modelunit_dollar;
                
                % Calculate gini
                GiniTable = MomentsGenerator(scenario,DIST.DIST,Market,OPTs).giniTable;
                gini      = GiniTable.model(GiniTable.Gini=='wealth');

                % Save and display elasticities
                save(fullfile(save_dir, 'paramsTargets.mat') ...
                    , 'captoout', 'labelas', 'savelas', 'outperHH', 'gini' ...
                    , 'beta', 'gamma', 'sigma', 'modelunit_dollar', 'bequest_phi_1' );
                
                fprintf( '\n' );
                fprintf( 'Finished at: %s\n', datetime );
                for label = { {'Beta'          , beta              } , ...
                              {'Gamma'         , gamma             } , ...
                              {'Sigma'         , sigma             } , ...
                              {'Model$'        , modelunit_dollar  } , ...
                              {'phi_1'         , bequest_phi_1     } }
                    fprintf('\t%-25s= % 7.8f\n', label{1}{:})
                end
                fprintf( '--------------\n' );
                for label = { {'Capital/Output'        , captoout   } , ...
                              {'Labor elasticity'      , labelas    } , ...
                              {'Savings elasticity'    , savelas    } , ...
                              {'Output/HH'             , outperHH   } , ...
                              {'Wealth Gini'           , gini       } }
                    fprintf('\t%-25s= % 7.4f\n', label{1}{:})
                end
                fprintf('\n');
                
        end
        
        % Release MEX file to avoid locks.
        clear mex;
    end % solve
    
end % methods

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

