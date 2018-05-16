%%
% Dynamic model solver.
%
%%
classdef ModelSolver

methods (Static)

    % Solve dynamic model
    function [save_dir] = solve(scenario, callertag, resolve) %#ok<*FXUP>
        
        % Identify working directory
        save_dir = PathFinder.getWorkingDir(scenario);
        
        % Return if scenario already solved and resolve flag not set
        if scenario.isSolved() && (~exist('resolve', 'var') || isempty(resolve) || ~resolve)
            fprintf('Scenario already solved.\n');
            return
        end
        
        if ~exist('callertag' , 'var'), callertag  = ''; end
        
        % Append caller tag to save directory name and generate calling tag
        %   Obviates conflicts between parallel solver calls
        save_dir   = [save_dir                                                      , callertag];
        callingtag = [sprintf('^%s_%s', scenario.counterdeftag, scenario.economytag), callertag];
        
        % Clear or create save directory
        if exist(save_dir, 'dir'), rmdir(save_dir, 's'), end, mkdir(save_dir)
        
        
        
        %% PARAMETERS
        
        economy             = scenario.economy;
                
        % Unpack parameters from baseline definition
        beta                = scenario.beta ;
        gamma               = scenario.gamma;
        sigma               = scenario.sigma;
        modelunit_dollar    = scenario.modelunit_dollar;
        
        closure_year        = scenario.ClosureYear - scenario.TransitionFirstYear;
        
        % Identify baseline run 
        isbase = scenario.isCurrentPolicy();
        
        % Immigration policies
        legal_scale         = scenario.legal_scale      ;
        prem_legal          = scenario.prem_legal       ;
        amnesty             = scenario.amnesty          ;
        deportation         = scenario.deportation      ;

        % Define time constants
        s = ParamGenerator.timing(scenario);
        first_transition_year   = s.TransitionFirstYear;    % map model inputs (or outputs) to actual years
        T_life                  = s.T_life;                 % Total life years
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
        A               = s.A;              % Total factor productivity
        alpha           = s.alpha;          % Capital share of output
        depreciation    = s.depreciation;   % Depreciation rate
        
        % Load population growth parameters
        % Load age-dependent parameters
        s = ParamGenerator.demographics( scenario );
        birth_rate      = s.birth_rate;                 % Annual birth rate
        legal_rate      = s.legal_rate * legal_scale;   % Annual legal immigration rate
        illegal_rate    = s.illegal_rate;               % Annual illegal immigration rate
        surv            = s.surv;                       % Survival probabilities by age
        imm_age         = s.imm_age;                    % Immigrants' age distribution
        
        % Load Social Security parameters, 
        %   rem: all USdollars have been converted to modelunit dollars
        socialsecurity  = ParamGenerator.social_security( scenario );
        T_works         = socialsecurity.T_works     ;    % Total working years
        sstax_brackets  = socialsecurity.taxbrackets ;    % Payroll tax brackets (currentlaw is 0 to taxmax)
        sstax_rates     = socialsecurity.taxrates    ;    % Payroll tax rates (currentlaw is 12.4%)
        sstax_indices   = socialsecurity.taxindices  ;    % Type of index to use on tax bracket
        sstaxcredit     = socialsecurity.taxcredit   ;    % Benefit tax credit percentage
        ssincmaxs       = socialsecurity.ssincmaxs   ;    % Maximum wage allowed for benefit calculation
        ssincmins       = socialsecurity.ssincmins   ;    % Minimum wage allowed for benefit calculation
        
        % NOTE: Social Security benefits are calculated per cohort
        
        %%  Budget: CBO interest rates, expenditures, and debt
        budget = ParamGenerator.budget( scenario );
        
        %% Tax parameters
        %    rem: all US dollars have been converted to modelunit_dollars
        tax = ParamGenerator.tax( scenario );
        
        %  TBD: Revise this to clean up. Have separate structs for ordinary
        %  and preferred HH rates.
        
        %  The following brackets, burdens, rates vary by year
        pit.brackets      = tax.brackets;       % Tax func is linearized, these are income thresholds 
        pit.burdens       = tax.burdens;        % Tax burden (cumulative tax) at thresholds
        pit.rates         = tax.rates;          % Effective marginal tax rate between thresholds
        pit.prefbrackets  = tax.prefbrackets;
        pit.prefburdens   = tax.prefburdens;
        pit.prefrates     = tax.prefrates;
        pit.rateCapGain   = tax.rateCapGain;    % Capital gains tax
        
        captaxshares      = tax.captaxshare;    % Portion of capital income taxed at preferred rates

        % Tax parameters for current policy, steady-state
        sstax = ParamGenerator.tax( scenario.steady().currentPolicy() );  
        
        % Define parameters on residual value of bequest function.
        s = ParamGenerator.bequest_motive( scenario );
        bequest_phi_1 = s.phi1;                 % phi1 reflects parent's concern about leaving bequests to her children (THIS IS THE ONE WE WANT TO CALIBRATE FOR LATER!)
        bequest_phi_2 = s.phi2;                 % phi2 measures the extent to which bequests are a luxury good
        bequest_phi_3 = s.phi3;                 % phi3 is the relative risk aversion coefficient

        % Instantiate a Firm (SingleFirm)
        theFirm       = Firm( scenario, Firm.SINGLEFIRM );
        
        % Instantiate the Pass-Through Firm
        thePassThrough  = Firm( scenario, Firm.PASSTHROUGH );
        
        
        %% Aggregate generation function
        
        function [Aggregate, LABs, savings, DIST, OPTs, DIST_trans] = generate_aggregates(Market, DIST_steady, LABs_static, savings_static, DIST_static)
            
            % Define dynamic aggregate generation flag
            isdynamic = isempty(DIST_static) || isempty(LABs_static) || isempty(savings_static);
            
            % Set static optimal decision values to empty values for dynamic aggregate generation
            if isdynamic, LABs_static = cell(nstartyears); savings_static = cell(nstartyears); end
            
            % Initialize optimal decision value arrays
            os = {'V', 'LABOR', 'SAVINGS', 'CONSUMPTION', 'AVG_EARNINGS', 'TAXABLE_INC', 'OASI_BENEFITS', 'ORD_LIABILITY', 'PAYROLL_LIABILITY', 'PREF_LIABILITY'};
            for o = os, OPTs.(o{1}) = zeros(nz,nk,nb,T_life,T_model); end
            
            % Initialize array of cohort optimal labor values
            LABs    = cell(nstartyears);
            savings = cell(nstartyears);
            
            % Initialize population distribution array
            DIST = zeros(nz,nk,nb,T_life,ng,T_model);
            if (strcmp(scenario.economy, 'steady'))
                DIST_trans = zeros(nz,nk,nb,T_life,ng,T_model);
            else
                DIST_trans = {};
            end
            
            % Calculate indexing vectors
            Market.priceindices = ModelSolver.generate_index(scenario, Market.wages, nstartyears, T_model, T_life);
            
            % Calculate indexed policy variables
            ssincmins_indexed   = (ssincmins .* Market.priceindices.wage_inflations)';
            ssincmaxs_indexed   = (ssincmaxs .* Market.priceindices.wage_inflations)';
                
            [sstax_brackets_indexed, sstax_rates_indexed, sstax_burdens_indexed] ...
                    = ModelSolver.indexSSTax( sstax_brackets         ...
                                            , sstax_indices          ...
                                            , sstax_rates            ...
                                            , Market.priceindices );            

            % Package fixed dynamic optimization arguments into anonymous function
            solve_cohort_ = @(V0, LAB_static, saving_static, T_past, T_shift, T_active, T_works, ssbenefits, cohort_wageindexes) ...
                solve_cohort( ...
                    V0, LAB_static, saving_static, isdynamic, ...
                    nz, nk, nb, T_past, T_shift, T_active, T_works, T_model, ...
                    zs, transz, kv, bv, scenario.beta, scenario.gamma, scenario.sigma, surv, ...
                    bequest_phi_1, bequest_phi_2, bequest_phi_3, ...
                    ssbenefits, ssincmins_indexed, ssincmaxs_indexed, cohort_wageindexes, ...
                    sstax_brackets_indexed, sstax_burdens_indexed, sstax_rates_indexed, ...
                    sstaxcredit, pit.brackets, pit.burdens, pit.rates, ... 
                    captaxshares, ...
                    pit.prefbrackets, pit.prefburdens, pit.prefrates, ... 
                    pit.rateCapGain, ...
                    Market.beqs, ...
                    Market.wages, ...
                    Market.capshares_0, zeros(1,T_model), ... % portfolio allocations
                    zeros(1,T_model), zeros(1,T_model), 0,  ... % pass-through prices
                    Market.equityFundDividends, Market.equityFundPrices, Market.equityFundPrice0, ... % Equity returns and prices
                    Market.bondFundDividends, Market.bondFundPrices, Market.bondFundPrice0  ... % Bond returns and prices
                    ); 

            % Initialize series of terminal utility values
            V0s = zeros(nz,nk,nb,T_life);

            % Solve steady state / post-transition path cohort
            if isdynamic

                % Note: retire_year = 1 so that ssbenefits is
                % calculated for T_model = 1
                ssbenefits = ModelSolver.calculateSSBenefitForCohort(   ...
                                            socialsecurity.startyear_benefitbrackets(end, :)   ...
                                        ,   socialsecurity.startyear_benefitrates   (end, :)   ...
                                        ,   Market.priceindices                                ...
                                        ,   1                                                  ...
                                        ,   bv, T_model );

                % Solve dynamic optimization
                % (Note that active time is set to full lifetime)
                OPT = solve_cohort_(V0s(:,:,:,T_life), [], [], T_pasts(end), T_shifts(end), T_life, T_works(end), ssbenefits, Market.priceindices.cohort_wages(:,end));

                % Define series of terminal utility values
                V0s(:,:,:,1:T_life-1) = OPT.V(:,:,:,2:T_life);

            end


            switch economy

                case 'steady'

                    % Store optimal decision values
                    for o = os, OPTs.(o{1})(:,:,:,:,1) = OPT.(o{1}); end
                    LABs{1}    = OPT.LABOR;
                    savings{1} = OPT.SAVINGS;

                case {'open', 'closed'}

                    % Solve transition path cohorts
                    OPTs_cohort = cell(1,nstartyears);
                    parfor i = 1:nstartyears

                        % Extract terminal utility values
                        V0 = V0s(:,:,:,T_ends(i)); %#ok<PFBNS>

                        % Calculate cohort-based year-varying benefits policy
                        ssbenefits = ModelSolver.calculateSSBenefitForCohort(   ...
                                            socialsecurity.startyear_benefitbrackets(i, :)  ...
                                        ,   socialsecurity.startyear_benefitrates   (i, :)  ...
                                        ,   Market.priceindices                             ...  
                                        ,   socialsecurity.retire_years(i)                  ...
                                        ,   bv, T_model );

                        % Solve dynamic optimization
                        OPTs_cohort{i} = solve_cohort_(V0, LABs_static{i}, savings_static{i}, T_pasts(i), T_shifts(i), T_actives(i), T_works(i), ssbenefits, Market.priceindices.cohort_wages(:,i));

                        LABs{i}    = OPTs_cohort{i}.LABOR;
                        savings{i} = OPTs_cohort{i}.SAVINGS;

                    end

                    % Construct optimal decision value arrays
                    for i = 1:nstartyears
                        for t = 1:T_actives(i)
                            age  = t + T_pasts (i);
                            year = t + T_shifts(i);
                            for o = os, OPTs.(o{1})(:,:,:,age,year) = OPTs_cohort{i}.(o{1})(:,:,:,t); end
                        end
                    end

            end


            if isdynamic

                switch economy

                    case 'steady'
                        
                        % Determine steady state age distribution without immigration
                        [v, ~] = eigs([birth_rate*ones(1,T_life); [diag(surv(1:T_life-1)), zeros(T_life-1,1)]], 1, 'largestreal');
                        DISTage0 = v / sum(v);
                        
                        % Define initial population distribution using steady state distributions over age and productivity without immigration
                        DIST_next = zeros(nz,nk,nb,T_life,ng);
                        DIST_next(:,:,:,:,g.citizen) = repmat(reshape(repmat(reshape(DISTage0, [1,T_life]), [nz,1]) .* DISTz(:,:,g.citizen), [nz,1,1,T_life]), [1,nk,nb,1]) / (nk*nb);
                        
                        % Specify number of distribution generation years as number of years required to achieve steady state
                        nyears = 153;

                    case {'open', 'closed'}
                        
                        % Define initial population distribution as steady state distribution
                        DIST_next = DIST_steady(:,:,:,:,:,1);
                        
                        % Specify number of distribution generation years as number of model years
                        nyears = T_model;

                end

                for year = 1:nyears

                    % Store population distribution for current year
                    DIST_year = DIST_next;
                    DIST(:,:,:,:,:,min(year, T_model)) = DIST_year;

                    % Extract optimal decision values for current year
                    K = OPTs.SAVINGS     (:,:,:,:,min(year, T_model));
                    B = OPTs.AVG_EARNINGS(:,:,:,:,min(year, T_model));

                    % Define population growth distribution
                    DIST_grow = zeros(nz,nk,nb,T_life,ng);
                    P = sum(DIST_year(:));

                    DIST_grow(:,1,1,1,g.citizen) = reshape(DISTz(:,1,g.citizen), [nz,1,1,1     ,1]) * P * birth_rate   ;
                    DIST_grow(:,1,1,:,g.legal  ) = reshape(DISTz(:,:,g.legal  ), [nz,1,1,T_life,1]) * P * legal_rate   .* repmat(reshape(imm_age, [1,1,1,T_life,1]), [nz,1,1,1,1]);
                    DIST_grow(:,1,1,:,g.illegal) = reshape(DISTz(:,:,g.illegal), [nz,1,1,T_life,1]) * P * illegal_rate .* repmat(reshape(imm_age, [1,1,1,T_life,1]), [nz,1,1,1,1]);

                    % Generate population distribution for next year
                    DIST_next = generate_distribution(DIST_year, DIST_grow, K, B, nz, nk, nb, T_life, ng, transz, kv, bv, surv);
                    assert(all(DIST_next(:)>=0), 'Negative mass of people at DIST_next.')

                    % Increase legal immigrant population for amnesty, maintaining distributions over productivity
                    DISTz_legal = DIST_next(:,:,:,:,g.legal) ./ repmat(sum(DIST_next(:,:,:,:,g.legal), 1), [nz,1,1,1,1]);
                    DISTz_legal(isnan(DISTz_legal)) = 1/nz;

                    DIST_next(:,:,:,:,g.legal) = DIST_next(:,:,:,:,g.legal) + repmat(sum(amnesty*DIST_next(:,:,:,:,g.illegal), 1), [nz,1,1,1,1]).*DISTz_legal;

                    % Reduce illegal immigrant population for amnesty and deportation
                    DIST_next(:,:,:,:,g.illegal) = (1-amnesty-deportation)*DIST_next(:,:,:,:,g.illegal);

                    switch economy, case 'steady'
                        DIST_trans(:,:,:,:,:,1) = DIST_next;
                    end

                end

            else
                DIST = DIST_static;
            end
            
            
            % Normalize steady state population distribution
            switch economy, case 'steady'
                DIST_trans = DIST_trans / sum(DIST(:));
                DIST = DIST / sum(DIST(:));
            end
            
            % Generate aggregates
            assert(all(DIST(:)>=0),'WARNING! Negative mass of people at DIST.')
            DIST_gs = reshape(sum(DIST, 5), [nz,nk,nb,T_life,T_model]);
            f = @(F) sum(sum(reshape(DIST_gs .* F, [], T_model), 1), 3);
            
            Aggregate.pops     = f(1);                                                                                   % Population
            Aggregate.bequests = f(OPTs.SAVINGS .* repmat(reshape(1-surv, [1,1,1,T_life,1]), [nz,nk,nb,1,T_model]));     % Bequests
            Aggregate.labs     = f(OPTs.LABOR);                                                                          % Labor
            Aggregate.labeffs  = f(OPTs.LABOR .* repmat(reshape(zs, [nz,1,1,T_life,1]), [1,nk,nb,1,T_model]));           % Effective labor
            Aggregate.lfprs    = f(OPTs.LABOR > 0.01) ./ f(1);                                                           % Labor force participation rate
            Aggregate.incs     = f(OPTs.TAXABLE_INC);                                                                    % Income
            Aggregate.pits     = f(OPTs.ORD_LIABILITY);                                                                  % Personal income tax
            Aggregate.ssts     = f(OPTs.PAYROLL_LIABILITY);                                                              % Social Security tax
            Aggregate.cits     = f(OPTs.PREF_LIABILITY);                                                                 % Capital income tax
            Aggregate.bens     = f(OPTs.OASI_BENEFITS);                                                                  % Social Security benefits
            Aggregate.cons     = f(OPTs.CONSUMPTION);                                                                    % Consumption
            Aggregate.assets_0 = f(repmat(reshape(kv, [1,nk,1,1,1]), [nz, 1,nb,T_life,T_model]));                        % Assets before re-pricing
            Aggregate.assets_1 = Aggregate.assets_0 .* (ones(1,T_model) + Market.capgains') ...                          % Assets after re-pricing            
                                    .* (Market.capshares_0./Market.capshares_1);                                         % Note: The definition of assets_1 corresponds to beginning of period assets at new policy prices, that is, accounting for eventual capital gains.
            
        end
        
        
        
        
        %% Helper function to calculate debt
        %     Inputs: (1) Aggregate (Static or Dynamic), 
        %             (2) Current iteration residuals: Gtilde, Ctilde, Ttilde
        %             (3) Initial debt -- i.e. t=1 debt
        %             (4) Rates on debt
        %             (5) T_model
        %     Outputs: None, but side-effect --> revised Aggregate which includes debts series
        function [debts] = calculate_debts(Aggregate, Gtilde, Ctilde, Ttilde, debt_1, debtrates, T_model) 
            debts = [debt_1 zeros(1,T_model-1)];
            for t = 1:T_model-1
                debts(t+1) = Gtilde(t) - Ctilde(t) + Aggregate.bens(t)  ...
                             - Ttilde(t) - Aggregate.revs(t)               ...
                             + debts(t)*(1 + debtrates(t))                 ...
                             ;
            end
        end

        
        %% Static aggregate generation
        
        if ~isbase && ~strcmp(economy, 'steady')
            
            baselineScenario = scenario.currentPolicy();
            base_generator = @() ModelSolver.solve(baselineScenario, callingtag);
            base_dir = PathFinder.getWorkingDir(baselineScenario);
            
            % Load baseline market conditions, optimal labor values, and population distribution
            Market = hardyload('market.mat'      , base_generator, base_dir);
            
            s      = hardyload('decisions.mat'   , base_generator, base_dir);
            LABs_static    = s.LABs;
            savings_static = s.savings;
            
            s      = hardyload('distribution.mat', base_generator, base_dir);
            DIST_static = s.DIST;
            
            
            % Generate static aggregates
            % (Intermediary structure used to filter out extraneous fields)
            [Static, ~, ~, Static_DIST, Static_OPTs, ~] = ...
                generate_aggregates(Market, {}, LABs_static, savings_static, DIST_static);
            
            % Copy additional static aggregates from baseline aggregates
            Dynamic_base = hardyload('dynamics.mat', base_generator, base_dir);
            
            for series = {'caps', 'caps_domestic', 'caps_foreign', 'capincs', 'labincs', 'outs', 'investment', 'debts_domestic', 'debts_foreign', 'Gtilde', 'Ttilde', 'Ctilde' }
                Static.(series{1}) = Dynamic_base.(series{1});
            end
            
            % Calculate static budgetary aggregate variables
            [corpDividends, corpTaxs]  = theFirm.dividends( Static.caps', Static.investment', ...
                                                            Market.rhos', Market.wages'       ...
                                                           );
            Static.corpTaxs  = corpTaxs';
            Static.dividends = corpDividends';
            Static.revs = Static.pits + Static.ssts + Static.cits + Static.corpTaxs;            
            
            % In general, we have:
            % Static.debts = Static.debts_domestic + Static.debts_foreign;
            % But this is not true in the static economies.
            % Since static decisions are fixed at the baseline, domestic and
            % foreign debts do not change. But total debt changes due to new
            % tax policy, which implies the equality above no longer holds.
            % Notice that debts is a combination of the actual static debts and
            % the residual mismatch from markets not clearing
            %    Rem: Dynamic_base(1) is supposed to be D' debt carried
            %    from steady state (before policy change)
            Static.debts = calculate_debts( Static, Static.Gtilde, Static.Ctilde, Static.Ttilde, Dynamic_base.debts(1), budget.debtrates, T_model); 

            % Total assets
            % Note: tot_assets is a sum of choice variables, those are constant at baseline values
            Static.tot_assets_1 = theFirm.priceCapital' .* Static.caps + ...
                                  Static.debts_domestic + Static.debts_foreign;
                        
            % Save static aggregates
            save(fullfile(save_dir, 'statics.mat')              , '-struct', 'Static')
            save(fullfile(save_dir, 'Static_decisions.mat')     , '-struct', 'Static_OPTs')
            save(fullfile(save_dir, 'Static_distribution.mat')  , 'Static_DIST')        
            
        end
        
        
        
        %% Dynamic aggregate generation
        
        switch economy
            
            case 'steady'
                
                % Load initial guesses
                Dynamic0.outs       = 3.35;
                Dynamic0.caps       = 12.20; 
                captoout            = Dynamic0.caps / Dynamic0.outs;
                
                Market0.beqs        = 0.1662;                                   % beqs are results from previous runs.
                Market0.capshares_0 = captoout / (captoout + budget.debttoout); % capshare = (K/Y / (K/Y + D/Y)), where K/Y = captoout = 3 and D/Y = debttoout.
                Market0.capshares_1 = Market0.capshares_0;                      % capshare = (K/Y / (K/Y + D/Y)), where K/Y = captoout = 3 and D/Y = debttoout.
                Market0.rhos        = 7.0652;                                   % rhos are results from previous runs.
                Market0.invtocaps   = 0.0078 + depreciation;                    % I/K = pop growth rate 0.0078 + depreciation

                Dynamic0.debts      = Dynamic0.outs * budget.debttoout;
                Dynamic0.assets_0   = theFirm.priceCapital*Dynamic0.caps + Dynamic0.debts;
                Dynamic0.assets_1   = Dynamic0.assets_0;
                Dynamic0.labeffs    = Dynamic0.caps / Market0.rhos; 
                Dynamic0.investment = Dynamic0.caps * Market0.invtocaps;
                
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
                s = hardyload('distribution.mat', steady_generator, steady_dir);
                DIST_steady = s.DIST_trans;
                
        end
        
        
        switch economy
            
            case 'open'
                
                if ~isbase
                    
                    % Calculate government expenditure adjustments
                    Dynamic_base = hardyload('dynamics.mat', base_generator, base_dir);
                    
                    Gtilde = budget.outlays_by_GDP     .* Dynamic_base.outs            ...
                             - Static.bens;
                    Ttilde = budget.tax_revenue_by_GDP .* Dynamic_base.outs            ...
                             - Static.revs; 
                    Ctilde = zeros(1,T_model);
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
                Ctilde = Dynamic_open.Ctilde;
        end
        
        
        
        
        %%
        % Define marketing clearing tolerance and initialize error term
        tolerance.rhos      = 1e-5;
        tolerance.beqs      = 1e-4;
        tolerance.invtocaps = 5e-4;
        isConverged         = false;
        
        % Set damper to update guesses
        %    0 = not dampened, i.e., completely updated to new value
        % 	 1 = fully dampened, i.e., stays the same
        switch economy
            case 'steady'
                damper.rhos      = 0.5;
                damper.beqs      = 0.5;
                damper.capshares = 0.5;
            case 'open'
                damper.rhos      = 1.0;      % In open economy, it's set by theFirm.calculateKLRatio function.
                damper.beqs      = 0.0;
                damper.capshares = 0.0;
            case 'closed'
                damper.rhos      = 0.0;
                damper.beqs      = 0.0;
                damper.capshares = 0.0; 
        end
        
        % Initialize iteration count and set maximum number of iterations
        iter    =   0;
        itermax = 100;
        
        % Create file for logging iterations
        % (Note that this should be done after any parallel pool is started to avoid file access issues)
        iterlog = fopen(fullfile(save_dir, 'iterations.csv'), 'w');
        
        % Display header
        fprintf( 'Started at: %s \n', datetime );
        fprintf( '%s\n', scenario.shortDescription );
        
        while (~isConverged && iter < itermax)
            
            % Increment iteration count
            iter = iter + 1;
            fprintf(' Iteration %2d  ...  ', iter)
            isinitial = iter == 1;
            
            % Capital prices
            Market.capgains(1,1) = (theFirm.priceCapital(1) - theFirm.priceCapital0)/theFirm.priceCapital0;
            for t = 2:T_model
               Market.capgains(t,1) = (theFirm.priceCapital(t) - theFirm.priceCapital(t-1))/theFirm.priceCapital(t-1);
            end

            % Define market conditions in the first iteration
            if isinitial
                Market.beqs        = Market0.beqs       *ones(1,T_model);
                Market.capshares_0 = Market0.capshares_0*ones(1,T_model);
                Market.capshares_1 = Market0.capshares_1*ones(1,T_model);
                Market.invtocaps   = Market0.invtocaps  *ones(1,T_model);
                
                Dynamic.assets_0    = Dynamic0.assets_0     *ones(1,T_model);
                Dynamic.assets_1    = Dynamic0.assets_1     *ones(1,T_model);
                Dynamic.debts       = Dynamic0.debts        *ones(1,T_model);
                Dynamic.caps        = Dynamic0.caps         *ones(1,T_model); 
                Dynamic.labeffs     = Dynamic0.labeffs      *ones(1,T_model);
                Dynamic.investment  = Dynamic0.investment   *ones(1,T_model);
                
                switch economy

                    case {'steady', 'closed'}

                        Market.rhos      = Market0.rhos*ones(1,T_model);
                        outs             = Dynamic0.outs;
                        
                    case 'open'

                        % Rem: Returns are fixed to match steady-state in
                        % open economy. That is, after-tax returns for
                        % capital are fixed, including cap gains.
                        % 
                        % First period capital (inherited from steady state)
                        Dynamic.caps(1) = Market0.capshares_0 * Dynamic0.assets_0;
                        % Define the pre-tax returns necessary to return
                        % the world rate from steady-state.
                        effectiveDividendRate = ( Market0.equityFundDividends*(1 - sstax.rateForeignCorpIncome) ...
                                                  - Market.capgains ) ...
                                                ./ (1 - tax.rateForeignCorpIncome);
                        klRatio = theFirm.calculateKLRatio( effectiveDividendRate   , ... 
                                                            Dynamic.caps'           , ...
                                                            Dynamic.labeffs'        , ...
                                                            Market0.invtocaps );
                        
                        rhos        = Market0.rhos*ones(1,T_model);
                        Market.rhos = klRatio';
                        %% TBD: Do we need to set caps to have foreign investment here?
                        %       Seems like yes.
                end
                
                Market.debtrates = budget.debtrates;
                Market.MPKs      = A*alpha*( Market.rhos.^(alpha-1) ); % This is just for reporting
            
            % end initial loop iteration
            else  
                
                Market.beqs        = damper.beqs*Market.beqs + (1 - damper.beqs)*beqs;
                Market.rhos        = damper.rhos*Market.rhos + (1-damper.rhos)*rhos;
                Market.capshares_0 = damper.capshares*Market.capshares_0 + (1-damper.capshares)*capshares_0;
                Market.capshares_1 = damper.capshares*Market.capshares_1 + (1-damper.capshares)*capshares_1;
                
                Market.MPKs      = A*alpha*( Market.rhos.^(alpha-1) );
                Market.invtocaps = invtocaps;
                
                switch economy

                    case 'open'
                        % NOTE: For open economy, capshares will NOT converge
                        % because the portfolio allocation is fixed by steady-state
                        % and we do not allow it to change even as the economy's 
                        % mix of capital vs. debt changes.
                        % Overwrite the first period capital
                        Dynamic.caps(1) = Market.capshares_0(1) * Dynamic.assets_0(1);
                        klRatio     = theFirm.calculateKLRatio( effectiveDividendRate   , ...
                                                                Dynamic.caps'           , ...
                                                                Dynamic.labeffs'        , ...
                                                                Market0.invtocaps );

                        rhos        = Market.rhos;
                        Market.rhos = klRatio';
                        
                end
            end
            
            
            % Compute prices
            Market.wages               = A*(1-alpha)*(Market.rhos.^alpha);
            [corpDividends, corpTaxs]  = theFirm.dividends( Dynamic.caps'           ...
                                                        ,   Dynamic.investment'     ...
                                                        ,   Market.rhos'            ...
                                                        ,   Market.wages'           ...
                                                        );
            
            % 'Price' of assets -- HH own equal shares of both bond & equity funds
            % (equityFund/bondFund)Dividends are actually dividend rates
            Market.equityFundPrice0     = theFirm.priceCapital0;
            Market.equityFundPrices     = theFirm.priceCapital';  
            Market.equityFundDividends  = (corpDividends ./ (Dynamic.caps' .* theFirm.priceCapital))';
            
            Market.bondFundPrice0       = 1;
            Market.bondFundPrices       = ones(1,T_model);
            Market.bondFundDividends    = budget.debtrates; %rem: dividendrate is per $ of assets
            
            
            % Generate dynamic aggregates
            [Dynamic, LABs, savings, DIST, OPTs, DIST_trans] = generate_aggregates(Market, DIST_steady, {}, {}, {});
            

            % Calculate additional dynamic aggregates
            % (Note that open economy requires capital calculation before debt calculation 
            % while closed economy requires the reverse)
            
            % Taxes
            Dynamic.corpDividends = corpDividends';
            Dynamic.corpTaxs      = corpTaxs';
            Dynamic.revs          = Dynamic.pits + Dynamic.ssts      ...
                                  + Dynamic.cits + Dynamic.corpTaxs  ...
                                   ;
            
            switch economy
                
                case 'steady'
                    
                    % Calculate debt, capital, and output
                    % (Numerical solver used due to absence of closed form solution)
                    f_debts = @(outs ) budget.debttoout*outs;
                    f_caps  = @(debts) (Dynamic.assets_1 - debts) ./ theFirm.priceCapital;
                    f_outs  = @(caps ) A*(max(caps, 0).^alpha).*(Dynamic.labeffs.^(1-alpha));
                    x_ = fsolve(@(x) x - [f_debts(x(3)); f_caps(x(1)); f_outs(x(2))], zeros(3,1), optimoptions('fsolve', 'Display', 'none'));
                    Dynamic.debts = x_(1);
                    Dynamic.caps  = x_(2);
                    Dynamic.outs  = x_(3);

                    Dynamic.debts_domestic = Dynamic.debts;
                    Dynamic.debts_foreign  = zeros(1,T_model);
                    Dynamic.caps_domestic  = Dynamic.caps;
                    Dynamic.caps_foreign   = zeros(1,T_model);
                    Dynamic.invest_foreign = zeros(1,T_model);
                    Dynamic.tot_assets_0   = Dynamic.assets_0;
                    Dynamic.tot_assets_1   = Dynamic.assets_1;
                    
                    % Calculate income
                    Dynamic.labincs = Dynamic.labeffs .* Market.wages;
                    Dynamic.capincs = Market.MPKs .* Dynamic.caps;
                    capincs_foreign = Market.equityFundDividends .* Dynamic.caps_foreign;
                    Dynamic.GNP     = Dynamic.outs - capincs_foreign;
                    
                    % Proxy for gross investment in physical capital
                    DIST_gs            = reshape(sum(DIST, 5), [nz,nk,nb,T_life,T_model]);
                    assets_tomorrow    = sum(sum(reshape(DIST_gs .* OPTs.SAVINGS, [], T_model), 1), 3);
                    Dynamic.investment = (Market.capshares_1 * (assets_tomorrow - Dynamic.bequests))./ theFirm.priceCapital' ...
                                         - (1 - depreciation) * Dynamic.caps;
                                     
                    % Update guesses
                    rhos      = Dynamic.caps / Dynamic.labeffs;
                    beqs      = Dynamic.bequests / sum(DIST_trans(:));          % Note: capgains is zero in steady state, so bequests don't need to be changed
                    invtocaps = Dynamic.investment ./ Dynamic.caps;
                    capshares_0 = (Dynamic.assets_0 - Dynamic.debts) ./ Dynamic.assets_0;
                    capshares_1 = (Dynamic.assets_1 - Dynamic.debts) ./ Dynamic.assets_1;

                case 'open'
                    
                    % Calculate capital and output
                    Dynamic.caps = Market.rhos .* Dynamic.labeffs;
                    Dynamic.outs = A*(max(Dynamic.caps, 0).^alpha).*(Dynamic.labeffs.^(1-alpha));
                    
                    % Note: Dynamic.assets_0 represents current assets at old prices.
                    Dynamic.caps_domestic  =  (Market.capshares_1 .* Dynamic.assets_1) ./ theFirm.priceCapital';
                    Dynamic.caps_foreign   = Dynamic.caps - Dynamic.caps_domestic;
                    Dynamic.invest_foreign = [Dynamic.caps_foreign(2:T_model) Dynamic.caps_foreign(T_model)] ...
                                              - (1 - depreciation) * [Dynamic.caps_foreign(1:T_model-1) Dynamic.caps_foreign(T_model-1)];
                    
                    if isbase
                        Gtilde = budget.outlays_by_GDP     .*Dynamic.outs              ...
                                 - Dynamic.bens;
                        Ttilde = budget.tax_revenue_by_GDP .*Dynamic.outs              ...
                                 - Dynamic.revs;
                        Ctilde = zeros(1,T_model);
                    end
                    
                    % Rem: debts(1) is an exogenous calculation: (D/Y) * Y
                    % for t=1 where D/Y is from data.
                    % Debt should be D' from steady state, but that is not right to
                    % calibrate to real world
                    Dynamic.debts = calculate_debts( Dynamic, Gtilde, Ctilde, Ttilde, budget.debttoout * Dynamic.outs(1), budget.debtrates, T_model);

                    Dynamic.Gtilde = Gtilde;
                    Dynamic.Ttilde = Ttilde;
                    Dynamic.Ctilde = Ctilde;
                    
                    Dynamic.debts_domestic = (1 - Market.capshares_1) .* Dynamic.assets_1;
                    Dynamic.debts_foreign  = Dynamic.debts - Dynamic.debts_domestic;
                    
                    Dynamic.tot_assets_0   = [theFirm.priceCapital0 theFirm.priceCapital(1:T_model-1)'] ...
                                               .* Dynamic.caps + Dynamic.debts;
                    Dynamic.tot_assets_1   = theFirm.priceCapital' .* Dynamic.caps + Dynamic.debts;
                    
                    % Calculate income
                    Dynamic.labincs = Dynamic.labeffs .* Market.wages;
                    Dynamic.capincs = Market.MPKs .* Dynamic.caps;
                    capincs_foreign = Market.equityFundDividends .* Dynamic.caps_foreign;
                    Dynamic.GNP     = Dynamic.outs - capincs_foreign;
                    
                    % Gross investment in physical capital
                    %   T_model investment converges to final steady
                    %   state investment, K_ss*(pop_growth + depreciation).
                    %   In an open economy, capital gets to its optimal level
                    %   instantaneously, so it's ok to use the steady state
                    %   rate of capital replacement here.
                    Dynamic.investment = [Dynamic.caps(2:T_model)   0] - (1 - depreciation) * ...
                                         [Dynamic.caps(1:T_model-1) 0];
                    Dynamic.investment(T_model) = Market0.invtocaps * Dynamic.caps(T_model);
                    
                    % Update guesses
                    % Note: Bequests should be priced according to the new policy because it
                    %       corresponds to yesterday's assets that were collected and sold by the
                    %       government yesterday after some people died, but redistributed today
                    %       after the new policy took place.
                    %       So we apply today's prices to yesterday's bequests and capshares.
                    beqs      = [Dynamic0.bequests * (1 + Market0.capshares_1 * Market.capgains(1)), ...
                                 Dynamic.bequests(1:T_model-1) .* (1 + Market.capshares_1(1:T_model-1) .* Market.capgains(2:T_model)') ...
                                ] ./ Dynamic.pops;
                    invtocaps = Dynamic.investment ./ Dynamic.caps;
                    capshares_0 = (Dynamic.assets_0 - Dynamic.debts) ./ Dynamic.assets_0;
                    capshares_1 = (Dynamic.assets_1 - Dynamic.debts) ./ Dynamic.assets_1;

                case 'closed'
                                    
                    % Calculate Ctilde to close debt growth
                    %closure_debttoout = Dynamic.debts(closure_year)/Dynamic.outs(closure_year);
                    %cont_debttoout    = Dynamic.debts(closure_year:T_model) ./ Dynamic.outs(closure_year:T_model);
                    %cont_Ctilde       = (cont_debttoout - closure_debttoout) .* Dynamic.outs(closure_year:T_model);
                    %tmpCtilde         = [zeros(1:closure_year) cont_Ctilde]

                    % Rem: debts(1) is an exogenous calculation: (D/Y) * Y
                    % for t=1 where D/Y is from data.
                    % Debt should be D' from steady state, but that is not right to
                    % calibrate to real world
                    Dynamic.debts = calculate_debts( Dynamic, Gtilde, Ctilde, Ttilde, budget.debttoout * outs(1), budget.debtrates, T_model);

                    Dynamic.Gtilde = Gtilde;
                    Dynamic.Ttilde = Ttilde;
                    Dynamic.Ctilde = Ctilde;
                    
                    Dynamic.debts_domestic = Dynamic.debts;
                    Dynamic.debts_foreign  = zeros(1,T_model);
                    
                    % Calculate capital and output
                    % Note: Dynamic.assets represents current assets at new prices.
                    Dynamic.caps = (Dynamic.assets_1 - Dynamic.debts) ./ theFirm.priceCapital';
                    outs         = A*(max(Dynamic.caps, 0).^alpha).*(Dynamic.labeffs.^(1-alpha));
                    Dynamic.outs = outs;  % outs var is used to keep last iteration values
                    
                    Dynamic.caps_domestic  = Dynamic.caps;
                    Dynamic.caps_foreign   = zeros(1,T_model);
                    Dynamic.invest_foreign = zeros(1,T_model);
                    Dynamic.tot_assets_0   = Dynamic.assets_0;
                    Dynamic.tot_assets_1   = Dynamic.assets_1;
                    
                    % Calculate income
                    Dynamic.labincs = Dynamic.labeffs .* Market.wages;
                    Dynamic.capincs = Market.MPKs .* Dynamic.caps;
                    capincs_foreign = Market.equityFundDividends .* Dynamic.caps_foreign;
                    Dynamic.GNP     = Dynamic.outs - capincs_foreign;
                    
                    % Gross investment in physical capital
                    %   T_model investment eventually converges to final steady
                    %   state investment, K_ss*(pop_growth + depreciation),
                    %   but closed economy capital has not yet converged so
                    %   the last period is a more smooth guess
                    Dynamic.investment = [Dynamic.caps(2:T_model)   0] - (1 - depreciation) * ...
                                         [Dynamic.caps(1:T_model-1) 0];
                    Dynamic.investment(T_model) = Market.invtocaps(T_model-1) * Dynamic.caps(T_model);

                    % Update guesses
                    % Note: Dynamic.assets represents current assets at new prices.
                    %       Bequests should also be priced according to the new policy.
                    %       So we apply today's prices to yesterday's bequests and capshares.
                    rhos      = Dynamic.caps ./ Dynamic.labeffs;
                    beqs      = [Dynamic0.bequests * (1 + Market0.capshares_1 * Market.capgains(1)), ...
                                 Dynamic.bequests(1:T_model-1) .* (1 + Market.capshares_1(1:T_model-1) .* Market.capgains(2:T_model)') ...
                                ] ./ Dynamic.pops;
                    invtocaps = Dynamic.investment ./ Dynamic.caps;
                    capshares_0 = (Dynamic.assets_0 - Dynamic.debts) ./ Dynamic.assets_0;
                    capshares_1 = (Dynamic.assets_1 - Dynamic.debts) ./ Dynamic.assets_1;

            end
            
            % Calculate market clearing series
            clearing.rhos      = max(abs((Market.rhos      - rhos)      ./ rhos     ));
            clearing.beqs      = max(abs((Market.beqs      - beqs)      ./ beqs     ));
            clearing.invtocaps = max(abs((Market.invtocaps - invtocaps) ./ invtocaps));
            clearing.capshares = max(abs((Market.capshares_1 - capshares_1) ./ capshares_1));
                    
            % Check convergence
            isConverged = (clearing.rhos      < tolerance.rhos     ) && ...
                          (clearing.beqs      < tolerance.beqs     ) && ...
                          (clearing.invtocaps < tolerance.invtocaps);
            
            fprintf('Errors: K/L = %7.6f beqs = %7.6f I/K = %7.6f capshares = %7.6f\n', ...
                    clearing.rhos, clearing.beqs, clearing.invtocaps, clearing.capshares)
            fprintf(iterlog, '%u,%0.6f,%0.6f,%0.6f,%0.6f\n', iter, ...
                    clearing.rhos, clearing.beqs, clearing.invtocaps, clearing.capshares);
            
        end % GE loop
        
        fprintf('\n')
        fclose(iterlog);
        
        Dynamic.is_converged = isConverged;
        % Issue warning if did not converge
        if (~Dynamic.is_converged)
            warning('Model did not converge.')
        end
        
        % Save market conditions, HH policies, and dynamic aggregates
        save(fullfile(save_dir, 'market.mat'   )    , '-struct', 'Market' )
        save(fullfile(save_dir, 'dynamics.mat' )    , '-struct', 'Dynamic')
        save(fullfile(save_dir, 'decisions.mat')    , 'OPTs', 'LABs', 'savings')
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
                captoout = (Dynamic.assets_1 - Dynamic.debts) / Dynamic.outs;
                
                
                % Calculate labor elasticity
                LAB_  = LABs{1};
                DIST_ = sum(DIST.DIST(:,:,:,:,:,1), 5);

                workind = (LAB_ > 0.01);

                workmass = sum(DIST_(workind));
                frisch   = sum(DIST_(workind) .* (1 - LAB_(workind)) ./ LAB_(workind)) * (1 - gamma*(1-sigma))/sigma;
                
                labelas = frisch / workmass;
                
                
                % Calculate savings elasticity
                ratedev = 0.01;
                Market_dev = Market;
                
                % Increase rates of return to HH by ratedev
                %   Note: Cap gains is zero in steady state, so 
                %         return to HH is only on equity + debt.
                Market_dev.equityFundDividends = Market.equityFundDividends * (1 + ratedev);
                Market_dev.bondFundDividends   = Market.bondFundDividends   * (1 + ratedev);
                
                [Dynamic_dev] = generate_aggregates(Market_dev, {}, {}, {}, {});
                
                savelas = (Dynamic_dev.assets_1 - Dynamic.assets_1) / (Dynamic.assets_1 * ratedev);
                
                % Calculate $GDP/HH
                outperHH = (Dynamic.outs./Dynamic.pops)./scenario.modelunit_dollar;
                
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
        
        
        %%
        
        % Create completion indicator file
        fclose(fopen(fullfile(save_dir, 'solved'), 'w'));
        
        % Release MEX file to avoid locks
        clear mex; %#ok<CLMEX>
        
    end % solve
    
end % methods



methods (Static, Access = private ) 
    
    % Calculate SS benefits for a cohort
    % Inputs:   
    %       brackets: Initial year bracket thresholds. These 
    %                 grow by wage_inflations index until cohort retires
    %                 and then freeze.
    %       rates   : Initial year replacement rates for each cohort. These
    %                 grow by given policy adjustments (which can include
    %                 any net adjustments due to COLA (REM: No adjustment
    %                 implies a COLA = CPI growth.)
    %       adjustments     : Adjustments (as ratio) to the rates.
    %       priceindices    : Price index structure with <wage_inflations> index
    %                       used to grow the brackets.
    %       retire_year     : Model year of first year of retirement for this
    %                       cohort.
    %       bv, T_model     : Benefits discretized grid, and time to end of
    %                       model run.
    function ssbenefits = calculateSSBenefitForCohort(  brackets                    ...
                                                     ,  rates                       ...
                                                     ,  priceindices                ...
                                                     ,  retire_year                 ...
                                                     ,  bv, T_model )
        
        % Fetch index used to grow brackets
        wage_index = priceindices.wage_inflations;
                                                 
        % Cohort based benefits by year
        %   REM: Every household is at a grid point, so can 
        %        calculate benefits directly here (outside of
        %        solve_cohort).
        %    NOTE: We set ssbenefits to -Inf otherwise -- it should not be
        %    used in solve_cohort (this will blow it up just in case)
        ssbenefits  = ones(T_model, size(bv,1)) .* -Inf;
        
        % Build indexed bracket cutoffs for each year (for the cohort)
        bracket_idx = [     wage_index( 1:min(retire_year, T_model) )                           ...
                        ,   ones(1, T_model - retire_year) ];
        bracket_idx = bracket_idx(1:min(end,T_model))';
        adjbrackets = repmat(brackets, [T_model, 1])                            ... 
                        .* repmat( bracket_idx, [1, size(brackets,1)] );
        
        % Build benefit replacement rates adjusted by policy
        adjrates    = repmat(rates, [T_model, 1]);
       
        % Build cumulative benefits matrix to aid benefit calculation
        adjtotben   = cumsum(diff(adjbrackets, 1, 2).*adjrates(:, 1:end-1), 2); 
        adjtotben   = [zeros(size(adjbrackets, 1), 1), adjtotben];  % rem: first column is zero
        
        % Calculate benefits for each year of retirement until end of time.
        for t = max(retire_year,1):T_model
            for ib = 1:size(bv,1)
                thebracket       = find(adjbrackets(t,:) <= bv(ib), 1, 'last');
                ssbenefits(t,ib) = adjtotben(t,thebracket)            ...
                                    + adjrates(t,thebracket)*(bv(ib) - adjbrackets(t,thebracket));
            end
        end
        
    end % calculateSSBenefitForCohort
    
    
    %% Calculate SS tax brackets using required indexing
    %    If brackets overlap because of indexing, reorder the brackets 
    %   Inputs:
    %       brackets        = bracket cutoffs in nominal $, 
    %                       (rem: first bracket cutoff is zero)
    %       rates           = rates for brackets
    %       bracketindices  = index to use for cutoff: (e.g. 'reals',
    %                       'nominals', 'wage_inflations', but not 'cohort_wages' ) 
    %       priceindices    = struct with indices
    function [ssbrackets, ssrates, ssburdens] = indexSSTax(     brackets        ...
                                                            ,   bracketindices  ...
                                                            ,   rates           ...
                                                            ,   priceindices ) 
        if( any(strcmp(bracketindices, 'cohort_wages')) )
            throw MException('calculateSSTaxBrackets:INDEX', 'Cannot use index type <cohort_wages>.' );
        end
        
        indices    = zeros(size(brackets));
        for i = 1:size(indices, 2)
            indices(:,i) = priceindices.(bracketindices{i});
        end
        ssbrackets = brackets .* indices;

        % Sort brackets and rates just in case of overlap
        [ssbrackets, sortindex] = sort(ssbrackets, 2);
        
        ssrates = zeros(size(rates));
        for r = 1:size(sortindex, 1)
            ssrates(r, :) = rates(r, sortindex(r, :) );
        end

        % Calculate tax burdens
        ssburdens = cumsum(diff(ssbrackets, 1, 2).*ssrates(:, 1:end-1), 2); 
        ssburdens = [zeros(size(ssbrackets, 1), 1), ssburdens];  % rem: first burden is zero

    end % indexSSTax
    

    %%
    % Create indexes for benefits and taxmax calculations and import CPI index
    %
    %   Inputs:
    %       scenario     = used to import a couple of variables, 
    %       Market_wages = T_model-dimension vector, 
    %       nstartyears  = number of cohorts, 
    %       startyears   = period of birth of a cohort (first period of transition is t = 0), 
    %       T_model      = number of periods, 
    function index = generate_index(scenario, Market_wages, nstartyears, T_model, T_life)

        index.wage_inflations = Market_wages./Market_wages(1);            % Time-varying indexes
        index.cohort_wages    = ones(T_model, nstartyears);               % Time- and cohort-varying indexes
        index.nominals        = 1./ParamGenerator.budget( scenario ).CPI; % Time-varying reciprocal CPI indexes from CBO
        index.reals           = ones(size(Market_wages));                 % Vector to 'index' real variables
        
        realage_entry = ParamGenerator.timing(scenario).realage_entry;
        % Indexes for the boundary cohorts
        cohortage60_at_1      = T_life + 1 - (60 - realage_entry);
        cohortage60_at_Tmodel = cohortage60_at_1 + T_model - 1;
        
        for i = cohortage60_at_1:cohortage60_at_Tmodel
            year_turn60 = i - (60 - realage_entry);
            index.cohort_wages(:,i) = Market_wages(year_turn60)./Market_wages(:);
        end

        % TBD: Erase these couple lines when the time comes to implement the policies above
        index.wage_inflations = ones(size(Market_wages));              % Time-varying indexes
        index.cohort_wages    = ones(T_model, nstartyears);            % Time- and cohort-varying indexes
    
    end % generate_index

    
end % private methods

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

