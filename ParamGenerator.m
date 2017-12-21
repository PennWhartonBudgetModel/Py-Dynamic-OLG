%%
% Reader and generator of parameters for dynamic model.
%
%%
classdef ParamGenerator

methods (Static)
    
    
    %% TIMING
    %      Includes policy-specific SS.NRA
    function s = timing(scenario)
        
        s.realage_entry         = 20    ;   % Real age of model age=0   
        s.first_transition_year = 2018  ;   % TBD: This will come from Scenario
        s.T_life                = 80    ;   % Death year
        s.Tmax_work             = 52    ;   % Last possible year of working age     
        
        switch scenario.economy
            case 'steady'
                s.T_model    = 1;                           % Steady state total modeling years
                s.startyears = 0;                           % Steady state cohort start year
            case {'open', 'closed'}
                s.T_model    = 25;                          % Transition path total modeling years
                s.startyears = (-s.T_life+1):(s.T_model-1); % Transition path cohort start years
        end
        
    end % timing
       
    
    %% DICRETIZATION GRIDS
    %   prem_legal : productivity premium of avg. (by year) legal
    %          immigrant. Note: This is a ratio.
    function s = grids(scenario)
        
        % (Shock process definitions from Storesletten, Telmer, and Yaron, 2004)
        % (Method of discretizing shock process derived from Adda and Cooper, 2003)
        
        % Define function to compute cutoff points for equally weighted slices of a normal distribution
        f = @(n, sig) n*sig*-diff(normpdf(norminv(linspace(0, 1, n+1))));
        
        % Determine permanent and transitory shocks
        nperm  = 2; zperm  = f(nperm , sqrt(0.2105));
        ntrans = 2; ztrans = f(ntrans, sqrt(0.0630));
        
        % Determine persistent shocks
        npers = 2;
        
        pep = 0.990;                        % Lagged productivity coefficient
        sigep = sqrt(0.018);                % Conditional standard deviation of productivity
        sigpers = sigep/(sqrt(1-pep^2));
        
        zpers = f(npers, sigpers);
        
        % Construct Markov transition matrix for persistent shocks by approximating an AR(1) process
        persv = sigpers*norminv(linspace(0, 1, npers+1));
        transpers = zeros(npers,npers);
        for ipers = 1:npers
            integrand = @(x) exp(-x^2/(2*sigpers^2)) * diff(normcdf((persv(ipers:ipers+1) - pep*x)/sigep));
            for jpers = 1:npers
                transpers(ipers,jpers) = (npers/(sqrt(2*pi*(sigpers^2)))) * quadgk(@(v) arrayfun(integrand, v), persv(jpers), persv(jpers+1));
            end
        end
        
        % Determine initial distribution over persistent shocks
        DISTpers = diff(normcdf(persv/sqrt(0.124)));
        
        % Define deterministic lifecycle productivities
        timing    = ParamGenerator.timing(scenario);
        T_life    = timing.T_life;
        T_workMax = timing.Tmax_work;
        % Life-cycle productivity from Conesa et al. 2017 - average for healthy workers
        zage      = read_series('ConesaEtAl_WageAgeProfile.csv', 1 + timing.realage_entry, PathFinder.getMicrosimInputDir());
        
        % Calculate total productivity
        ndem = nperm; nz = ntrans*npers;
        
        zs = max(0, repmat(reshape(zage                       , [1 ,T_life,1   ]), [nz,1     ,ndem]) ...
                  + repmat(reshape(zperm                      , [1 ,1     ,ndem]), [nz,T_life,1   ]) ...
                  + repmat(reshape(kron(ones(1,npers), ztrans) ...
                                 + kron(zpers, ones(1,ntrans)), [nz,1     ,1   ]), [1 ,T_life,ndem]));

        zs_check = zs(:,1:T_workMax,:);
        assert(all(zs_check(:) > 0), 'WARNING! Productivity shock grid contains zero.')
        
        % Construct Markov transition matrix for all shocks / productivities
        transz = kron(transpers, (1/ntrans)*ones(ntrans,ntrans));
        
        % Determine initial distribution over all shocks / productivities
        DISTz0 = kron(DISTpers , (1/ntrans)*ones(1     ,ntrans));

        % Include a fifth super large and rare shock
        nz = 5;
        zs(5,:,:) = zs(4,:,:) * 15;
        transz = [transz(1,1) transz(1,2) transz(1,3) transz(1,4) 0.00;
                  transz(2,1) transz(2,2) 0.02        0.02        (1-(transz(2,1)+transz(2,2))-2*0.02);
                  transz(3,1) transz(3,2) 0.47        0.47        (1-(transz(3,2)+transz(3,2))-2*0.47) ;
                  0.03        0.03        0.47        0.46         0.01                    ;
                  0.15        0.05        0.05        0.25         0.50];
        DISTz0 = [DISTz0 0];
        
        % Determine productivity distributions for ages
        DISTz_g      = zeros(nz,T_life);
        DISTz_g(:,1) = reshape(DISTz0, [nz,1]);
        
        for age = 2:T_life
            DISTz_g(:,age) = transz' * DISTz_g(:,age-1);
        end
        
        % Define population group index mapping
        groups = {'citizen', 'legal', 'illegal'};
        ng = length(groups);
        for ig = 1:ng, g.(groups{ig}) = ig; end
        
        % Replicate age productivity distributions for population groups
        DISTz = repmat(DISTz_g, [1,1,ng]);
        
        % Shift productivity distributions for legal and illegal immigrants towards highest and lowest productivity levels respectively
        %   prem_legal   : productivity premium of avg. (by year) legal
        %          immigrant. Note: This is a ratio.
        %          This value is a policy param on the Scenario.
        %   prem_illegal : productivity premium of avg. (by year) illegal
        %          immigrant. Note: This is a ratio.
        %          This value is calculated by PWBM
        %          at \\SHARED_DRIVE\PWBM_MicroSIM\Outputs\StructuralModelInputs
        prem_legal   = scenario.prem_legal;
        prem_illegal = 0.9; %  TODO: The real value is 0.62043 -- lowest shock is above this;
        for age = 1:T_life
            
            zmean = mean(zs(:,age,:), 3);
            
            if ( zmean .* DISTz(:,age,g.citizen) ~= zmean .* DISTz(:,age,g.legal) * prem_legal )
                zlegal   = sum(zmean .* DISTz(:,age,g.legal  )) * prem_legal  ; 
                plegal   = (zmean(nz) - zlegal  ) / (zmean(nz)*(nz-1) - sum(zmean(1:nz-1)));
                DISTz(:,age,g.legal  ) = [plegal*ones(nz-1,1); 1 - plegal*(nz-1)    ];
            end

            zillegal = sum(zmean .* DISTz(:,age,g.illegal)) * prem_illegal;
            pillegal = (zmean(1)  - zillegal) / (zmean(1) *(nz-1) - sum(zmean(2:nz  )));
            DISTz(:,age,g.illegal) = [1 - pillegal*(nz-1); pillegal*ones(nz-1,1)];
            
        end
        
        % Checks
        assert(all(transz(:) >= 0), 'WARNING! Negative transition probabilities.')
        assert(all(DISTz (:) >= 0), 'WARNING! Negative initial distribution of people DISTz.')
        if scenario.prem_legal==1
            citizen_legal = abs(DISTz(:,:,g.citizen)-DISTz(:,:,g.legal));
            assert(all(citizen_legal(:) < 1e-14), ...
                   'WARNING! Legal immigrants distribution does not match natives distribution although prem_legal = %f.\n', scenario.prem_legal )
        end
        
        % Define savings and average earnings discretization vectors
        f  = @(lb, ub, n, curv) lb + (ub-lb)*((0:n-1)/(n-1))'.^curv;
        nk = 12; kv = f(1e-3, 1/(275*scenario.modelunit_dollar), nk-4, 4);   % savings vector --- the grid is built to range from approx 10 to 1.5 million dollars --- 8.7230e-05 corresponds to the last modelunit_dollar value in steady state
        scale = 1;                                                           % scaling parameter to continue building the capital grid 
        for ik = nk-3:nk                                                     % this loop builds the top of capital grid such that:
            scale = 3.5*scale;                                               % i.  it re-starts at around 3.5 million dollars
            kv(ik) = scale*1e+6*scenario.modelunit_dollar;                   % ii. and its last point is around 150 million dollars 
        end
        ssincmax = 1.185e5;                                                  % 118,500 is the maximum annual labor income for Social Security tax purposes
        nb =  5; bv = f(   0, ssincmax*scenario.modelunit_dollar, nb  , 2);  % average earnings vector --- Upper bound of average earnings defined as maximum possible Social Security benefit
        bv = [bv; 15*bv(end)];                                               % max point arbitrarily set to 15x the second largest one
        nb = nb + 1;

        s.ndem     = ndem;
        s.g        = g;
        s.ng       = ng;
        s.nz       = nz;
        s.transz   = transz;
        s.DISTz    = DISTz;
        s.zs       = zs;
        s.nk       = nk;
        s.kv       = kv;
        s.nb       = nb;
        s.bv       = bv;
        
    end % grids
    
    
    %% TAXES
    % 
    % Generate tax policy parameters according to predefined plans.
    % 
    function s = tax( scenario )
        
        timing                  = ParamGenerator.timing(scenario);
        % Find tax plan ID corresponding to scenario tax parameters
        taxplanid               = find_taxplanid(scenario);
        
        T_model                 = timing.T_model;
        switch scenario.economy
            case 'steady'
                first_year      = timing.first_transition_year - 1;
            otherwise
                first_year      = timing.first_transition_year;
        end    
        
        bracketsfile    = strcat('PIT_', taxplanid, '.csv' );
        bracketsfile    = fullfile( PathFinder.getTaxPlanInputDir(), bracketsfile );      

        [brackets, rates, burdens] = read_brackets_rates( bracketsfile, first_year, T_model );                               ...
        
        % Convert US dollar amounts into modelunit_dollars
        s.burdens       = burdens  .*scenario.modelunit_dollar;   % Cumulative tax burden
        s.brackets      = brackets .*scenario.modelunit_dollar;   % PIT tax brackets, rem: first one is zero
        s.rates         = rates;                                  % Rate for above each bracket threshold

        % Get the capital and tax treatment allocation params and store them.
        %  Input files CIT_<taxplanid>.CSV expected to have the
        %  following structure:
        %      <Tax variable> as header, <Value> under that header
        %  The tax variable names are defined below.
        filename = strcat('CIT_', taxplanid, '.csv');
        tax_vars = read_tax_vars( filename );
        % Calculate combined tax rate and share for Capital
        s.captaxshare           = tax_vars.shareCapitalCorporate + tax_vars.shareCapitalPreferred; 
        s.taucap                = (   tax_vars.rateCorporate * tax_vars.shareCapitalCorporate   ...
                                    + tax_vars.ratePreferred * tax_vars.shareCapitalPreferred   ...
                                  ) ...
                                  / (tax_vars.shareCapitalCorporate + tax_vars.shareCapitalPreferred );
        s.taucapgain            = 0;
        
        % Pass along all parameters as well
        for f = fieldnames(tax_vars)'
            s.(f{1}) = tax_vars.(f{1});
        end
        
        % Warn if parameters are outside expectations
        if( (s.captaxshare < 0) || (s.captaxshare > 1) )
            fprintf( 'WARNING! captaxshare=%f outside expecations.\n', s.captaxshare );
        end
        if( (s.shareCapitalExpensing < 0) || (s.shareCapitalExpensing > 1) )
            fprintf( 'WARNING! shareCapitalExpensing=%f outside expectations.\n', s.shareCapitalExpensing );
        end
        if( (s.taucap < 0) || (s.taucap > 1) )
            fprintf( 'WARNING! taucap=%f outside expectations.\n', s.taucap );
        end        
        if( (s.taucapgain < 0) || (s.taucapgain > 1) )
            fprintf( 'WARNING! taucapgain=%f outside expectations.\n', s.taucapgain );
        end  
        
        % Catch fatal errors
        if( (s.shareCapitalOrdinary + s.shareCapitalPreferred + s.shareCapitalCorporate) ~= 1 ) 
            error( 'Capital tax shares must sum to 1.' );
        end
        if( (s.shareLaborOrdinary + s.shareLaborPreferred + s.shareLaborCorporate) ~= 1 )
            error( 'Labor tax shares must sum to 1.' );
        end
        if( s.shareLaborOrdinary ~= 1 )
            error( 'Current model does not handle allocating taxable labor income outside Ordinary rates.' );
        end
    end  % tax()

    
    %% DEMOGRAPHICS
    %     Includes:
    %          survival probabilities
    %        , immigrant age distribution
    %        , birth rate
    %        , legal immigration rate
    %        , illegal immigration rate
    function s = demographics( scenario ) 
        
        timing    = ParamGenerator.timing( scenario );
        
        param_dir = PathFinder.getMicrosimInputDir();
        survival  = read_series('SIMSurvivalProbability.csv'        , 1 + timing.realage_entry, param_dir);
        imm_age   = read_series('SIMImmigrantAgeDistribution.csv'   , 1 + timing.realage_entry, param_dir);
        s.surv    = survival';
        s.imm_age = imm_age';

        s.birth_rate   = 0.018923919;    % Annual birth rate
        s.legal_rate   = 0.002371966;    % Annual legal immigration rate
        s.illegal_rate = 0.000836294;    % Annual illegal immigration rate

    end % demographics
    
    
    %% PRODUCTION
    %     Includes:
    %          TFP
    %        , depreciation
    %        , capital share
    function s = production( scenario )
        
        s.A             = 1;                    % Total factor productivity
        s.alpha         = 0.34;                 % Capital share of output
        s.d             = 0.056;                % Actual depreciation rate
        if( scenario.is_low_return )
            s.d = 0.08;                         % "Depreciation rate" to generate r=risk-free rate         
        end

    end % production
        
    
    %% SOCIAL SECURITY
    %
    function s = social_security( scenario )
        
        timing                  = ParamGenerator.timing(scenario);
        T_model                 = timing.T_model;
        T_life                  = timing.T_life;
        first_transition_year   = timing.first_transition_year;
        nstartyears             = length(timing.startyears);
        realage_entry           = timing.realage_entry;
        
        %  OLD STUFF: TBD Revisit and revise
        s.taxcredit     = 0.15;     % Benefit tax credit percentage
        s.ssincmaxs     = repmat(1.185e5*scenario.modelunit_dollar, [1,T_model]); % Maximum income subject to benefit calculation
        s.ssincmins     = zeros(1,T_model);                                       % Minimum income subject to benefit calculation

        % Get T_works (retirement ages)
        nrafile     = strcat(   'NRA_'                                                                          ...
                            ,   find_policy_id( scenario                                                        ...
                                       ,   {'SSNRAPolicy'}                                                      ...
                                       ,   fullfile( PathFinder.getSocialSecurityNRAInputDir(), 'Map.csv' ) )   ...                    ...
                            ,   '.csv' );
        switch scenario.economy
            case 'steady'
                first_year   = first_transition_year - 1;
                survivalprob = ParamGenerator.demographics(scenario).surv;
                T_works      = read_series(nrafile, first_transition_year - (T_life + realage_entry + 1), PathFinder.getSocialSecurityNRAInputDir());
                mass         = ones(T_life,1); for i = 2:T_life; mass(i) = mass(i-1)*survivalprob(i-1); end;
                T_works      = round(sum((mass.*T_works(1:T_life))/sum(mass))) - realage_entry;
            case {'open', 'closed'}
                first_year   = first_transition_year;
                T_works      = read_series(nrafile, first_transition_year - (T_life + realage_entry), PathFinder.getSocialSecurityNRAInputDir());
                T_works      = T_works(1:nstartyears) - realage_entry;
        end
        s.T_works           = T_works;
        
        % Read tax brackets and rates on payroll 
        %   Pad if file years do not go to T_model, truncate if too long
        %   Calculate cumulative liability to speed up calculation 
        %   Convert from US dollars to modelunit dollars
        matchparams     = {'SSTaxPolicy'};
        mapfile         = fullfile( PathFinder.getSocialSecurityTaxInputDir(), 'Map.csv' );
        sstaxid         = find_policy_id( scenario, matchparams, mapfile );
        bracketsfile    = fullfile( PathFinder.getSocialSecurityTaxInputDir()      ...
                                ,   strcat('PayrollTax_'    , sstaxid, '.csv' ) );
        indexingfile    = fullfile( PathFinder.getSocialSecurityTaxInputDir()      ...
                                ,   strcat('BracketsIndex_' , sstaxid, '.csv' ) );

        [brackets, rates, burdens] = read_brackets_rates  ( bracketsfile, first_year, T_model );                               ...
        [indices]                  = read_brackets_indices( indexingfile );
        
        if( size(indices,2) ~= size(brackets,2) )
            throw(MException('social_security:TAXBRACKETS','SSTaxBrackets and BracketsIndexes must have same number of brackets.'));
        end    
    
        s.taxburdens    = burdens  .*scenario.modelunit_dollar;     % Cumulative tax burden
        s.taxbrackets   = brackets .*scenario.modelunit_dollar;     % Payroll tax brackets, rem: first one is zero
        s.taxrates      = rates                               ;     % Rate for above each bracket threshold
        s.taxindices    = indices                             ;     % Type of index to use for the bracket change
        
        % Fetch initial benefits for each cohort 
        %   REM: Benefits are per month in US dollars 
        %        in year = first_transition_year - 1
        first_birthyear = first_year - (T_life + realage_entry);
        matchparams     = {'SSBenefitsPolicy'};
        mapfile         = fullfile( PathFinder.getSocialSecurityBenefitsInputDir(), 'Map.csv' );
        policy_id       = find_policy_id( scenario, matchparams, mapfile );
        
        bracketsfile    = fullfile( PathFinder.getSocialSecurityBenefitsInputDir() ...
                                ,   strcat('InitialBenefits_', policy_id , '.csv' ) );      
        
        [brackets, rates, ~]    = read_brackets_rates( bracketsfile, first_birthyear, nstartyears );                               ...
        
        s.startyear_benefitbrackets   = 12*scenario.modelunit_dollar*brackets;    
        s.startyear_benefitrates      = rates;
        
        % Year-based policy for benefit rates adjustments
        filename = strcat('BenefitsAdjustment_', policy_id , '.csv' );
        
        adjrates = read_series_withpad( filename                                            ...
                                    ,   first_year                                          ...
                                    ,   PathFinder.getSocialSecurityBenefitsInputDir()      ...
                                    ,   first_year + T_model );      
    
        % Verify that adj has same number of brackets
        if( size(adjrates, 2) ~= size(brackets, 2) )
            throw(MException('ParamGenerator.social_security:SIZE','SSBenefits and BenefitsAdjustments must have same number of brackets.'));
        end
        
        s.benefits_adjustment = adjrates;
    end % social_security
    
    
    %% BUDGET AND INTEREST RATES
    %
    function s = budget( scenario )
        
        s = ParamGenerator.timing( scenario );
        first_transition_year   = s.first_transition_year;
        T_model                 = s.T_model;
        

        %  CBO interest rates, expenditures, and debt
        % Input: InterestRates.csv -- interest rate (as pct) 
        %           Format is 
        %           Year	DebtHeldByPublic	EffectiveInterestRateOnDebt	NonInterestDeficit	CPI
        % Input: SIMGDP.csv -- nominal GDP from baseline SIM
        %           Format is (Year), (NominalGDP) w/ header row.
        % Input: SIMRevenues.csv -- nominal tax revenues from baseline SIM
        %           Format is (Year), (Revenues) w/ header row.
        % Input: SIMExpenditures.csv -- nominal non-interest expenditures from baseline SIM
        %           Format is (Year), (Expenditures) w/ header row.
        % Input: CBONonInterestSpending.csv -- NIS as pct GDP
        %           Format is (Year), (PctGDP) w/ header row.
        % Input: CBOSocialSecuritySpending.csv -- SS spending as pct GDP
        %           Format is (Year), (PctGDP) w/ header row.
        % Input: CBOMedicareSpending.csv -- Medicare spending as pct GDP
        %           Format is (Year), (PctGDP) w/ header row.
        % Output: 
        %       debttoout, fedgovtnis, cborates, GEXP_by_GDP
        cbo_param   = PathFinder.getCboInputDir     ();
        sim_param   = PathFinder.getMicrosimInputDir();
        
        first_year                  = first_transition_year - 1;    % first year from which to read series
        SIMGDP                      = read_series( 'SIMGDP.csv'                     , first_year, sim_param );
        SIMRevenues                 = read_series( 'SIMRevenues.csv'                , first_year, sim_param );
        SIMExpenditures             = read_series( 'SIMExpenditures.csv'            , first_year, sim_param );
        SIMGDPPriceIndex            = read_series( 'SIMGDPPriceIndex.csv'           , first_year, sim_param );
        CBONonInterestSpending      = read_series( 'CBONonInterestSpending.csv'     , first_year, cbo_param );
        CBOSocialSecuritySpending   = read_series( 'CBOSocialSecuritySpending.csv'  , first_year, cbo_param );
        CBOMedicareSpending         = read_series( 'CBOMedicareSpending.csv'        , first_year, cbo_param );

        cbo_series      = read_series( 'InterestRates.csv', first_year, cbo_param );
        CBODebt         = cbo_series(:, 1);
        CBORates        = cbo_series(:, 2);
        CBOCPI          = cbo_series(2:end, 4);  % rem: CPI is only for transition path
        
        if( size(SIMGDP) ~= size(SIMRevenues) ...
            | size(SIMGDP) ~= size(SIMExpenditures) ...
            | size(SIMGDP) ~= size(CBORates) ) ...
            throw(MException('read_cbo_parameters:SIZE','Inputs are different sizes.'));
        end;
        if( size(CBONonInterestSpending) ~= size(CBOSocialSecuritySpending) ...
            | size(CBONonInterestSpending) ~= size(CBOMedicareSpending) ...
           )
            throw(MException('read_cbo_spending:SIZE','Inputs are different sizes.'));
        end;
        
        GEXP_by_GDP     = CBONonInterestSpending./100 ...
                        - CBOSocialSecuritySpending./100 ...
                        - CBOMedicareSpending./100 ...
                    ;
        GEXP_by_GDP     = GEXP_by_GDP';
                    
        deficit_nis     = SIMRevenues - SIMExpenditures;
        debt            = zeros(size(deficit_nis));
        debt(1)         = CBODebt(1);
        for i = 2:size(deficit_nis)
            debt(i) = debt(i-1)*(1+CBORates(i)/100.0) - deficit_nis(i);
        end;
        
        growth_deflator      = zeros(size(SIMGDPPriceIndex), 'double');
        growth_deflator(1)   = 1.0;
        for i = 2:size(growth_deflator)
            growth_deflator(i) = SIMGDPPriceIndex(i)/SIMGDPPriceIndex(i-1);
        end;
        
        CBO_rates_growth_adjusted   = ((100.0+CBORates)./growth_deflator)./100.0 - 1.0;    
        deficit_nis_fraction_GDP    = deficit_nis./SIMGDP;                
        debt_percent_GDP            = debt./SIMGDP;                       
        
        % Calculate CPI index -- normalize to 1 for first year
        if( strcmp(scenario.economy, 'steady' ) )
            CBOCPI = 1;
        else
            CBOCPI = CBOCPI ./ CBOCPI(1);
        end
        
        % Name, transpose, truncate vars to correspond to currently used variables 
        % NOTE: The series must go out to T_model or import will break.
        %       Breaking is good here. (Though gracefully would be better)
        s.GEXP_by_GDP       = GEXP_by_GDP;                                % Expenditures as pct gdp
        s.debt              = debt;                                       % Debt as pct gdp
        s.debttoout         = debt_percent_GDP(1);                        % Debt-to-output ratio of initial year
        s.debttoout_trans1  = debt_percent_GDP(2);                        % Debt-to-output ratio of first transition year
        s.fedgovtnis        = deficit_nis_fraction_GDP(2:T_model + 1)';   % starts from first transition path year
        s.cborates          = CBO_rates_growth_adjusted(2:T_model + 1)';  % starts from first transition path year
        s.cbomeanrate       = nanmean(CBO_rates_growth_adjusted(2:end));  % Mean growth-adjusted interest rate over modeling period
        
        % Warn if parameters are outside expectations
        if( (s.cbomeanrate < -0.03) || (s.cbomeanrate > 0.08) )
            fprintf( 'WARNING! cbomeanrate=%f outside expecations.\n', cbomeanrate );
        end
        if( (s.debttoout < 0.5) || (s.debttoout > 1.5) )
            fprintf( 'WARNING! debttoout=%f ouside expecations.\n', debttoout );
        end
        if( max(abs(s.cborates)) > 0.08 )
            fprintf( 'WARNING! cborates outside expectations.\n' );
        end
        if( max(abs(s.fedgovtnis)) > 0.1 )
            fprintf( 'WARNING! fedgovtnis outside expectations.\n' );
        end

        % TAX REVENUE TARGETS (if given)                 
        % Tax revenues as fraction of GDP are loaded from
        % single-series CSV files which contain data by
        % tax plan.
        % Input: Revenues_<taxplanid>.csv -- Estimate of tax
        %           revenues as percent GDP
        %           Format is (Year), (PctRevenues) w/ header row.
        taxplaninputdir = PathFinder.getTaxPlanInputDir();
        taxplanid = find_taxplanid( scenario );
        filename = strcat('Revenues_', taxplanid, '.csv');
        tax_revenue_by_GDP = read_series(filename, first_transition_year, taxplaninputdir );
        tax_revenue_by_GDP = tax_revenue_by_GDP'; 
        if( T_model - length(tax_revenue_by_GDP) < 0 )
            tax_revenue_by_GDP = tax_revenue_by_GDP(1:T_model);
        else
            tax_revenue_by_GDP  = [tax_revenue_by_GDP, ...
                tax_revenue_by_GDP(end)*ones(1, T_model-length(tax_revenue_by_GDP))];
        end
        s.tax_revenue_by_GDP = tax_revenue_by_GDP; 
        s.CPI = CBOCPI';
        
    end % budget
    
    
    %% BEQUEST MOTIVE
    %
    function s = bequest_motive(scenario)
        % phi1 reflects parent's concern about leaving bequests to her children (-9.5 in De Nardi's calibration)
        % phi2 measures the extent to which bequests are a luxury good
        % phi3 is the relative risk aversion coefficient
        
        s.phi1 = scenario.bequest_phi_1;   
        s.phi2 = 11.6;                     
        s.phi3 = 1 - (1 - scenario.sigma)*scenario.gamma;    

    end % bequest_motive
        
    

    
    %% ELASTICITY INVERSION
    % 
    % Invert target elasticities using calibration points
    %   Reusable inverter constructed in the process
    % 
    function [inverse] = invert(targets)
        
        % TODO -- revisit the Calibration process
        %   for now just pick one of two hardcoded targets
        % BEGIN TEMP
        inverse = struct(                             ...
            'beta'              , 0.98600000            , ...
            'gamma'             , 0.75000000            , ...
            'sigma'             , 1.24000000            , ...
            'bequest_phi_1'     , 0                     , ...
            'modelunit_dollar'  , 4.359874681178362e-05   ...
            );

        if( isfield( targets, 'is_low_return' ) )
            if ( targets.is_low_return ) 
                inverse = struct(                               ...
                'beta'              , 1.003341000000000     ,   ...
                'gamma'             , 0.680000000000000     ,   ...
                'sigma'             , 1.500000000000000     ,   ...
                'bequest_phi_1'     , 0                     ,   ...
                'modelunit_dollar'  , 4.135682750000000e-05     ...
                );
            end
        end
        
        return;
        %  END TEMP
        
        % Load calibration points from calibration input directory
        s = load(fullfile(PathFinder.getCalibrationInputDir(), 'calibration.mat'));
        paramv  = s.paramv ;
        targetv = s.targetv;
        solved  = s.solved ;
        
        % Determine list of calibration parameters
        paramlist = fieldnames(paramv)';
        
        % Construct inverse interpolants that map individual target values to calibration parameters
        %   Target values include capital-to-output ratio in addition to elasticities
        for p = paramlist
            interp.(p{1}) = scatteredInterpolant(...
                targetv.captoout(solved)', ...
                targetv.labelas (solved)', ...
                targetv.savelas (solved)', ...
                paramv.(p{1})(solved)', 'nearest');
        end
        
        % Construct elasticity inverter by consolidating inverse interpolants
        %   Capital-to-output ratio target fixed at 3
        function [inverse] = f_(targets)
            captoout = 3;
            for p_ = paramlist
                inverse.(p_{1}) = interp.(p_{1})(captoout, targets.labelas, targets.savelas);
            end
        end
        f = @f_;
        
        % Invert target elasticities if provided
        if exist('targets', 'var'), inverse = f(targets); else, inverse = struct(); end
        
    end
    
    
    % TAXPLAN ID
    function id = getTaxPlanID( scenario )
        id = find_taxplanid( scenario );
    end
    
end % methods

end % class ParamGenerator


%%
%  Helper function to find tax plan ID corresponding to scenario tax parameter values
function [taxplanid] = find_taxplanid( scenario )
    
    % Load tax plan ID map from tax plan input directory
    taxplanidmap = table2struct(readtable(fullfile(PathFinder.getTaxPlanInputDir(), 'Map.csv')));
    
    % Define mapping from dynamic model tax parameter names to tax plan parameter names
    parammap = struct(...
        'base_brackets'                 , 'BaseBrackets'                , ...
        'has_buffet_rule'               , 'HasBuffetRule'               , ...
        'has_double_standard_deduction' , 'HasDoubleStandardDeduction'  , ...
        'has_limit_deductions'          , 'HasLimitDeductions'          , ...
        'has_expand_child_credit'       , 'HasExpandChildCredit'        , ...
        'no_aca_income_tax'             , 'NoACAIncomeTax'              , ...
        'corporate_tax_rate'            , 'CorporateTaxRate'            , ...
        'has_special_pass_through_rate' , 'HasSpecialPassThroughRate'   , ...
        'has_immediate_expensing'       , 'HasImmediateExpensing'       , ...
        'has_repeal_corporate_expensing', 'HasRepealCorporateExpensing' );
    
    % Identify tax plans with parameter values matching scenario tax parameter values
    %   Default values specified for tax plan parameters not represented in dynamic model
    match = arrayfun(@(row) ...
        row.('HasAGISurcharge_5m'         ) == false && ...
        row.('HasPITRateOnCarriedInterest') == false && ...
        all(cellfun(@(param) ...
            isequal(scenario.(param), row.(parammap.(param))), fieldnames(parammap) ...
        )), ...
        taxplanidmap ...
    );
    
    % Check for singular match
    assert(sum(match) > 0, 'No tax plan found with parameter values matching scenario tax parameter values.'           );
    assert(sum(match) < 2, 'More than one tax plan found with parameter values matching scenario tax parameter values.');
    
    % Extract ID of matching tax plan
    taxplanid = num2str(taxplanidmap(match).ID);
    
end


%%
%  Helper function to find a Policy ID corresponding to scenario parameter values
%        matchparams : cell array of param names to match
%        mapfile     : fullfile name of map file with format 
%           (ID) (Param1) (Param2) ... (ParamN)
function [id] = find_policy_id( scenario, matchparams, mapfile )
    
    % Load plan ID map from input directory
    map = table2struct(readtable(mapfile));
    
    % Identify policies with parameter values matching scenario parameter values
    match = arrayfun(@(row) ...
        all(cellfun(@(param) ...
            isequal(scenario.(param), row.(param)), matchparams) ...
        ), ...
        map ...
    );
    
    % Check for singular match
    assert(sum(match) > 0, 'No ID found with parameter values matching scenario parameter values.'           );
    assert(sum(match) < 2, 'More than one ID found with parameter values matching scenario parameter values.');
    
    % Extract ID of matching plan
    id = num2str(map(match).ID);
    
end %find_policy_id


%%
%  Helper function to read CSV file with format:
%     Header has variable names, then one row of values
function [tax_vars] = read_tax_vars( filename )

    warning( 'off', 'MATLAB:table:ModifiedVarnames' );          % for 2016b
    warning( 'off', 'MATLAB:table:ModifiedAndSavedVarnames' );  % for 2017a
    
    filepath = fullfile(PathFinder.getTaxPlanInputDir(), filename);
    if ~exist(filepath, 'file')
        err_msg = strcat('Cannot find file = ', strrep(filepath, '\', '\\'));
        throw(MException('read_tax_vars:FILENAME', err_msg ));
    end;
        
    % Expected format is
    %    rateCorporate	ratePreferred	rateOrdinary	expensingShare	shareCapitalCorporate	shareCapitalPreferred	shareCapitalOrdinary	shareLaborCorporate	shareLaborPreferred	shareLaborOrdinary

    T           = readtable( filepath );
    tax_vars    = table2struct(T);

end % read_tax_vars()


%%
%  Helper function to read CSV files with format: 
%    (Year), (Bracket1), ... (BracketN), (Rate1), ... (RateN)
%    filename     : fullfile of CSV to read
%    first_year   : don't read years before this param
%    T_years      : read this many years 
function [brackets, rates, burdens] = read_brackets_rates( ...
                                        filename, first_year, T_years )

    warning( 'off', 'MATLAB:table:ModifiedVarnames' );          % for 2016b
    warning( 'off', 'MATLAB:table:ModifiedAndSavedVarnames' );  % for 2017a

    if ~exist(filename, 'file')
        err_msg = strcat('Cannot find file = ', strrep(filename, '\', '\\'));
        throw(MException('read_brackets_rates:FILENAME', err_msg ));
    end;
        
    T           = readtable(filename);
    T_arr       = table2array(T);
    
    % Find first year
    years       = T_arr(:,1);
    year_start  = find( years == first_year, 1);
    if( isempty(year_start) )
        throw(MException('read_brackets_rates:FIRSTINDEX','Cannot find first index in file.'));
    end    
    
    % Find all brackets
    % Enforce that first bracket must be zero.
    brackets    = T_arr(year_start:end,  contains(T.Properties.VariableNames, 'Bracket' ) );
    if( all(brackets(:,1)) > 0 )
        err_msg = strcat('First bracket must be 0 in file ', strrep(filepath, '\', '\\'));
        throw(MException('read_brackets_rates:BRACKET0', err_msg ));
    end 
    
    % Find all rates
    rates       = T_arr(year_start:end,  contains(T.Properties.VariableNames, 'Rate' ) );

    % Pad brackets and rates if not long enough, truncate if too long
    num_years = size(brackets,1);
    if( T_years - num_years <= 0 )
        brackets    = brackets(1:T_years,:);
        rates       = rates(1:T_years,:);
    else
        brackets    = [brackets; ...
            repmat(brackets(end,:)  , [T_years-num_years, 1])   ];
        rates       = [rates; ...
            repmat(rates(end,:)     , [T_years-num_years, 1])   ];
    end
    
    % Calculate cumulative tax burdens along brackets dimension
    burdens         = cumsum(diff(brackets, 1, 2).*rates(:, 1:end-1), 2); 
    burdens         = [zeros(size(brackets, 1), 1), burdens];  % rem: first burden is zero
end % read_brackets_rates()


%%
% Read a CSV file in format (Bracket1Index),...(BracketNIndex)
%    This file defines the indexing methodology to use for the brackets.
function [indices] = read_brackets_indices( filename )

    warning( 'off', 'MATLAB:table:ModifiedVarnames' );          % for 2016b
    warning( 'off', 'MATLAB:table:ModifiedAndSavedVarnames' );  % for 2017a

    if ~exist(filename, 'file')
        err_msg = strcat('Cannot find file = ', strrep(filename, '\', '\\'));
        throw(MException('read_brackets_indices:FILENAME', err_msg ));
    end;
        
    T       = readtable(filename);
    indices = table2array(T(1, :));
    
    % Validated that indices are in the allowed set:
    %    reals, nominals, wage_inflations, cohort_wages
    
    
end % read_brackets_indices


%%
% Read a CSV file in format (Index), (Value1), (Value2), ... (ValueN)
%       first_index     : first index to read (previous part of file is ignored
%       last_index      : (optional) If given, copy the last available
%                       value so that series goes from first_index ...
%                       last_index. If the original series is too long,
%                       truncate.
%    For time series, (Index) is (Year), 
%    For other series (e.g. age_survival_probability, (Index) is (Age)
function [series] = read_series_withpad(filename, first_index, param_dir, last_index )

    warning( 'off', 'MATLAB:table:ModifiedVarnames' );          % for 2016b
    warning( 'off', 'MATLAB:table:ModifiedAndSavedVarnames' );  % for 2017a
 
    % Check if file exists 
    filepath    = fullfile(param_dir, filename);
    if ~exist(filepath, 'file')
        err_msg = strcat('Cannot find file = ', strrep(filepath, '\', '\\'));
        throw(MException('read_series:FILENAME', err_msg ));
    end;
        
    T           = readtable(filepath);
    indices     = table2array(T(:,1));
    vals        = table2array(T(:,2:end));
    
    idx_start   = find( indices == first_index, 1);
    if( isempty(idx_start) )
        throw(MException('read_series:FIRSTINDEX','Cannot find first index in file.'));
    end;
    
    series      = vals(idx_start:end, : );

    if( ~isempty(last_index) )
        % Pad or truncate
        num_indices     = size(series,1);
        num_required    = last_index - first_index;

        if( num_required - num_indices <= 0 )
            series    = series(1:num_required,:);
        else
            series    = [series; ...
                repmat(series(end,:)  , [num_required-num_indices, 1])   ];
        end
    end;

end % read_series_withpad


%  Easier signature for read_series
function [series] = read_series(filename, first_index, param_dir )
    series = read_series_withpad( filename, first_index, param_dir, [] );
end % read_series





%%  END FILE

