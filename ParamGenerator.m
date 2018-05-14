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
        s.T_life                = 80    ;   % Death year
        s.Tmax_work             = 52    ;   % Last possible year of working age     
        
        % Time range for transition path
        s.TransitionFirstYear = scenario.TransitionFirstYear;
        s.TransitionLastYear  = scenario.TransitionLastYear;
        
        switch scenario.economy
            case 'steady'
                s.T_model    = 1;                           % Steady state total modeling years
                s.startyears = 0;                           % Steady state cohort start year
            case {'open', 'closed'}
                s.T_model    = s.TransitionLastYear - s.TransitionFirstYear;
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
        filename  = fullfile(PathFinder.getMicrosimInputDir(), 'ConesaEtAl_WageAgeProfile.csv' );
        series    = read_series(filename, 'Age', 1 + timing.realage_entry, []);
        zage      = series.Wage;
        
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
        
        %
        % TODO: We are working toward replacing the shock generation with
        %       a calibrated (from PWBMsim data) transition.
        %       As an incremental step, we move to a larger grid structure
        %       and copy the current transz for all years/types.
        %
        
        % Make the larger sized transz grid and replicate the current transz grid.
        transz      = [transz, zeros(nz); zeros(nz), transz];
        transitions = repmat(transz, [1, 1, T_life]);
        zs          = [zs(:,:,1); zs(:,:,2)];
        DISTz       = [DISTz; DISTz] / 2;
        
        s.ndem     = ndem;
        s.g        = g;
        s.ng       = ng;
        s.nz       = ndem*nz;
        s.transz   = transitions;
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
        taxplanmapfile = fullfile(PathFinder.getTaxCalculatorInputDir(), 'Map.csv');
        taxplanid   = find_policy_id(taxplanmapfile, scenario                         , {'TaxCode'}, 'ID');
        sstaxplanid = find_policy_id(taxplanmapfile, scenario.steady().currentPolicy(), {'TaxCode'}, 'ID');
                
        T_model                 = timing.T_model;
        switch scenario.economy
            case 'steady'
                first_year      = timing.TransitionFirstYear - 1;
                last_year       = first_year;
            otherwise
                first_year      = timing.TransitionFirstYear;
                last_year       = timing.TransitionLastYear - 1;
        end    
        
        file = fullfile(PathFinder.getTaxCalculatorInputDir(), strcat('OrdinaryRates_', taxplanid, '.csv' ));
        [brackets, rates] = read_brackets_rates_indices(file, first_year, T_model);
        
        % TBD: This should be done in ModelSolver as for SocialSecurity
        % Calculate cumulative tax burdens along brackets dimension
        burdens         = cumsum(diff(brackets, 1, 2).*rates(:, 1:end-1), 2); 
        burdens         = [zeros(size(brackets, 1), 1), burdens];  % rem: first burden is zero
        
        % Convert US dollar amounts into modelunit_dollars
        s.burdens       = burdens  .*scenario.modelunit_dollar;   % Cumulative tax burden
        s.brackets      = brackets .*scenario.modelunit_dollar;   % PIT tax brackets, rem: first one is zero
        s.rates         = rates;                                  % Rate for above each bracket threshold

        % Get the capital and tax treatment allocation params. 
        %    Rem: These are time-varying, so read results are vectors.
        file = fullfile(PathFinder.getTaxCalculatorInputDir(), strcat('CapitalTaxes_', taxplanid, '.csv'));
        tax_vars = read_series(file, 'year', first_year, last_year);
        
        % Portion of corporate income taxed at preferred rates
        %   TBD: Fix this when John R. gives new inputs
        s.captaxshare   = tax_vars.shareCapitalPreferred ...
                         ./ (tax_vars.shareCapitalPreferred + tax_vars.shareCapitalOrdinary);
        
        % Capital gains tax == zero for now
        s.rateCapGain   = zeros(T_model,1);
        
        % Pass along all CIT parameters 
        for f = fieldnames(tax_vars)'
            s.(f{1}) = tax_vars.(f{1});
        end
        
        % Read Preferred Rates
        file = fullfile(PathFinder.getTaxCalculatorInputDir(), strcat('PreferredRates_', taxplanid, '.csv'));
        [brackets, rates] = read_brackets_rates_indices(file, first_year, T_model);
        
        % TBD: This should be done in ModelSolver as for SocialSecurity
        % Calculate cumulative tax burdens along brackets dimension
        burdens         = cumsum(diff(brackets, 1, 2).*rates(:, 1:end-1), 2); 
        burdens         = [zeros(size(brackets, 1), 1), burdens];  % rem: first burden is zero
        
        % Convert US dollar amounts into modelunit_dollars
        s.prefburdens   = burdens  .*scenario.modelunit_dollar;   % Cumulative tax burden
        s.prefbrackets  = brackets .*scenario.modelunit_dollar;   % Preferred Rate tax brackets, rem: first one is zero
        s.prefrates     = rates;                                  % Rate for above each bracket threshold
        
        % Foreign tax withholding -- TEMP
        %   TBD: Read these from file
        s.rateForeignCorpIncome         = 0.179 ;   % Tax rate on corp. distributions to foreigners
        s.rateForeignPassThroughIncome  = 0.10  ;   % Tax rate on pass-through distributions to foreigners
        
        % Allocation of capital income between Corp & Pass-Through -- TEMP
        %   TBD: Read these from file
        s.shareIncomeCorp               = 1.00;  % Amount of capital to corps
        s.shareIncomePassThrough        = 1 - s.shareIncomeCorp;
        
        % Calculate Q-Tobin
        switch scenario.economy
            case 'steady'
                s.qtobin0   = 1 - tax_vars.shareCapitalExpensing(1) * s.rateCorporate(1);
                s.qtobin    = s.qtobin0;
            otherwise
                % Read tax params from steady-state, current-policy to make t=0 values
                file = fullfile(PathFinder.getTaxCalculatorInputDir(), strcat('CapitalTaxes_', sstaxplanid, '.csv'));
                sstax_vars = read_series(file, 'year', first_year - 1, first_year - 1);

                s.qtobin0   = 1 - sstax_vars.shareCapitalExpensing(1) * sstax_vars.rateCorporate(1);
                s.qtobin    = 1 - tax_vars.shareCapitalExpensing .* s.rateCorporate;
        end  
               
        
        % Warn if parameters are outside expectations
        if( any(s.captaxshare < 0) || any(s.captaxshare > 1) )
            fprintf( 'WARNING! captaxshare outside expecations.\n', s.captaxshare );
        end
        if( any(s.shareCapitalExpensing < 0) || any(s.shareCapitalExpensing > 1) )
            fprintf( 'WARNING! shareCapitalExpensing=%f outside expectations.\n', s.shareCapitalExpensing );
        end
        
        %% TBD: Enlarge the space of error checks
        
        % Catch fatal errors
        if( any(s.shareCapitalOrdinary + s.shareCapitalPreferred + s.shareCapitalCorporate) ~= 1 ) 
            error( 'Capital tax shares must sum to 1.' );
        end
        if( any(s.shareLaborOrdinary + s.shareLaborPreferred + s.shareLaborCorporate) ~= 1 )
            error( 'Labor tax shares must sum to 1.' );
        end
        if( any(s.shareLaborOrdinary ~= 1) )
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
        
        filename  = fullfile( PathFinder.getMicrosimInputDir(), 'SurvivalProbability.csv' );
        series    = read_series( filename, 'Age', 1 + timing.realage_entry, [] );
        survival  = series.SurvivalProbability;
        if( length(survival) ~= timing.T_life )
            throw(MException('timing:survival','Survival probabilities must exist for every age.'));
        end
        
        filename  = fullfile( PathFinder.getMicrosimInputDir(), 'ImmigrantAgeDistribution.csv' );
        series    = read_series( filename, 'Age', 1 + timing.realage_entry, [] );
        imm_age   = series.ImmigrantAgeDistribution;
        if( length(imm_age) ~= timing.T_life )
            throw(MException('timing:imm_age','Immigrant age distributions must exist for every age.'));
        end
        
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
        
        % TBD: These are hard coded for now, until we resolve sourcing
        s.A             = 1;                    % Total factor productivity
        s.alpha         = 0.34;                 % Capital share of output
        s.risk_premium  = 0;
        s.depreciation  = 0.056;                % Actual depreciation rate
        if( scenario.IsLowReturn )
            s.risk_premium  = 0.08 - s.depreciation; % "Depreciation rate" to generate r=risk-free rate         ;
        end
        
    end % production
        
    
    %% SOCIAL SECURITY
    %
    function s = social_security( scenario )
        
        timing                  = ParamGenerator.timing(scenario);
        T_model                 = timing.T_model;
        T_life                  = timing.T_life;
        first_transition_year   = timing.TransitionFirstYear;
        nstartyears             = length(timing.startyears);
        realage_entry           = timing.realage_entry;
        
        %  OLD STUFF: TBD Revisit and revise
        s.taxcredit = 0.15;     % Benefit tax credit percentage

        % Get OASIcalculator scenario ID
        id = find_policy_id(fullfile(PathFinder.getOASIcalculatorInputDir(), 'map.csv'), scenario, {
            'TaxRate', 'TaxMax', 'DonutHole', 'COLA', 'PIA', 'NRA', 'CreditEarningsAboveTaxMax', 'FirstYear', 'GradualChange'
        }, 'id_OASIcalculator');
        
        % Get T_works (retirement ages)
        nrafile = fullfile(PathFinder.getOASIcalculatorInputDir(), strcat('retirementAges_', id, '.csv'));
        
        switch scenario.economy
            case 'steady'
                first_year   = first_transition_year - 1;
                survivalprob = ParamGenerator.demographics(scenario).surv;
                series       = read_series( nrafile, 'birthYear', first_year - (T_life + realage_entry), [] );
                T_works      = round(series.NRA);
                mass         = ones(T_life,1); for i = 2:T_life; mass(i) = mass(i-1)*survivalprob(i-1); end
                T_works      = round(sum((mass.*T_works(1:T_life))/sum(mass))) - realage_entry;
            case {'open', 'closed'}
                first_year   = first_transition_year;
                series       = read_series( nrafile, 'birthYear', first_year - (T_life + realage_entry), [] );
                T_works      = round(series.NRA);
                T_works      = T_works(1:nstartyears) - realage_entry;
        end
        s.T_works           = T_works;
        retire_years        = zeros(nstartyears, 1);
        for i = 1:nstartyears
            retire_years(i)  =  (i - T_life) + T_works(i);
        end
        s.retire_years      = retire_years;
        
        %  TAXATION
        
        % Read tax brackets and rates on payroll 
        %   Pad if file years do not go to T_model, truncate if too long
        %   Calculate cumulative liability to speed up calculation 
        %   Convert from US dollars to modelunit dollars
        file = fullfile(PathFinder.getOASIcalculatorInputDir(), strcat('TaxParameters_', id, '.csv' ));
        [brackets, rates, indices] = read_brackets_rates_indices(file, first_year, T_model);
        
        if( size(indices,2) ~= size(brackets,2) )
            throw(MException('social_security:TAXBRACKETS','SSTaxBrackets and BracketsIndexes must have same number of brackets.'));
        end    
    
        s.taxbrackets   = brackets .*scenario.modelunit_dollar;     % Payroll tax brackets, rem: first one is zero
        s.taxrates      = rates                               ;     % Rate for above each bracket threshold
        s.taxindices    = indices                             ;     % Type of index to use for the bracket change
        
        % BENEFITS
        
        % Get range of income which is credited toward benefit calculation
        file = fullfile(PathFinder.getOASIcalculatorInputDir(), strcat('BenefitParameters_', id, '.csv' ));
        [brackets] = read_brackets_rates_indices(file, first_year, T_model);
        
        s.ssincmins = (brackets(:, 1) * scenario.modelunit_dollar)';
        s.ssincmaxs = (brackets(:, 2) * scenario.modelunit_dollar)';
        
        % Fetch initial benefits for each cohort 
        %   REM: Benefits are per month in US dollars 
        %        in year = first_transition_year - 1
        first_birthyear = first_year - (T_life + realage_entry);
        
        file = fullfile(PathFinder.getOASIcalculatorInputDir(), strcat('PIAParameters_', id, '.csv' ));
        [brackets, rates] = read_brackets_rates_indices(file, first_birthyear, nstartyears);
        
        s.startyear_benefitbrackets   = 12*scenario.modelunit_dollar*brackets;
        s.startyear_benefitrates      = rates;
        
    end % social_security
    
    
    %% BUDGET AND INTEREST RATES
    %
    function s = budget( scenario )
        
        timing                  = ParamGenerator.timing( scenario );
        first_transition_year   = timing.TransitionFirstYear;
        T_model                 = timing.T_model;
        
        switch scenario.economy
            case 'steady'
                first_year  = first_transition_year - 1;    
                last_year   = first_year;
            otherwise
                first_year  = first_transition_year;
                last_year   = first_year + T_model - 1;
        end
        
        
        projections_file = fullfile(PathFinder.getProjectionsInputDir(), 'Projections.csv');
        projections_series      = read_series(projections_file, 'Year', first_year               , last_year                );
        projections_past_series = read_series(projections_file, 'Year', first_transition_year - 1, first_year               );
        projections_debt_series = read_series(projections_file, 'Year', first_transition_year - 1, first_transition_year - 1);
        projections_full_series = read_series(projections_file, 'Year', first_year               , []                       );
        
        taxcalculator_id = find_policy_id(fullfile(PathFinder.getTaxCalculatorInputDir(), 'Map.csv'), scenario, {'TaxCode'}, 'ID');
        taxcalculator_file = fullfile(PathFinder.getTaxCalculatorInputDir(), strcat('Aggregates_', taxcalculator_id, '.csv'));
        taxcalculator_series      = read_series(taxcalculator_file, 'Year', first_year               , last_year );
        taxcalculator_past_series = read_series(taxcalculator_file, 'Year', first_transition_year - 1, first_year);
        
        oasicalculator_id = find_policy_id(fullfile(PathFinder.getOASIcalculatorInputDir(), 'map.csv'), scenario, {
            'TaxRate', 'TaxMax', 'DonutHole', 'COLA', 'PIA', 'NRA', 'CreditEarningsAboveTaxMax', 'FirstYear', 'GradualChange'
        }, 'id_OASIcalculator');
        oasicalculator_file = fullfile(PathFinder.getOASIcalculatorInputDir(), strcat('aggregate_', oasicalculator_id, '.csv'));
        oasicalculator_series      = read_series(oasicalculator_file, 'year', first_year               , last_year );
        oasicalculator_past_series = read_series(oasicalculator_file, 'year', first_transition_year - 1, first_year);
        
        f_revenues = @(t, o) ...
            ...
            + t.RevenuesIndividualIncomeTax         ...
            + t.RevenuesPayrollTaxExSocialSecurity  ...
            + t.RevenuesCorporateIncomeTax          ...
            + t.RevenuesEstateAndGiftTaxes          ...
            + t.RevenuesExciseTaxes                 ...
            + t.RevenuesCustomsDuties               ...
            + t.RevenuesMiscellaneousReceipts       ...
            ...
            + o.RevenuesPayrollTaxSocialSecurity    ...
            ...
        ;
        revenues      = f_revenues(taxcalculator_series     , oasicalculator_series     );
        past_revenues = f_revenues(taxcalculator_past_series, oasicalculator_past_series);
        
        f_noninterestspending = @(p, t, o) ...
            ...
            + p.OutlaysDiscretionary        ...
            + p.OutlaysMedicare             ...
            + p.OutlaysMedicaid             ...
            + p.OutlaysFederalRetirement    ...
            + p.OutlaysVeteransPrograms     ...
            + p.OutlaysOtherPrograms        ...
            + p.OutlaysOffsettingReceipts   ...
            ...
            + t.OutlaysIncomeSecurity       ...
            ...
            + o.OutlaysSocialSecurity       ...
            ...
        ;
        noninterestspending      = f_noninterestspending(projections_series     , taxcalculator_series     , oasicalculator_series     );
        past_noninterestspending = f_noninterestspending(projections_past_series, taxcalculator_past_series, oasicalculator_past_series);
        
        % Calculate debt to get DEBT/GDP for first_year
        deficit_nis     = past_revenues - past_noninterestspending;
        debt            = zeros(size(deficit_nis));
        debt(1)         = projections_debt_series.DebtHeldByThePublic(1);
        for i = 2:size(deficit_nis)
            debt(i) = debt(i-1)*(1+projections_past_series.AverageInterestRateOnDebt(i)/100.0) - deficit_nis(i);
        end
        s.debttoout     = debt(end) / projections_series.GDP_FY(1);
        
        % Rates
        %    For interest rate in steady-state, we use avg. rate across all data
        %    NOTE: EffectiveInterestRateOnDebt is in NOMINAL terms and we
        %    deflate by GDPPriceIndex
        if( strcmp(scenario.economy, 'steady') )
            gdpPriceIndex   = projections_full_series.ChainedCPIU;
            interest_rate   = projections_full_series.AverageInterestRateOnDebt / 100;
        else
            gdpPriceIndex   = projections_series.ChainedCPIU;
            interest_rate   = projections_series.AverageInterestRateOnDebt / 100;
        end
        
        deflator        = zeros(size(gdpPriceIndex));
        deflator(1)     = 1.0;
        for i = 2:size(deflator)
            deflator(i) = gdpPriceIndex(i)/gdpPriceIndex(i-1);
        end
        rates_adjusted  = ((1 + interest_rate)./deflator) - 1.0;    
        
        if( strcmp(scenario.economy, 'steady') )
            s.debtrates = nanmean( rates_adjusted );
        else
            s.debtrates = rates_adjusted';
        end
       
        % Consumption good price index
        s.CPI = (projections_series.ChainedCPIU / projections_series.ChainedCPIU(1))'; % normalize to 1 for first_year

        % TAX REVENUE AND EXPENDITURE TARGETS           
        s.tax_revenue_by_GDP = (revenues            ./ projections_series.GDP_FY)';
        s.outlays_by_GDP     = (noninterestspending ./ projections_series.GDP_FY)';
        
        
        
        % WARNINGS if parameters are outside expectations
        if( any( abs(s.debtrates) > 0.05 ) )
            fprintf( 'WARNING! debtrates outside expectations.\n' );
        end
        if( (s.debttoout < 0.6) || (s.debttoout > 1.0) )
            fprintf( 'WARNING! debttoout=%f outside expectations.\n', debttoout );
        end
        for t = 2:T_model
            if( abs((s.CPI(t)/s.CPI(t-1))-1 > 0.05 ) )
                fprintf( 'WARNING! cpi outside expectations. \n' );
            end
        end
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
        inverse = struct(                                 ...
            'beta'              , 0.976000000000000     , ...
            'gamma'             , 0.779000000000000     , ...
            'sigma'             , 1.600000000000000     , ...
            'bequest_phi_1'     , 0                     , ...
            'modelunit_dollar'  , 3.952983000000000e-05   ...
            );

        if( isfield( targets, 'is_low_return' ) )
            if ( targets.is_low_return ) 
                inverse = struct(                               ...
                'beta'              , 0.99400000     ,   ...
                'gamma'             , 0.72600000     ,   ...
                'sigma'             , 1.50000000     ,   ...
                'bequest_phi_1'     , 0                     ,   ...
                'modelunit_dollar'  , 3.85000000e-05     ...
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
    
    
end % methods

end % class ParamGenerator



%%
%  Helper function to find a policy ID corresponding to a given scenario's parameter values
%        mapfile     : Name of map file mapping parameter values to policy IDs
%        scenario    : Scenario of interest
%        matchparams : Cell array of parameter names identifying parameters to use for matching
%        idcolumn    : Name of map file column containing policy IDs
function [id] = find_policy_id(mapfile, scenario, matchparams, idcolumn)
    
    % Load plan ID map from input directory
    map = table2struct(readtable(mapfile));
    
    % Identify policies with parameter values matching scenario parameter values
    match = arrayfun( ...
        @(row) all(cellfun( ...
            @(param) isequal(scenario.(param), row.(param)), ...
            matchparams) ...
        ), ...
        map ...
    );
    
    % Check for singular match
    assert(sum(match) > 0, 'No ID found with parameter values matching scenario parameter values.'           );
    assert(sum(match) < 2, 'More than one ID found with parameter values matching scenario parameter values.');
    
    % Extract ID of matching plan
    id = num2str(map(match).(idcolumn));
    
end %find_policy_id


%%
%  Helper function to read CSV file with format:
%     Header has variable names, then one row of values
function [vars] = read_vars( filename )

    warning( 'off', 'MATLAB:table:ModifiedVarnames' );          % for 2016b
    warning( 'off', 'MATLAB:table:ModifiedAndSavedVarnames' );  % for 2017a
    
    if ~exist(filename, 'file')
        err_msg = strcat('Cannot find file = ', strrep(filename, '\', '\\'));
        throw(MException('read_vars:FILENAME', err_msg ));
    end;
        
    T           = readtable( filename );
    vars        = table2struct(T);

end % read_vars()


%%
%  Helper function to read CSV files with format: 
%    (Year), (Bracket1), ... (BracketN), (Rate1), ... (RateN)
%    filename     : fullfile of CSV to read
%    first_year   : don't read years before this param
%    T_years      : read this many years 
function [brackets, rates, indices] = read_brackets_rates_indices(filename, first_year, T_years)

    warning( 'off', 'MATLAB:table:ModifiedVarnames' );          % for 2016b
    warning( 'off', 'MATLAB:table:ModifiedAndSavedVarnames' );  % for 2017a

    if ~exist(filename, 'file')
        err_msg = strcat('Cannot find file = ', strrep(filename, '\', '\\'));
        throw(MException('read_brackets_rates:FILENAME', err_msg ));
    end
        
    T = readtable(filename);
    
    % Find first year
    years       = table2array(T(:,1));
    year_start  = find( years == first_year, 1);
    if( isempty(year_start) )
        throw(MException('read_brackets_rates:FIRSTINDEX','Cannot find first index (%u) in file (%s).', first_year, filename));
    end    
    
    % Find all brackets, rates, and indices
    brackets = table2array(T(year_start:end, contains(T.Properties.VariableNames, 'Bracket')));
    rates    = table2array(T(year_start:end, contains(T.Properties.VariableNames, 'Rate'   )));
    indices  = table2array(T(1             , contains(T.Properties.VariableNames, 'Index'  )));
    
    % Enforce that first bracket must be zero.
    % If there are no brackets, make all zeros brackets
    if( isempty(brackets) )
        brackets = zeros(size(rates));
    end
    if( all(brackets(:,1)) > 0 )
        err_msg = strcat('First bracket must be 0 in file ', strrep(filepath, '\', '\\'));
        throw(MException('read_brackets_rates:BRACKET0', err_msg ));
    end
    
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
    
end % read_brackets_rates_indices()


%%
% Read a CSV file in format
%    header: IndexName, VarName1, VarName2, ... VarNameN
%    data  : (Index)  , (Value1), (Value2), ... (ValueN)
%       index_name      : string name of index variable -- i.e. IndexName
%       first_index     : first index to read (previous part of file is ignored
%       last_index      : (optional) If given, copy the last available
%                       value so that series goes from first_index ...
%                       last_index. If the original series is too long,
%                       truncate.
%    For time series, (Index) is (Year), 
function [series] = read_series(filename, index_name, first_index, last_index )

    warning( 'off', 'MATLAB:table:ModifiedVarnames' );          % for 2016b
    warning( 'off', 'MATLAB:table:ModifiedAndSavedVarnames' );  % for 2017a
 
    % Check if file exists 
    if ~exist(filename, 'file')
        err_msg = strcat('Cannot find file = ', strrep(filename, '\', '\\'));
        throw(MException('read_series:FILENAME', err_msg ));
    end
        
    T = readtable(filename, 'TreatAsEmpty', {'NA'});
    
    idx_drop    = ( T.(index_name) < first_index );
    if( all(idx_drop) )
        throw(MException('read_series:FIRSTINDEX','Cannot find first index in file.'));
    end
    
    % Remove unused table rows
    T( idx_drop, : ) = [];
    
    if( ~isempty(last_index) )
        % Truncate if needed
        idx_drop = ( T.(index_name) > last_index );
        T( idx_drop, : ) = [];
        
        % Pad if needed 
        num_add  = last_index - T.(index_name)(end);
        if( num_add > 0 )
            T    = [T; repmat(T(end,:), [num_add, 1])];
        end
    end
    
    series = table2struct( T, 'ToScalar', true );
    
    % Do not return index column
    series = rmfield( series, index_name );
 
end % read_series



%%  END FILE

