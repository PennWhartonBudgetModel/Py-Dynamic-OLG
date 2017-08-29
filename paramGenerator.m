%%
% Reader and generator of parameters for dynamic model
% 
%%


classdef paramGenerator

properties (Constant)
end % properties

methods (Static)

    %% TIMING
    % 
    function s = timing()
        s.T_life    = 80;
        s.T_work    = 47;
        s.T_model   = 25;
    end % timing
       
    %% DICRETIZATION GRIDS
    %   pass in 
    %     T_life     : years of life 
    %     prem_legal : productivity premium of avg. (by year) legal
    %          immigrant. Note: This is a ratio.
    function s = grids( T_life, prem_legal )
        % (Shock process definitions from Storesletten, Telmer, and Yaron, 2004)
        % (Method of discretizing shock process derived from Adda and Cooper, 2003)
        
        % Define function to compute cutoff points for equally weighted slices of a normal distribution
        f = @(n, sig) n*sig*-diff(normpdf(norminv(linspace(0, 1, n+1))));
        
        % Determine permanent and transitory shocks
        nperm  = 2; zperm  = f(nperm , sqrt(0.2105));
        ntrans = 2; ztrans = f(ntrans, sqrt(0.0630)); 
        
        % Determine persistent shocks
        npers = 2;
        
        pep = 0.973;                        % Lagged productivity coefficient
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
        % (Estimated coefficients from Barro and Barnes Medicaid working paper)
        zage = polyval([-5.25e-7, 1.05e-4, -8.1467e-3, 0.2891379, -1.203521], 19+(1:T_life));
        
        % Calculate total productivity
        ndem = nperm; nz = ntrans*npers;
        
        zs = max(0, repmat(reshape(zage                       , [1 ,T_life,1   ]), [nz,1     ,ndem]) ...
                  + repmat(reshape(zperm                      , [1 ,1     ,ndem]), [nz,T_life,1   ]) ...
                  + repmat(reshape(kron(ones(1,npers), ztrans) ...
                                 + kron(zpers, ones(1,ntrans)), [nz,1     ,1   ]), [1 ,T_life,ndem]));
        
        % Construct Markov transition matrix for all shocks / productivities
        transz = kron(transpers, (1/ntrans)*ones(ntrans,ntrans));
        
        % Determine initial distribution over all shocks / productivities
        DISTz0 = kron(DISTpers , (1/ntrans)*ones(1     ,ntrans));
        
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
        for age = 1:T_life
            
            zmean = mean(zs(:,age,:), 3);
            
            zlegal   = sum(zmean .* DISTz(:,age,g.legal  )) * prem_legal;
            zillegal = sum(zmean .* DISTz(:,age,g.illegal)) * 0.9       ;
            
            plegal   = (zmean(nz) - zlegal  ) / (zmean(nz)*(nz-1) - sum(zmean(1:nz-1)));
            pillegal = (zmean(1)  - zillegal) / (zmean(1) *(nz-1) - sum(zmean(2:nz  )));
            
            DISTz(:,age,g.legal  ) = [plegal*ones(nz-1,1); 1 - plegal*(nz-1)    ];
            DISTz(:,age,g.illegal) = [1 - pillegal*(nz-1); pillegal*ones(nz-1,1)];
            
        end
        
        % Define savings and average earnings discretization vectors
        % (Upper bound of average earnings defined as maximum possible Social Security benefit)
        f = @(lb, ub, n) lb + (ub-lb)*((0:n-1)/(n-1))'.^2;
        nk = 10; kv = f(1e-3, 120           , nk);
        nb =  5; bv = f(0   , 1.5*max(zs(:)), nb);
        
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
    function s = tax( taxplan )
        % Get PIT tax rates.
        %  Input files TPCPIT_<taxplan>.CSV expected to have the
        %  following structure:
        %      <Income $>, <Marginal Tax Rate>
        %  The marginal tax rate represents the rate charged 
        %  above the listed income threshold (and implicitly
        %  until the next income threshold, if extant).
        %  First threshold must be zero to make explicit the
        %  fact that tax is charged on income from zero to some $ amount.
        filename    = strcat('TPCPIT_', taxplan, '.csv');
        [income_thresholds, tax_rates] = read_tax_rates( filename );

        % Calculate tax liability at each threshold.
        %  REM: Last tax rate is for last threshold to Inf, so do not
        %  apply.
        tax_        = cumsum(diff(income_thresholds).*tax_rates(1:end-1)); 
        income_tax  = [0; tax_];

        % No deductions, re-calculate effective tax if existent
        %  NOTE: If the deduction thresholds are different, we need to make
        %        a single vector of the union of the two threshold vectors
        %        and calculate effective income tax at those points

        % Store income tax and deduction functions.
        %   REM: Transpose into row vectors
        s.tax_thresholds      = income_thresholds';
        s.tax_burden          = income_tax';
        s.tax_rates           = tax_rates';

        % Get the capital tax params and store them.
        %  Input files TPCCIT_<taxplan>.CSV expected to have the
        %  following structure:
        %      <Tax variable> as header, <Value> under that header
        %  The tax variable names are defined below.
        filename                = strcat('TPCCIT_', taxplan, '.csv');
        tax_vars                = read_tax_vars( filename );
        s.captaxshare = tax_vars.CapitalTaxShare; 
        s.expshare    = tax_vars.ExpensingShare; 
        s.taucap      = tax_vars.CapitalTaxRate; 
        s.taucapgain  = 0           ;
        
        % Warn if parameters are outside expectations
        if( (s.captaxshare < 0) || (s.captaxshare > 1) )
            fprintf( 'WARNING! captaxshare=%f outside expecations.\n', captaxshare );
        end
        if( (s.expshare < 0) || (s.expshare > 1) )
            fprintf( 'WARNING! expshare=%f outside expectations.\n', expshare );
        end
        if( (s.taucap < 0) || (s.taucap > 1) )
            fprintf( 'WARNING! taucap=%f outside expectations.\n', taucap );
        end        
        if( (s.taucapgain < 0) || (s.taucapgain > 1) )
            fprintf( 'WARNING! taucapgain=%f outside expectations.\n', taucapgain );
        end  
    end  % tax()

    
    %% DEMOGRAPHICS
    %     Includes:
    %          survival probabilities
    %        , immigrant age distribution
    %        , birth rate
    %        , legal immigration rate
    %        , illegal immigration rate
    function s = demographics()
        
        survival  = read_series('XXXSurvivalProbability.csv', [], dirFinder.param);
        imm_age   = read_series('XXXImmigrantAgeDistribution.csv', [], dirFinder.param);
        s.surv    = survival';
        s.imm_age = imm_age';

        s.birth_rate   = 0.0200;    % Annual birth rate
        s.legal_rate   = 0.0016;    % Annual legal immigration rate
        s.illegal_rate = 0.0024;    % Annual illegal immigration rate

    end % demographics
    
    
    %% PRODUCTION
    %     Includes:
    %          TFP
    %        , depreciation
    %        , capital share
    function s = production()
        
        s.A     = 1;      % Total factor productivity
        s.alpha = 0.45;   % Capital share of output
        s.d     = 0.085;  % Depreciation rate

    end % production
        
    
    %% SOCIAL SECURITY
    %
    function s = social_security( modelunit_dollars, bv, T_model )
        
        ssthresholds = [856, 5157]*12*modelunit_dollars;    % Thresholds for earnings brackets
        ssrates      = [0.9, 0.32, 0.15];                   % Marginal benefit rates for earnings brackets
        ss_scale     = 1.6;                                 % Benefit scaling factor used to match total outlays as a percentage of GDP
        
        ssbenefit = [ max(min(bv, ssthresholds(1)) - 0              , 0) , ...
                      max(min(bv, ssthresholds(2)) - ssthresholds(1), 0) , ...
                      max(min(bv, Inf            ) - ssthresholds(2), 0) ] * ssrates';
        
        s.ssbenefits  = repmat(ssbenefit                , [1,T_model]);  % Benefits
        s.sstaxs      = repmat(0.124                    , [1,T_model]);  % Tax rates
        s.ssincmaxs   = repmat(1.185e5*modelunit_dollars, [1,T_model]);  % Maximum taxable earnings
        
        s.sstaxcredit = 0.15;     % Benefit tax credit percentage

    end % social_security
    
    
    %% BUDGET AND INTEREST RATES
    %
    function s = budget( first_transition_year, T_model )
        %  CBO interest rates, expenditures, and debt
        % Input: CBOInterestRate.csv -- interest rate (as pct) 
        %           Format is (Year), (PctRate) w/ header row.
        %           NOTE: Initial year is unused.
        % Input: SIMGDP.csv -- nominal GDP from baseline SIM
        %           Format is (Year), (NominalGDP) w/ header row.
        % Input: SIMRevenues.csv -- nominal tax revenues from baseline SIM
        %           Format is (Year), (Revenues) w/ header row.
        % Input: SIMExpenditures.csv -- nominal non-interest expenditures from baseline SIM
        %           Format is (Year), (Expenditures) w/ header row.
        % Input: CBOPublicDebt.csv  -- dollar debt of the gvt
        %           Format is (Year), (Debt) w/ header row.
        %           NOTE: Only the year corresponding to first_year is
        %           used. The other debt amounts are calculated from SIM
        %           series.
        % Input: CBONonInterestSpending.csv -- NIS as pct GDP
        %           Format is (Year), (PctGDP) w/ header row.
        % Input: CBOSocialSecuritySpending.csv -- SS spending as pct GDP
        %           Format is (Year), (PctGDP) w/ header row.
        % Input: CBOMedicareSpending.csv -- Medicare spending as pct GDP
        %           Format is (Year), (PctGDP) w/ header row.
        % Output: 
        %       debttoout, fedgovtnis, cborates, GEXP_by_GDP
        first_year                  = first_transition_year - 1;    % first year from which to read series
        CBODebt                     = read_series( 'CBOPublicDebt.csv', first_year, dirFinder.param );
        CBORates                    = read_series( 'CBOInterestRate.csv', first_year, dirFinder.param );
        SIMGDP                      = read_series( 'SIMGDP.csv', first_year, dirFinder.param );
        SIMRevenues                 = read_series( 'SIMRevenues.csv', first_year, dirFinder.param );
        SIMExpenditures             = read_series( 'SIMExpenditures.csv', first_year, dirFinder.param );
        CBONonInterestSpending      = read_series( 'CBONonInterestSpending.csv', first_year, dirFinder.param );
        CBOSocialSecuritySpending   = read_series( 'CBOSocialSecuritySpending.csv', first_year, dirFinder.param );
        CBOMedicareSpending         = read_series( 'CBOMedicareSpending.csv', first_year, dirFinder.param );

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
        
        growth_GDP      = zeros(size(SIMGDP), 'double');
        growth_GDP(1)   = 1.0;
        for i = 2:size(growth_GDP)
            growth_GDP(i) = SIMGDP(i)/SIMGDP(i-1);
        end;
        
        CBO_rates_growth_adjusted   = ((100.0+CBORates)./growth_GDP)./100.0 - 1.0;    
        deficit_nis_fraction_GDP    = deficit_nis./SIMGDP;                
        debt_percent_GDP            = debt./SIMGDP;                       
        
        % Name, transpose, truncate vars to correspond to currently used variables 
        % NOTE: The series must go out to T_model or import will break.
        %       Breaking is good here. (Though gracefully would be better)
        s.GEXP_by_GDP = GEXP_by_GDP;                                % Expenditures as pct gdp
        s.debt        = debt;                                       % Debt as pct gdp
        s.debttoout   = debt_percent_GDP(2);                        % Debt-to-output ratio of initial year
        s.fedgovtnis  = deficit_nis_fraction_GDP(2:T_model + 1)';   % starts from first transition path year
        s.cborates    = CBO_rates_growth_adjusted(2:T_model + 1)';  % starts from first transition path year
        s.cbomeanrate = nanmean(CBO_rates_growth_adjusted(2:end));  % Mean growth-adjusted interest rate over modeling period
        
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

    
    end % budget
    
end % methods

end % class paramGenerator


%%
%  Helper function to read CSV files with format: (Income), (Rate)
function [incomes, taxrates] = read_tax_rates( filename )

    warning( 'off', 'MATLAB:table:ModifiedVarnames' );          % for 2016b
    warning( 'off', 'MATLAB:table:ModifiedAndSavedVarnames' );  % for 2017a

    % Check if file exists and generate if necessary
    filepath    = fullfile(dirFinder.param(), filename);
    if ~exist(filepath, 'file')
        err_msg = strcat('Cannot find file = ', strrep(filepath, '\', '\\'));
        throw(MException('read_tax_table:FILENAME', err_msg ));
    end;
        
    T           = readtable(filepath, 'Format', '%f%f');
    incomes     = table2array(T(:,1));
    taxrates    = table2array(T(:,2));
    
    % Enforce that first threshold must be zero.
    if( incomes(1) > 0 )
        err_msg = strcat('First income threshold must be 0 in file ', strrep(filepath, '\', '\\'));
        throw(MException('read_tax_table:INCOME_THRESHOLD', err_msg ));
    end 
        

end % read_tax_rates()


%%
%  Helper function to read CSV file with format:
%     Header has variable names, then one row of values
function [tax_vars] = read_tax_vars( filename )

    warning( 'off', 'MATLAB:table:ModifiedVarnames' );          % for 2016b
    warning( 'off', 'MATLAB:table:ModifiedAndSavedVarnames' );  % for 2017a
    
    % Check if file exists and generate if necessary
    filepath    = fullfile(dirFinder.param(), filename);
    if ~exist(filepath, 'file')
        err_msg = strcat('Cannot find file = ', strrep(filepath, '\', '\\'));
        throw(MException('read_tax_vars:FILENAME', err_msg ));
    end;
        
    T           = readtable(filepath, 'Format', '%f%f%f%f');
    tax_vars    = table2struct(T);

end % read_tax_vars()


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

    warning( 'off', 'MATLAB:table:ModifiedVarnames' );          % for 2016b
    warning( 'off', 'MATLAB:table:ModifiedAndSavedVarnames' );  % for 2017a
 
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

end % read_series

%%  END FILE
