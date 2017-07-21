%%
% Reader and generator of parameters for dynamic model
% 
%%


classdef paramGenerator

properties (Constant)
end % properties

methods (Static)

    %%
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
        s.tax_income          = income_tax';
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
    
    end  % tax()

    
    %% Demographic params
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

