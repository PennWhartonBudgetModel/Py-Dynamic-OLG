%%
% Generate tax policy parameters according to predefined plans.
% 
function [] = generate_tax()

    % To have access to dirFinder
    addpath( '..' );
    
    for taxplans = {'base', 'ryan'}, taxplan = taxplans{1};

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
        s.(taxplan).tax_thresholds      = income_thresholds';
        s.(taxplan).tax_income          = income_tax';
        s.(taxplan).tax_rates           = tax_rates';

        % Get the capital tax params and store them.
        %  Input files TPCCIT_<taxplan>.CSV expected to have the
        %  following structure:
        %      <Tax variable> as header, <Value> under that header
        %  The tax variable names are defined below.
        filename                = strcat('TPCCIT_', taxplan, '.csv');
        tax_vars                = read_tax_vars( filename );
        s.(taxplan).captaxshare = tax_vars.CapitalTaxShare; 
        s.(taxplan).expshare    = tax_vars.ExpensingShare; 
        s.(taxplan).taucap      = tax_vars.CapitalTaxRate; 
        s.(taxplan).taucapgain  = 0           ;
    end  % for

    % Save parameters
    save(fullfile(dirFinder.param(), 'tax.mat'), '-struct', 's');

    % DEBUG: show_tax_rates(s.('base').pit_coefs, s.('base').inc_tax, incv);

end  % generate_tax()

%%
%  Helper function to read CSV files with format: (Income), (Rate)
function [incomes, taxrates] = read_tax_rates( filename )

    warning( 'off', 'MATLAB:table:ModifiedVarnames' );

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

    warning( 'off', 'MATLAB:table:ModifiedVarnames' );

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
%   Utility to plot the tax liability
function [] = show_tax_rates( pit_coefs, actual_tax, incomes )

    incv = incomes;
    gouveiastrauss  = @(pit_coefs) pit_coefs(1).* ...
                (1 - (pit_coefs(2).*(incv.^pit_coefs(3)) + 1).^-(1./pit_coefs(3)));
    taxes           = gouveiastrauss(pit_coefs).*incv;

    y(:,1)          = taxes./incv;
    y(:,2)          = max(actual_tax,0)./incv;
    y(:,3)          = actual_tax;
    
    plot (incv(1:50), y(1:50,3));

end % show_tax_rates

%%
