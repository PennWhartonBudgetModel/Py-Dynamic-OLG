%%
% Generate tax policy parameters according to predefined plans.
% 
function [] = generate_tax()

    % Identify parameter directory
    param_dir = dirFinder.param();

    for taxplans = {'base', 'ryan'}, taxplan = taxplans{1};

        % Get tax rates
        filename    = strcat('TPCPIT_', taxplan, '.csv');
        [income_thresholds, tax_rates] = read_tax_rates( filename, param_dir );

        % calculate tax liability at each threshold
        if( income_thresholds(1) > 0 )
            inc       = [0; income_thresholds];
        else
            inc       = income_thresholds;
        end;
        tax_        = cumsum(diff(inc).*tax_rates); 
        income_tax  = [0; tax_];

        % no deductions, re-calculate effective tax if existent
        %  TODO: If these thresholds are different, we need to make
        %        a single vector of the union of the two thresholds
        %        and calculate effective income tax at those points
        

        % Store income tax and deduction functions
        %   REM: Transpose into row vectors
        s.(taxplan).tax_thresholds      = inc';
        s.(taxplan).tax_income          = income_tax';
        s.(taxplan).tax_rates           = [tax_rates; tax_rates(end)]';

        % Get the capital tax params and store them
        filename                = strcat('TPCCIT_', taxplan, '.csv');
        tax_vars                = read_tax_vars( filename, param_dir );
        s.(taxplan).captaxshare = tax_vars.CapitalTaxShare; 
        s.(taxplan).expshare    = tax_vars.ExpensingShare; 
        s.(taxplan).taucap      = tax_vars.CapitalTaxRate; 
        s.(taxplan).taucapgain  = 0           ;
    end  % for

    % Save parameters
    save(fullfile(param_dir, 'tax.mat'), '-struct', 's');

    % show_tax_rates(s.('base').pit_coefs, s.('base').inc_tax, incv);

        %  OLD CODE: FIT TO GS FUNCTION
        % Define vector of sample incomes and find corresponding tax rates
        %    Rem: taxrates are marginal rates
        % topincome   = 1e6; % incthresholds(end);
        % incv        = linspace(1, topincome, topincome/1e3)'; % put point every $1000
        % inct        = arrayfun(@(inc) taxrates(find(inc <= incthresholds, 1)), incv);
        % inc_tax     = cumsum(inct.*1e3);

        % Fit income tax function using least squares
        %   Fit is on tax liability.
        %   GS is  ATR = b(1) - b(1)*( b(2)* y^b(3)+1)^(-1/b(3))
        %   we subtract b(4)/y  as deduction
        
%         gouveiastrauss  = @(pit_coefs) pit_coefs(1).* ...
%                 (1 - (pit_coefs(2).*(incv.^pit_coefs(3)) + 1).^-(1./pit_coefs(3)));
%                 % - pit_coefs(4);    % deduction addtion
%         starting_guess = [0.25, 0.03, 0.77];
%         pit_coefs       = lsqnonlin(@(pit_coefs) gouveiastrauss(pit_coefs).*incv - inc_tax, starting_guess, [], [], optimoptions(@lsqnonlin, 'Display', 'off'));
% END OLD CODE

end  % generate_tax()

%%
%  Helper function to read CSV files with format: (Income), (Rate)
function [incomes, taxrates] = read_tax_rates(filename, param_dir )

    warning( 'off', 'MATLAB:table:ModifiedVarnames' );

    % Check if file exists and generate if necessary
    filepath    = fullfile(param_dir, filename);
    if ~exist(filepath, 'file')
        err_msg = strcat('Cannot find file = ', strrep(filepath, '\', '\\'));
        throw(MException('read_tax_table:FILENAME', err_msg ));
    end;
        
    T           = readtable(filepath, 'Format', '%f%f');
    incomes     = table2array(T(:,1));
    taxrates    = table2array(T(:,2));

end % read_tax_rates()


%%
%  Helper function to read CSV file with format:
%     Header has variable names, then one row of values
function [tax_vars] = read_tax_vars(filename, param_dir )

    warning( 'off', 'MATLAB:table:ModifiedVarnames' );

    % Check if file exists and generate if necessary
    filepath    = fullfile(param_dir, filename);
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
