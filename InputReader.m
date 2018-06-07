classdef InputReader

methods (Static)
    
    %%
    %  Helper function to find input scenario ID corresponding to dynamic model scenario
    function [id] = find_input_scenario_id(mapfile, dynamic_model_scenario)
        
        % Load mapping table from file
        map = readtable(mapfile, 'ReadRowNames', true);
        
        % Specify input parameters to use for matching
        %   Assumes a 1-to-1 relationship between input parameters and identically named dynamic model scenario properties
        matchparams = { ...
            'TaxCode', ...
            'TaxRate', 'TaxMax', 'DonutHole', 'COLA', 'PIA', 'NRA', 'CreditEarningsAboveTaxMax', 'FirstYear', 'GradualChange', 'ProgressiveBenefitCut' ...
        };
        
        % Filter input parameters for those present in mapping table
        matchparams = matchparams(cellfun(@(param) any(strcmp(param, map.Properties.VariableNames)), matchparams));
        
        % Identify matching input scenarios
        match = arrayfun( ...
            @(input_scenario) all(cellfun( ...
                @(param) isequal(dynamic_model_scenario.(param), input_scenario.(param)), ...
                matchparams ...
            )), ...
            table2struct(map) ...
        );
        
        % Check for singular match
        if( sum(match) < 1 )
            error( 'No input scenario found corresponding to dynamic model scenario. Mapfile %s', mapfile );
        else
            if( sum(match) > 1)
                fprintf( '\nDuplicate ids: ');
                for i = 1:length(match)
                    if( match(i) )
                        s = map.Properties.RowNames(i);
                        fprintf( ' %s ', char(s) );
                    end
                end
                fprintf('\n');
                error( 'More than one input scenario found corresponding to dynamic model scenario. Mapfile %s', mapfile);
            end
        end
        
        % Extract ID of matching input scenario
        id = map.Properties.RowNames{match};
        
    end
    
    
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
            err_msg = strcat('First bracket must be 0 in file ', strrep(filename, '\', '\\'));
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
    
    
end

end
