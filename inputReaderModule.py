import numpy as np
import pandas as pd
import numbers
import os
import warnings

class InputReader:

    #  Helper function to read CSV files with format: 
    #    (Year), (Bracket1), ... (BracketN), (Rate1), ... (RateN)
    #    filename     : fullfile of CSV to read
    #    firstYear    : don't read years before this param
    #    Tyears       : read this many years 
    @staticmethod
    def read_brackets_rates_indices(filename, firstYear, Tyears):
        
        if not os.path.exists(filename):
            raise Exception('Cannot find file %s' % filename)
        
        T = pd.read_csv(filename)
        
        try:
            (brackets, rates, indices) = InputReader.parse_brackets_rates_indices(T, firstYear, Tyears)
        except:
            raise Exception('Problem in file %s ' % filename)
        return (brackets, rates, indices)


    #  Helper function to read CSV files with format: 
    #    (BirthYear), (Year), (Bracket1), ... (BracketN), (Rate1), ... (RateN)
    #    filename     : fullfile of CSV to read
    #    firstYear    : don't read years before this param
    #    Tyears       : read this many years 
    #    firstCohort  : birthyear of oldest cohort
    #    lastCohort   : birthyear of youngest cohort
    #    
    #    Output: a struct with fields named 'b<Year>' where <Year> is the
    #            birthYear of the cohort. Each field is a struct with 
    #            fields 'brackets', 'rates', 'indices' and each of these is
    #            a time-indexed matrix.
    @staticmethod
    def read_cohort_brackets_rates_indices(filename, firstYear, Tyears, firstCohort, lastCohort):
        
        if not os.path.exists(filename):
            raise Exception('Cannot find file %s' % filename)
        
        T = pd.read_csv(filename)
        cohorts = {}
        
        # Check that all needed cohorts are present, 
        #   if missing at end, we pad by copying last cohort
        lastAvailableCohort = max(T['BirthYear'])
        
        # Extract sub-table for each cohort and parse
        for birthYear in range(firstCohort, lastAvailableCohort+1):
            T_cohort = T[T.BirthYear == birthYear]
            try:
                (brackets, rates, indices) = InputReader.parse_brackets_rates_indices(T_cohort, firstYear, Tyears)
            except:
                raise Exception('Problem in file %s with cohort %d' % (filename, birthYear) )
            cohortField = 'b%u' % birthYear
            cohorts[cohortField] = {}
            cohorts[cohortField]['brackets']  = brackets
            cohorts[cohortField]['rates']     = rates
            cohorts[cohortField]['indices']   = indices
        
        # Pad if needed by copying last cohort
        endCohort   = cohorts['b%u' % lastAvailableCohort]
        numBrackets = endCohort['brackets'].shape[1]
        
        # Helper function to shift and fill with NaN's
        @staticmethod
        def shiftAndFill(x, offset):
            n = Tyears - offset
            if not isinstance(x, numbers.Number):
                # TBD: Until we switch to time-varying indices,
                #   the first line needs to be a real index and not <blank>
                f = np.tile(x[1,:], [offset, 1])
            else:
                f = np.tile(None, [offset, numBrackets])
            y = np.vstack(f, x[1:n, :])
            return y
        
        for birthYear in range(lastAvailableCohort,lastCohort):
            offset      = birthYear - lastAvailableCohort
            cohortField = 'b%u' % birthYear
            cohorts[cohortField]['brackets']  = shiftAndFill(endCohort['brackets'], offset)
            cohorts[cohortField]['rates']     = shiftAndFill(endCohort['rates'], offset)
            cohorts[cohortField]['indices']   = shiftAndFill(endCohort['indices'], offset)
        
        return cohorts
    
    ##
    #  Helper function to read CSV files with format: 
    #    (Year), (Age), (LegalStatus), (Emigrated), (Immigrated),
    #    (Born), (Died), (Population)
    #    filename     : fullfile of CSV to read
    #    firstYear    : don't read years before this param
    #    T_model      : read this many years 
    #    realage_entry: real age of entry in the economy
    #    T_life       : maximum length of years alive
    #    g            : population group index mapping
    #    
    #    Output: a struct with fields named 'g'. Each 'g' is a struct with 
    #            fields "Emigrated Immigrated Born Died Population", each 
    #            which is a (T_model x T_life) array where Year = firstYear
    #            is the first year to be considered.
    @staticmethod
    def read_demographics(filename, firstYear, T_model, realage_entry, T_life, g, microSimStatus, statusMap, scenario):
        
        if not os.path.exists(filename):
            raise Exception('Cannot find file %s' % filename)
        
        T = pd.read_csv(filename)
                
        groups    = {}
        groupList = list(g.keys())
        microsimGroups = {}
        microsimGroupList = list(microSimStatus.keys())
        
        varlist   = ['EmigratedThisYear', 'ImmigratedThisYear', 'BornThisYear', 'DiedThisYear', 'Population']

        for group in range(len(microsimGroupList)):
            # Extract sub-table for each group
            
            T_group = T.loc[T.LegalStatus == group]
            
            microsimGroups[microsimGroupList[group]] = {
            'EmigratedThisYear': np.zeros((T_model, T_life)),
            'ImmigratedThisYear': np.zeros((T_model, T_life)),
            'BornThisYear': np.zeros((T_model, T_life)),
            'DiedThisYear': np.zeros((T_model, T_life)),
            'Population': np.zeros((T_model, T_life))
            }
            
            for age in range(realage_entry + 1, realage_entry + T_life + 1):

                # Find age
                if sum(T_group.Age.values == age) == 0:
                    raise Exception('Cannot find age %u for group %s in input table.' % (age, group))
                    
                # Extract sub-sub-table for each group & age
                T_group_age = T_group.loc[T_group.Age == age]

                # Find first year & check for missing years
                if not firstYear in T_group_age.Year.values:
                    raise Exception('Cannot find first index (%u) for age %u of group %s in input table.' % (firstYear, age, microsimGroupList[group]))
                years = (T_group_age.Year >= firstYear) 
                gap = np.diff(T_group_age.loc[years, 'Year'].values)
                if any(gap>1):
                    missing_years = '%d ' % (T_group_age.Year[gap > 1] + 1)
                    raise Exception('Missing year(s) %sfor age %u of group %s in input table.' % (missing_years, age, microsimGroupList[group]))
                
                # Truncating the data if microsim is longer than OLG.
                # Throwing an error if OLG is longer than microsim.
                
                for var in varlist:
                    
                    temp_array = np.array(T_group_age.loc[years,var].values)
                
                    # Pad series if not long enough, truncate if too long
                    numYears = len(temp_array)
                    if T_model - numYears <= 0:
                        temp_array = temp_array[0:T_model]
                    else:
                        raise Exception('read_demographics:PADDING Microsim does not extend far enough into the future.')

                    microsimGroups[microsimGroupList[group]][var][:,age-realage_entry-1] = temp_array
                    
        # Initializing the structure for the demographics.     
        for group in range(len(groupList)):
            groups[groupList[group]] = {}
            for var in varlist:
                groups[groupList[group]][var] = np.zeros(microsimGroups[microsimGroupList[0]][var].shape)
        
        # Taking each of the microsim modules and adding them the OLG
        # groups using the mapping.
        assert len(statusMap[:,0]) == len(np.unique(statusMap[:,0])), 'Microsim demographics do not uniquely map to OLG demographics.'
        assert len(statusMap[:,0]) <= len(microsimGroups.keys()),'Mapping of demographics has indices that do not exist in microsim inputs.'
        if len(statusMap[:,0]) < len(microsimGroups.keys()):
            print ('Immigration status exists microsim model is not used by OLG model.')
        for groupIndex in range(len(statusMap[:,0])):
            for var in varlist:
                groups[statusMap[groupIndex,1]][var] = groups[statusMap[groupIndex,1]][var] + microsimGroups[statusMap[groupIndex,0]][var]
        
        
        # This is the error checking.  It makes sure that the law of motion
        # holds for all groups, over all ages, and for all years.  If not,
        # it sends a warning. This only reports when we're actually using
        # the new demographics.
        if scenario.UseNewDemographics:
            for group in range(len(g.keys())):
                for iterYear in range(1, T_model):
                    for age in range(1,T_life):
                        if (groups[groupList[group]]['Population'][iterYear,age] > 10000):
                            a = 10
                        
                        popError = (groups[groupList[group]]['Population'][iterYear,age]         + groups[groupList[group]]['EmigratedThisYear'][iterYear-1,age-1] + 
                                   groups[groupList[group]]['DiedThisYear'][iterYear-1,age-1]   - groups[groupList[group]]['Population'][iterYear-1,age-1] -
                                   groups[groupList[group]]['ImmigratedThisYear'][iterYear,age] - groups[groupList[group]]['BornThisYear'][iterYear,age])
                        if abs(popError) > 0:
                                print('Inputs are not consistent for age %s in year %s in group %s. Error is: %s.' % (str(age), str(iterYear), str(group), str(popError)))
            
        return groups
       
    ##
    # Read a CSV file in format
    #    header: IndexName, VarName1, VarName2, ... VarNameN
    #    data  : (Index)  , (Value1), (Value2), ... (ValueN)
    #       index_name      : string name of index variable -- i.e. IndexName
    #       first_index     : first index to read (previous part of file is ignored
    #       last_index      : (optional) If given, copy the last available
    #                       value so that series goes from first_index ...
    #                       last_index. If the original series is too long,
    #                       truncate.
    #    For time series, (Index) is (Year), 
    @staticmethod
    def read_series(filename, index_name, first_index, last_index=None):
        
        # Check if file exists 
        if not os.path.exists(filename):
            raise Exception('Cannot find file %s' % filename)
        
        T = pd.read_csv(filename)
        
        # Remove unused table rows
        drop = T[index_name] < first_index
        if sum(drop) == len(T):
            raise Exception('Cannot find first index in file.')
        
        T = T.drop(T.loc[drop,].index)
        
        # Pad if needed 
        num_add  = T[index_name].iloc[0] - first_index
        if num_add > 0:
            print('WARNING! File %s begins %u periods after %u.\n' % (filename, num_add, first_index))
            columns = T.columns.values
            N = T.to_numpy()
            N = np.vstack((np.tile(N[1,:], [num_add, 1]), N))
            T = pd.DataFrame(N, columns = columns)
            
        if last_index != None:
            # Truncate if needed
            drop = T[index_name] > last_index
            T = T.drop(T.loc[drop,].index)
            
            # Pad if needed 
            num_add  = last_index - T[index_name].iloc[-1]
            if num_add > 0:
                print('WARNING! Padding file %s since it ends %u periods before last index %u.\n' % (filename, num_add, last_index))
                columns = T.columns.values
                N = T.to_numpy()
                N = np.vstack((N, np.tile(N[-1,:], [num_add, 1])))
                T = pd.DataFrame(N, columns = columns)
        
        series = T.to_dict('list')
        
        #Turn list elements to numpy arrays
        for k in series.keys():
            series[k] = np.array(series[k])
        
        # Do not return index column
        series.pop(index_name, None)
        
        return series
    
    ##
    #  Helper function to read CSV files with format: 
    #    (Year), (Age), (Value1), ... (ValueN)
    #    filename     : fullfile of CSV to read
    #    firstYear    : don't read years before this param
    #    Tyears       : read this many years 
    def read_shocks(filename, firstYear, Tyears, T_life):
        
        if not os.path.exists(filename):
            raise Exception('Cannot find file %s' % filename)
        
        T = pd.read_csv(filename)
        
        # Find first year
        if not firstYear in T.Year.values:
            raise Exception('Cannot find first index (%u) in input table.' % firstYear)
        
        # Find ranges of ages and states
        firstAge    = min(T.Age.values)
        lastAge     = max(T.Age.values)
        T_year_columns = [True if 'Value' in x else False for x in T.columns.values]
        numStates   = np.sum(T_year_columns)        
        if lastAge - firstAge + 1 != T_life :
            raise Exception('Input table has different life span.')
        
        # Create the right size arrays
        values      = np.zeros((Tyears, T_life, numStates))
        
        for year in range(firstYear,firstYear+Tyears):
            
            y = year - firstYear            
            T_year_rows = (T.Year.values == year)

            # Check if need to pad, and pad if so
            if( np.sum(T_year_rows) == 0 ):
                # At least firstYear must have been read (or previous
                # error), so take year-1
                values[y, :, :]         = values[y-1, :, :]
                continue
            
            values[y, :, :] = T.values[T_year_rows,:][:,T_year_columns]
        
        return values
            
    
    ##
    #  Helper function to read CSV files with format: 
    #    (Year), (Age), (Value1), ... (ValueN), (Transition_1_1), ... (Transition_N_N)
    #    filename     : fullfile of CSV to read
    #    firstYear    : don't read years before this param
    #    Tyears       : read this many years 
    @staticmethod
    def read_transitions(filename, firstYear, Tyears, T_life):
        
        if not os.path.exists(filename):
            raise Exception('read_transitions: Cannot find file %s' % filename)
        
        T = pd.read_csv(filename)
        
        # Find first year
        if not firstYear in T.Year.values:
            raise Exception('read_transitions: Cannot find first index (%u) in input table.' % firstYear)
        
        # Find ranges of ages and states
        firstAge    = min(T.Age.values)
        lastAge     = max(T.Age.values)
        numStates   = np.sum([True if 'Transition_1_' in x else False for x in T.columns.values])
        
        if( lastAge - firstAge + 1 != T_life ):
            raise Exception('read_transitions: Input table has different life span.')
        
        # Create the right size arrays
        transitions = np.zeros((Tyears, T_life, numStates, numStates))
        
        for year in range(firstYear,firstYear+Tyears):
            
            y = year - firstYear           
            T_year_rows = (T.Year.values == year)

            # Check if need to pad, and pad if so
            if( np.sum(T_year_rows) == 0 ):
                # At least firstYear must have been read (or previous
                # error), so take year-1
                transitions[y, :, :, :] = transitions[y-1, :, :, :]
                continue
            
            for initState in range(1,numStates+1):
                # Read all transitions out of initState
                T_columns = [True if ('Transition_%u_' % initState) in x else False for x in T.columns.values]
                trans = T.values[T_year_rows,:][:, T_columns]
                
                # Check that transitions out sum to 100% probability,
                s = np.sum(trans, axis = 1) - 1
                if any(s > 1e-15) or any(s < -1e-15) :
                    raise Exception('read_transitions: Transition probabilities out of state do not sum to 1.')
                else:
                    # Fix-up, by adjusting first valid column if within numerical deviation
                    if any(s) :
                        for i in range(trans.shape[0]):
                            for j in range(trans.shape[1]):
                                if( (trans[i,j] > 1e-15) and (trans[i,j] < 1 - 1e-15) ):
                                    trans[i,j] = trans[i,j] - s[i]
                                    break
                
                transitions[ y, :, initState-1, :] = trans

        return transitions
        
    ##
    #  Helper function to read CSV files with format: 
    #    (Year), (Age), (Mass_1_1), ... (Mass_ng_nz)
    #    filename     : fullfile of CSV to read
    #    initialYear  : don't read years before this param
    #    Tlife        : life spam
    @staticmethod
    def read_initDIST(filename, initialYear, T_life):
        
        if not os.path.exists(filename):
            raise Exception('read_initDIST: Cannot find file %s' % filename)
        
        T = pd.read_csv(filename)       
        T_year_rows = (T.Year.values == initialYear)
        
        if sum(T_year_rows) == 0:
            raise Exception('Cannot find index (%u) in input table.' % initialYear)
        
        # Select sub-table with the appropriate year
        T = T.iloc[T_year_rows,:]
        
        # Find ranges of ages and states
        firstAge    = min(T.Age.values)
        lastAge     = max(T.Age.values)
        numStates   = sum(['Mass_' in s for s in T.columns.values])
        numGroups   = sum(T.Age.values == firstAge)
        
        if lastAge - firstAge + 1 != T_life:
            raise Exception('read_initDIST: Input table has different life span.')
        
        # Create the distribution array
        initDIST = np.zeros((numStates, T_life, numGroups))
        
        for groupState in range(numGroups):
            rows_group = (T.Group.values == groupState + 1)
            dist_group = T.iloc[rows_group,:]
            dist_group = dist_group.drop(dist_group.columns.values[0:3], axis = 1)
            dist_group = dist_group.values

            # Check that distributions sum to 100% of mass,
            s = np.sum(dist_group, axis = 1)
            if any(s - 1 > 5e-15) or any(s - 1 < -5e-15):
                raise Exception('Mass of households across states do not sum to 1.')
            else:
                # Fix-up by re-scaling
                if any(s):
                    dist_group = dist_group / np.tile(s[:,np.newaxis], [1, numStates])

            initDIST[:, :, groupState] = np.transpose(dist_group)
        
        return initDIST

    ##
    #  Helper function to read Matlab table with format: 
    #    (Year), (Bracket1), ... (BracketN), (Rate1), ... (RateN)
    #    tableData    : data in Matlab table format
    #    firstYear    : don't read years before this param
    #    Tyears       : read this many years 
    @staticmethod
    def parse_brackets_rates_indices(tableData, firstYear, Tyears):
        
        # Find first year
        if not firstYear in tableData.Year.values:
            raise Exception('parse_brackets_rates: Cannot find first index (%u) in input table.' % firstYear)   
        yearStart = tableData[tableData.Year == firstYear].index[0]
        
        # Find all brackets, rates, and indices
        brackets = tableData[tableData.columns[['Bracket' in s for s in tableData.columns.values]]]
        brackets = brackets.loc[yearStart:,].to_numpy()
        
        rates    = tableData[tableData.columns[['Rate' in s for s in tableData.columns.values]]]
        rates = rates.loc[yearStart:,].to_numpy()
        
        indices    = tableData[tableData.columns[['Index' in s for s in tableData.columns.values]]]
        indices = indices.loc[yearStart:,].to_numpy()
        
        # Enforce that first bracket must be zero.
        # If there are no brackets, make all zeros brackets
        if len(brackets) == 0:
            brackets = np.zeros(rates.shape)
        
        #check = [1 if v > 0 else 0 for v in brackets[:,0]]
        #if 0 in check:
        if all(brackets[:,0]>0):
            raise Exception('First bracket must be 0 in input table.')
        
        # Pad inputs if not long enough, truncate if too long
        numYears = brackets.shape[0]
        if Tyears - numYears <= 0:
            brackets    = brackets  [0:Tyears,:]
            rates       = rates     [0:Tyears,:]
            indices     = indices   [0:Tyears,:]
        else:
            brackets    = np.vstack(brackets, np.tile(brackets[-1,:], [Tyears-numYears, 1]))
            rates       = np.vstack(rates, np.tile(rates[-1,:], [Tyears-numYears, 1]))
            indices     = np.vstack(indices, np.tile(indices[-1,:], [Tyears-numYears, 1]))
        
        return (brackets, rates, indices)

'''
    ##
    #  Helper function to read Matlab table with format: 
    #    (Year), (Age), (Value1), ... (ValueN),(Transition_1_1), ... (Transition_N_N)  
    #   where we assume an NxN transition matrix
    #    tableData    : data in Matlab table format
    #    firstYear    : don't read years before this param
    #    Tyears       : read this many years 
    @staticmethod
    def parse_transitions(tableData, firstYear, Tyears, Tlife):
        
        # Find first year
        if not firstYear in tableData.Year:
            raise Exception('Cannot find first index (%u) in input table.' % firstYear)  
        
        # Find ranges of ages and states
        firstAge    = min(tableData.Age)
        lastAge     = max(tableData.Age)
        numStates   = sum('Value' in s for s in tableData.columns.values)
        
        if lastAge - firstAge + 1 != Tlife:
            raise Exception('Input table has different life span.')
        
        # Create the right size arrays
        values      = np.zeros(Tyears, Tlife, numStates)
        transitions = np.zeros(Tyears, Tlife, numStates, numStates)
    
        for year in range(firstYear-1,firstYear+Tyears-1):
            
            y = year - firstYear + 1
            
            T_year_rows = tableData.Year == year

            # Check if need to pad, and pad if so
            if sum(T_year_rows) == 0:
                # At least firstYear must have been read (or previous
                # error), so take year-1
                values[y, :, :]         = values[y-1, :, :]
                transitions[y, :, :, :] = transitions[y-1, :, :, :]
            
            selection = tableData[tableData.columns[['Value' in s for s in tableData.columns.values]]]
            selection = selection.loc[T_year_rows,]
            values[y, :, :] = selection.to_numpy()
            for initState in range(numStates):
                # Read all transitions out of of initState
                trans = tableData[tableData.columns[[('Transition_%u_' % initState) in s for s in tableData.columns.values]]]
                trans = tableData.loc[T_year_rows,]
                trans = trans.to_numpy()
                
                # Check that transitions out sum to 100% probability,
                s = sum(trans, 2) - 1
                check = [1 if v > 1e-15 or v < -1e-15 else 0 for v in s]
                if 1 in check:
                    raise Exception('Transition probabilities out of state do not sum to 1.')
                else:
                    # Fix-up, by adjusting first valid column if within numerical deviation
                    check2 = [1 if v != 0 else 0 for v in s]
                    if 1 in check2:
                        for i in range(trans.shape[0]):
                            for j in range(trans.shape[1]):
                                if (trans[i,j] > 1e-15) and (trans[i,j] < 1 - 1e-15):
                                    trans[i,j] = trans[i,j] - s[i]
                                    break
                
                transitions[y, :, initState, :] = trans
        
        return (values, transitions)
'''