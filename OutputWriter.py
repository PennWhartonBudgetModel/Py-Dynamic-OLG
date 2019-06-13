import os
import numpy as np
import pandas as pd
import scipy.io as sio
import PathFinder
import InputReader
import Scenario

class OutputWriter:
            
    # Define mapping from series names to dynamic model variable names and series sources
    #      for "sources", if source_name is [], use the series_name (for convenience)
    series_names = {
       'GDP'                                   : {'var_name', 'outs'                               , 'source', 'projections'    , 'source_name', []                    },   
       'GDP_FY'                                : {'var_name', 'outs'                               , 'source', 'projections'    , 'source_name', []                    },   
       'GNP'                                   : {'var_name', 'GNI'                                , 'source', 'projections'    , 'source_name', []                    },   
       'CapitalServices'                       : {'var_name', 'caps'                               , 'source', 'projections'    , 'source_name', []                    },   
       'LaborInput'                            : {'var_name', 'labeffs'                            , 'source', 'projections'    , 'source_name', []                    },   
       'WagesAndSalaries'                      : {'var_name', 'labincs'                            , 'source', 'projections'    , 'source_name', []                    },   
       'CompensationOfEmployees'               : {'var_name', 'labincs'                            , 'source', 'projections'    , 'source_name', []                    },   
       'Employment'                            : {'var_name', 'labs'                               , 'source', 'projections'    , 'source_name', []                    },   
       'AverageInterestRateOnDebt'             : {'var_name', '_add_inflation_to_interestrate'     , 'source', 'Market'         , 'source_name', 'bondDividendRates'   },   
     # BUDGET
       'OutlaysDiscretionary'                  : {'var_name', '_nonindexed'                        , 'source', 'projections'    , 'source_name', []                    },   
       'OutlaysMedicare'                       : {'var_name', '_nonindexed'                        , 'source', 'projections'    , 'source_name', []                    },   
       'OutlaysMedicaid'                       : {'var_name', '_nonindexed'                        , 'source', 'projections'    , 'source_name', []                    },   
       'OutlaysFederalRetirement'              : {'var_name', '_nonindexed'                        , 'source', 'projections'    , 'source_name', []                    },   
       'OutlaysVeteransPrograms'               : {'var_name', '_nonindexed'                        , 'source', 'projections'    , 'source_name', []                    },   
       'OutlaysOtherPrograms'                  : {'var_name', '_nonindexed'                        , 'source', 'projections'    , 'source_name', []                    },   
       'OutlaysOffsettingReceipts'             : {'var_name', '_nonindexed'                        , 'source', 'projections'    , 'source_name', []                    },   
       'OutlaysIncomeSecurity'                 : {'var_name', '_nonindexed'                        , 'source', 'taxcalculator'  , 'source_name', []                    },   
       'OutlaysSocialSecurity'                 : {'var_name', 'bens'                               , 'source', 'oasicalculator' , 'source_name', []                    },   
       'RevenuesPayrollTaxExSocialSecurity'    : {'var_name', 'labincs'                            , 'source', 'taxcalculator'  , 'source_name', []                    },   
       'RevenuesPayrollTaxSocialSecurity'      : {'var_name', 'ssts'                               , 'source', 'oasicalculator' , 'source_name', []                    },   
       'RevenuesIndividualIncomeTax'           : {'var_name', 'pits'                               , 'source', 'taxcalculator'  , 'source_name', []                    },   
       'RevenuesCorporateIncomeTax'            : {'var_name', 'corpTaxs'                           , 'source', 'taxcalculator'  , 'source_name', []                    },   
       'RevenuesEstateAndGiftTaxes'            : {'var_name', 'caps'                               , 'source', 'taxcalculator'  , 'source_name', []                    },   
       'RevenuesExciseTaxes'                   : {'var_name', 'outs'                               , 'source', 'taxcalculator'  , 'source_name', []                    },   
       'RevenuesCustomsDuties'                 : {'var_name', 'GNI'                                , 'source', 'projections'    , 'source_name', []                    },   
       'RevenuesMiscellaneousReceipts'         : {'var_name', '_nonindexed'                        , 'source', 'taxcalculator'  , 'source_name', []                    },   
       'SS_Cost'                               : {'var_name', 'bens'                               , 'source', 'oasicalculator' , 'source_name', []                    },   
       'SS_NonInterestIncome'                  : {'var_name', 'ssts'                               , 'source', 'oasicalculator' , 'source_name', []                    },   
       'SS_TaxablePayroll'                     : {'var_name', 'laborIncomeSubjectToSocialSecuritys', 'source', 'oasicalculator' , 'source_name', []                    },   
       'SS_TaxesOnBenefits'                    : {'var_name', 'averageEffectivePITRatesProductSocialSecurityBenefits', 'source', 'oasicalculator' , 'source_name', []  },   
       'SS_PayrollTaxes'                       : {'var_name', 'ssts'                               , 'source', 'oasicalculator' , 'source_name', []                    }   
    }

    dynamic_outvars = {
        'ssts'              : 'PayrollTax'              , 
        'constax'           : 'ConsumptionTax'          , 
        'corpTaxs'          : 'CorpTax'                 , 
        'pits'              : 'HH_IncomeTax'            , 
        'caps_domestic'     : 'Capital_Domestic'        , 
        'caps_foreign'      : 'Capital_Foreign'         , 
        'debts_domestic'    : 'Debt_Domestic'           , 
        'debts_foreign'     : 'Debt_Foreign'            , 
        'debts'             : 'GovernmentDebt'          , 
        'outs'              : 'Output'                  , 
        'GNI'               : 'GrossNationalIncome'     , 
        'bens'              : 'SocialSecurityBenefits'  , 
        'caps'              : 'Capital'                 , 
        'labeffs'           : 'EfficientLabor'          , 
        'labs'              : 'Labor'                   , 
        'labincs'           : 'LaborIncome'             , 
        'capincs'           : 'CapitalIncome'           , 
        'investment'        : 'Investment'              , 
        'corpDebts'         : 'CorpDebt'                , 
        'passDebts'         : 'PassThroughDebt'         , 
        'cons'              : 'Consumption'             , 
        'pops'              : 'PopulationHouseholds'      
    } 
    
    market_outvars = { 
        'MPKs'               : 'MPK'                     , 
        'capsharesPM'        : 'CapitalShares'           , 
        'equityPrices'       : 'EquityPrices'            , 
        'equityDividendRates': 'EquityDividendRates'     , 
        'bondDividendRates'  : 'BondDividendRates'       , 
        'capgains'           : 'Capital_Gains'           , 
        'wages'              : 'MPL'                       
    } 
    

    #  Find scenarios not existant in 'series' interface 
    #  and write them.
    #     Inputs:   scenarios:  array of Scenario objects
    @staticmethod
    def writeScenarios(scenarios):
        
        numWritten      = 0
        numScenarios    = len(scenarios)
        if numScenarios == 0:
            return numWritten
        
        # Load mapping table from file
        mapfile = os.path.join(PathFinder.getSeriesOutputDir(), 'map.csv')
        if not os.path.exists(mapfile):
            map = pd.DataFrame()
        else:
            map = pd.read_csv(mapfile)
        
        # Identify any scenarios not in current map file and add them.
        paramNames      = scenarios[0].getParams().keys()
        numMapScenarios = len(map.index)
        firstScenario   = 0
        
        # In case of new map file, add first viable scenario and then
        #  continue with existant map
        if numMapScenarios == 0:
           for firstScenario in range(numScenarios):
               scenario = scenarios[firstScenario]
               if scenario.isSolved():
                   id = 1
                   success          = OutputWriter.exportSeries(scenario, id)
                   numWritten       = numWritten + success
                   firstScenario    = firstScenario + 1
                   
                   # Create table and get out
                   map              = pd.DataFrame(scenario.getParams())
                   map.ID           = id
                   numMapScenarios  = 1
                   break
               else:
                   print('WARNING! Worklist scenario %u is not solved. Cannot export.\n' % firstScenario)

           if numMapScenarios == 0:
               return numWritten
        
        # Convert map to struct array, but without the ID column
        #   We assume that the IDs are going to go from 1:N,
        #   so we can renumber afterwards and the IDs for existing
        #   mapfile entries should remain the same.
        mapScenarios    = map.to_dict(orient = 'list')
        mapScenarios.pop('ID')
        
        # Iteratate over worklist
        for i in range(firstScenario,numScenarios + 1):
            scenario        = scenarios[i]
            scenarioParams  = scenario.getParams()
            
            # Search the map for any matching scenarios
            f = lambda mapEntry: scenarioParams == mapEntry
            match = np.array([f(x) for x in mapScenarios])

            # Check for singular match
            numMatch = sum(match)
            if numMatch == 0:
                if scenario.isSolved():
                    numMapScenarios = numMapScenarios + 1 
                    id              = numMapScenarios
                        
                    success    = OutputWriter.exportSeries(scenario, id)
                    numWritten = numWritten + success
                        
                    # Add to the struct array
                    mapScenarios['numMapScenarios'] = scenarioParams
                else:
                    print('WARNING! Worklist scenario %u is not solved. \n' % i)
                        
            elif numMatch == 1:
                    # Already in mapfile, overwrite (exportSeries() will warn)
                    id = np.where(match > 0)
                    if scenario.isSolved():
                        success    = OutputWriter.exportSeries( scenario, id )
                        numWritten = numWritten + success
                    else:
                        print( 'WARNING! Worklist scenario %u is not solved. \n' % i )
                        
            else:
                raise Exception( 'Non-unique scenarios in mapfile.')   
        
        # Recreate the IDs and write the mapfile
        map = pd.DataFrame(mapScenarios)
        map.to_csv(mapfile)
        
        # Generate interface dependencies file
        with open(os.join.path(PathFinder.getSeriesOutputDir(), 'dependencies.csv'), 'w') as fid:
            fid.write('Component,Interface,Version\n')
            for r in PathFinder.getInputSet():
                fid.write('%s,%s,%s\n' % (r[0,0], r[0,1], r[0,2]))
        
        return numWritten
    
    
        
    ## Write out the scenario for the 'series' interface
    #    Inputs:    
    #               scenario:   Scenario to write 
    #               ID:         ID of scenario in the mapfile
    #    Outputs:
    #               Aggregates.csv  file with delta'd real series
    #               Dynamics.csv    file with raw model outputs
    #    Returns: 1 if wrote, 0 otherwise
    @staticmethod
    def exportSeries( scenario, ID ):
        
        success     = 0
        outputDir   = OutputWriter.getOutputDir( ID )
        if not os.path.isifile('dir'):
            os.mkdir(outputDir)
        else:
            print( 'WARNING! OutputWriter.exportSeries() is overwriting scenario ID %s\n' % int(ID))
        
        # Load dynamic and static variables
        cacheDir    = PathFinder.getCacheDir(scenario)
        Dynamic     = sio.loadmat(os.path.join(cacheDir, 'dynamics.mat'))
        Market      = sio.loadmat(os.path.join(cacheDir, 'market.mat'  ))
        if scenario.isCurrentPolicy() or scenario.postShock().isCurrentPolicy():
            Static       = Dynamic
            StaticMarket = Market
        else:
            Static       = sio.loadmat(os.path.join(cacheDir, 'statics.mat'))
            StaticMarket = sio.loadmat(os.path.join(PathFinder.getCacheDir(scenario.currentPolicy()), 'market.mat' ))

        # Load Dynamics for Microsim baseline add factor fix-up   
        # NOTE: Fix-up factor is going to be Dynamic_base/Dynamic_open_base
        try:
            Dynamic_open_base   = sio.loadmat(os.join.path(PathFinder.getCacheDir(scenario.currentPolicy().open()), 'dynamics.mat'))
            Dynamic_base        = sio.loadmat(os.join.path(PathFinder.getCacheDir(scenario.baseline()), 'dynamics.mat'))
        except:
            raise Exception('WARNING! Cannot read files to make "Dynamic baseline". Skipping...\n' )
            return success
        
        ## Write AGGREGATES.csv
        
        ##
        # Build the source series.
        #   For static series, read from input interfaces
        firstYear = scenario.TransitionFirstYear
        lastYear  = scenario.TransitionLastYear - 1
        numYears  = lastYear - firstYear + 1

        pathFinder = PathFinder(scenario)
        source_series = {}
        
        projections_file = pathFinder.getProjectionsInputPath( 'Projections' )
        source_series['projections'] = InputReader.read_series(projections_file, 'Year', firstYear, lastYear)
               
        taxcalculator_file = pathFinder.getTaxCalculatorInputPath( 'Aggregates' )
        source_series['taxcalculator'] = InputReader.read_series(taxcalculator_file, 'Year', firstYear, lastYear )
        
        oasicalculator_file = pathFinder.getOASIcalculatorInputPath( 'aggregates' );
        source_series['oasicalculator'] = InputReader.read_series(oasicalculator_file, 'Year', firstYear, lastYear );
             
        source_series['Market'] = Market
        
        # Add GDP deflator changes series
        p_series        = InputReader.read_series(projections_file, 'Year', firstYear-1, lastYear)
        gdp_deflator    = p_series['GDPDeflator']
        inflation_rate  = np.ones((1,numYears))
        for i in range(numYears):
            inflation_rate[i] = gdp_deflator[i+1]/gdp_deflator[i]
        
        # Construct dynamic scaling series

        ##  Helper function 
        @staticmethod
        def makeDeltaSeries( Source1, Source2, var_name ):
            
            if var_name == '_add_inflation_to_interestrate':
                delta   = inflation_rate
            elif var_name == '_asis':
                delta   = np.ones((1, numYears))
            elif var_name == '_nonindexed':
                series1 = [np.ones((1,10)), Source1['outs'][10:-1]]
                series2 = [np.ones((1,10)), Source2['outs'][10:-1]]
                delta   = series1 / series2
                delta   = delta[0:numYears]     # truncate in case too long
            else:
                series1 = Source1[var_name]
                series2 = Source2[var_name]
                delta   = series1 / series2
                    
            return delta

        dynamic_series = {}
                
        # Iterate over series names
        for o in OutputWriter.series_names.keys():
            series_name = o
            var_name    = OutputWriter.series_names[series_name]['var_name']
            source      = OutputWriter.series_names[series_name]['source']
            source_name = OutputWriter.series_names[series_name]['source_name']
            if len(source_name) == 0: 
                source_name = series_name
            
            # Calculate Dynamic/Static delta
            delta = makeDeltaSeries( Dynamic, Static, var_name )
            
            # Calculate fixup factor for microsim (as delta open_baseline/baseline)
            fixup = makeDeltaSeries( Dynamic_base, Dynamic_open_base, var_name )

            # Calculate scaling series
            v_scale = delta * fixup   # Adjust by fix-up factor
            v_scale[np.isnan(v_scale)] = 1
            
            # Apply to source series
            v_source = source_series[OutputWriter.series_names[series_name]['source']]['source_name']
            if np.size(v_source,2) == numYears:
                v_source = np.transpose(v_source)   # frickin' matlab and its vector direction
              
            if var_name == '_add_inflation_to_interestrate':
                dynamic_series[series_name] = ((1+v_source) * delta) - 1 
            else:
                dynamic_series[series_name] = v_source * v_scale
            
        # Write series to file
        series_table = pd.DataFrame(dynamic_series)
        
        series_table.to_csv((os.path.join(outputDir, 'Aggregates.csv')))
                
        ## Write DYNAMICS.CSV
            
        Dynamic['outvars']      = OutputWriter.dynamic_outvars
        Market['outvars']       = OutputWriter.market_outvars
        Static['outvars']       = Dynamic.outvars
        StaticMarket['outvars'] = Market.outvars
        
        # Create new column names for Static
        for o in Static['outvars'].keys():
            p = Static['outvars'][o]
            Static['outvars'][o] = 'STATIC_' + p
        
        for o in StaticMarket['outvars'].keys():
            p = StaticMarket['outvars'][o]
            StaticMarket['outvars'][o] = 'STATIC_' + p
        
        # Load steady state variables
        Dynamic_steady = sio.loadmat(os.join.path(PathFinder.getCacheDir(scenario.currentPolicy().steady()), 'dynamics.mat'))
        Market_steady  = sio.loadmat(os.join.path(PathFinder.getCacheDir(scenario.currentPolicy().steady()), 'market.mat'))
                
        # Append steady state variables and reset the first year
        firstYear = firstYear - 1
        for o in Dynamic['outvars'].keys():
            Dynamic[o] = np.hstack((Dynamic_steady[o], Dynamic[o]))
            Static[o]  = np.hstack((Dynamic_steady[o], Static[o]))
        
        for o in Market['outvars'].keys():
            Market[o]       = np.hstack((Market_steady[o], Market[o]))
            StaticMarket[o] = np.hstack((Market_steady[o], StaticMarket[o]))
                
        # Concatenate structs
        output_series = {}

        for M in [Dynamic, Market, Static, StaticMarket]:
            for o in M['outvars'].keys():
                p = M['outvars'][o]
                output_series[p] = M[o]
        
        # Write series to file
        series_table = pd.DataFrame(output_series)
        series_table.to_csv(os.join.path(outputDir, 'Dynamics.csv' ))
        
        success = 1
        
        return success
    
    ## Identify output location
    @staticmethod
    def getOutputDir( ID ):
        seriesDir   = PathFinder.getSeriesOutputDir()
        outputDir   = os.join.path(seriesDir, str(ID))
    
        return outputDir
