##
#
# Dynamic Model Generate Report Output

from paramGeneratorModule import ParamGenerator
from pathFinderModule import PathFinder

import scipy.io as sio
from scipy import interpolate
import warnings
import os
import numpy as np
import math
import pandas as pd
import shutil
import openpyxl
import numbers

class GenerateReports: 
    
    #External dependencies: Pathfinder.m, ParamGenerator.m
    
    # If there is no structure passed with the information on how to
    # structure the output, the program uses the one in the properties by
    # default. 
    
    # Group 1: (main elements, always required)
    
    # variable: This is a list of variables separated by commas. The variable
    # names can be in markets or aggregates, or it can be a new variable if
    # it's a function of the household distributions.
    
    # type: One of these for each variable in. It's either level
    # (to report in level), per_diff (percent difference), delta (divided),
    # or lev_diff (level difference).
    
    # relation: If type is level, this is "b", "s", or "d" for
    # base, static, or counter respectively. Any other choice in Line 5 is
    # "x_y" where x and y are "b", "s", or "d" respectively. For example,
    # if it were "c_b" and "per_diff", it would be the percent difference
    # counter results compared to baseline. Only "b", "s", or "d" are
    # available for the "moments" reports.
    
    # outputNames: these are the names displayed on the graph or on the
    # table. Note that on tables, the spaces are replaced with underscores
    # because spaces are not part of valid table variable names. 
    
    
    # Group 2: (file i/o, optional)
    
    # filename: The name of the file to which the output is saved. If left
    # as an empty string, the program generates the file name from the
    # fields in the scenario structure as it generates working folder
    # directories.
    
    # outputFormat: CSV or XLS (note: no xlsx format, that's too slow and
    # unreliable)
    
    # outputSheet: The name of the sheet to which the output is saved, if
    # outputFormat is XLS. Otherwise this is ignored.
    
    
    # Group 3: (graphing, optional)
    
    # axis: If this exists, the program assumes you're making
    # graphs instead of text output. The values are 'l' and 'r', and should
    # be the same number as the number of variables.  Note: the
    # "outputNames" are used for the legend.  The 'l' corresponds to left
    # axis, and the 'r' corresponds to the right axis.  
    
    # xlabel, llabel, rlabel (optional): These are the labels for the x-axis, the left
    # y-axis, and the right y-axis respectively. The x-label will be
    # replaced with the "years" for any variable that is indexed by the
    # years.  These can be ommitted.  
    
    # title (optional): This is the title of the graph.  If ommitted, the title is
    # blank.
    
    
    # Group 4: (percentiles, optional)
    
    # firstVar: If this exists, the program assumes you're
    # making a new variable from the household distributions. This is the
    # variable for the percentile you're reporting. For example, "TAX"
    # would report the tax paid by the household.
    
    # secondVar: This is the
    # variable against which you're building the quantile. So if you had
    # TAX for firstVar and TOTINC for the second var, it would be the tax
    # paid by a particular quantile of total income.  Using the same
    # variable for the second means that the quantile is of itself.  So TAX
    # would be a particular quantile of total tax paid. 
    
    # percentTop: The
    # upper bound of the quantile being constructed. Must be greater than 0
    # and less than or equal to 100. Must be greater than percentBottom (below).
    
    # percentBottom: The
    # lower bound of the quantile being constructed. Must be greater than
    # or equal to 0, and less than 100.  Must also be less than
    # percentTop. If you want to construct a CDF, just set this to 0.
    
    # threshold: instead of
    # returning the quantile, it returns the value at the top of the
    # quantile. Valid inputs are 'y' and 'n'.
    
    
    # Group 5: (other options, all are optional except for 'year' in
    # certain circumstances):
    
    # startYear: This is the start year for the series on the
    # graph or table. This is only available for variables indexed by time.
    
    # endYear: This is the end year for the series on the graph
    # or table. This is only available for variables indexed by time. Must
    # be greater than startYear.
    
    # year: This is the year that is being reported
    # for various moments. This is required for variables not indexed by
    # time.
    
    # scaleFactor: A scale factor that's multiplied by the
    # variable. Note, although not disabled for "deltas" and "per_diffs",
    # it shouldn't do anything if those options are chosen.
    
    # dollar: (optional, but required if scaleFactor exists): If you want
    # to use the modelunit_dollar, will be multiplied by scaleFactor,
    # choose 'y'.  If you do not want to use this number, choose 'n'.
                           
    # cohort: This is the cohort
    # that you want to follow. The variables you choose must be in OPTs if
    # you choose this option. Average values by cohort will be selected for
    # each year of the cohort you select. 
    
    # Group 6: (summary statistics): Variable name must be from the
    # household problem or one of the states (AGE, LABOR, RETIRED,
    # WORKING, EARNINGS_HISTORY, EARNINGS_SHOCK, CITIZEN).  'year' is a required
    # option.
    
    # firstCond: If creating summary statistics, this one is requred. This
    # is the first conditional.  For example, choosing AGE will create a
    # summary statistic based on the range of ages selected. 
    
    # firstCondUpper (Lower) (optional, but required if firstCond is there: 
    # This is the upper and lower bound for the
    # first conditional, inclusive.  For example, upper 30 and lower 25
    # leads you to report the variable for everyone ages 25 to 30
    # inclusive. If selecting a varibale denominated in model units such as
    # taxable income or consumption, choose the value in USD, and the
    # program will automatically translate to model units.
    
    # secondCond: (optional):  This
    # is the first conditional.  For example, choosing AGE will create a
    # summary statistic based on the range of ages selected.
    
    # secondCondUpper (Lower) (optional, but required if secondCond is there: 
    # This is the upper and lower bound for the
    # second conditional, inclusive.  For example, upper 30 and lower 25
    # leads you to report the variable for everyone ages 25 to 30
    # inclusive. If selecting a varibale denominated in model units such as
    # taxable income or consumption, choose the value in USD, and the
    # program will automatically translate to model units.
    

    # This is the dictionary of variables that are from the markets.mat
    # file. The variable on the left is the name in the reporting program,
    # and the variable on the right of each pair is the name of the
    # variable in markets.mat. If this variable does not exist in the
    # markets.mat, a warning will be thrown.
    
    dictMarkets     = [     ['beqs', 'beqs'],
                            ['capgains', 'capgains'],
                            ['capsharesPM', 'capsharesPM'],
                            ['invtocaps', 'invtocaps'],
                            ['equityDividendRates', 'equityDividendRates'],
                            ['rhos', 'rhos'],
                            ['MPKs', 'MPKs'],
                            ['wages', 'wages'],
                            ['priceindices', 'priceindices'],
                            ['equityPrices', 'equityPrices'],
                            ['bondPrices', 'bondPrices'],
                            ['corpDividendRates','corpDividendRates'],
                            ['bondDividendRates', 'bondDividendRates'],
                            ['passDividendRates','passDividendRates']]
                       
    # This is the dictionary of variables that are from the <simType>.mat
    # file. The variable on the left is the name in the reporting program,
    # and the variable on the right of each pair is the name of the
    # variable in <simType>.mat. If this variable does not exist in the
    # <simType>.mat, a warning will be thrown.
    
    dictAggregates = [      ['lumpSumTaxes', 'lumpSumTaxes'],
                            ['infraSpending', 'infraSpending'],
                            ['constax', 'constax'],
                            ['pops', 'pops'],
                            ['bequests', 'bequests'],
                            ['labs', 'labs'],
                            ['labeffs', 'labeffs'],
                            ['lfprs', 'lfprs'],
                            ['incs', 'incs'],
                            ['pits', 'pits'],
                            ['ssts', 'ssts'],
                            ['bens', 'bens'],
                            ['cons', 'cons'],
                            ['assets', 'assetsPM'],
                            ['corpDividends', 'corpDividends'],
                            ['corpTaxs', 'corpTaxs'],
                            ['revs', 'revs'],
                            ['corpDebts', 'corpDebts'],
                            ['caps', 'caps'],
                            ['outs', 'outs'],
                            ['caps_domestic', 'caps_domestic'],
                            ['caps_foreign', 'caps_foreign'],
                            ['invest_foreign', 'invest_foreign'],
                            ['debts', 'debts'],
                            ['Gtilde', 'Gtilde'],
                            ['Ttilde', 'Ttilde'],
                            ['Ctilde', 'Ctilde'],
                            ['debts_domestic', 'debts_domestic'],
                            ['debts_foreign', 'debts_foreign'],
                            ['tot_assets', 'tot_assetsPM'],
                            ['labincs', 'labincs'],
                            ['capincs', 'capincs'],
                            ['GNI', 'GNI'],
                            ['investment', 'investment'],
                            ['laborIncomes', 'laborIncomes'],
                            ['laborIncomeSubjectToSocialSecuritys','laborIncomeSubjectToSocialSecuritys'],
                            ['capitalIncomes','capitalIncomes'],
                            ['capitalIncomeSubjectToPITs','capitalIncomeSubjectToPITs'],
                            ['averageEffectivePITRates','averageEffectivePITRates'],
                            ['corpForeignWithholding','corpForeignWithholding'],
                            ['passForeignWithholding','passForeignWithholding'],
                            ['averageEffectivePITRatesProductSocialSecurityBenefits','averageEffectivePITRatesProductSocialSecurityBenefits']]   
                       
    # This is the dictionary of variables that are from the decisions.mat
    # file. The variable on the left is the name in the reporting program,
    # and the variable on the right of each pair is the name of the
    # variable in decisions.mat. If this variable does not exist in the
    # decisions.mat, a warning will be thrown. These are the only variables
    # that can be used to construct the percentiles.
    
    dictDecisions  = [      ['V', 'V'],
                            ['CONSUMPTION_TAX', 'CONSUMPTION_TAX'],
                            ['LABOR', 'LABOR'],
                            ['SAVINGS', 'SAVINGS'],
                            ['CONSUMPTION', 'CONSUMPTION'],
                            ['AVG_EARNINGS', 'AVG_EARNINGS'],
                            ['TAXABLE_INC', 'TAXABLE_INC'],
                            ['OASI_BENEFITS', 'OASI_BENEFITS'],
                            ['ORD_LIABILITY', 'ORD_LIABILITY'],
                            ['PAYROLL_LIABILITY', 'PAYROLL_LIABILITY'],
                            ['PREF_LIABILITY', 'PREF_LIABILITY']]  
    
    @staticmethod
    def sampleReport(structName):
        
        # This has a bunch of sample reports.  Bundling this in with the
        # Generate Reports module.  Pass the name of the report that you
        # want, and the report structure will be returned.
        
        reports = {}
    
        # CBO reports for Jagadeesh
        
        reports['CBOage'] = {'filename': 'cboResults', # If empty, then will default to scenario name.
               'outputFormat': 'xls', # Choose "csv" or "xls".
               'outputSheet': 'test', # Sheetname to write output to if format is "xls".
               'outputNames': 'Consumption by Age',
               'variable': 'cons_age',
               'type': 'level',
               'relation': 'c',
               'year': '2048'}
        
        # This was a report for the 2018 CBO presentation on Dynamic
        # models. These are the aggregate values we used to generate the
        # graphs.
        
        reports['cboAggBase'] = {'filename': 'cboAGG', # If empty, then will default to scenario name.
               'outputFormat': 'xls', # Choose "csv" or "xls".
               'outputSheet': 'base', # Sheetname to write output to if format is "xls".
               'outputNames': 'GDP,GNI,Consumption,Investment,Capital,Hours,Average Hours per Worker,Wage Rate,Total Equity Rates,Govt Interest Rate,Social Security Taxes,Personal Income Taxes,Corporate Income Taxes,Social Security Benefits,Government Spending,Debt,Debt to GDP,Interest Paid,Savings,employment,Gtilde,Ttilde,Ctilde,revenues,Total Rates',
               'variable': 'outs,GNI,cons,investment,caps,labs,hoursWorker,wages,equityDividendRates,bondDividendRates,ssts,pits,corpTaxs,bens,spending,debts,debtToGDP,interestPaid,savings,employment,Gtilde,Ttilde,Ctilde,revs,totrates',
               'type': 'level,level,level,level,level,level,level,level,level,level,level,level,level,level,level,level,level,level,level,level,level,level,level,level,level',
               'relation': 'b,b,b,b,b,b,b,b,b,b,b,b,b,b,b,b,b,b,b,b,b,b,b,b,b'}
        
        # This is a report used for testing the fiscal variables; it
        # outputs base and counter for all of the determinants of the
        # deficits and debts.
                                            
        reports['fiscal'] = {'filename': 'fiscal', # If empty, then will default to scenario name.
               'outputFormat': 'xls', # Choose "csv" or "xls".
               'outputSheet': 'fiscal', # Sheetname to write output to if format is "xls".
               'variable': 'constax, pits,ssts,corpTaxs,corpForeignWithholding,passForeignWithholding,Gtilde,bens,Ttilde,Ctilde,bondDividendRates,revs,debts,caps,outs,constax,pits,ssts,corpTaxs,corpForeignWithholding,passForeignWithholding,Gtilde,bens,Ttilde,Ctilde,bondDividendRates,revs,debts,caps,outs',
               'outputNames': 'b_constax,b_pits,b_ssts,b_corpTaxs,b_corpForeignWithholding,b_passForeignWithholding,b_Gtilde,b_bens,b_Ttilde,b_Ctilde,b_bondDividendRates,b_revs,b_debts,b_caps,b_outs,c_constax,c_pits,c_ssts,c_corpTaxs,c_corpForeignWithholding,c_passForeignWithholding,c_Gtilde,c_bens,c_Ttilde,c_Ctilde,c_bondDividendRates,c_revs,c_debts,c_caps,c_outs',
               'type': 'level, level,level,level,level,level,level,level,level,level,level,level,level,level,level,level,level,level,level,level,level,level,level,level,level,level,level,level,level,level',
               'relation': 'b,b,b,b,b,b,b,b,b,b,b,b,b,b,b,c,c,c,c,c,c,c,c,c,c,c,c,c,c,c'}
        
        # This is a report used for testing the fiscal variables; it
        # outputs base and counter for all of the determinants of the
        # deficits and debts.
                                            
        reports['infra'] = {'filename': 'infra', # If empty, then will default to scenario name.
               'outputFormat': 'xls', # Choose "csv" or "xls".
               'outputSheet': 'tmp', # Sheetname to write output to if format is "xls".
               'variable': 'debts,outs,labs,wages,caps,labeffs,lumpSumTaxes,infraSpending,constax,pits,ssts,corpTaxs,corpForeignWithholding,passForeignWithholding,Gtilde,bens,Ttilde,Ctilde,bondDividendRates,revs,debts,caps,labeffs,outs',
               'outputNames': 'ddebts,douts,dlabs,dwages,dcaps,dlabeffs,lumpSumTaxes,infraSpending,constax,pits,ssts,corpTaxs,corpForeignWithholding,passForeignWithholding,Gtilde,bens,Ttilde,Ctilde,bondDividendRates,revs,debts,caps,labeffs,outs',
               'type': 'per_diff,per_diff,per_diff,per_diff,per_diff,per_diff,level       ,level        ,level  ,level,level, level,level                 ,level                 ,level ,level,level,level, level            ,level,level,level,level,level',
               'relation': 'c_b,c_b,c_b,c_b,c_b,c_b,c,c,c,c,c,c,c,c,c,c,c,c,c,c,c,c,c,c'}
                                            
                                            
        # This was a report for the 2018 CBO presentation on Dynamic
        # models. These are the densities at each of the grid points for
        # capital.
        
        reports['cboAggWealth'] = {'filename': 'cboResults', # If empty, then will default to scenario name.
               'outputFormat': 'xls', # Choose "csv" or "xls".
               'outputSheet': 'baseWealth', # Sheetname to write output to if format is "xls".
               'outputNames': 'Density of Capital for Age 65 in 2031,Density of Capital for Age 65 in 2050',
               'variable': 'capDensity,capDensity',
               'type': 'level,level',
               'relation': 'c,c',
               'year': '2031,2050'}
                                            
        # This was a report for the 2018 CBO presentation on Dynamic
        # models. This is the average consumption for each decile of
        # consumption for retirees.
                                            
        reports['cboAggCons'] = {'filename': 'cboResults', # If empty, then will default to scenario name.
               'outputFormat': 'xls', # Choose "csv" or "xls".
               'outputSheet': 'baseCons', # Sheetname to write output to if format is "xls".
               'outputNames': 'Consumption for Retirees in 2031,Consumption for Retirees in 2050',
               'variable': 'consDensity,consDensity',
               'type': 'level,level',
               'relation': 'c,c',
               'year': '2031,2050'}
                                            
                                            
        # This is a test of EV.  Use the field "firstCond" and you will get
        # a summary statistic.  You can condition on one or two variables,
        # secondCond and its affiliates are optional.  In this example, we
        # use EV (which is computed explicitly by the program) or you can
        # use any of the other household variables (including AGE, RETIRED,
        # TAXABLE_INC, CONSUMPTION, etc.).  
                                            
        reports['EVTest'] = {'outputNames': 'EV',
               'variable': 'EV',
               'year': '2019',
               'firstCond': 'AGE',
               'firstCondUpper': '60',
               'firstCondLower': '45',
               'secondCond': 'PRODUCTIVITY_INDEX',
               'secondCondUpper': '1',
               'secondCondLower': '1',
               'relation': 'c'}
                                            
        reports['summaryTest'] = {'outputNames': 'TAXABLE_INC',
               'variable': 'TAXABLE_INC',
               'year': '2018',
               'firstCond': 'AGE',
               'firstCondUpper': '30',
               'firstCondLower': '30',
               'secondCond': 'TAXABLE_INC',
               'secondCondUpper': '100000',
               'secondCondLower': '50000',
               'relation': 'c'}
        
        # This was a report for the 2018 CBO presentation on Dynamic
        # models. This is the lorenz curve for consumption in the baseline
        # and in the anticipated policy (as well as the 45-degree line) in
        # 2031. Change "year" to get the 2050 value.
        
        reports['cboLorenz'] = {'variable': 'c_lorenz, c_lorenz, c_45',
               'outputNames': 'Baseline, Anticipated Policy, Line 45 Degrees',
               'type': 'level,level,level',
               'relation': 'b, c, b',
               'axis': 'l, l, l',
               'title': 'Lorenz Curves 2031',
               'year': '2031, 2031, 2031'}
        
        # This was a report for the 2018 CBO presentation on Dynamic
        # models. This is the lorenz curve for consumption for retirees in
        # the unanticipated shock model. This graph is layered on top of
        # the cboLorenz graph to get one graph with the baseline,
        # anticipated, unanticipated, and 45-degree lines.  
        
        reports['cboLorenz2'] = {'variable': 'c_lorenz',
               'outputNames': 'Unanticipated Policy',
               'type': 'level',
               'relation': 'c',
               'axis': 'l',
               'title': 'Lorenz Curves 2031',
               'year': '2031'}
        
        # This is a reporting structure that outputs the values from
        # report_deltasValues.m.
                                            
        reports['deltas'] = {'filename': 'deltas',      # If empty, then will default to scenario name.
               'outputFormat': 'xls',      # Choose "csv" or "xls".
               'outputSheet': 'deltas', # Sheetname to write output to if format is "xls".
               'variable': 'MPKs,bondDividendRates,equityDividendRates,equityPrices,wages,rhos,outs,caps,caps_foreign,labs,labeffs,assets,cons,investment,debts,revs,tax,ssts,pits,corpTaxs,totinc,totincwss,totrates,MPKs,bondDividendRates,equityDividendRates,equityPrices,wages,rhos,totrates,tax,pits,ssts,corpTaxs',
               'outputNames': 'c_MPKs,c_bondDividendRates,c_equityDividendRates,c_equityPrices,c_wages,c_rhos,delta_outs,delta_caps,delta_caps_foreign,delta_labs,delta_labeffs,delta_assets,delta_cons,delta_investment,delta_debts,delta_revs,delta_tax,delta_ssts,delta_pits,delta_corpTaxs,delta_totinc,delta_totincwss,c_totrates,b_MPKs,b_bondDividendRates,b_equityDividendRates,b_equityPrices,b_wages,b_rhos,b_totrates,c_tax,c_pits,c_ssts,c_corpTaxs',
               'type': 'level,level,level,level,level,level,delta,delta,delta,delta,delta,delta,delta,delta,delta,delta,delta,delta,delta,delta,delta,delta,level,level,level,level,level,level,level,level,level,level,level,level',
               'relation': 'c,c,c,c,c,c,c_s,c_s,c_s,c_s,c_s,c_s,c_s,c_s,c_s,c_s,c_s,c_s,c_s,c_s,c_s,c_s,c,b,b,b,b,b,b,b,c,c,c,c'}

        # This is a reporting structure that is an example of how to use
        # the graphing functions of the reporting programs. It shows
        # percent differences between baseline / counter economies for
        # output and capital, and shows them for 2021 through 2030.
                                            
        reports['exampleGraph'] = {'variable': 'outs, caps',
               'outputNames': 'Output per change,Capital per change',
               'type': 'per_diff,per_diff',
               'relation': 'c_b,c_b',
               'axis': 'l,r',
               'title': 'Output, Capital in Closed Economy',
               'startYear': '2021',
               'endyear': '2025'}
        
        # This report breaks down savings by income; it was used to test
        # the Larson bill, which had large effects on savings in the top 5
        # percent of the income distribution.
                                            
        reports['savingsIncome'] = {'filename': 'savings',
               'outputFormat': 'xls',
               'outputSheet': 'savings',
               'variable': 'sav05, sav10, sav15, sav20, sav25, sav30, sav35, sav40, sav45, sav50, sav55, sav60, sav65, sav70, sav75, sav80, sav85, sav90, sav95, sav00',
               'outputNames': 'sav05, sav10, sav15, sav20, sav25, sav30, sav35, sav40, sav45, sav50, sav55, sav60, sav65, sav70, sav75, sav80, sav85, sav90, sav95, sav00',
               'type': 'level,level,level,level,level,level,level,level,level,level,level,level,level,level,level,level,level,level,level,level',
               'relation': 'c,c,c,c,c,c,c,c,c,c,c,c,c,c,c,c,c,c,c,c',
               'firstVar': 'SAVINGS,SAVINGS,SAVINGS,SAVINGS,SAVINGS,SAVINGS,SAVINGS,SAVINGS,SAVINGS,SAVINGS,SAVINGS,SAVINGS,SAVINGS,SAVINGS,SAVINGS,SAVINGS,SAVINGS,SAVINGS,SAVINGS,SAVINGS',
               'secondVar': 'TAXABLE_INC,TAXABLE_INC,TAXABLE_INC,TAXABLE_INC,TAXABLE_INC,TAXABLE_INC,TAXABLE_INC,TAXABLE_INC,TAXABLE_INC,TAXABLE_INC,TAXABLE_INC,TAXABLE_INC,TAXABLE_INC,TAXABLE_INC,TAXABLE_INC,TAXABLE_INC,TAXABLE_INC,TAXABLE_INC,TAXABLE_INC,TAXABLE_INC',
               'percentTop': '5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100',
               'percentBottom': '0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95',
               'threshold': 'n,n,n,n,n,n,n,n,n,n,n,n,n,n,n,n,n,n,n,n'}
                                            
        # This is a reporting structure that is an example of how to the
        # reporting plot moments from a particular year. In this case, it
        # plots average consumption and tax by age for the year 2021.
                        
        reports['exampleMomentsGraph'] = {'variable': 'cons_age, tax_age',
               'outputNames': 'Consumption by Age, Tax by Age',
               'type': 'level,level',       
               'relation': 'c, c',
               'axis': 'l, l',       
               'title': 'Moments by Age',
               'xlabel': 'Age',
               'llabel': 'Model Units',
               'rlabel': 'Y2-Axis Label',
               'year': '2021, 2021'}

        # This is a reporting structure that is an example for how to save
        # data to a table without actually saving it to an excel file.
        # createOutput will just return this table directly to whatever
        # program calls it.
                                        
        reports['exampleMomentsTable'] = {'title': 'Distributions of Assets',
               'variable': 'a_distmodel_thresh,a_distmodel_cumul,a_distdata_thresh,a_distdata_cumul',
               'type': 'level,level,level,level',
               'relation': 'c,c,c,c',
               'outputNames': 'Model in 2016 Dollars, Model CDF, Data in 2016 Dollars, Data CDF',
               'year': '2018, 2018, 2018, 2018'}
        
        # This is a reporting structure that is an example for how to use
        # the cohort reporting features. The "cohort" value is the cohort
        # value in the steady state.  It follows the average value of a
        # "decisions" variable for that cohort through the range of years
        # covered. Also note that these values can be paired with aggregate
        # values like output.  Under "cohort", the value is "n" to indicate
        # that it's another aggregate value.
                                            
        reports['exampleMomentsGraphCohort']  = {'filename': 'oasi',      # If empty, then will default to scenario name.
               'outputFormat': 'xls',      # Choose "csv" or "xls".
               'outputSheet': 'OASI of Retired People',      # Sheetname to write output to if format is "xls".
               'cohort': '0,1,60,n',
               'variable': 'OASI_BENEFITS,OASI_BENEFITS,OASI_BENEFITS,outs',
               'outputNames': 'OASI_base,OASI_counter,OASI_old,outs',
               'type': 'level,level,level,per_diff',
               'relation': 'b,c,c,c_b',
               'axis': 'l,l,l,r',
               'startYear': '2019',
               'endYear': '2025'}
    
        reportTemplate = reports[structName]
        
        return reportTemplate

    @staticmethod
    def createOutput(scenario, reportTemplate):
        
        # Used to create the output sheet.
        
        # scenario (required): This can be either a scenario structure
        # that is passed directly, or it can be a character array with the
        # filename (with path) for a .mat file with a scenario structure
        # stored in it.
        
        # reportTemplate (required): This can be either an reportTemplate
        # structure (see above for the format) passed directly, or it can
        # be a character array with the filename (with path) for a .mat
        # file with a reportTemplate structure stored in it. This is an
        # OPTIONAL argument, if none is provided, the program defaults to
        # the reportTemplate structure class property.
        
        # We're checking to see if "scenario" is a filename, if so, we'll
        # load it (must be a .mat) file. Otherwise we'll just use the
        # structure that's passed (after 'de-nesting' it, if necessary.
        # Checking to see if it's a steady-state economy, in which case the
        # transition is skipped.
        
        if isinstance(scenario, str):
            scenario = sio.loadmat(scenario)
        scenario = GenerateReports.deNestStructure(scenario)
        
        steadyState = False
        if scenario.isSteady():
            steadyState = True 
            print('This is a STEADY STATE report, only one year is reported.')

        # If the report template is a file, get the file, otherwise use the
        # structure we're passed.
        
        if isinstance(reportTemplate, str):
            reportTemplate = GenerateReports.getOutputStructure(reportTemplate)
        reportTemplate = GenerateReports.deNestStructure(reportTemplate)
        
        # Creating the output directory if it does not exist
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            os.mkdir(os.path.join(os.getcwd(), 'Output'))
        
        # if the filename is there but an empty string, the filename will
        # be created as the scenario pathname.  If there is no filename
        # field, this will not be created and there will be no output to an
        # Excel or CSV file.
        
        if 'filename' in reportTemplate.keys(): 
            if reportTemplate['filename'] == '':
                reportTemplate['filename'] = np.array([scenario.basedeftag, '_', scenario.counterdeftag])
                          
        # Parsing the structure that defines the output table, making sure
        # that it conforms to specification, returns errors if there are
        # mistakes. Graphs is a boolean if this is meant to
        # create graphs instead of tables for output. Scales is a boolean
        # if there are one or more variables in there that are being
        # scaled.
        # Creates VARS, which is a structure that keeps track of the
        # variables and their attributes we'll need to construct for the
        # output of the table.
        (vars, graphs, scale, experiments) = GenerateReports.parseVars(reportTemplate, scenario, steadyState)
        
        # Load the data into two or four structures. outputDataStructure is
        # structure that will hold all reportable data.
        # Markets are typically market prices. (markets and aggregates get
        # merged into aggregates later).  Dist is the distribution of the
        # population. Decisions are the household decisions. Dist and
        # decisions are only loaded if "distributions" is true. "T" is the
        # structure generated by ParamGenerators with the timings in the
        # model.
        (outputDataStructure, markets, T, dist, decisions) = GenerateReports.openFiles(scenario, steadyState, experiments)

        # Copy the steady state into the first observation
        (outputDataStructure, dist, decisions) = GenerateReports.copySteadyState(outputDataStructure, markets, dist, decisions, steadyState, T, experiments) 

        # Create new aggregate variables for reporting
        outputDataStructure = GenerateReports.createAggVars(outputDataStructure, experiments)
        
        # We construct quantiles with this section.
        (outputDataStructure, vars) = GenerateReports.generateDistributions(outputDataStructure, dist, decisions, vars, scenario, T, reportTemplate.variable, experiments)

        # Generate and export the table and create plots: the table is used
        # for both output to excel or csv, and it is used for the graphs. 
        # outputTable is returned 
        outputTable = GenerateReports.createTable(outputDataStructure, T, vars, scale, graphs)

        # If it's a graph, plot the series
        # If it's not a graph and there's a filename field,
        # save the table to excel or csv
       
        if graphs:
            GenerateReports.plotSeries(outputTable, vars)
        if 'filename' in reportTemplate.keys():
            GenerateReports.saveTable(outputTable, reportTemplate)
        
        return outputTable
      
    
    def MomentTable(scenario, outputVar, varName1, list1, bound1, varName2, list2, bound2, outputYear, economyType):

        # This function generates a large table of conditional moments
        # Parameters:
        # scenario (the scenario being pulled)
        # outputVar: the output variable, typically from OPTs
        # varName1: conditional variable #1
        # list1: the list of varName1's you use
        # bound1: 0 if we're using the actual variable, 1 if we're using
        #   ranges. For example, if you want 20 year olds, list1 = [20] and
        #   this is set to 0.  If we want ranges of width 5 years, we set
        #   list1 = [5,10,15,20] and set this to 0.  Then we'll get 6-10,
        #   11-15, 16-20 year olds.  
        # varName2, list2, bound2: same, but for the second conditional
        #   variable
        # outputYear: the year in which the output is generated. This is a
        #   string.
        # economyType: 'b' for base, 's' for static, 'c' for counter. This
        #   is a string.
    
        # The program starts by creating a report with a different variable for
        # each element of the table.

        T = ParamGenerator.timing(scenario)
    
        # This is the base year to be reported. 
        baseyear = outputYear
        baseage = T['realage_entry']

        # These are going to be offsets depending on whether want a list of
        # items or the range of the items.
        firstOffset = 1*bound1
        secondOffset = 1*bound2

        # Initializing variables such as list lengths and the output table
        firstLen = len(list1)
        secondLen = len(list2)
        # outputTable = np.zeros((firstLen-firstOffset, secondLen-secondOffset))

        # Initializing the fields of the report
        report = {}
        report['year'] = '';
        report['relation'] = '';
        report['outputNames'] = '';
        report['variable'] = '';
        report['firstCondUpper'] = '';
        report['secondCondUpper'] = '';
        report['firstCondLower'] = '';
        report['secondCondLower'] = '';
        report['firstCond'] = '';
        report['secondCond'] = '';
        report['type'] = '';
    
        # Out loop for the first conditional variable  
        for indexOne in range(firstOffset, firstLen):

            for indexTwo in range(secondOffset, secondLen):
            
                # Getting the top and bottom values for the range of the y-axis
                # variable.  Note, if firstOffset = 0, then they're the same
                # value.  We put .00001*firstOffset so that we're getting the range down to
                # but not including the lower value of the range.  
                theYear = baseyear
                valueOne = list1[indexOne]
                valueOneBottom = list1[indexOne - firstOffset] + .00001*firstOffset

                # We don't want ranges that include unborn people, so we're
                # checking to see if that's an issue.  
                if (valueOneBottom < baseage) and varName1 == 'AGE':
                    assert firstOffset == 0, 'Ranges for AGE variable are not compatible with unborn cohorts.'
                    theYear = theYear + ( baseage - valueOne)
                    valueOne = baseage
                    valueOneBottom = baseage
            
                # Getting the top and bottom values for the range of the x-axis
                # variable. Note, if secondOffset = 0, then they're the same
                # value.
                valueTwo = list2[indexTwo]
                valueTwoBottom = list2(indexTwo - secondOffset) + .00001*secondOffset
            
                # We don't want ranges that include unborn people, so we're
                # checking to see if that's an issue.  
                if (valueTwo < baseage) and varName2 == 'AGE':
                    assert secondOffset == 0, 'Ranges for AGE variable are not compatible with unborn cohorts.'
                    theYear = theYear + ( baseage - valueTwo)
                    valueTwo = baseage
                    valueTwoBottom = baseage
            
                firstCondUpper = format(valueOne,'25.10f').strip()
                firstCondLower = format(valueOneBottom,'25.10f').strip()
                secondCondUpper = format(valueTwo,'25.10f').strip()
                secondCondLower = format(valueTwoBottom,'25.10f').strip()
               
                # Creating the elements of the report structure that will be
                # passed to "createOutput" to generate the actual values. We're
                # creating one large report structure so that we do not have to
                # load the files over and over again; they get read once and
                # then all of the values we need are generated in one pass.
                if report['year'] == '':
                    joinStr = ''
                else:
                    joinStr = ', '
            
                report['type'] += joinStr + 'level'
                report['year'] += joinStr + str(theYear)
                report['relation'] += joinStr + economyType
                report['outputNames'] += joinStr + outputVar + '_' + str(indexOne-firstOffset) + '_' + str(indexTwo-secondOffset)
                report['variable'] += joinStr + outputVar + '_' + str(indexOne-firstOffset) + '_' + str(indexTwo-secondOffset)
                report['firstCondUpper'] += joinStr + firstCondUpper
                report['secondCondUpper'] += joinStr + secondCondUpper
                report['firstCondLower'] += joinStr + firstCondLower
                report['secondCondLower'] += joinStr + secondCondLower
                report['firstCond'] += joinStr + varName1
                report['secondCond'] += joinStr + varName2
    
    
        tmpTable = GenerateReports.createOutput(scenario, report)
    
        # Deconstructing the output from createOuput and arranging it into a
        # nice table format for easy consumption.
        rowNames = np.zeros(firstLen-firstOffset)
        columnNames = np.zeros(secondLen-secondOffset)
        outputTable = np.zeros((firstLen-firstOffset, secondLen-secondOffset))
        
        for indexOne in range(firstLen-firstOffset):
            if firstOffset > 0:
                rowNames[indexOne] = format(list1[indexOne] + .00001, '15.5f').strip() + ' to ' + str(list1[indexOne+firstOffset])
            else:
                rowNames[indexOne] = list1[indexOne]
        
            for indexTwo in range(secondLen-secondOffset):
                tableIndex = outputVar + '_' + str(indexOne) + '_' + str(indexTwo)
                outputTable[indexOne, indexTwo] = tmpTable[tableIndex]
            
                if secondOffset > 0:
                    columnNames[indexTwo] = format(list2[indexTwo] + .00001,'15.5f').strip() + ' to ' + str(list2[indexTwo+secondOffset])
                else:
                    columnNames[indexTwo] = list2[indexTwo]

    
        # Filling in the labels on the table.  The first line fills in the
        # y-axis values (note, if it's a range, it only reports the upper value
        # for the range).  The second line fills in the x-axis values, and the
        # top-left value is "NaN" or the bottom boundary of the range.
        
        outputTable = pd.DataFrame(outputTable, columns = columnNames)
        outputTable.set_index(pd.Index(np.hstack((varName1 + ' / ' + varName2, rowNames))))
        
        return outputTable

    @staticmethod
    def createAggVars(outputDataStructure, experiments):
        
        # This creates variables that are a function of existing aggregate
        # variables. These variables are economically interesting, but not
        # generated by the main dynamic model. Make sure to create both the
        # .data and .xaxis components if you're adding new variables to
        # this list.
        
        # Contributions to tax and growth
        for p in experiments:
            for q in ['revs', 'pits', 'ssts', 'corpTaxs']:
                newvarname = 'c_revs_' + q
                assert q in outputDataStructure[p].keys(), 'Creating new aggregate variables error: ' + q + ' is not a recognized variable.'
                outputDataStructure[p][newvarname]['data'] = (outputDataStructure[p][q]['data'] - outputDataStructure['base'][q]['data'])  / outputDataStructure['base']['revs']['data'] * 100
                outputDataStructure[p][newvarname]['xaxis'] = outputDataStructure['base'][q]['xaxis']
        
        # Contributions to GDP
        for p in experiments:
            for q in ['outs', 'cons', 'investment', 'Gtilde']:
                newvarname = 'c_outs_' + q
                assert q in outputDataStructure[p].keys(), 'Creating new aggregate variables error: ' + q + ' is not a recognized variable.'
                outputDataStructure[p][newvarname]['data'] = (outputDataStructure[p][q]['data'] - outputDataStructure['base'][q]['data'] ) / outputDataStructure['base']['outs']['data'] * 100
                outputDataStructure[p][newvarname]['xaxis'] = outputDataStructure[p][q]['xaxis']
            outputDataStructure[p]['c_outs_NX']['data'] = outputDataStructure[p]['c_outs_outs']['data'] - outputDataStructure[p]['c_outs_cons']['data'] - outputDataStructure[p]['c_outs_investment']['data'] - outputDataStructure[p]['c_outs_Gtilde']['data']   # This one negatively contributes to revenues
            outputDataStructure[p]['c_outs_NX']['xaxis'] = outputDataStructure[p]['c_outs_outs']['xaxis']
        
        for p in experiments:
            outputDataStructure[p]['totrates']['data'] = outputDataStructure[p]['capsharesPM']['data'] * outputDataStructure[p]['equityDividendRates']['data'] + (1 - outputDataStructure[p]['capsharesPM']['data']) * outputDataStructure[p]['bondDividendRates']['data']
            outputDataStructure[p]['labeffspops']['data'] = outputDataStructure[p]['labeffs']['data'] / outputDataStructure[p]['pops']['data']
            outputDataStructure[p]['assetspop']['data'] =  outputDataStructure[p]['assets']['data'] / outputDataStructure[p]['pops']['data']
            outputDataStructure[p]['debtToGDP']['data'] = outputDataStructure[p]['debts']['data'] / outputDataStructure[p]['outs']['data'] * 100
            outputDataStructure[p]['spending']['data'] = outputDataStructure[p]['bens']['data'] + outputDataStructure[p]['Gtilde']['data']
            outputDataStructure[p]['interestPaid']['data'] = outputDataStructure[p]['bondDividendRates']['data'] * outputDataStructure[p]['debts']['data']
            outputDataStructure[p]['totrates']['xaxis'] = outputDataStructure[p]['capsharesPM']['xaxis']
            outputDataStructure[p]['labeffspops']['xaxis'] = outputDataStructure[p]['labeffs']['xaxis']
            outputDataStructure[p]['assetspop']['xaxis'] =  outputDataStructure[p]['assets']['xaxis']
            outputDataStructure[p]['debtToGDP']['xaxis'] = outputDataStructure[p]['debts']['xaxis']
            outputDataStructure[p]['spending']['xaxis'] = outputDataStructure[p]['bens']['xaxis']
            outputDataStructure[p]['interestPaid']['xaxis'] = outputDataStructure[p]['bondDividendRates']['xaxis']
        
        return outputDataStructure
    
    @staticmethod
    def generateDistributions(outputDataStructure, dist, decisions, vars, scenario, T, variableList, experiments):
        
        pathFinder = PathFinder(scenario)
        
        # This is the function that creates a time series of the distributions of indicators that
        # charcaterize the households. For example, we will have all of the
        # household taxes paid by state in here; same with labor income, same with
        # corporate income tax, etc.
        households = {}
    
        T_life  = T['T_life']         # Total life years
        T_model = T['T_model'] + 1    # Transition path model years, plus one for steady state      
        yearRange = list(range((T['TransitionFirstYear']-1),T['TransitionLastYear']))
        
        # Define grids
        
        s = ParamGenerator.grids(scenario)
        sSteady = ParamGenerator.grids(scenario.currentPolicy().steady())
        
        nz     = s['nz']         # num labor productivity shock
        nk     = s['nk']         # num asset points
        nb     = s['nb']         # num avg. earnings points
        ng     = s['ng']
        zs     = np.concatenate((sSteady['zs'], s['zs']))         # shocks grid (by demographic type and age)
        kv     = s['kv']         # capital grid
        bv     = s['bv']         # earnings grid

        # Reshaping the distribution so that it will match the decision
        # matricies below.
        distrib = {}
        
        if 'base' in experiments:
            distrib['base'] = np.reshape(dist['base']['DIST'], (-1, T_model))
        if 'counter' in experiments:
            distrib['counter'] = np.reshape(dist['counter']['DIST'], (-1, T_model))
        if 'static' in experiments:
            distrib['static'] = np.reshape(dist['static']['DIST'], (-1, T_model))
        
        # Making sure that all of our distribution is positive.
        
        if 'base' in experiments:
            assert all(x >= 0 for x in distrib['base']),'WARNING! Negative mass of people at DIST (baseline).'
        if 'counter' in experiments:
            assert all(x >= 0 for x in distrib['counter']),'WARNING! Negative mass of people at DIST (counter).'
        if 'static' in experiments:
            assert all(x >= 0 for x in distrib['static']),'WARNING! Negative mass of people at DIST (static).'

        # Creating a matrix of the variables with a distribution. The size
        # is (nz*nk*nb*T_life*ng) x (T_model). 
        
        for economy in experiments:
            temp = np.zeros((nz, nk, nb, T_life, ng, T_model))
            for var_name in decisions[economy]:
                for q in range(ng):
                    temp[:,:,:,:,q,:] = decisions[economy][var_name]               
                households[economy][var_name] = np.reshape(temp, (-1, T_model))
       
        # Creating variables for the grids that are will be two-dimensional
        # and with size (nz*nk*nb*T_life) x (T_model), so that each array
        # is one year of the simulation.
        zLocGrid = np.outer(np.array(range(nz)), np.ones(80))
        zLocArray = np.reshape(np.tile(np.reshape(zLocGrid, (nz,1,1,T_life,1,1)), [1,nk,nb,1,ng,T_model]),(-1, T_model))
        
        z = np.reshape(np.tile(np.reshape(np.transpose(zs, axes = [3,2,1]), (nz,1,1,T_life,1,T_model)), [1,nk,nb,1,ng,1]), (-1, T_model))
        k = np.reshape(np.tile(np.reshape(kv, (1,nk,1,1,1,1)), [nz,1,nb,T_life,ng,T_model]), (-1, T_model))
        age = np.reshape(np.tile(np.reshape(np.array(range(T_life)), (1,1,1,T_life,1,1)), [nz,nk,nb,1,ng,T_model]), (-1, T_model))
        citizen = np.reshape(np.tile(np.reshape(np.array(range(ng)), (1,1,1,1,ng,1)), [nz,nk,nb,T_life,1,T_model]), (-1, T_model))
        b = np.reshape(np.tile(np.reshape(bv, (1,1,nb,1,1,1)), [nz,nk,1,T_life,ng,T_model]), (-1, T_model))
        nossLocation = (b == min(b)) 
        
        # Creating big matricies of the aggregates variables that are
        # applied to the household decision making process. These matrices
        # have the same size (nz*nk*nb*T_life) x (T_model) so that we can
        # do element-by-element multiplication. 
            
        # Creating additional variables such as taxes (TAX) and labor
        # income (LABINC) from existing variables above.
        
        for economy in experiments:

            wages = np.reshape(np.tile(np.reshape(outputDataStructure[economy]['wages']['data'], (1, 1, 1, 1, 1, T_model)), [nz,nk,nb,T_life,ng,1]), (-1, T_model))
            capsharesPM = np.reshape(np.tile(np.reshape(outputDataStructure[economy]['capsharesPM']['data'], (1, 1, 1, 1, 1, T_model)), [nz,nk,nb,T_life,ng,1]), (-1, T_model))
            bondDividendRates = np.reshape(np.tile(np.reshape(outputDataStructure[economy]['bondDividendRates']['data'], (1, 1, 1, 1, 1, T_model)), [nz,nk,nb,T_life,ng,1]), (-1, T_model))
            equityrates = np.reshape(np.repmat(np.reshape(outputDataStructure[economy]['equityDividendRates']['data'], (1, 1, 1, 1, 1, T_model)), [nz,nk,nb,T_life,ng,1]), (-1, T_model))
            
            retiredLocation = households[economy]['OASI_BENEFITS']  > 0
            
            households[economy]['LABOR_INCOME'] =  wages * z * households[economy]['LABOR'] 
            households[economy]['WORKING'] = households[economy]['LABOR'] > 0
            households[economy]['ASSETS'] = k
            households[economy]['PRODUCTIVITY_INDEX'] = zLocArray
            households[economy]['TAX']    = households[economy]['PREF_LIABILITY'] + households[economy]['ORD_LIABILITY'] + households[economy]['PAYROLL_LIABILITY']
            households[economy]['TOTINC'] = households[economy]['LABOR_INCOME'] + households[economy]['ASSETS'] * (np.ones(capsharesPM.shape) - capsharesPM) * bondDividendRates + households[economy]['ASSETS'] * (capsharesPM) * equityrates
            households[economy]['TOTINCWSS'] = households[economy]['TOTINC'] + households[economy]['OASI_BENEFITS']
            households[economy]['CAPITAL_INCOME'] = households[economy]['ASSETS'] * (np.ones(capsharesPM) - capsharesPM) * bondDividendRates + households[economy]['ASSETS'] * capsharesPM * equityrates
            households[economy]['AGE'] = age
            households[economy]['EARNINGS_SHOCK'] = z
            households[economy]['EARNINGS_HISTORY'] = b
            households[economy]['RETIRED'] = retiredLocation
            households[economy]['CITIZEN'] = citizen
        
        # Section for computing summary statistics
        if (vars[1]['summary'] == 1):
            
            for varNum in range(len(vars)):
            
                economy = experiments(max(experiments.shape))
            
                # Checking to see if the first boundaries are for something
                # that's not in model units (age, earnings shocks, working or
                # retired indicators, and amoun to labor.  If not, convert the
                # bounds to model units.
            
                if not (vars[varNum]['firstcond'] =='PRODUCTIVITY_INDEX' or vars[varNum]['firstcond'] == 'AGE' or vars[varNum]['firstcond'] == 'CITIZEN' or vars[varNum]['firstcond'] == 'EARNINGS_SHOCK' or vars[varNum]['firstcond'] == 'RETIRED' or vars[varNum]['firstcond'] == 'WORKING' or vars[varNum]['firstcond'] == 'LABOR'):
                    vars[varNum]['firstcondtop'] =  vars[varNum]['firstcondtop'] * scenario.modelunit_dollar
                    vars[varNum]['firstcondbottom'] =  vars[varNum]['firstcondbottom'] * scenario.modelunit_dollar
                else:
                    if vars[varNum]['firstcond'] == 'AGE':
                        vars[varNum]['firstcondtop'] = vars[varNum]['firstcondtop'] - T['realage_entry'] + 1
                        vars[varNum]['firstcondbottom'] =  vars[varNum]['firstcondbottom'] - T['realage_entry'] + 1
            
                # Same thing for the second set of boundaries, if they are
                # being used.
            
                if 'secondcond' in vars[varNum].keys():    
                    if not (vars[varNum]['secondcond'] == 'PRODUCTIVITY_INDEX' or vars[varNum]['secondcond'] == 'AGE' or vars[varNum]['secondcond'] == 'CITIZEN' or vars[varNum]['secondcond'] == 'EARNINGS_SHOCK' or vars[varNum]['secondcond'] == 'RETIRED' or vars[varNum]['secondcond'] == 'WORKING' or vars[varNum]['secondcond'] == 'LABOR'):
                        vars[varNum]['secondcondtop'] =  vars[varNum]['secondcondtop'] * scenario.modelunit_dollar
                        vars[varNum]['secondcondbottom'] =  vars[varNum]['secondcondbottom'] * scenario.modelunit_dollar
                    else:      
                        if vars[varNum]['secondcond'] == 'AGE':
                            vars[varNum]['secondcondtop'] = vars[varNum]['secondcondtop'] - T['realage_entry']
                            vars[varNum]['secondcondbottom'] =  vars[varNum]['secondcondbottom'] - T['realage_entry']
            
                theYear = vars[varNum]['year'] - min(yearRange) + 1
                
                # We're initializing a matrix of indicators, 1 if it's
                # a summary statistic calculations, 0 if it is not.
                # Checks to see if the first condition is satisfied, sets
                # conditionalMatrix to 1 if it is.  Then it checks to see
                # if the second condition is also satisfied, if not, the
                # value is reduced back to 0.

                conditionalMatrix = np.zeros((households['base']['LABOR']).shape)
                conditionalMatrix = conditionalMatrix + math.trunc(((households['base'][vars[varNum]['firstcond']][:,theYear] <= vars[varNum]['firstcondtop']) + (households['base'][vars[varNum]['firstcond']][:,theYear] >= vars[varNum]['firstcondbottom']))/2)
                if 'secondcond' in vars[varNum].keys():    
                    conditionalMatrix = conditionalMatrix + math.trunc(((households['base'][vars[varNum]['secondcond']][:,theYear] <= vars[varNum]['secondcondtop']) + (households['base'][vars[varNum]['secondcond']][:,theYear] >= vars[varNum]['secondcondbottom']))/2)
                    conditionalMatrix = math.trunc(conditionalMatrix / 2)

                # Initializing the value matrix and the density
                # matrix to zeros.  

                sampleMass = np.zeros((max(conditionalMatrix.shape),1))
                sampleSummary = np.zeros((max(conditionalMatrix.shape),1))

                # Looping over all nz*nb*T_life*ng*nk data points.  We'll
                # evaluate the EV only at the points we need.

                for numZ, numB, numAge, numG, numK in zip(iter(nz), iter(nb), iter(T_life), iter(ng), iter(nk)):

                    # Translating the index in the multidimensional array into the
                    # index for the single-dimensional array
                    currentStateIndexNum = (nk*nb*T_life*nz*(numG-1) + nk*nb*nz*(numAge - 1) + nk*nz*(numB-1) + nz*(numK - 1) + (numZ))  

                    # If these are people we're interested in (as defined
                    # by the reporting program), we'll compute their EVs.
                    # Otherwise they get skipped and their EV and mass is
                    # 0.
                    if conditionalMatrix[currentStateIndexNum]:
                    
                        if vars[varNum]['variable'] == 'EV':
                            # Creates an array of indices that reference the
                            # value function in the baseline across different assets.
                            # Find the capital in the counter that gives the
                            # same utility.
                            # if distrib.base(currentStateIndexNum,theYear) > 0
                            # a = 1;
                            # end
                            baseValueFunctionIndexArray = np.ones((1,nk))*(nk*nb*T_life*nz*(numG-1) + nk*nb*nz*(numAge - 1) + nk*nz*(numB-1) + (numZ)) + np.arange(0, nz*nk-1, nz)   
                            
                            # This is the standard cubic spline, it does
                            # not work very well in our experience, but I
                            # left it here in case someone wants to try it
                            # again when we have more grid points
                        
                            #kSpline = spline(households.base.V(baseValueFunctionIndexArray,theYear), kv, households.(economy{1}).V(currentStateIndexNum,theYear));
                            
                            # This is the shape-preserving cubic spline, it does
                            # not work very well at the bottom end of the grid (extrapolation) in our experience, but I
                            # left it here in case someone wants to try it
                            # again when we have more grid points
                            
                            #kSpline = interp1(households.base.V(baseValueFunctionIndexArray,theYear), kv, households.(economy{1}).V(currentStateIndexNum,theYear),'pchip');
                            
                            # This is the linear interpolation /
                            # extrapolation. This seems to work the best
                            # of the three methods currently.
                            
                            interp_f = interpolate.interp1d(households['base']['V'][baseValueFunctionIndexArray,theYear], kv)
                            
                            kSpline = interp_f(households[economy]['V'][currentStateIndexNum,theYear])
                            
                            sampleSummary[currentStateIndexNum] = (kSpline - k[currentStateIndexNum, theYear])/ scenario.modelunit_dollar
                        
                        elif vars['variable'] == 'CE':
                        
                            # From MACROECONOMIC EFFECTS OF MEDICARE
                            # https://www.nber.org/papers/w23389.pdf
                            # Page 15.
                            
                            sampleSummary[currentStateIndexNum] = ((households[economy]['V'][currentStateIndexNum,theYear] / households['base']['V'][currentStateIndexNum,theYear])**(1/(scenario.gamma * (1- scenario.sigma))) - 1) *100
                            
                        else:
          
                            # If it's not a welfare measure, locate the values
                            # at each of the relevent states.  If it's
                            # something in model units, convert to dollars.
                            
                            # When we're constructing the table, we use
                            # variable names with suffix of _ROW_COL,
                            # and we're deleting those because those
                            # arne't actual variable names. This will
                            # work with a table or if we're generating
                            # a single number, as long as there are no
                            # additional "_" characters in the
                            # variable.
                            varName = vars[varNum]['variable']
                            for charNum in range(len(varName)-2, -1, -1):
                                if varName[charNum] == '_' and np.all(np.in1d(varName[min(charNum+1,len(varName))-1:len(varName)], list('0123456789'))):
                                    varName = varName[0:max(1,charNum-1)]


                            if vars[varNum]['variable'] == 'AGE' or vars[varNum]['variable'] == 'EARNINGS_SHOCK' or vars[varNum]['variable'] == 'RETIRED' or vars[varNum]['variable'] == 'WORKING' or vars[varNum]['variable'] == 'LABOR':
                                sampleSummary[currentStateIndexNum] = households[economy][vars[varNum]['variable']][currentStateIndexNum, theYear]
                            else:
                                sampleSummary[currentStateIndexNum] = households[economy][vars[varNum]['variable']][currentStateIndexNum, theYear]/ scenario.modelunit_dollar
                        
                    sampleMass[currentStateIndexNum] = distrib['base'][currentStateIndexNum,theYear]

            outputDataStructure['base'][vars['variable']]['data'][1:(max(yearRange)-min(yearRange)+1)] = None
            outputDataStructure[economy][vars['variable']]['data'][1:(max(yearRange)-min(yearRange)+1)] = None
            outputDataStructure['base'][vars['variable']]['data'][theYear] = None
            if np.sum(sampleMass) == 0:
                outputDataStructure[economy][vars['variable']]['data'][theYear] = None
            else:
                outputDataStructure[economy][vars['variable']]['data'][theYear] = np.sum(sampleSummary * sampleMass) / np.sum(sampleMass)

            outputDataStructure['base'][vars['variable']]['xaxis'] = yearRange
            outputDataStructure[economy][vars['variable']]['xaxis'] = yearRange
        
        for economy in experiments:
                                            
            # Labor income, labor income with social security, and taxes
            # paid by households, total by year.
                                            
            outputDataStructure[economy]['totinc']['data'] = np.sum(households[economy]['TOTINC'] * distrib[economy])
            outputDataStructure[economy]['totinc']['xaxis'] = yearRange
            outputDataStructure[economy]['savings']['data'] = np.sum(households[economy]['SAVINGS'] * distrib[economy])
            outputDataStructure[economy]['savings']['xaxis'] = yearRange
            outputDataStructure[economy]['totincwss']['data'] = np.sum(households[economy]['TOTINCWSS'] * distrib[economy])
            outputDataStructure[economy]['totincwss']['xaxis'] = yearRange
            outputDataStructure[economy]['tax']['data'] = np.sum(households[economy]['TAX'] * distrib[economy])
            outputDataStructure[economy]['tax']['xaxis'] = yearRange
            outputDataStructure[economy]['hoursWorker']['data'] = np.sum(households[economy]['LABOR'] * distrib[economy]) / np.sum(households[economy]['WORKING'] * distrib[economy])
            outputDataStructure[economy]['hoursWorker']['xaxis'] = yearRange
            outputDataStructure[economy]['employment']['data'] = np.sum(households[economy]['WORKING'] * distrib[economy])
            outputDataStructure[economy]['employment']['xaxis'] = yearRange
            
            # Consumption and savings for people who have no social
            # security income and are past the maximum working age

            outputDataStructure[economy]['assets_noss']['data'] = np.sum(k * distrib[economy] * nossLocation * retiredLocation) / np.sum(distrib[economy] * nossLocation * retiredLocation) / scenario.modelunit_dollar
            outputDataStructure[economy]['assets_noss']['xaxis'] = yearRange
            outputDataStructure[economy]['cons_noss']['data'] = np.sum(households[economy]['CONSUMPTION'] * distrib[economy] * nossLocation * retiredLocation) / np.sum(distrib[economy] * nossLocation * retiredLocation) / scenario.modelunit_dollar
            outputDataStructure[economy]['cons_noss']['xaxis'] = yearRange
            
            # Consumption, savings, and benefits for people past the
            # maximum working age
            
            outputDataStructure[economy]['assets_retired']['data'] = np.sum(k * distrib[economy] * retiredLocation) / sum(distrib[economy] * retiredLocation) / scenario.modelunit_dollar
            outputDataStructure[economy]['assets_retired']['xaxis'] = yearRange
            outputDataStructure[economy]['cons_retired']['data'] = np.sum(households[economy]['CONSUMPTION'] * distrib[economy] * retiredLocation) / np.sum(distrib[economy] * retiredLocation) / scenario.modelunit_dollar
            outputDataStructure[economy]['cons_retired']['xaxis'] = yearRange
            outputDataStructure[economy]['bens_retired']['data'] = np.sum(households[economy]['OASI_BENEFITS'] * distrib[economy] * retiredLocation) / np.sum(distrib[economy] * retiredLocation) / scenario.modelunit_dollar
            outputDataStructure[economy]['bens_retired']['xaxis'] = yearRange

        # This loop goes through all of the variables, looks to see whether
        # or not they are functions of the distributions, and constucturs
        # the variables that are stored in the "aggregate" structure but
        # created from the distributions of the household variables
        
        for p in range(max(vars.shape)):
            if vars[p]['dist']:                      
                for economy in experiments:
                    outputDataStructure[economy][vars[p]['variable']]['data'] = GenerateReports.getPercentiles(households[economy][vars[p]['firstVar']], households[economy][vars[p]['secondVar']], distrib[economy], T_model, vars[p]['percentTop'], vars[p]['percentBottom'], vars[p]['threshold'])
                    outputDataStructure[economy][vars[p]['variable']]['xaxis'] = yearRange 

        # This loop goes through the economies and constructs a number of
        # moments that we might find interesting. The variables are saved in
        # the "aggregates" variable.
        
        for economy in experiments:
            
            # Actual data for the distribution of assets loaded from the
            # shared drive
                  
            file = pathFinder.getMicrosimInputPath('SIM_NetPersonalWealth_distribution')

            tempRead = pd.read_csv(file)
            tempRead.append([100, None, 1])     # Append last point for graph
            for q in range(T_model):
                outputDataStructure[economy]['percentileList']['data'][:,q] = np.array(tempRead.iloc[:,0])
                outputDataStructure[economy]['percentileList']['xaxis'][:,q] = np.array(tempRead.iloc[:,0])
                outputDataStructure[economy]['a_distdata_thresh']['data'][:,q] = np.array(tempRead.iloc[:,1])
                outputDataStructure[economy]['a_distdata_thresh']['xaxis'][:,q] = np.array(tempRead.iloc[:,1])
                outputDataStructure[economy]['a_distdata_cumul']['data'][:,q] = np.array(tempRead.iloc[:,2])      
                outputDataStructure[economy]['a_distdata_cumul']['xaxis'][:,q] = np.array(tempRead.iloc[:,0])  
                    
            # Actual data for the distribution of labor income loaded from
            # the shared drive           
            
            file = pathFinder.getMicrosimInputPath('SIM_PreTaxLaborInc_distribution')
            tempRead = pd.read_csv(file)
            tempRead.append([100, None, 1])    # Append last point for graph
            for q in range(T_model):
                outputDataStructure[economy]['l_distdata_thresh']['data'][:,q] = np.array(tempRead.iloc[:,1])
                outputDataStructure[economy]['l_distdata_thresh']['xaxis'][:,q] = np.array(tempRead.iloc[:,0])
                outputDataStructure[economy]['l_distdata_cumul']['data'][:,q] = np.array(tempRead.iloc[:,2]) 
                outputDataStructure[economy]['l_distdata_cumul']['xaxis'][:,q] = np.array(tempRead.iloc[:,0])

            # simulated gini coefficients for labor and assets
            outputDataStructure[economy]['a_ginidata']['data'][0,0:T_model] = 0.857                         # Number from SIM
            outputDataStructure[economy]['a_ginidata']['xaxis'][0,0:T_model] = 1                             # Number from SIM
            outputDataStructure[economy]['l_ginidata']['data'][0,0:T_model] = 0.4858                        # Number from SIM
            outputDataStructure[economy]['l_ginidata']['xaxis'][0,0:T_model] = 1                             # Number from SIM
            
            # Distribution of assets for age 65 people
            
            ageIndicator = (age == (65 - T['realage_entry']))
            if 'capDensity' in variableList:
                for t in range(T_model):
                    capDensity = np.zeros(nk)
                    totDensity = 0
                    for q in range(nk):
                        capIndicator = (k == kv[q])
                        capDensity[q] = sum(distrib[economy][:,t] * capIndicator[:,t])   # .* ageIndicator(:,t) .* capIndicator(:,t) );
                        totDensity = totDensity + capDensity[q]
                    capDensity = capDensity / totDensity
                    outputDataStructure[economy]['capDensity']['data'][:,t]= capDensity
                    outputDataStructure[economy]['capDensity']['xaxis'][:,t] = kv / scenario.modelunit_dollar
            
            # Distribution of consumption for retirees 
            
            if 'consDensity' in variableList:
                retCons = households[economy]['CONSUMPTION'] * retiredLocation /(scenario.modelunit_dollar)
                retDist = distrib[economy] * retiredLocation
                numObs = max(retDist.shape)
                retDist = retDist / (np.ones(numObs) * np.sum(retDist, axis=0))
                for t in range(T_model):
                    outputDataStructure[economy]['consDensity']['xaxis'][:,t] = np.arange(10,100,10)
                for q in range(max(outputDataStructure[economy]['consDensity']['xaxis'][:,t].shape)):
                    outputDataStructure[economy]['consDensity']['data'][q,:] = GenerateReports.getPercentiles(retCons, retCons, retDist, T_model, outputDataStructure[economy]['consDensity']['xaxis'][q,0], 0, 0)
            
            if 'c_lorenz' in variableList:  
                retCons = households[economy]['CONSUMPTION'] * retiredLocation /(scenario.modelunit_dollar)
                retDist = distrib[economy] * retiredLocation
                numObs = max(retDist.shape)
                retDist = retDist / (np.ones(numObs) * np.sum(retDist, axis=0))
                (outputDataStructure[economy]['c_gini']['data'], outputDataStructure[economy]['c_lorenz']['xaxis'], outputDataStructure[economy]['c_lorenz']['data']) = GenerateReports.gini(retDist, retCons, T_model)
                outputDataStructure[economy]['c_45']['xaxis'] = outputDataStructure['economy']['c_lorenz']['xaxis']
                outputDataStructure[economy]['c_45']['data'] = outputDataStructure[economy]['c_lorenz']['xaxis']
            
            # Average consumption by age
            
            if 'cons_age' in variableList:                           # Only run this code if the variable is needed because it's slow
                (outputDataStructure[economy]['cons_age']['data'], outputDataStructure[economy]['cons_age']['xaxis']) = GenerateReports.getConditionalMoment(households[economy]['CONSUMPTION'], age, distrib[economy],1/(scenario.modelunit_dollar))
            
            # Average labor by age
            
            if 'labor_age' in variableList:                           # Only run this code if the variable is needed because it's slow
                (outputDataStructure[economy]['labor_age']['data'], outputDataStructure[economy]['labor_age']['xaxis']) = GenerateReports.getConditionalMoment(households[economy]['LABOR'], age, distrib[economy],1)
            
            # Average assets by age
            
            if 'assets_age' in variableList:                         # Only run this code if the variable is needed because it's slow
                (outputDataStructure[economy]['assets_age']['data'], outputDataStructure[economy]['assets_age']['xaxis']) = GenerateReports.getConditionalMoment(households[economy]['ASSETS'], age, distrib[economy],1/(scenario.modelunit_dollar))

            # Average taxes paid by age
            
            if 'tax_age' in variableList:                            # Only run this code if the variable is needed because it's slow
                (outputDataStructure[economy]['tax_age']['data'], outputDataStructure[economy]['tax_age']['xaxis']) = GenerateReports.getConditionalMoment(households[economy]['TAX'], age, distrib[economy],1/(scenario.modelunit_dollar))
            
            # Distribution of assets
            
            if 'dist_k' in variableList:                             # Only run this code if the variable is needed because it's slow
                for ik in range(nk):
                    outputDataStructure[economy]['dist_k']['data'][ik,:] = np.sum(distrib[economy]) * (kv[ik] == k)
                outputDataStructure[economy]['dist_k']['xaxis'] = math.round(kv/scenario.modelunit_dollar)

            # Consumption by distribution of assets
            
            if 'cons_wealth' in variableList:                           # Only run this code if the variable is needed because it's slow
                (outputDataStructure[economy]['cons_wealth']['data'], outputDataStructure[economy]['cons_wealth']['xaxis']) = GenerateReports.getConditionalMoment(households[economy]['CONSUMPTION'], k, distrib['economy'],1/(scenario.modelunit_dollar))
                outputDataStructure[economy]['cons_wealth']['xaxis'] = outputDataStructure[economy]['cons_wealth']['xaxis'] / scenario.modelunit_dollar
            
            # Labor hours by distribution of assets
            
            if 'labor_wealth' in variableList:                           # Only run this code if the variable is needed because it's slow
                (outputDataStructure[economy]['labor_wealth']['data'], outputDataStructure[economy]['labor_wealth']['xaxis']) = GenerateReports.getConditionalMoment(households[economy]['LABOR'], k, distrib[economy],1)
                outputDataStructure[economy]['labor_wealth']['xaxis'] = outputDataStructure[economy]['labor_wealth']['xaxis'] / scenario.modelunit_dollar
            
            # Creating the gini coefficients and lorenz curves for assets
            # and labor income.
            
            if 'l_lorenz' in variableList or 'l_gini' in variableList:  
                (outputDataStructure[economy]['l_gini']['data'], outputDataStructure[economy]['l_lorenz']['xaxis'], outputDataStructure[economy]['l_lorenz']['data']) = GenerateReports.gini(distrib[economy], households[economy]['LABOR_INCOME'], T_model)
                outputDataStructure[economy]['l_gini']['xaxis'] = yearRange
                outputDataStructure[economy]['l_45']['data'] = outputDataStructure[economy]['l_lorenz']['xaxis']
                outputDataStructure[economy]['l_45']['xaxis'] = outputDataStructure[economy]['l_lorenz']['xaxis']
            if 'a_lorenz' in variableList or 'a_gini' in variableList:  
                (outputDataStructure[economy]['a_gini']['data'], outputDataStructure[economy]['a_lorenz']['xaxis'], outputDataStructure[economy]['a_lorenz']['data']) = GenerateReports.gini(distrib[economy], households[economy]['ASSETS'], T_model)
                outputDataStructure[economy]['a_gini']['xaxis'] = yearRange
                outputDataStructure[economy]['a_45']['data'] = outputDataStructure[economy]['a_lorenz']['xaxis']
                outputDataStructure[economy]['a_45']['xaxis'] = outputDataStructure[economy]['a_lorenz']['xaxis']
            
            # The tables of common percentiles
            
            if 'thresh' in variableList or 'cumul' in variableList:  
                for q in range(max(outputDataStructure[economy]['percentileList']['data'][:,0].shape)):
                    outputDataStructure[economy]['a_distmodel_thresh']['data'][q,:] = GenerateReports.getPercentiles(households[economy]['ASSETS']/scenario.modelunit_dollar, households[economy]['ASSETS']/scenario.modelunit_dollar, distrib[economy], T_model, outputDataStructure[economy]['percentileList']['data'][q,0], 0, 1)
                    outputDataStructure[economy]['a_distmodel_cumul']['data'][q,:] = GenerateReports.getPercentiles(households[economy]['ASSETS']/scenario.modelunit_dollar, households[economy]['ASSETS']/scenario.modelunit_dollar, distrib[economy], T_model, outputDataStructure[economy]['percentileList']['data'][q,0], 0, 0)
                    outputDataStructure[economy]['l_distmodel_thresh']['data'][q,:] = GenerateReports.getPercentiles(households[economy]['LABOR_INCOME']/scenario.modelunit_dollar, households[economy]['LABOR_INCOME']/scenario.modelunit_dollar, distrib[economy], T_model, outputDataStructure[economy]['percentileList']['data'][q,0], 0, 1)
                    outputDataStructure[economy]['l_distmodel_cumul']['data'][q,:] = GenerateReports.getPercentiles(households[economy]['LABOR_INCOME']/scenario.modelunit_dollar, households[economy]['LABOR_INCOME']/scenario.modelunit_dollar, distrib[economy], T_model, outputDataStructure[economy]['percentileList']['data'][q,0], 0, 0)
                for q in range(T_model):
                    outputDataStructure[economy]['a_distmodel_cumul']['data'][:,q] = outputDataStructure[economy]['a_distmodel_cumul']['data'][:,q] / outputDataStructure[economy]['a_distmodel_cumul']['data'][max(outputDataStructure[economy]['percentileList']['data'][:,0].shape),q]
                    outputDataStructure[economy]['l_distmodel_cumul']['data'][:,q] = outputDataStructure[economy]['l_distmodel_cumul']['data'][:,q] / outputDataStructure[economy]['l_distmodel_cumul']['data'][max(outputDataStructure[economy]['percentileList']['data'][:,0].shape),q]
                    outputDataStructure[economy]['a_distmodel_thresh']['xaxis'][:,q] = outputDataStructure[economy]['percentileList']['data'][:,0]
                    outputDataStructure[economy]['a_distmodel_cumul']['xaxis'][:,q] = outputDataStructure[economy]['percentileList']['data'][:,0]
                    outputDataStructure[economy]['l_distmodel_thresh']['xaxis'][:,q] = outputDataStructure[economy]['percentileList']['data'][:,0]
                    outputDataStructure[economy]['l_distmodel_cumul']['xaxis'][:,q] = outputDataStructure[economy]['percentileList']['data'][:,0]
            
            # Creating average variables for cohorts if we want them. These
            # variables are defined in the reportStructure structure.
            # Otherwise we don't create them.
            
            for q in range(max(vars.shape)):
                if 'cohort' in vars[q].keys():    
                    if vars[q]['cohort'] != 'n':
                        varName = vars[q]['variable']
                        cohortNum = vars[q]['cohort']
                        if cohortNum < 0:
                            cohortSuffix = 'm' + str(cohortNum)
                        else:
                            cohortSuffix = str(cohortNum)
                        newVarName = varName + '_' + vars[q]['relation'] + '_cohort' + cohortSuffix
                        assert varName in households[economy].keys(), varName + ' is not a field in household decisions (OPTs).'
                        vars[q]['variableCohort'] = newVarName
                        outputDataStructure[economy][newVarName]['data'] = GenerateReports.computeCohortAverage(households[economy]['varName'], distrib[economy], age, cohortNum, T, 1/(scenario.modelunit_dollar))
                        outputDataStructure[economy][newVarName]['xaxis'] = yearRange
        
        return (outputDataStructure, vars)
    
    @staticmethod
    def saveTable(outputTable, reportTemplate):
        
        # This takes the output table and saves it to the specified Excel or
        # comma-separated file.

        outputFile = os.path.join(os.getcwd(), 'Output', reportTemplate['filename'] + '.' + reportTemplate['outputFormat'])
        
        # Export to Excel
        if reportTemplate['outputFormat'] == 'xls':
            with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    
                    # If the file exists, we create a garbage sheet called
                    # TEMPSHEET. TEMPSHEET is deleted at the end. This is because
                    # we want to delete the old output sheet of the same name, but
                    # Excel won't do that if there's no other sheet left in the
                    # workbook.
                    if os.path.exists(outputFile):
                        
                        #use pandas for exporting table to excel
                        writer = pd.ExcelWriter(outputFile)
                        outputTable.to_excel(writer, 'TEMPSHEET')
                        writer.save()
                        
                        #use openpyxl to remove sheets
                        wb = openpyxl.load_workbook(outputFile)
                        wb.remove_sheet(wb.get_sheet_by_name(reportTemplate['outputSheet']))
                        
                        #use pandas again to add sheet with right name
                        outputTable.to_excel(writer, reportTemplate['outputSheet'])
                        writer.save()
                        
                        #use openpyxl to remove sheets
                        wb.remove_sheet(wb.get_sheet_by_name('TEMPSHEET'))
                        wb.remove_sheet(wb.get_sheet_by_name('Sheet1'))
                        wb.remove_sheet(wb.get_sheet_by_name('Sheet2'))
                        wb.remove_sheet(wb.get_sheet_by_name('Sheet3'))
                    
                    else:
                    
                        #use openpyxl to create excel file
                        wb = openpyxl.Workbook()
                        wb.save(outputFile)
                        
                        #use pandas for exporting table to excel
                        writer = pd.ExcelWriter(outputFile)
                        outputTable.to_excel(writer, reportTemplate['outputSheet'])
                        writer.save()
                        
                        #use openpyxl to remove sheets
                        wb.remove_sheet(wb.get_sheet_by_name('Sheet1'))
                        wb.remove_sheet(wb.get_sheet_by_name('Sheet2'))
                        wb.remove_sheet(wb.get_sheet_by_name('Sheet3'))

        # Export to CSV
        if reportTemplate['outputFormat'] == 'csv':
            if os.path.exists(outputFile):
                shutil.rmtree(outputFile)
            outputTable.to_csv(outputFile)   
    
    @staticmethod
    def createTable(outputDataStructure, T, vars, scale, graph):

        # This function goes through each of the variables and creates a
        # table that will eventually be exported to XLS or CSV. The table
        # is created with the variables in the "vars" structure.  
        # Moments is a boolean; if false, then we're
        # working with variables that have years. Otherwise we're working
        # with the xaxis defined in the moments variable.
        
        numVars = len(vars)
        outputTable = pd.DataFrame()
        yearRange = np.array(range(T['TransitionFirstYear'] - 2, T['TransitionLastYear']-1))
         
        # This loop goves through each variable and "constructs" it to
        # report on the output table. It locates the variable, differences
        # it (in the appropriate manner, if applicable) from the correct
        # counterfactual (static or basline.
        
        for p in range(numVars):
            
            if 'variableCohort' in vars[p].keys():
                vars[p]['variable'] = vars[p]['variableCohort']
                
            if not 'year' in vars[p]:
                try:
                    startindex = np.where((outputDataStructure['base'][vars[p]['variable']]['xaxis'] - max(min(outputDataStructure['base'][vars[p]['variable']]['xaxis']),vars[p]['startYear'])) in [0,1])
                except:
                    startindex = 1
                try:
                    endindex = np.where((outputDataStructure['base'][vars[p]['variable']]['xaxis'] - min(max(outputDataStructure['base'][vars[p]['variable']]['xaxis']),vars[p]['endYear'])) in [0,1])
                except:
                    endindex = max(outputDataStructure['base'][vars[p]['variable']]['xaxis'])   
            else:
                startindex = 1
                endindex = max(outputDataStructure['base'][vars[p]['variable']]['xaxis'][:,1])

            relation = vars[p]['relation']
            
            try:
                yearIndex = np.where(yearRange == vars[p]['year'] * np.ones(yearRange.shape))
            except:
                yearIndex = 0
                
            if 'year' in vars[p] and len(yearIndex) == 0:
                raise Exception('Year ' + str(vars[p]['year']) + ' does not exist in this simulation.')
            
            varEcon = []
            
            if relation[0] == 'b':
                varEcon[0] = ['base']
            elif relation[0] == 's':
                varEcon[0] = ['static']
            elif relation[0] == 'c':
                varEcon[0] = ['counter']
                
            if max(relation.shape) > 1:
                if relation[2] == 'b':
                    varEcon[1] = ['base']
                elif relation[2] == 's':
                    varEcon[1] = ['static']
                elif relation[2] == 'c':
                    varEcon[1] = ['counter']
            
            # Checks to make sure that the variable exists
            assert vars[p]['variable'] in outputDataStructure['base'], 'Error creating table: ' + vars[p]['variable'] + ' is not a recognized variable.'
            
            if yearIndex == 0:
                outputTable['year'] = outputDataStructure[varEcon[0]][vars[p]['variable']]['xaxis'][startindex:endindex]
                if p > 1 and graph:
                    indexName = 'index' + str(p)
                    outputTable[indexName] = outputDataStructure[varEcon[0]][vars[p]['variable']]['xaxis'][startindex:endindex]
                d1 = outputDataStructure[varEcon[0]][vars[p]['variable']]['data'][startindex:endindex]
                if max(relation.shape) > 1: 
                    d2 = outputDataStructure[varEcon[1]][vars[p]['variable']]['data'][startindex:endindex]
            else:
                indexName = 'index' + str(p)
                if p == 1:
                    outputTable['index'] = outputDataStructure[varEcon[0]][vars[p]['variable']]['xaxis'][startindex:endindex,yearIndex]
                else:
                    tmp = outputDataStructure[varEcon[0]][vars[p]['variable']]['xaxis'][startindex:endindex,yearIndex]
                    if not min(outputTable['index'] == tmp) or graph:
                        outputTable[indexName] = tmp
                d1 = outputDataStructure[varEcon[0]][vars[p]['variable']]['data'][startindex:endindex,yearIndex] 
                vars[p]['type'] = 'level'
            
            # If the variable is scaled or uses the dollar deflator, then
            # run this little section
            
            if scale:    
               d1 = d1 * vars[p]['dollar'] * vars[p]['scale']
               if max(relation.shape) > 1: 
                   d2 = d2 * vars[p]['dollar'] * vars[p]['scale']
            
            # This puts the actual values in the table. Note that in
            # per_diff and delta, there is a check to make sure that there
            # are nonzero values for some of the series. That's because
            # some of the open economy series are set to zero in the closed
            # economy and dividing by them will result in NaN. They are
            # instead set to zero.
            
            if vars[p]['type'] == 'level':
                outputTable[vars[p]['description']] = d1
            elif vars[p]['type'] == 'per_diff':
                outputTable[vars[p]['description']] = ((d1 - d2)/d2)*100
            elif vars[p]['type'] == 'lev_diff':
                outputTable[vars[p]['description']] = d1 - d2
            elif vars[p]['type'] == 'delta':
                outputTable[vars[p]['description']] = d1 / d2 

        return outputTable
    
    '''
    @staticmethod
    def plotSeries(outputTable, vars):
        
        # This creates variables that are a function of existing aggregate
        # variables. These variables are economically interesting, but not
        # generated by the main dynamic model.
        
        # If the start or end year is "out-of-bounds", the steady-state or
        # final transition year is used as the start or end year.
    
        a = findall(0,'type','figure');
        lineStyles = {'-','--',':','-.','-*','--o'};
        try numPlots = max(size(a(1).Children(1).String)); catch; numPlots = 0; end
        colorMatrix = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.4940, 0.1840, 0.5560; 0.3010, 0.7450, 0.9330; 0.6350, 0.0780, 0.1840; 0.4660, 0.6740, 0.1880];
        figureTitle = vars{1}.title;
        
        for p = 1:max(size(vars))
            xAxis = outputTable{:, (p-1)*2 + 1 };
            yAxis = outputTable{:, (p)*2 };
            hold on
            
            % Plotting on the left or right axis, as specified in the
            % output structure.
            if (vars{p}.axis == 'l')
                yyaxis left;
                plot(xAxis, yAxis, lineStyles{mod(numPlots+p-1,6)+1}, 'LineWidth',2, 'Color',colorMatrix(mod(numPlots+p-1,6)+1, :));
                legendSuffix = '';
            else
                yyaxis right;
                plot(xAxis, yAxis, lineStyles{mod(numPlots+p-1,6)+1}, 'LineWidth',2, 'Color',colorMatrix(mod(numPlots+p-1,6)+1, :));
                legendSuffix = ' (right axis)';
            end      
            title(figureTitle);
            
            % This looks to see if the legend exists. If it does exist, it
            % appends to the existing legend.  That way, we can add more
            % series to existing figures; if we run the reporting program
            % for multiple scenarios and it "holds" the other scenarios.
            % That way we can compare different scenarios on the same
            % graph.
            try
                tt = findall(gcf, 'Type', 'Legend'); bb = tt.String; bb{find(strncmp(bb, 'data', length(4)))} = [strrep(vars{p}.description,'_',' ') legendSuffix];
            catch
                bb{1} = [strrep(vars{p}.description,'_',' ') legendSuffix];
            end
            legend(bb);
            legend('Location','northwest')
            hold off
        end  % End of the loop for the variables graphed
    end  % End of function plotSeries
    '''
    
    @staticmethod
    def parseVars(reportTemplate, scenario, steadyState):
        
        # Parse vars is a function that goes through the input structure
        # that lays out the structure of the output sheet and parses it for
        # key information. If we
        
        graphs = 0
        scale = 0
        
        tmp_var = np.array([s.strip() for s in reportTemplate['variable'].split(',')])
        tmp_rel = np.array([s.strip() for s in reportTemplate['relation'].split(',')])
        
        validSummary = False
        secondCond = False
        if 'firstCond' in reportTemplate.keys():
            validSummary = True
            relations = ['b', 's', 'c']
            reportTemplate['type'] = 'level'
            assert not 'filename' in reportTemplate.keys(), 'Welfare is not compatible with filename option'
            assert not 'outputFormat' in reportTemplate.keys(), 'Welfare is not compatible with outputFormat option'
            assert not 'outputSheet' in reportTemplate.keys(), 'Welfare is not compatible with outputSheet option'
            assert not 'axis' in reportTemplate.keys(), 'Welfare is not compatible with axis option'
            assert not 'xlabel' in reportTemplate.keys(), 'Welfare is not compatible with xlabel option'
            assert not 'llabel' in reportTemplate.keys(), 'Welfare is not compatible with llabel option'
            assert not 'rlabel' in reportTemplate.keys(), 'Welfare is not compatible with rlabel option'
            assert not 'title' in reportTemplate.keys(), 'Welfare is not compatible with title option'
            assert not 'firstVar' in reportTemplate.keys(), 'Welfare is not compatible with firstVar option'
            assert not 'percentTop' in reportTemplate.keys(), 'Welfare is not compatible with percentTop option'
            assert not 'threshold' in reportTemplate.keys(), 'Welfare is not compatible with threshold option'
            assert not 'scaleFactor' in reportTemplate.keys(), 'Welfare is not compatible with scaleFactor option'
            assert not 'dollar' in reportTemplate.keys(), 'Welfare is not compatible with dollar option'
            assert not 'cohort' in reportTemplate.keys(), 'Welfare is not compatible with cohort option'
            assert not 'startYear' in reportTemplate.keys(), 'Welfare is not compatible with startYear option'
            assert not 'endYear' in reportTemplate.keys(), 'Welfare is not compatible with endYear option'
            
            assert 'firstCondUpper' in reportTemplate.keys(), 'First conditional needs an upper bound.'
            assert 'firstCondLower' in reportTemplate.keys(), 'First conditional needs a lower bound.'
            
            tmp_firstCond = np.array([s.strip() for s in reportTemplate['firstCond'].split(',')])
            tmp_firstCondTop = np.array([s.strip() for s in reportTemplate['firstCondUpper'].split(',')])
            tmp_firstCondBottom = np.array([s.strip() for s in reportTemplate['firstCondLower'].split(',')])
            tmp_year = np.array([s.strip() for s in reportTemplate['year'].split(',')])
            
            assert len(tmp_year) == len(tmp_var), 'Wrong number of years.'
            assert len(tmp_firstCondTop) == len(tmp_var), 'Wrong number of top values for the first conditional.'
            assert len(tmp_firstCondBottom) == len(tmp_var), 'Wrong number of bottom values for the first conditional.'
            assert len(tmp_var) == len(tmp_rel), 'The number of variables does not match the number of relations.'
            assert tmp_var != 'EV' and tmp_rel != 'c', 'For EV, counter is the only available report.'
            assert tmp_var != 'CE' and tmp_rel != 'c', 'For CE, counter is the only available report.'
            
            if 'secondCond' in reportTemplate:
                secondCond = True
                tmp_secondCond = np.array([s.strip() for s in reportTemplate['secondCond'].split(',')])
                assert 'secondCondUpper' in reportTemplate.keys(), 'Second conditional needs an upper bound.'
                assert 'secondCondLower' in reportTemplate.keys(), 'Second conditional needs a lower bound.'
                
                tmp_secondCondTop = np.array([s.strip() for s in reportTemplate['secondCondUpper'].split(',')])
                tmp_secondCondBottom = np.array([s.strip() for s in reportTemplate['secondCondLower'].split(',')])
                assert len(tmp_secondCondTop) == len(tmp_var), 'Wrong number of top values for the second conditional.'
                assert len(tmp_secondCondBottom) == len(tmp_var), 'Wrong number of bottom values for the second conditional.'
        
        # Loading, splitting the comma-delimited strings into an array of
        # variables, with leading and lagging spaces trimmed from the
        # strings. These are the required inputs.
        
        tmp_out = np.array([s.strip() for s in reportTemplate['outputNames'].split(',')])
        tmp_typ = np.array([s.strip() for s in reportTemplate['type'].split(',')])
       
        # These are the acceptable values for each field in the output
        # structure.
        types = np.array(['level', 'per_diff', 'lev_diff', 'delta'])
        relations = np.array(['c_b', 'c_s', 's_b', 'b', 's', 'c'])
        filetypes = np.array(['xls', 'csv'])
        graphaxes = np.array(['l','r'])
        dollar = np.array(['y','n'])
        threshold = dollar
        relationsIfYear = np.array(['b', 's', 'c'])
        
        # If doing distributions or scale factors, we're checking to make
        # sure that the number and names of the fields are correct.
        validDist = ('firstVar' in reportTemplate.keys()) + ('secondVar' in reportTemplate.keys()) + ('percentTop' in reportTemplate.keys()) + ('percentBottom' in reportTemplate.keys())  + ('threshold' in reportTemplate.keys())
        assert (validDist == 0) or (validDist == 5), 'You have the wrong number of valid fields for distributions, you need firstVar, secondVar, percentTop, and percentBottom fields.'
        validGraphs = 'axis' in reportTemplate.keys()
        validScale = ('scaleFactor' in reportTemplate.keys()) + ('dollar' in reportTemplate.keys())
        assert (validScale == 0) or (validScale == 2), 'You have the wrong number of valid fields for scale factors, you need scaleFactor and dollar fields.'
        validYear = 'year' in reportTemplate.keys()
        validCohort = 'cohort' in reportTemplate.keys()
        validStartYear = 'startyear' in reportTemplate.keys()
        validEndYear = 'endyear' in reportTemplate.keys()
        
        # Error checking basic syntax of output structure: needs the
        # correct number of variables, must be a valid output file type,
        # etc.
        if 'filename' in reportTemplate.keys(): 
            assert not reportTemplate['outputFormat'] in filetypes, 'File type is not valid. Must be xls or csv.'
            if reportTemplate['outputFormat'] == 'xls':
                assert 'outputSheet' in reportTemplate.keys(), 'There is no sheetname for the xls file.'

        assert len(tmp_var) == len(tmp_out), 'The number of variables does not match the number of output names.'
        assert len(tmp_var) == len(tmp_typ), 'The number of variables does not match the number of types.'
        assert len(tmp_var) == len(tmp_rel), 'The number of variables does not match the number of relations.'

        # This loads in the contructed variables, assuming that we're
        # constructing densities.If the graphing components are not in the structure,
        # this part is skipped.
        
        if (validDist > 0):
            tmp_firstVar = np.array([s.strip() for s in reportTemplate['firstVar'].split(',')])
            tmp_secondVar = np.array([s.strip() for s in reportTemplate['secondVar'].split(',')])
            tmp_pt = np.array([s.strip() for s in reportTemplate['percentTop'].split(',')])
            tmp_pb = np.array([s.strip() for s in reportTemplate['percentBottom'].split(',')])
            tmp_thresh = np.array([s.strip() for s in reportTemplate['threshold'].split(',')])
            
            assert len(tmp_var) == len(tmp_firstVar), 'The number of variables does not match the number of first variables for densities.'
            assert len(tmp_var) == len(tmp_secondVar), 'The number of variables does not match the number of second variables for densities.'
            assert len(tmp_var) == len(tmp_pt), 'The number of variables does not match the number of percentiles at the top.'
            assert len(tmp_var) == len(tmp_pb), 'The number of variables does not match the number of percentiles at the bottom.'
            assert len(tmp_var) == len(tmp_thresh), 'The number of variables does not match the number of thresholds.'
       
        
        # This loads cohorts, if necessary   
        if validCohort:
            tmp_cohort = np.array([s.strip() for s in reportTemplate['cohort'].split(',')]) 
            assert len(tmp_var) == len(tmp_cohort), 'The number of variables does not match the number of cohorts.'
            checkNum = np.zeros(len(tmp_cohort))
            for q in range(len(tmp_cohort)):
               if np.isnan(float(tmp_cohort[q])):
                   checkNum[q] = (tmp_cohort[q] == 'n') + (len(tmp_cohort[q]) == 0)
               else:
                   checkNum[q] = isinstance(float(tmp_cohort[q]), numbers.Number)

               assert min(checkNum) > 0, 'Cohorts can only be integers or the letter n or empty'
        
        if validGraphs > 0:
            tmp_axis = np.array([s.strip() for s in reportTemplate['axis'].split(',')])
            assert len(tmp_var) == len(tmp_axis), 'The number of variables does not match the number of graphs.'
            graphs = 1
        
        # This parses the start and end years by variable
        
        if validStartYear:
            tmp_start = float(reportTemplate['startyear'].strip())
        if validEndYear:
            tmp_end = float(reportTemplate['endyear'].strip())

        if validYear:
            tmp_year = np.array([s.strip() for s in reportTemplate['year'].split(',')])
            assert max(tmp_var) == max(tmp_year), 'The number of variables does not match the number of years.'
        
        # This loads in the contructed variables, assuming we have a scale factor 
        
        if (validScale > 0):
            tmp_scale = np.array([s.strip() for s in reportTemplate['scaleFactor'].split(',')])
            tmp_dollar = np.array([s.strip() for s in reportTemplate['dollar'].split(',')])
            assert len(tmp_var) == len(tmp_scale), 'The number of variables does not match the number of scale values.'
            assert len(tmp_var) == len(tmp_dollar), 'The number of variables does not match the number of dollar scale values.'
        
        # Checking to make sure that certain mutually exclusive options are
        # not selected
        
        if (validStartYear or validEndYear) and validYear:
            raise Exception('Selected a specific year to report a moment and a start/end year.')
        if (validCohort) and validYear:
            raise Exception('Selected a specific year to report a moment and construction a cohort.')
        if (validDist) and validYear:
            raise Exception('Selected a specific year to report a moment and distribution.')
        if (not validYear) and validGraphs and steadyState:
            raise Exception('Cannot create figure of this series in the steady state: there is only a one value to plot.')

        # This part constructs the "vars" cell matrix from the output
        # structure. The part of the program that makes the tables for
        # output uses the "vars" cell matrix to contruct the tables. This
        # also checks to make sure that a number of the field values are
        # valid, and throws errors if they are not.
        
        '''
        CONTINUE CHECKING FROM HERE
        '''
        
        for p in range(len(tmp_var)):
            vars[p]['variable'] = tmp_var[p]    
            vars[p]['description'] = tmp_out[p].replace(' ','_') 
            assert tmp_typ[p] in types, vars[p]['variable'] + ' type ' + tmp_typ[p] + ' is not valid.'
            vars[p]['type'] = tmp_typ[p]
            assert tmp_rel[p] in relations, vars[p]['variable'] + ' relation ' + tmp_rel[p] + ' is not valid.'
            assert not ((tmp_typ[p] == 'level') and (len(tmp_rel[p])>1)), 'Level output is not compatible with the type of relation.'
            if validYear:
                assert tmp_rel[p] in relationsIfYear, vars[p]['variable'] + ' relation ' + tmp_rel[p] + ' is not valid for this type of variable.'
            vars[p]['relation'] = tmp_rel[p] 
            vars[p]['dist'] = False
            vars[p]['startYear'] = -10000
            vars[p]['endYear'] = 10000
            if 'startyear' in reportTemplate.keys():
                vars[p]['startYear'] = tmp_start
            if 'endyear' in reportTemplate.keys():
                vars[p]['endYear'] = tmp_end
            if validYear:
                vars[p]['year'] = float(tmp_year[p])
            assert vars[p]['endYear'] >= vars[p]['startYear'], float(vars[p].endYear) + ' (end year) is before ' + float(vars[p].startYear) + ' (start year). Check the start date and end date.'
           
            if (validDist > 0):  
                if not np.isnan(float(tmp_pt[p])):
                    vars[p]['dist'] = True
                    vars[p]['firstVar'] = tmp_firstVar[p]
                    vars[p]['secondVar'] = tmp_secondVar[p]
                    vars[p]['percentTop'] = float(tmp_pt[p]) 
                    assert not np.isnan(float(tmp_pt[p])), vars[p]['variable'] + ' percentile bottom ' + tmp_pb[p] + ' is not valid.'
                    vars[p]['percentBottom'] = float(tmp_pb[p])
                    assert (vars[p]['percentTop'] > 0) and (vars[p]['percentBottom'] >= 0) and (vars[p]['percentTop'] <=100) and (vars[p]['percentBottom'] < 100), vars[p]['variable'] + ' has a percentile that is out of bounds. (Note: top must be greater than 0.)'
                    assert vars[p]['percentTop'] > vars[p]['percentBottom'], vars[p]['variable'] + ' bottom percentile is bigger than top percentile.'
                    assert tmp_thresh[p] in threshold, tmp_thresh[p] + ' threshold is not equal to y or n.'
                    vars[p]['threshold'] = (tmp_thresh[p] == 'y')
            
            if (validGraphs > 0):
                
                vars[p]['axis']  = tmp_axis[p]
                assert tmp_axis[p] in graphaxes, vars[p]['axis'] + ' axis side ' + tmp_typ[p] + ' is not valid.'
               
                try:
                    vars[p]['xlabel'] = reportTemplate['xlabel']
                except:
                    vars[p]['xlabel'] = ''
                try:
                    vars[p]['llabel'] = reportTemplate['llabel']
                except:
                    vars[p]['llabel'] = ''
                try:
                    vars[p]['rlabel'] = reportTemplate['rlabel']
                except:
                    vars[p]['rlabel'] = ''
                try:
                    vars[p]['title'] = reportTemplate['title']
                except:
                    vars[p]['title'] = '' 
                
                if not validYear:
                    vars[p]['xlabel'] = "Years"
            
            if (validScale > 0):                
                scale = True
                vars[p]['scale']  = float(tmp_scale[p])
                if len(vars[p]['scale']):
                    vars[p]['scale'] = 1
                assert tmp_dollar[p] in dollar, tmp_dollar[p] + ' scale is not equal to y or n.'
                vars[p]['dollar'] = 1
                if (tmp_dollar[p] == 'y'):
                    vars[p]['dollar'] = 1/scenario['modelunit_dollar']
            
            vars[p]['summary'] = 0
            if validSummary > 0:   
                vars[p]['summary'] = 1
                vars[p]['firstcond'] = tmp_firstCond[p]
                vars[p]['firstcondtop'] = float(tmp_firstCondTop[p])
                vars[p]['firstcondbottom'] = float(tmp_firstCondBottom[p])
                assert vars[p]['type'] == 'level', 'Summary statistics must be in levels.'
                assert vars[p]['firstcondtop'] >= vars[p]['firstcondbottom'], 'Lower bound for first conditional is larger than upper bound.'
               
                if 'secondCond' in reportTemplate.keys():
                    vars[p]['secondcond'] = tmp_secondCond[p]
                    vars[p]['secondcondtop'] = float(tmp_secondCondTop[p])
                    vars[p]['secondcondbottom'] = float(tmp_secondCondBottom[p])
                    assert vars[p]['secondcondtop'] >= vars[p]['secondcondbottom'], 'Lower bound for second conditional is larger than upper bound.'
            
            vars[p]['cohort'] = 'n'
            if 'cohort' in reportTemplate.keys():
               vars[p]['cohort'] = tmp_cohort[p] 
        
        # Extracting the experiments we'll need to open. Always pulls the
        # baseline experiment, but pulls static and counter only if
        # requested by the user.
        
        experiments = ['base']
        relTest = ['c_s', 's_b', 's']
        if (relTest[0] in tmp_rel) or (relTest[1] in tmp_rel) or (relTest[2] in tmp_rel):
            experiments.append('static')
        relTest = ['c_b', 'c_s', 'c']
        if (relTest[0] in tmp_rel) or (relTest[1] in tmp_rel) or (relTest[2] in tmp_rel):
            experiments.append('counter')

        return (vars, graphs, scale, experiments)
  
    @staticmethod
    def deNestStructure(orig_structure):  
        
        # If we're loading a variable from a ".mat" file, it loads it as a
        # structure within a structure, and we have to "unpack" it. This
        # isn't a problem if we're using the scenario that's defined in a
        # place like model_tester or the batch processor. We can fix this
        # by changing how the "scenario.mat" is saved. We can just save its
        # elements directly into the file, instead of the whole "scenario"
        # structure.
        
        field_names = orig_structure.keys()
        if len(field_names) == 1:
            output_structure = orig_structure[field_names[0]]
        else:
            output_structure = orig_structure
        
        return output_structure
    
    @staticmethod
    def getOutputStructure(outfile):  
        # Loading in the file that tells us what the output file is going
        # to look like.  It has to be a '.mat' file.
        
        # This part is if we're loading it as a '.mat' file
        if '.mat' in outfile:
            outputStructure = sio.loadmat(outfile)
         
        return outputStructure
  
    @staticmethod
    def openFiles(scenario, steadyState, experiments):
        # This is the function that loads the data associated with the
        # specified scenarios
        
        # If we're reporting steady state values 
        
        # Get the directories associated with either the open or closed
        # economies
        steadyDir = PathFinder.getCacheDir(scenario.currentPolicy().steady())
        baseDir = PathFinder.getCacheDir(scenario.currentPolicy())
        simDir = PathFinder.getCacheDir(scenario)

        T = ParamGenerator.timing(scenario)
        if steadyState:
            T['TransitionLastYear'] = T['TransitionFirstYear']
            T['T_model'] = 0
        
        decisions = {}
        dist = {}
        outputDataStructure = {}
        markets = {}
        
        # Load all of the aggregate values
        outputDataStructure['steady'] = sio.loadmat(os.join.path(steadyDir, 'dynamics.mat'))
        if (not steadyState) and ('base' in experiments):
            outputDataStructure['base']     = sio.loadmat(os.join.path(baseDir,   'dynamics.mat'))
        if (not steadyState) and ('counter' in experiments):
            outputDataStructure['counter']  = sio.loadmat(os.join.path(simDir,    'dynamics.mat')) 
        if (not steadyState) and ('static'in experiments):
            outputDataStructure['static']   = sio.loadmat(os.join.path(simDir,    'statics.mat')) 
        if steadyState:
            if 'base' in experiments:    
                outputDataStructure['base'] = outputDataStructure['steady']
            if 'counter' in experiments:
                outputDataStructure['counter'] = outputDataStructure['steady']
            if 'static' in experiments: 
                outputDataStructure['static'] = outputDataStructure['steady']

        # Load all of the market values
        markets['steady']  = sio.loadmat(os.join.path(steadyDir, 'market.mat'))
        if (not steadyState) and ('base' in experiments):
            markets['base'] = sio.loadmat(os.join.path(baseDir, 'market.mat'))
        if (not steadyState) and ('counter' in experiments): 
            markets['counter'] = sio.loadmat(os.join.path(simDir, 'market.mat')) 
        if (not steadyState) and ('static' in experiments):
            markets['static']  = sio.loadmat(os.join.path(baseDir, 'market.mat'))
        if steadyState:
            if 'base' in experiments:
                markets['base'] = markets['steady']
            if 'counter' in experiments:
                markets['counter'] = markets['steady']
            if 'static' in experiments:
                markets['static'] = markets['steady']

        # Load all of the distributions
        dist['steady'] = sio.loadmat(os.join.path(steadyDir, 'distribution.mat'))
        if (not steadyState) and ('base' in experiments):
            dist['base'] = sio.loadmat(os.join.path(baseDir, 'distribution.mat'))
        if (not steadyState) and ('counter' in experiments):
            dist['counter'] = sio.loadmat(os.join.path(simDir, 'distribution.mat')) 
        if (not steadyState) and ('static' in experiments):
            dist['static'] = sio.loadmat(os.join.path(simDir, 'Static_distribution.mat')) 
        if steadyState:
            if ('base' in experiments):
                dist['base'] = dist['steady']
            if ('counter' in experiments):
                dist['counter'] = dist['steady']
            if ('static' in experiments):
                dist['static'] = dist['steady']

        # Load all of the decisions
        decisions['steady'] = sio.loadmat(os.join.path(steadyDir, 'decisions.mat'))
        if (not steadyState) and ('base' in experiments):
            decisions['base'] = sio.loadmat(os.join.path(baseDir, 'decisions.mat'))
        if (not steadyState) and ('counter' in experiments):
            decisions['counter'] = sio.loadmat(os.join.path(simDir, 'decisions.mat')) 
        if (not steadyState) and ('static' in experiments):
            decisions['static'] = sio.loadmat(os.join.path(simDir, 'Static_decisions.mat')) 
        if steadyState:
            if 'base' in experiments:
                decisions['base'] = decisions['steady']
            if 'counter' in experiments:
                decisions['counter'] = decisions['steady']
            if 'static' in experiments:
                decisions['static'] = decisions['steady']

        return (outputDataStructure, markets, T, dist, decisions)
    
    
    @staticmethod
    def copySteadyState(outputDataStructure, markets, dist, decisions, steadyState, T, experiments):  
    # This function is going to copy all of the steady state values into
    # the first value for the aggregates and markets
        
        yearRange = np.array(range(T['TransitionFirstYear'] - 1, T['TransitionLastYear']))
        include_ss = lambda x, xss: np.hstack(xss, x)
        
        outputDataStructure2 = {}
        
        for p in experiments:
            fieldNames = GenerateReports.dictAggregates
            for o in range(len(fieldNames)):  # copying the aggregate values
                
                if fieldNames[o, 1] in outputDataStructure[p].keys():
                    # If the field exists in both, merge the two fields with steady state going in the first element, if not, set first element to NaN
                    thesize = max(outputDataStructure[p][fieldNames[o,1]].shape)
                    if fieldNames[o, 1] in outputDataStructure['steady']:
                        outputDataStructure2[p][fieldNames[o,1]]['data'] = include_ss(np.reshape(outputDataStructure[p][fieldNames[o,1]], (1,thesize)), outputDataStructure['steady'][fieldNames[o,1]][0])
                    else:
                        outputDataStructure2[p][fieldNames[o,1]]['data'] = include_ss(np.reshape(outputDataStructure[p][fieldNames[o,1]], (1,thesize)), np.isnan(1))

                    outputDataStructure2[p][fieldNames[o,0]]['xaxis'] = yearRange
                    if steadyState:
                        outputDataStructure2[p][fieldNames[o,0]]['data'][1] = []
                else:
                    print('Variable ' + fieldNames[o,1] + ' is not a recognized variable in the ' + p + ' output.\n')
        
        outputDataStructure = outputDataStructure2
        
        for p in experiments:
            fieldNames = GenerateReports.dictMarkets
            for o in range(len(fieldNames)):  # copying the aggregate values
                
                if fieldNames[o, 1] in markets[p]:
                    # If the field exists in both, merge the two fields with steady state going in the first element, if not, set first element to NaN
                    thesize = max(markets[p][fieldNames[o,1]].shape)
                    if fieldNames[o,1] in markets['steady']: 
                        outputDataStructure[p][fieldNames[o,0]]['data'] = include_ss(np.reshape(markets[p][fieldNames[o,1]],(1,thesize)), markets['steady'][fieldNames[o,1]][0])
                    else:
                        outputDataStructure[p][fieldNames[o,0]]['data'] = include_ss(np.reshape(markets[p][fieldNames[o,1]],(1,thesize)), np.nan(1))
                    outputDataStructure[p][fieldNames[o,0]]['xaxis'] = yearRange
                    if steadyState:
                        outputDataStructure[p][fieldNames[o,0]]['data'][1] = []
                else:
                    print('Variable ' + fieldNames[o,1] + ' is not a recognized variable in the ' + p + ' market output.\n')

        if 'static' in experiments:
            if not steadyState: 
                dist['static']['DIST'] = dist['static']['Static_DIST']
            else:
                dist['static']['DIST'] = dist['static']['DIST'] 

        decisions['steady'] = decisions['steady']['OPTs']
        if 'base' in experiments:
            decisions['base'] = decisions['base']['OPTs']    
        if 'counter' in experiments:
            decisions['counter'] = decisions['counter']['OPTs']
        if 'static' in experiments and steadyState:
            decisions['static'] = decisions['steady']

        decisions2 = {}
        
        # working with the distributions and not the aggregate variables
        for p in experiments:
            dimsDist = max(dist[p]['DIST'].shape)
            if steadyState:
                dist[p]['DIST'] = dist['steady']['DIST']
            else:
                dist[p]['DIST'] = np.concatenate(dist['steady']['DIST'], dist[p]['DIST'], axis = dimsDist - 1) 
            
            fieldNames = GenerateReports.dictDecisions
            for o in range(len(fieldNames)):  # copying the aggregate values

                # If the field exists in both, merge the two fields
                # with steady state going in the first element, if not,
                # it will throw an error
                if fieldNames[o,1] in decisions[p].keys():
                    dims_decs = max(decisions[p][fieldNames[o,1]].shape)
                    if fieldNames[o,1] in decisions['steady']:  
                        if steadyState:
                            decisions2[p][fieldNames[o,0]] = decisions['steady'][fieldNames[o,1]]
                        else:
                            decisions2[p][fieldNames[o,0]] = np.concatenate(decisions['steady'][fieldNames[o,1]], decisions[p][fieldNames[o,1]], axis = dims_decs)
                    else:
                        raise Exception('Missing variable in distribution: ' + fieldNames[o]) 
                    
                else:
                    print('Variable ' + fieldNames[o,1] + ' is not a recognized variable in the decisions ' + p + ' decisions output.\n')

        decisions = decisions2
        dist.pop('steady')
        if (not steadyState and 'static' in experiments):
            dist['static'].pop('Static_DIST')
          
        return (outputDataStructure, dist, decisions)
  
    @staticmethod
    def computeCohortAverage(x, distrib, age, cohort, T, scale):
        
        # This function computes the average of the variable "x" for cohort
        # "cohort", which the age in steady-state. Negative numbers are
        # allowed (people who haven't been born yet). If the people don't
        # exist, the value returned is 0.
        
        ageVector = [float(cohort)]
        for p in range(1,T['TransitionLastYear'] - T['TransitionFirstYear'] + 1):
            ageVector[p] = ageVector[p-1] + 1
        
        ageMatrix = (np.tile(ageVector, [max(distrib.shape),1]) == age)
        conditionalMoment = np.array([[]])
        conditionalMoment[:,0] = sum(x * distrib * ageMatrix) / sum(distrib * ageMatrix) * scale
        conditionalMoment[np.isnan(conditionalMoment)]=0
        
        return conditionalMoment
    
    @staticmethod
    def getConditionalMoment(x, index, distrib, scale):
        
        # This function gets and returns the values of 
        
        uniqueIndex = np.unique(index)
        conditionalMoment = np.zeros((max(uniqueIndex.shape), min(x.shape)))
        xaxis = conditionalMoment
        
        for p in range(max(uniqueIndex.shape)):
           conditionalMoment[p,:] = np.sum(x * distrib * (index == uniqueIndex[p]), axis = 0) / np.sum(distrib * (index == uniqueIndex[p]), axis = 0) * scale
           xaxis[p,:] = uniqueIndex[p]
        
        return (conditionalMoment, xaxis)
    
    @staticmethod
    def gini(distrib, x, T_model):
        
        # putting zeros up front for this and initializing other variables
        x = np.vstack(np.zeros((1,T_model)), x)
        ginicoeff = np.zeros((1,T_model))
        distrib = np.vstack(np.zeros((1,T_model)), distrib)
        lorenzPop = np.zeros(x.shape)
        lorenz_var = lorenzPop
        
        z = x * distrib
        for t in range(T_model):
           order = np.argsort(x[:,t], axis = 1)
           distrib[:,t] = distrib[order,t]
           z[:,t] = z[order,t]
           distrib[:,t] = np.cumsum(distrib[:,t], axis = 0)
           z[:,t] = np.cumsum(z[:,t])
           lorenzPop[:,t] = distrib[:,t] / distrib[-1,t] 
           lorenz_var[:,t] = z[:,t] / z[-1,t] 
           ginicoeff[t] = 1 - np.sum((lorenz_var[0:,t] + lorenz_var[1:,t]) * np.diff(lorenzPop[:,t]))
        
        return (ginicoeff, lorenzPop, lorenz_var)
    
    @staticmethod
    def getPercentiles(x, y, hhDist, T_model, percentTop, percentBottom, threshold):

        # Inputs:  x       = array with the base variable
        #          y       = array with the variable of interest
        #          Example: y is tax, and x is income, the outcome would be 
        #          the tax paid by quantile of income
        #          hh_dist    = measure of households at each state of x
        #          T_model = number of transition periods
        #          percentTop and percentBottom  = The range for the
        #          percentile desired
        # Outputs: series = the series that they want

        # Initializing variables
        cdf = np.zeros((1,T_model)) 
        sort_x = np.zeros(x[:,1].shape)
        x_t = sort_x
        dist_t = sort_x
        y_t = sort_x
    
        for t in range(T_model):

            # Sort variables
           
            x_t[:,t] = np.sort(x[:,t])
            sort_x[:,t] = np.sort(x[:,t])
            dist_t[:,t] = hhDist[sort_x[:,t],t]
            y_t[:,t] = x[sort_x[:,t],t]
            
            # Find cdf at the designated percentile at the top
            check = np.cumsum(dist_t[:,t]) >= percentTop/100
            i = np.where(check != 0)[0,0]
            if percentTop >= 100:
                i = max(x_t[:,t].shape)  # If percent is 100, then just choose the maximum value directly
            cdfTop = np.sum(y_t[0:i,t] * dist_t[0:i,t], axis = 1)
            thresh = y_t[i,t]
            
            # Find cdf at the designated percentile at the bottom
            check = np.cumsum(dist_t[:,t]) >= percentBottom/100
            i = np.where(check != 0)[0,0]
            if percentBottom >= 100:
                i = max(x_t[:,t].shape)    # If percent is 100, then just choose the maximum value directly
            cdfBottom = np.sum(y_t[0:i,t] * dist_t[0:i,t])
            if percentBottom == 0:
                cdfBottom = 0

            cdf[t] = cdfTop - cdfBottom      
            if threshold:
                cdf[t] = thresh

        return cdf

"""
% These are extraneous functions pulled from the Matlab community, they are
% used to clean up the excel spreadsheets as needed

function RemoveSheet(excelFileName,sheetName)
% RemoveSheet - removes the sheets that are automatically added to excel
% file. 
% When Matlab writes data to a new Excel file, the Excel software
% automatically creates 3 sheets (the names are depended on the user
% languade). This appears even if the user has defined the sheet name to be
% added. 
%
% Usage:
% RemoveSheet123(excelFileName) - remove "sheet1", "sheet2","sheet3" from
% the excel file. excelFileName is a string of the Excel file name.
% RemoveSheet123(excelFileName,sheetName) - enables the user to enter the
% sheet name when the language is other than English.
% sheetName is the default sheet name, without the number.
%
%
%                       Written by Noam Greenboim
%                       www.perigee.co.il
%
    %% check input arguments
    if nargin < 1 || isempty(excelFileName)
        error('Filename must be specified.');
    end

    if ~ischar(excelFileName)
        error('Filename must be a string.');
    end

    try
        excelFileName = validpath(excelFileName);
    catch 
        error('File not found.');
    end

    if ~ischar(sheetName)
        error('Default sheet name must be a string.');
    end

    % Open Excel file.
    objExcel = actxserver('Excel.Application');
    objExcel.DisplayAlerts = false;
    objExcel.Workbooks.Open(excelFileName); % Full path is necessary!
    % Delete sheets.
    try
          objExcel.ActiveWorkbook.Worksheets.Item([sheetName]).Delete;
    end
    % Save, close and clean up.
    objExcel.ActiveWorkbook.Save;
    objExcel.ActiveWorkbook.Close;
    objExcel.Quit;
    objExcel.delete;
end

function filenameOut = validpath(filename)
    % VALIDPATH builds a full path from a partial path specification
    %   FILENAME = VALIDPATH(FILENAME) returns a string vector containing full
    %   path to a file. FILENAME is string vector containing a partial path
    %   ending in a file or directory name. May contain ..\  or ../ or \\. The
    %   current directory (pwd) is prepended to create a full path if
    %   necessary. On UNIX, when the path starts with a tilde, '~', then the
    %   current directory is not prepended.
    %
    %   See also XLSREAD, XLSWRITE, XLSFINFO.
    
    %   Copyright 1984-2012 The MathWorks, Inc.
    
    %First check for wild cards, since that is not supported.
    if strfind(filename, '*') > 0
        error(message('MATLAB:xlsread:Wildcard', filename));
    end
    
    % break partial path in to file path parts.
    [Directory, file, ext] = fileparts(filename);

    if ~isempty(ext)
        filenameOut = getFullName(filename);
    else
        extIn = matlab.io.internal.xlsreadSupportedExtensions;
        for i=1:length(extIn)
            try                                                                %#ok<TRYNC>
                filenameOut = getFullName(fullfile(Directory, [file, extIn{i}]));
                return;
            end
        end
        error(message('MATLAB:xlsread:FileDoesNotExist', filename));    
    end
end

function absolutepath=abspath(partialpath)
    
    % parse partial path into path parts
    [pathname, filename, ext] = fileparts(partialpath);
    % no path qualification is present in partial path; assume parent is pwd, except
    % when path string starts with '~' or is identical to '~'.
    if isempty(pathname) && strncmp('~', partialpath, 1)
        Directory = pwd;
    elseif isempty(regexp(partialpath,'(.:|\\\\)', 'once')) && ...
            ~strncmp('/', partialpath, 1) && ...
            ~strncmp('~', partialpath, 1);
        % path did not start with any of drive name, UNC path or '~'.
        Directory = [pwd,filesep,pathname];
    else
        % path content present in partial path; assume relative to current directory,
        % or absolute.
        Directory = pathname;
    end
    
    % construct absolute filename
    absolutepath = fullfile(Directory,[filename,ext]);
end

function filename = getFullName(filename)
    FileOnPath = which(filename);
    if isempty(FileOnPath)
        % construct full path to source file
        filename = abspath(filename);
        if isempty(dir(filename)) && ~isdir(filename)
            % file does not exist. Terminate importation of file.
            error(message('MATLAB:xlsread:FileDoesNotExist', filename));
        end
    else
        filename = FileOnPath;
    end
end
"""
