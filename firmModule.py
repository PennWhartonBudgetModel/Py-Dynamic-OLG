##
# Firms sector
# -- currently contains only Single Firm production

import numpy as np

from helpers import makeIterable

class Firm:
    
    ## 
    # Calculate leverageCost from B/K ratio target
    # and calibration target interest rate and tax params
    @staticmethod
    def calculateLeverageCost(leverageRatio, shareInterestDeduction, taxRate, interestRate, leverageSensitivity):
        
        if interestRate <= 0:
            raise Exception('Interest rate must be positive for leverage cost calculation.')
        
        nu              = leverageSensitivity
        taxBenefitRate  = shareInterestDeduction * taxRate * interestRate
        
        x = (taxBenefitRate/(1 - shareInterestDeduction * taxRate))
        d = - (1/(nu-1))
        h = leverageRatio * x ** d
        
        return h #calculateLeverageCost static method

    ## Global variables
    priceCapital0 = None
    userCostCapital0 = None
    corpLeverageCost = None                 # nu_corp
    passLeverageCost = None                 # nu_pass
     
    TFP = None                              # A
    capitalShare = None                     # alpha
    depreciationRate = None                 # delta
    riskPremium = None                      # for high/low return
    interestRate = None                     # market return on capital (and debt)
    allowBusinessDebt = None                # whether to allow borrowing
    capitalAdjustmentCost = None            # eta
    leverageSensitivity = None              # leverage cost function h parameter
    corpInitLeverageRatio = None            # stored value of initial leverage ratio for corps
    passInitLeverageRatio = None            # stored value of initial leverage ratio for pass-throughs
                
    shareIncomeCorp = None                  # portion of total capital which is corp (vs. pass-through)
    shareIncomePass = None                  # portion of total capital which is pass-through (vs. corp)
    rateExpensingCorp = None                # share of expensed corporate investment
    rateExpensingPass = None                # share of expensed pass-through investment
    shareTaxbaseCorp = None                 # match taxbase to empirics
    shareTaxbasePass = None                 # match taxbase to empirics
    rateCreditsCorp = None                  # tax credits as rate on corp GDP
    shareInterestDeductionPass = None       # phi_int
    shareInterestDeductionCorp = None       # phi_int
        
    ratePass = None                         # implied rate on pass-through capital income
    rateCorpStatutory  = None               # statutory tax rate on corp profits
    rateOtherExpensesCorp = None            # unmodelled expenses as rate on corp GDP
    rateOtherExpensesPass = None            # unmodelled expenses as rate on pass GDP
    rateOtherDeductionsCorp  = None         # tax deductions as rate on corp GDP
    rateOtherDeductionsPass = None          # tax deductions as rate on pass GDP
              
    rateForeignerCorpIncome  = None         # tax on foreign corp capital income
    rateForeignerPassIncome  = None         # tax on foreign pass capital income
    rateForeignBusinessIncome  = None       # tax on foreign capital income (weighted average of corp & pass)
        
    capital  = None                         # Capital series
    caps_foreign   = None                   # Foreign capital series
    labor = None                            # Efficient labor series
    KLratio  = None                         # Capital to labor ratio, ( desired not capital/labor)
    invtocaps_end  = None                   # Investment/capital at time T
    invtocaps_1  = None                     # Investment/capital at time 1
        
    priceCapital_ = None                    # cached value at time of construction

    ##
    # Constructor
    #     INPUTS:   Aggregate = Dynamic or Static aggregates
    #               Market    = the prices of stuff
    #               paramsTax = ParamGenerator.tax()
    #               paramsProduction = ParamGenerator.production()
    def __init__(self, Aggregate, Market, paramsTax, paramsProduction):
        
        self.TFP                        = paramsProduction['TFP']
        self.capitalShare               = paramsProduction['alpha']
        self.depreciationRate           = paramsProduction['depreciation']
        self.riskPremium                = paramsProduction['risk_premium']
        self.allowBusinessDebt          = paramsProduction['allowBusinessDebt']
        self.leverageSensitivity        = paramsProduction['leverageSensitivity']
        self.capitalAdjustmentCost      = paramsProduction['capitalAdjustmentCost']
        self.corpInitLeverageRatio      = paramsProduction['leverageCorp0']
        self.passInitLeverageRatio      = paramsProduction['leveragePass0']
            
        self.shareIncomeCorp            = paramsTax['shareIncomeCorp']
        self.shareIncomePass            = paramsTax['shareIncomePass']
        self.rateExpensingCorp          = paramsTax['rateExpensingCorp']
        self.rateExpensingPass          = paramsTax['rateExpensingPass']
        self.shareTaxbaseCorp           = paramsTax['shareTaxbaseCorp']
        self.shareTaxbasePass           = paramsTax['shareTaxbasePass']
        self.rateCreditsCorp            = paramsTax['rateCreditsCorp']
        self.shareInterestDeductionPass = paramsTax['shareInterestDeductionPass']
        self.shareInterestDeductionCorp = paramsTax['shareInterestDeductionCorp']

        self.ratePass                   = paramsTax['ratePass']
        self.rateCorpStatutory          = paramsTax['rateCorpStatutory']
        self.rateOtherExpensesCorp      = paramsTax['rateOtherExpensesCorp']
        self.rateOtherExpensesPass      = paramsTax['rateOtherExpensesPass']
        self.rateOtherDeductionsCorp    = paramsTax['rateOtherDeductionsCorp']
        self.rateOtherDeductionsPass    = paramsTax['rateOtherDeductionsPass']
            
        self.rateForeignerCorpIncome    = paramsTax['rateForeignerCorpIncome']
        self.rateForeignerPassIncome    = paramsTax['rateForeignerPassIncome']
        self.rateForeignBusinessIncome  = paramsTax['rateForeignBusinessIncome']
            
        self.capital                    = makeIterable(Aggregate['caps'])
        self.caps_foreign               = makeIterable(Aggregate['caps_foreign'])
        self.labor                      = makeIterable(Aggregate['labeffs'])
        self.KLratio                    = Market['rhos']
        try:
            self.invtocaps_end              = Market['invtocaps'][-1]
            self.invtocaps_1                = Market['invtocaps'][0]
        except:
            self.invtocaps_end              = Market['invtocaps']
            self.invtocaps_1                = Market['invtocaps']
        self.interestRate               = Market['equityDividendRates']
        if not isinstance(self.interestRate, np.ndarray):
            self.interestRate = np.array([self.interestRate])
            
        self.corpLeverageCost           = Market['corpLeverageCost']
        self.passLeverageCost           = Market['passLeverageCost']
        
        # Calculate the price of capital (p_K, see docs)
        # If userCostCapital0 already available use that, otherwise
        # calculate. The calculation should only occur in steady
        # state.
        if 'userCostCapital0' in Market.keys():
            self.userCostCapital0 = Market['userCostCapital0']
        else:
            self.resetPriceCapital() # rem: sets self.userCostCapital0
            
        if 'priceCapital0' in Market.keys():
            self.priceCapital0 = Market['priceCapital0']
        else:
            self.priceCapital0 = 1
            
        # Initialize local instance price of capital
        self.priceCapital_  = self.priceCapital()
        
    ##
    # Recalculate capital prices, user cost
    # This only makes sense for 'steady' economy
    def resetPriceCapital(self):
        
        #overload float/int and unit-length iterables
        try:
            shareIncomeCorp = self.shareIncomeCorp[0]
        except:
            shareIncomeCorp = self.shareIncomeCorp
        try:
            rateExpensingCorp = self.rateExpensingCorp[0]
        except:
            rateExpensingCorp = self.rateExpensingCorp
        try:
            rateCorpStatutory = self.rateCorpStatutory[0]
        except:
            rateCorpStatutory = self.rateCorpStatutory
        try:
            shareIncomePass = self.shareIncomePass[0]
        except:
            shareIncomePass = self.shareIncomePass
        try:
            rateExpensingPass = self.rateExpensingPass[0]
        except:
            rateExpensingPass = self.rateExpensingPass
        try:
            ratePass = self.ratePass[0]
        except:
            ratePass = self.ratePass
        
        self.userCostCapital0 = (shareIncomeCorp * (1 - rateExpensingCorp * rateCorpStatutory) 
                                    + shareIncomePass * (1 - rateExpensingPass * ratePass)          
                                    + self.capitalAdjustmentCost * self.invtocaps_1 )
                       
    ##
    # Recalculate leverage cost
    # This only makes sense for 'steady' economy
    # but, we take t=1 element (just in case)
    def resetLeverageCost(self):
       
        #overload float/int and unit-length iterables
        try:
            shareInterestDeductionCorp = self.shareInterestDeductionCorp[0]
        except:
            shareInterestDeductionCorp = self.shareInterestDeductionCorp
        try:
            rateCorpStatutory = self.rateCorpStatutory[0]
        except:
            rateCorpStatutory = self.rateCorpStatutory
        try:
            interestRate = self.interestRate[0]
        except:
            interestRate = self.interestRate
        try:
            shareInterestDeductionPass = self.shareInterestDeductionPass[0]
        except:
            shareInterestDeductionPass = self.shareInterestDeductionPass
        try:
            ratePass = self.ratePass[0]
        except:
            ratePass = self.ratePass
        try:
            interestRate = self.interestRate[0]
        except:
            interestRate = self.interestRate
        
        
        # Calculate or use leverage cost
        # Rem: the leverage cost is size invariant, so set capital=1
        #     also, capital from initLeverageRatio is in $, so already
        #     scaled
        self.corpLeverageCost = Firm.calculateLeverageCost(
                self.corpInitLeverageRatio,
                shareInterestDeductionCorp,
                rateCorpStatutory,
                interestRate,
                self.leverageSensitivity)
        self.passLeverageCost = Firm.calculateLeverageCost(
                self.passInitLeverageRatio,
                shareInterestDeductionPass,
                ratePass,
                interestRate,
                self.leverageSensitivity)
        
    ##
    # Wage level the firm is willing to pay
    # under competitive market assumptions, this is the MP_L
    def wageRequired(self):
        alpha = self.capitalShare
        wages = self.TFP * (1 - alpha) * (self.KLratio ** alpha)
        return wages
    
    ##
    # MPK -- Just for reporting and indexing, not a received payment
    def MPK(self):
        alpha = self.capitalShare
        r     = self.TFP * alpha * (self.KLratio ** (alpha-1))
        return r
    
    ##
    # Output -- Just for reporting and indexing, not a received payment
    # Used with 'capital' input to solve steady state
    # NOTE: If output is imaginary (i.e. economy has crashed), force
    # it to be zero.
    def output(self, capital = None, labor = None):
        
        if capital is None:
            capital = self.capital
        if labor is None:
            labor = self.labor
            
        capital = makeIterable(capital)
        labor = makeIterable(labor)
        
        alpha = self.capitalShare
        out = self.TFP * (capital ** alpha) * (labor **(1 - alpha)) 
        
        for i in range(len(out)):
            if not np.isreal(out[i]):
                out[i] = 0

        return out
    
    ##
    # Calculate various business payments and taxes
    def distributions(self):
        
        wage = self.wageRequired()
        (divs, cits, debts, corp_aux_vars) = self.corpDividends(self.shareIncomeCorp * self.capital, self.KLratio, wage)
        
        x = {}
        x['corpDividends'] = divs
        x['corpTaxs'] = cits
        x['corpDebts'] = debts
        x['corpDividendRates'] = divs / (self.shareIncomeCorp * self.capital * self.priceCapital())
        
        x.update(corp_aux_vars)
            
        if self.shareIncomePass == 0:
            x['passDividends'] = np.zeros(shape = self.capital.shape)
            x['passTaxIncRates'] = np.zeros(shape = self.capital.shape)
            x['passDebts'] = np.zeros(shape = self.capital.shape)
            x['passDividendRates'] = np.zeros(shape = self.capital.shape)
            x['pass_aux_vars'] = np.zeros(shape = self.capital.shape)
            x['passTaxableIncome'] = np.zeros(shape = self.capital.shape)
            x['passDebtInterest'] = np.zeros(shape = self.capital.shape)
            x['passDepreciation'] = np.zeros(shape = self.capital.shape)
            x['passExpensing'] = np.zeros(shape = self.capital.shape)
            x['passDebtTaxBenefit'] = np.zeros(shape = self.capital.shape)
        else:
            (divs, taxableIncome, debts, pass_aux_vars) = self.passDividends(self.shareIncomePass * self.capital, self.KLratio, wage)
            x['passDividends'] = divs
            x['passTaxIncRates'] = taxableIncome / (self.shareIncomePass * self.capital * self.priceCapital())
            x['passDebts'] = debts
            x['passDividendRates'] = divs / (self.shareIncomePass * self.capital * self.priceCapital())
            x.update(pass_aux_vars)
            
        x['equityDividendRates'] = ((self.shareIncomeCorp * x['corpDividendRates'] +
                            self.shareIncomePass * x['passDividendRates']) / (
                            self.shareIncomeCorp + self.shareIncomePass))
            
        x['corpDividendsForeign'] = x['corpDividendRates'] * (self.shareIncomeCorp * self.caps_foreign)
        x['passDividendsForeign'] = x['passDividendRates'] * (self.shareIncomePass * self.caps_foreign)
        x['corpForeignWithholding'] = np.maximum(np.zeros((x['corpDividendsForeign'] * self.rateForeignerCorpIncome).shape), x['corpDividendsForeign'] * self.rateForeignerCorpIncome)
        x['passForeignWithholding'] = np.maximum(np.zeros((x['passDividendsForeign'] * self.rateForeignerPassIncome).shape), x['passDividendsForeign'] * self.rateForeignerPassIncome)
        return x
    
    ##
    def dividends(self, capital, klRatio, wage):
        
        capital = makeIterable(capital)
        
        (corpDivs, _, _, _) = self.corpDividends(self.shareIncomeCorp * capital, klRatio, wage)
        if self.shareIncomePass == 0:
            passDivs = np.zeros(shape = self.capital.shape)
        else:
            (passDivs, _, _, _) = self.passDividends(self.shareIncomePass * capital, klRatio, wage)
        
        divs = corpDivs + passDivs
        return divs
    
    ##
    def corpDividends(self, capital, klRatio, wage):
        # Inputs : capital
        #          klRatio & wage 
        # Outputs: divs = dividend is the return of the corporation net of tax.
        #          cits = corporate income taxes paid by the firms
        
        capital = makeIterable(capital)
        
        # Labor and capital
        labor = capital / klRatio
        
        # Total revenues
        y = self.TFP * ((capital**self.capitalShare) * (labor**(1-self.capitalShare)))
        
        # Wage payments
        wagesPaid = wage * labor
           
        # Replace depreciated capital (rem: at cost to buy it)
        depreciation = self.depreciationRate * capital * self.priceCapital_
        
        # Risk premium (rem: at cost to buy it)
        risk = self.riskPremium * capital * self.priceCapital_
        
        # Investment
        investment = capital[1:] - (1 - self.depreciationRate)*capital[0:-1]
        investment = np.append(investment, capital[-1] * self.invtocaps_end)

        # Investment expensing 'subsidy'
        expensing = np.maximum(self.rateExpensingCorp * investment * self.priceCapital_, np.zeros((self.rateExpensingCorp * investment * self.priceCapital_).shape))
        #expensing = max(self.rateExpensingCorp * investment * self.priceCapital_, 0)
        
        # Capital adjustment cost = eta/2 * (I/K)**2 * K
        adjustmentCost = (self.capitalAdjustmentCost/2) * (investment * investment) / capital
        
        # Find optimal debt, interest tax benefit, and leverage cost
        iscorp = True
        (debts, debtCost, debtTaxBenefit) = self.calculateDebt(capital, iscorp)
        
        # Combine to get net tax
        taxbase = ((y # Output
              - wagesPaid # Labor costs
              - self.interestRate * debts # Interest paid
              - depreciation # Depreciated K
              - self.rateOtherExpensesCorp * y # Other expenses
              - adjustmentCost # K adjustment cost
              ) * self.shareTaxbaseCorp # Scaling factor to match data
              - self.rateOtherDeductionsCorp * y) # Deductions

        cits = (taxbase * self.rateCorpStatutory # Tax net of other deductions
          - expensing * self.rateCorpStatutory # Investment expensing
          - debtTaxBenefit # Interest deduction
          - y * self.rateCreditsCorp) # Tax credits
                
        # Calculate returns to owners
        divs = (y # revenues
            - wagesPaid # labor costs
            - depreciation # replace depreciated capital
            - risk # discount money lost due to risk
            - adjustmentCost # cap adjustment cost
            - self.rateOtherExpensesCorp * y # Other expenses
            - debtCost # Cost of leverage
            - cits) # net taxes
                  
        # Auxiliary variables
        corp_aux_vars = {}
        corp_aux_vars['corpTaxbase'] = taxbase
        corp_aux_vars['corpCits'] = cits
        corp_aux_vars['corpDebtInterest'] = (self.interestRate * debts)
        corp_aux_vars['corpDepreciation'] = depreciation
        corp_aux_vars['corpExpensing'] = (expensing * self.rateCorpStatutory)
        corp_aux_vars['corpDebtTaxBenefit'] = debtTaxBenefit
        corp_aux_vars['corpDebtCost'] = debtCost
            
        return (divs, cits, debts, corp_aux_vars)
    
    ##
    def passDividends(self, capital, klRatio, wage):
        # Inputs : capital
        #          klRatio & wage 
        # Outputs: divs = dividend is the return of the corporation net of tax.
        #          cits = corporate income taxes paid by the firms
        
        capital = makeIterable(capital)
        
        # Labor and capital
        labor = (capital / klRatio)
        
        # Total revenues
        y = self.TFP * ( (capital**self.capitalShare) * (labor**(1-self.capitalShare)))
        
        # Wage payments
        wagesPaid = wage * labor
        
        # Replace depreciated capital (rem: at cost to buy it)
        depreciation = self.depreciationRate * capital * self.priceCapital_
        
        # Risk premium (rem: at cost to buy it)
        risk = self.riskPremium * capital * self.priceCapital_
        
        # Investment
        investment = capital[1:] - (1 - self.depreciationRate)*capital[0:-1]
        investment = np.append(investment, capital[-1] * self.invtocaps_end)
            
        # Investment expensing 'subsidy'
        #expensing = max(self.rateExpensingPass * investment * self.priceCapital_, 0)
        expensing = np.maximum(self.rateExpensingPass * investment * self.priceCapital_, np.zeros((self.rateExpensingPass * investment * self.priceCapital_).shape))
        
        # Capital adjustment cost = eta/2 * (I/K)**2 * K
        adjustmentCost = (self.capitalAdjustmentCost/2) * (investment * investment) / capital
        
        # Find optimal debt, interest tax benefit, and leverage cost
        iscorp = False
        (debts, debtCost, debtTaxBenefit) = self.calculateDebt(capital, iscorp)
        
        # Calculate returns to owners
        divs  = (y # revenues
            - wagesPaid # labor costs
            - depreciation # replace depreciated capital
            - risk # discount money lost due to risk
            - adjustmentCost # cap adjustment cost
            - debtCost # Cost of leverage
            - self.rateOtherExpensesPass * y) # Other expenses
            
        # Taxable income
        taxableIncome = ((y # revenues
                   - wagesPaid # labor costs
                   - depreciation # replace depreciated capital
                   - risk # discount money lost due to risk
                   - adjustmentCost # cap adjustment cost
                   - self.rateOtherExpensesPass * y # Other expenses
                   - expensing # investment expensing
                   - debtTaxBenefit # interest deduction
                   - self.rateOtherDeductionsPass * y # other deductions
                   ) * self.shareTaxbasePass) # Scaling factor to match data
                        
        # Auxiliary variables
        pass_aux_vars = {}
        pass_aux_vars['passTaxableIncome']      = taxableIncome
        pass_aux_vars['passDebtInterest']       = (self.interestRate * debts)
        pass_aux_vars['passDepreciation']       = depreciation
        pass_aux_vars['passExpensing']          = expensing
        pass_aux_vars['passDebtTaxBenefit']     = debtTaxBenefit
              
        return (divs, taxableIncome, debts, pass_aux_vars)

    ##
    # Calculate capital gains rates from value of capital.
    def capitalGains(self, capital = None):
        if capital is None:
            priceCapital = self.priceCapital_
        else:
            priceCapital = self.priceCapital(capital)
        capgains = np.zeros(shape = priceCapital.shape)
        capgains[0] = (priceCapital[0] - self.priceCapital0)/self.priceCapital0
        for t in range(1, len(capgains)):
            capgains[t] = (priceCapital[t] - priceCapital[t-1])/priceCapital[t-1]
        return capgains
    
    ##
    # Calculate the price of capital
    def priceCapital(self, capital = None):
        if capital is None:
            if self.priceCapital_ != None:
                price = self.priceCapital_
                return price
            capital = self.capital
            
        capital = makeIterable(capital)
            
        # Gross Investment
        investment = capital[1:] - (1 - self.depreciationRate)*capital[0:-1]
        investment = np.append(investment, capital[-1] * self.invtocaps_end)
        invtocaps = investment / capital
        
        userCostCapital = (self.shareIncomeCorp * (1 - self.rateExpensingCorp * self.rateCorpStatutory)
            + self.shareIncomePass * (1 - self.rateExpensingPass * self.ratePass)
            + self.capitalAdjustmentCost * invtocaps)
        
        #limit to the same number of decimals as in Matlab to mantain same final result
        userCostCapital = np.around(userCostCapital, decimals = 15)
        
        # TEMP implies priceCapital0 = 1
        price = userCostCapital/self.userCostCapital0
        return price
    
    ##
    # Calculate K/L ratio from dividend rate
    def calculateKLRatio(self, fixedAfterTaxReturn, init_caps, labor):
        # Inputs : fixedAfterTaxReturn = dollars received as dividend per
        #                        dollar owned of equity
        #          init_caps = capital initial guess
        #          labor = efficient units of labor (from last iteration)
        # Outputs: KLratio that generates the dividendRate of inputs
        
        # Cast
        fixedAfterTaxReturn = makeIterable(fixedAfterTaxReturn)
        init_caps = makeIterable(init_caps)
        labor = makeIterable(labor)
        
        # Find dividends rate needed to fix returns to target
        dividendRate = ((fixedAfterTaxReturn - self.capitalGains() )
            / (1 - self.rateForeignBusinessIncome))

        # Initialize variables
        caps    = init_caps
        divRate = dividendRate
        
        tolerance = 1e-12
        err_div   = float('inf')
        
        while err_div > tolerance:
        
            # Update capital 
            #  if divRate > dividendRate -> caps gets bigger, and divRate gets smaller
            caps[1:] = caps[1:] * ((1+divRate[1:]) / (1+dividendRate[1:]))
            
            # To prevent nonsensical results with non-positive capital.
            # TBD: This may keep loop from converging under some circumstances.
            caps = np.maximum(np.full(caps.shape,1e-5), caps)
            
            K_by_L = caps / labor
            wage   = self.TFP * (1-self.capitalShare) * (K_by_L ** self.capitalShare)
            
            divs = self.dividends(caps, K_by_L, wage)
            divRate = divs / (caps * self.priceCapital(caps))
            
            try:
                err_div = max(abs((divRate[1:] - dividendRate[1:]) / dividendRate[1:]))
            except:
                err_div = -float('inf')
                
            # Update dividends rate since capitalGains change with caps
            dividendRate = (fixedAfterTaxReturn - self.capitalGains(caps)) / (1 - self.rateForeignBusinessIncome)
            
        # Calculate capital-labor ratio
        KLratio  = caps / labor
        return (KLratio, caps)
    
    ##
    # Calculate optimal debt
    # nu() = 1/nu (B/K_s)**nu , where K_s = K h p_K; 
    def calculateDebt(self, capital, isCorp):
        
        # iscorp = Boolean variable indicating type of firm
        # If not using debt, then jump out.
        if not self.allowBusinessDebt: 
            debt        = np.zeros(shape = capital.shape)
            debtCost    = np.zeros(shape = capital.shape)
            taxBenefit  = np.zeros(shape = capital.shape)
            return (debt, debtCost, taxBenefit)
        
        # Calculate optimal debt
        if isCorp:
            taxBenefitRate = self.shareInterestDeductionCorp * self.rateCorpStatutory * self.interestRate
            h = self.corpLeverageCost
            d = 1 - self.shareInterestDeductionCorp * self.rateCorpStatutory
        else:
            taxBenefitRate = self.shareInterestDeductionPass * self.ratePass * self.interestRate
            h = self.passLeverageCost
            d = 1 - self.shareInterestDeductionPass * self.ratePass
            
        nu = self.leverageSensitivity
        x = 1/(nu-1)
        capital_scaled = capital * self.priceCapital_ * h
        debt = ((taxBenefitRate / d ) ** x) * capital_scaled
        debtCost = ((1/nu) * (debt/capital_scaled ) ** nu) * capital_scaled
        taxBenefit = debt * taxBenefitRate
        
        return (debt, debtCost, taxBenefit)
    

