##
# Firms sector
# -- currently contains only Single Firm production

import numpy as np

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
        h = leverageRatio * x ^ d
        
        return h #calculateLeverageCost static method

    ## Global variables
    priceCapital0
    userCostCapital0
    corpLeverageCost                 # nu_corp
    passLeverageCost                 # nu_pass
     
    TFP                              # A
    capitalShare                     # alpha
    depreciationRate                 # delta
    riskPremium                      # for high/low return
    interestRate                     # market return on capital (and debt)
    allowBusinessDebt                # whether to allow borrowing
    capitalAdjustmentCost            # eta
    leverageSensitivity              # leverage cost function h parameter
    corpInitLeverageRatio            # stored value of initial leverage ratio for corps
    passInitLeverageRatio            # stored value of initial leverage ratio for pass-throughs
                
    shareIncomeCorp                  # portion of total capital which is corp (vs. pass-through)
    shareIncomePass                  # portion of total capital which is pass-through (vs. corp)
    rateExpensingCorp                # share of expensed corporate investment
    rateExpensingPass                # share of expensed pass-through investment
    shareTaxbaseCorp                 # match taxbase to empirics
    shareTaxbasePass                 # match taxbase to empirics
    rateCreditsCorp                  # tax credits as rate on corp GDP
    shareInterestDeductionPass       # phi_int
    shareInterestDeductionCorp       # phi_int
        
    ratePass                         # implied rate on pass-through capital income
    rateCorpStatutory                # statutory tax rate on corp profits
    rateOtherExpensesCorp            # unmodelled expenses as rate on corp GDP
    rateOtherExpensesPass            # unmodelled expenses as rate on pass GDP
    rateOtherDeductionsCorp          # tax deductions as rate on corp GDP
    rateOtherDeductionsPass          # tax deductions as rate on pass GDP
              
    rateForeignerCorpIncome          # tax on foreign corp capital income
    rateForeignerPassIncome          # tax on foreign pass capital income
    rateForeignBusinessIncome        # tax on foreign capital income (weighted average of corp & pass)
        
    capital                          # Capital series
    caps_foreign                     # Foreign capital series
    labor                            # Efficient labor series
    KLratio                          # Capital to labor ratio, ( desired not capital/labor)
    invtocaps_end                    # Investment/capital at time T
    invtocaps_1                      # Investment/capital at time 1
        
    priceCapital_                    # cached value at time of construction

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
            
        self.capital                    = Aggregate['caps']
        self.caps_foreign               = Aggregate['caps_foreign']
        self.labor                      = Aggregate['labeffs']
        self.KLratio                    = Market['rhos']
        self.invtocaps_end              = Market['invtocaps'][-1]
        self.invtocaps_1                = Market['invtocaps'][0]
        self.interestRate               = Market['equityDividendRates']
            
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
        
        self.userCostCapital0 = self.shareIncomeCorp * (1 - self.rateExpensingCorp * 
                          self.rateCorpStatutory + self.shareIncomePass * 
                          (1 - self.rateExpensingPass * self.ratePass) + 
                          self.capitalAdjustmentCost * self.invtocaps_1)
                       
    ##
    # Recalculate leverage cost
    # This only makes sense for 'steady' economy
    # but, we take t=1 element (just in case)
    def resetLeverageCost(self):
       
        # Calculate or use leverage cost
        # Rem: the leverage cost is size invariant, so set capital=1
        #     also, capital from initLeverageRatio is in $, so already
        #     scaled
        self.corpLeverageCost = Firm.calculateLeverageCost(
                self.corpInitLeverageRatio, self.shareInterestDeductionCorp,
                self.rateCorpStatutory, self.interestRate, self.leverageSensitivity)
        self.passLeverageCost = Firm.calculateLeverageCost(
                self.passInitLeverageRatio, self.shareInterestDeductionPass,
                self.ratePass, self.interestRate, self.leverageSensitivity)
        
    ##
    # Wage level the firm is willing to pay
    # under competitive market assumptions, this is the MP_L
    def wageRequired(self):
        alpha = self.capitalShare
        wages = self.TFP * (1 - alpha) * (self.KLratio ^ alpha)
        return wages
    
    ##
    # MPK -- Just for reporting and indexing, not a received payment
    def MPK(self):
        alpha = self.capitalShare
        r     = self.TFP * alpha * (self.KLratio ^ (alpha-1))
        return r
    
    ##
    # Output -- Just for reporting and indexing, not a received payment
    # Used with 'capital' input to solve steady state
    # NOTE: If output is imaginary (i.e. economy has crashed), force
    # it to be zero.
    def output(self, capital, labor):
        
        if nargin == 1:
            capital = self.capital
            labor   = self.labor
        
        alpha = self.capitalShare
        out = self.TFP * (capital ^ alpha) * (labor ^(1 - alpha)) 
        for i in range(len(out)):
            if not np.isreal(out(i)):
                out(i) = 0
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
        x['corpDividendRates'] = divs / (self.shareIncomeCorp * self.capital * self.priceCapital)
        
        x.update(corps_aux_vars)
            
        if self.shareIncomePass == 0:
            x['passDividends'] = np.zeros(shape = self.capital.shape())
            x['passTaxIncRates'] = np.zeros(shape = self.capital.shape())
            x['passDebts'] = np.zeros(shape = self.capital.shape())
            x['passDividendRates'] = np.zeros(shape = self.capital.shape())
            x['pass_aux_vars'] = np.zeros(shape = self.capital.shape())
            x['passTaxableIncome'] = np.zeros(shape = self.capital.shape())
            x['passDebtInterest'] = np.zeros(shape = self.capital.shape())
            x['passDepreciation'] = np.zeros(shape = self.capital.shape())
            x['passExpensing'] = np.zeros(shape = self.capital.shape())
            x['passDebtTaxBenefit'] = np.zeros(shape = self.capital.shape())
        else:
            (divs, taxableIncome, debts, pass_aux_vars) = self.passDividends(self.shareIncomePass * self.capital, self.KLratio, wage)
            x['passDividends'] = divs
            x['passTaxIncRates'] = taxableIncome / (self.shareIncomePass * self.capital * self.priceCapital)
            x['passDebts'] = debts
            x['passDividendRates'] = divs / (self.shareIncomePass * self.capital * self.priceCapital)
            x.update(pass_aux_vars)
            
        x['equityDividendRates'] = ((self.shareIncomeCorp * x['corpDividendRates'] +
                            self.shareIncomePass * x['passDividendRates']) / (
                            self.shareIncomeCorp + self.shareIncomePass))
            
        x['corpDividendsForeign'] = x['corpDividendRates'] * (self.shareIncomeCorp * self.caps_foreign)
        x['passDividendsForeign'] = x['passDividendRates'] * (self.shareIncomePass * self.caps_foreign)
        x['corpForeignWithholding'] = max(0, x['corpDividendsForeign'] * self.rateForeignerCorpIncome)
        x['passForeignWithholding'] = max(0, x['passDividendsForeign'] * self.rateForeignerPassIncome)
        return x
    
    ##
    def dividends(self, capital, klRatio, wage):
        
        corpDivs = self.corpDividends(self.shareIncomeCorp * capital, klRatio, wage)
        if self.shareIncomePass == 0:
            passDivs = np.zeros(shape = self.capital.shape())
        else:
            passDivs = self.passDividends(self.shareIncomePass * capital, klRatio, wage)
        
        divs = corpDivs + passDivs
        return divs
    
    ##
    def corpDividends(self, capital, klRatio, wage):
        # Inputs : capital
        #          klRatio & wage 
        # Outputs: divs = dividend is the return of the corporation net of tax.
        #          cits = corporate income taxes paid by the firms
        
        # Labor and capital
        labor = capital / klRatio
        
        # Total revenues
        y = self.TFP * ((capital^self.capitalShare) * (labor^(1-self.capitalShare)))
        
        # Wage payments
        wagesPaid = wage * labor
           
        # Replace depreciated capital (rem: at cost to buy it)
        depreciation = self.depreciationRate * capital * self.priceCapital_
        
        # Risk premium (rem: at cost to buy it)
        risk = self.riskPremium * capital * self.priceCapital_
        
        # Investment
        investment = (capital[1:-1] - (1 - self.depreciationRate)*capital[0:-2],
               capital[-1] * self.invtocaps_end)

        # Investment expensing 'subsidy'
        expensing = max(self.rateExpensingCorp * investment * self.priceCapital_, 0)
        
        # Capital adjustment cost = eta/2 * (I/K)^2 * K
        adjustmentCost = (self.capitalAdjustmentCost/2) * (investment * investment) / capital
        
        # Find optimal debt, interest tax benefit, and leverage cost
        iscorp = true
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
        
        # Labor and capital
        labor = (capital / klRatio)
        
        # Total revenues
        y = self.TFP * ( (capital^self.capitalShare) * (labor^(1-self.capitalShare)))
        
        # Wage payments
        wagesPaid = wage * labor
        
        # Replace depreciated capital (rem: at cost to buy it)
        depreciation = self.depreciationRate * capital * self.priceCapital_;
        
        # Risk premium (rem: at cost to buy it)
        risk = self.riskPremium * capital * self.priceCapital_;
        
        # Investment
        investment = (capital[1:-1] - (1 - self.depreciationRate)*capital[0:-2],
                capital[-1] * self.invtocaps_end)
            
        # Investment expensing 'subsidy'
        expensing = max(self.rateExpensingPass * investment * self.priceCapital_, 0)
        
        # Capital adjustment cost = eta/2 * (I/K)^2 * K
        adjustmentCost = (self.capitalAdjustmentCost/2) * (investment * investment) / capital
        
        # Find optimal debt, interest tax benefit, and leverage cost
        iscorp = false
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
        pass_aux_vars['passTaxableIncome']      = taxableIncome
        pass_aux_vars['passDebtInterest']       = (self.interestRate * debts)
        pass_aux_vars['passDepreciation']       = depreciation
        pass_aux_vars['passExpensing']          = expensing
        pass_aux_vars['passDebtTaxBenefit']     = debtTaxBenefit
              
        return (divs, taxableIncome, debts, pass_aux_vars)

    ##
    # Calculate capital gains rates from value of capital.
    def capitalGains(self, capital):
        if nargin == 1:
            priceCapital = self.priceCapital_
        else:
            priceCapital = self.priceCapital(capital)
        capgains = np.zeros(shape = priceCapital.shape())
        capgains[0] = (priceCapital[0] - self.priceCapital0)/self.priceCapital0
        for t in range(1, len(capgains) + 1):
            capgains[t] = (priceCapital(t) - priceCapital(t-1))/priceCapital(t-1)
        return capgains
    
    ##
    # Calculate the price of capital
    def priceCapital(self, capital):
        if nargin == 1:
            if len(self.priceCapital_) == 0:
                price = self.priceCapital_
                return price
            capital = self.capital
            
        # Gross Investment
        investment = (capital[1:-1] - (1 - self.depreciationRate)*capital[1:-2],
                capital[-1] * self.invtocaps_end)
        invtocaps = investment / capital
        
        userCostCapital = (self.shareIncomeCorp * (1 - self.rateExpensingCorp * self.rateCorpStatutory)
            + self.shareIncomePass * (1 - self.rateExpensingPass * self.ratePass)
            + self.capitalAdjustmentCost * invtocaps)
        
        # TEMP implies priceCapital0 = 1
        price = (1/self.userCostCapital0) * userCostCapital
        return price
    
    ##
    # Calculate K/L ratio from dividend rate
    def calculateKLRatio(self, fixedAfterTaxReturn, init_caps, labor):
        # Inputs : fixedAfterTaxReturn = dollars received as dividend per
        #                        dollar owned of equity
        #          init_caps = capital initial guess
        #          labor = efficient units of labor (from last iteration)
        # Outputs: KLratio that generates the dividendRate of inputs
        
        # Find dividends rate needed to fix returns to target
        dividendRate = ((fixedAfterTaxReturn - self.capitalGains() )
            / (1 - self.rateForeignBusinessIncome))

        # Initialize variables
        caps    = init_caps
        divRate = dividendRate
        
        tolerance = 1e-12
        err_div   = Inf
        
        while err_div > tolerance:
        
            # Update capital 
            #  if divRate > dividendRate -> caps gets bigger, and divRate gets smaller
            caps[1:-1] = caps[1:-1] * ((1+divRate[1:-1]) / (1+dividendRate[1:-1]))
            
            # To prevent nonsensical results with non-positive capital.
            # TBD: This may keep loop from converging under some circumstances.
            caps = max(1e-5, caps)
            
            K_by_L = caps / labor
            wage   = self.TFP * (1-self.capitalShare) * (K_by_L ^ self.capitalShare)
            
            divs = self.dividends(caps, K_by_L, wage)
            divRate = divs / (caps * self.priceCapital(caps))
            
            err_div = max(abs((divRate[1:-1] - dividendRate[1:-1]) / dividendRate[1:-1]))
            
            # Update dividends rate since capitalGains change with caps
            dividendRate = (fixedAfterTaxReturn - self.capitalGains[caps]) / (1 - self.rateForeignBusinessIncome)
            
            # Calculate capital-labor ratio
            KLratio  = caps / labor
        return (KLratio, caps)
    
    ##
    # Calculate optimal debt
    # nu() = 1/nu (B/K_s)^nu , where K_s = K h p_K; 
    def calculateDebt(self, capital, isCorp):
        
        # iscorp = Boolean variable indicating type of firm
        # If not using debt, then jump out.
        if not self.allowBusinessDebt: 
            debt        = np.zeros(shape = capital.shape())
            debtCost    = np.zeros(shape = capital.shape())
            taxBenefit  = np.zeros(shape = capital.shape())
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
        debt = ((taxBenefitRate / d ) ^ x) * capital_scaled
        debtCost = ((1/nu) * (debt/capital_scaled ) ^ nu) * capital_scaled
        taxBenefit = debt * taxBenefitRate
        
        return (debt, debtCost, taxBenefit)
    
    # Test
    def testMe(self):
        
        T_model = 25
        
        # Make sample inputs
        fixedAfterTaxReturn = np.ones(T_model) * 0.05
        foreignTaxRates     = np.ones(T_model) * 0.10
        labor = np.ones(T_model)
        init_caps = np.ones(T_model) * 12
        invtocap_Tmodel = 0.07
        
        (klRatio, tCaps) = self.calculateKLRatio(fixedAfterTaxReturn, foreignTaxRates,
                            init_caps, labor)
        
        tempCaps   = klRatio * labor
        tempI      = (np.diff(tempCaps), invtocap_Tmodel * tempCaps(T_model))
        tempWages  = self.TFP*(1-self.capitalShare)*(klRatio^self.capitalShare)
        
        (tempDividends, _) = self.dividends(tempCaps, tempI, klRatio, tempWages)  
        tempDivRates = (tempDividends / (tempCaps * self.priceCapital))                                
        
        print('diff in dividends')
        tempDivRates - effectiveDividendRate

