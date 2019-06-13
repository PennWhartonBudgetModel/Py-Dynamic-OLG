##
# Social Security Support 
#

class SocialSecurity:
	
	params         # Struct of params from ParamGenerator
	bv             # Grid of discretized average earnings
	priceIndices   # Indices of various prices for this iteration (wages, CPI, etc)
	startyears     # Offsets to t=1 of birthyears of cohorts
	realage_entry  # Actual age of birth in model
	T_model        # time of model
	
	incomeFloor            # wage indexed floor    for accepting earnings to count for benefits
	incomeCeiling          # wage indexed ceiling  for accepting earnings to count for benefits
	
	payrollTaxBrackets     # indexed payroll tax brackets
	payrollTaxRates        # indexed payroll tax rates
	payrollTaxBurdens      # cumulative payroll tax liability
	
	
	# Constructor
	#    Inputs:
	#       params          : struct from ParamGenerator.social_security
	#       bv              : discretized grid of average earnings
	#       priceIndices    : struct of price indicies
	def __init__(params, bv, priceIndices, startyears, realage_entry, T_model):
		
		self.params         = params
		self.bv             = bv
		self.priceIndices   = priceIndices
		self.startyears     = startyears
		self.realage_entry  = realage_entry
		self.T_model        = T_model
		
		# Calculate indexed policy variables, use only model periodindices
		t_1                 = priceIndices['T_modelStart']
		t_end               = priceIndices['T_modelEnd']
		wage_inflations     = priceIndices['wage_inflations'][t_1:t_end]
		self.incomeFloor    = (params['ssincmins'] * wage_inflations)
		self.incomeCeiling  = (params['ssincmaxs'] * wage_inflations)
		
		# Index the payroll tax structures and put in local vars
		self.indexPayrollTax()
		
		return self
	
	# Get time-varying, indexed Social Security benefits by cohort
	#   Inputs:
	#       i   : cohort number, if -1, then "steady" cohort
	def getBenefitsForCohort(self, i):
		
		if i == -1:
			# Note: retire_year = 1 so that ssbenefits is calculated 
			# for T_model = 1 and indexed for first startyear
			startYear   = self.startyears[0]
			retireYear  = max(self.T_model - 1, 1)
			brackets    = self.params['benefitParamsByCohort'][-1]['brackets']
			rates       = self.params['benefitParamsByCohort'][-1]['rates']
		else:
			startYear   = self.startyears[i]
			retireYear  = self.params['retire_years'][i]
			brackets    = self.params['benefitParamsByCohort'][i]['brackets']
			rates       = self.params['benefitParamsByCohort'][i]['rates']
			
		# Fetch index used to grow brackets
		wage_index = self.priceIndices['wage_inflations']  # full long index (has all cohorts)
		
		# Year cohort turns 62 in current model. 
		# Index benefits by wage index of that year.
		year62 = startYear + (62 - self.realage_entry)
		
		# Cohort based benefits by year
		#   REM: Every household is at a grid point, so can
		#        calculate benefits directly here (outside of
		#        solve_cohort).
		#   NOTE: We set ssbenefits to -Inf otherwise -- it should not be
		#   used in solve_cohort (this will blow it up just in case)
		benefits  = np.ones(self.T_model, size(this.bv,1)) * -Inf
		
		#Build cohort-specific bracket cutoffs for each year indexed by
		# wage_index of when it turns 62.
		w_year62    = year62 + self.priceIndices['T_modelStart'] # year in index vector when cohort turns 62
		bracket_idx = wage_index[w_year62 - 1] * np.ones(self.T_model,1)
		adjbrackets = brackets * np.tile(bracket_idx, (1, size(brackets,2)))
		
		# A possible step to implement here would be to deflate the adjusted
		# brackets by another index, e.g. CPI old people,  but currently
		# it's only deflated by CPI, which is already done in ParamGenerator
		
		# Build cumulative benefits matrix to aid benefit calculation
		? adjtotben   = cumsum(diff(adjbrackets, 1, 2)*rates(:, 1:end-1), 2); 
		? adjtotben   = [zeros(size(adjbrackets, 1), 1), adjtotben];  % rem: first column is zero
		
		# Calculate benefits for each year of retirement until end of time.
            for t = max(retireYear,1):this.T_model
                for ib = 1:size(this.bv,1)
                    thebracket       = find(adjbrackets(t,:) <= this.bv(ib), 1, 'last');
                    benefits(t,ib) = adjtotben(t,thebracket)            ...
                                        + rates(t,thebracket)*(this.bv(ib) - adjbrackets(t,thebracket));
                end
            end

        return benefits
  
        
   end % public instance methods
    
    
   methods (Access=private)
   
        %% Calculate SS tax brackets using required indexing
        %    If brackets overlap because of indexing, reorder the brackets 
        function [] = indexPayrollTax( this )  
            
            % For easier typing, write as local vars
            brackets        = this.params.taxbrackets;
            bracketindices  = this.params.taxindices(1,:);  % TBD: Make indices time-varying
            rates           = this.params.taxrates;  
            
            if( any(strcmp(bracketindices, 'cohort_wages')) )
                throw MException('SocialSecurity.indexSSTax:INDEX', 'Cannot use index type <cohort_wages>.' );
            end

            indices    = zeros(size(brackets));
            for i = 1:size(indices, 2)
                indexname = bracketindices{i};
                % TBD: This should be for all long indices
                if( strcmp(indexname, 'wage_inflations' ) )
                    t_1     = this.priceIndices.T_modelStart;
                    t_end   = this.priceIndices.T_modelEnd;
                else
                    t_1     = 1;
                    t_end   = size(indices, 1);
                end
                % Only use model period range of long index
                indices(:,i) = this.priceIndices.(indexname)(t_1:t_end);
            end
            ssbrackets = brackets .* indices;

            % Take initial brackets and see where they are after indexing and
            % sorting. Then, apply rates to new topology with rates applying
            % between moved brackets. Resolve any overlaps by taking higher rate.
            [ssbrackets, new_old_idx]   = sort(ssbrackets, 2); % rem: sort index from sort() is map from new to old location
            old_new_idx                 = zeros(size(new_old_idx));
            n_brackets                  = size(old_new_idx, 2);
            for i = 1:n_brackets          % map sort index from old to new locations 
                old_new_idx(:, new_old_idx(:, i)) = i; 
            end

            ssrates = -Inf * zeros(size(rates));
            for year = 1:size(rates,1) % TBD: Can this be vectorized?
                for r_old = 1:n_brackets
                    r_rate      = rates      (year, r_old);
                    r_new       = old_new_idx(year, r_old);
                    if( r_old == n_brackets ) 
                        r_newtop = n_brackets + 1; % extends all the way up
                    else
                        r_newtop = old_new_idx(year, r_old+1);
                    end
                    for rr = r_new:r_newtop-1
                        ssrates(year, rr)  = max(r_rate, ssrates(year, rr));
                    end
                end
            end

            % Calculate tax burdens
            ssburdens = cumsum(diff(ssbrackets, 1, 2).*ssrates(:, 1:end-1), 2); 
            ssburdens = [zeros(size(ssbrackets, 1), 1), ssburdens];  % rem: first burden is zero

            this.payrollTaxBrackets = ssbrackets;
            this.payrollTaxRates    = ssrates;
            this.payrollTaxBurdens  = ssburdens;
        end % indexPayrollTax

    end % private instance methods
    
end

