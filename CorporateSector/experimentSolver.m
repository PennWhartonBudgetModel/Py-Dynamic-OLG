classdef experimentSolver
    
    methods (Static)
        
        function [] = solve_experiment1()
            %  Dividend tax versus dividend income share: subjects a fraction of dividend income to the personal income tax
            
            
            
            % Identify and create results directories
            [~, save_dir] = identify_dirs(1);
            if exist(save_dir,'dir')
                rmdir(save_dir,'s')
            end
            mkdir(save_dir)
            

            %% Find the level of government expenditures for the given taxes.

            s_tax = load(fullfile('Parameters','tax_parameters.mat'));
            taxes = s_tax.taxes;
            [~, s_quantities] = solve_ss_equilibrium_global([], taxes);
            benchmark_expenditures = s_quantities.government.government_expenditures;
            clear('s_quantities')


            % For a fixed policy instrument, find the tax that clears the government budget constraint.

            %% Solving optimal income vs dividend tax
            div_tax_ub = taxes.hh.dividend;
            div_tax_lb = 0;
            n_divs = 20;
            dividend_taxes = linspace(div_tax_lb,div_tax_ub,n_divs);

            maxiter = 40;

            dividend_income_shares = zeros(1,n_divs);
            welfares = zeros(1,n_divs);

            parfor id = 1:n_divs

                exp_taxes = taxes;
                exp_taxes.hh.dividend = dividend_taxes(id);

                share_ub = 1;
                share_lb = 0;
                iter = 0;
                while true
                    iter = iter+1;
                    % Determine tax revenue
                    exp_taxes.hh.div_inc_share = (share_lb + share_ub)/2;
                    [~, quantities] = solve_ss_equilibrium_global([], exp_taxes);
                    rev = quantities.government.tax_revenue;

                    % Determine error and break if complete
                    error = abs(rev - benchmark_expenditures);
                    if (error<0.001)||(iter>maxiter) , break, end 

                    % Update criteria
                    if rev > benchmark_expenditures
                        share_ub = exp_taxes.hh.div_inc_share;
                    elseif rev<= benchmark_expenditures
                        share_lb = exp_taxes.hh.div_inc_share;
                    end

                end
                dividend_income_shares(id) = exp_taxes.hh.div_inc_share;
                welfares(id) = quantities.hh.welfare;

            end

            dividend_experiment_results.dividend_taxes         = dividend_taxes;
            dividend_experiment_results.dividend_income_shares = dividend_income_shares;
            dividend_experiment_results.welfares               = welfares;


            % Solve optimal solution (using grid-search)
            [~,max_location] = max(welfares);
            taxes.hh.dividend = dividend_taxes(max_location);
            taxes.hh.div_inc_share = dividend_income_shares(max_location);
            [prices, quantities] = solve_ss_equilibrium_global([], taxes);
            dividend_experiment_results.prices     = prices;
            dividend_experiment_results.quantities = quantities; %#ok<STRNU>

%             s_plot = plot(dividend_income_shares,welfares,'LineWidth',4);

            save(fullfile(save_dir,'experiment_results.mat'),'dividend_experiment_results')      





        end
        
        
        
        function [] = solve_experiment2()
            % This experiment solves the optimal dividend and corporate income tax when expensing is different from 1.
            
            
            % Identify and create results directories
            [~, save_dir] = identify_dirs(2);
            if exist(save_dir,'dir')
                rmdir(save_dir,'s')
            end
            mkdir(save_dir)



            %% From the benchmark tax structure, set CIT and expense share, solve equilibrium.
            taxes = paramGenerator.tax;

            [~, s_quantities] = solve_ss_equilibrium_global([], taxes);
            benchmark_expenditures = s_quantities.government.government_expenditures;
            benchmark_welfare = s_quantities.hh.welfare;
            clear('s_quantities')

            %% Find the dividend tax that clears market when CIT declines towards zero.
            
            tolerance = .00001;
            
            % Create grid of CIT's
            cit_ub = taxes.firm.income;
            cit_lb = 0;
            n_cits = 20;
            cits = linspace(cit_lb, cit_ub, n_cits);
            welfares = zeros(1,n_cits);
            dits     = zeros(1,n_cits);
            errors   = zeros(1,n_cits);
            
            % Setting expense share to 1 just to debug.
            taxes.firm.exp_share = 1;
            maxiter = 20;
            parfor ic = 1:n_cits
                experiment_taxes = taxes;
                experiment_taxes.firm.income = cits(ic);
                dub = .5;
                dlb = 0;
                iter = 0;
                while true
                    iter = iter+1;
                    experiment_taxes.hh.dividend = (dub+dlb)/2;
                    [~, quantities] = solve_ss_equilibrium_global([], experiment_taxes);
                    rev = quantities.government.tax_revenue;
                    error = abs(rev - benchmark_expenditures);
                    if (error<tolerance)||(iter>maxiter), break, end
                    
                    % Only print to screen if code is running serially.
                    pool = gcp('nocreate');
                    if ~isempty(pool)
                        fprintf('\nDividend tax upper bound = %0.6f\n', dub)
                        fprintf('\nDividend tax lower bound = %0.6f\n', dlb)
                    end
                    fprintf('\nError = %0.6f\n', error)
                    fprintf('\nIteration = %d\n', iter)
                    if rev>benchmark_expenditures
                        dub = experiment_taxes.hh.dividend;
                    elseif rev<=benchmark_expenditures
                        dlb = experiment_taxes.hh.dividend;
                    end
                end
                welfares(ic) = quantities.hh.welfare;
                dits(ic)     = experiment_taxes.hh.dividend;
                errors(ic)   = error;
            end

            
            experiment2.benchmark_welfare       = benchmark_welfare;
            experiment2.counterfactual_welfares = welfares;
            experiment2.cits                    = cits;
            experiment2.dits                    = dits; 
            experiment2.errors                  = errors;  %#ok<STRNU>


            save(fullfile(save_dir,'results.mat'),'experiment2')




        end
        
        
        
        function [] = solve_experiment3()
            % This experiment determines whether to tax corporate income with partial expensing or non-linear dividend taxation.
            
            
            % Identify and create results directories
            [~, save_dir] = identify_dirs(3);
            if exist(save_dir,'dir')
                rmdir(save_dir,'s')
            end
            mkdir(save_dir)
            
            

            % Generates benchmark tax parameter values.
            taxes = paramGenerator.tax;
            taxes.hh.dividend = 0;
            
            [~, s_quantities] = solve_ss_equilibrium_global([], taxes);
            benchmark_expenditures = s_quantities.government.government_expenditures;
            benchmark_welfare = s_quantities.hh.welfare;
            clear('s_quantities')


            %% Solving optimal income vs corporate income tax
            cit_ub = taxes.firm.income;
            cit_lb = 0;
            n_cits = 20;
            cits = linspace(cit_lb, cit_ub, n_cits);

            maxiter = 60;

            dividend_income_shares = zeros(1, n_cits);
            welfares = zeros(1, n_cits);
            errors   = zeros(1, n_cits);

            parfor ic = 1:n_cits

                experiment_taxes = taxes;
                experiment_taxes.firm.income = cits(ic);

                share_ub = 1;
                share_lb = 0;

                iter = 0;
                while true
                    iter = iter + 1;
                    % Determine tax revenue
                    experiment_taxes.hh.div_inc_share = (share_lb + share_ub)/2;
                    [~, quantities] = solve_ss_equilibrium_global([], experiment_taxes); 
                    rev = quantities.government.tax_revenue;

                    % Determine error and break if complete
                    error = abs(rev - benchmark_expenditures);
                    if (error<0.00001)||(iter>maxiter) , break, end 

                    % Update criteria
                    if rev > benchmark_expenditures
                        share_ub = experiment_taxes.hh.div_inc_share;
                    elseif rev<= benchmark_expenditures
                        share_lb = experiment_taxes.hh.div_inc_share;
                    end

                end

                dividend_income_shares(ic) = experiment_taxes.hh.div_inc_share;
                welfares(ic) = quantities.hh.welfare;
                errors(ic)   = error;

            end


            experiment3.benchmark_welfare     = benchmark_welfare;
            experiment3.counterfactual_welfares = welfares; 
            experiment3.cits                    = cits;
            experiment3.dividend_income_share   = dividend_income_shares;
            experiment3.errors                  = errors; %#ok<STRNU>
            

            save(fullfile(save_dir,'results.mat'),'experiment3')






        end

        
        
        
        
    end
    
end
            