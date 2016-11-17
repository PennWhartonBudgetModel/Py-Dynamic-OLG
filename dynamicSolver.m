%%
% Dynamic model solver.
% 
% Methods:
% 
%   open(basedef, counterdef, partag)
%       Solve open economy transition path.
% 
%   closed(basedef, counterdef, partag)
%       Solve closed economy transition path.
% 
%%


classdef dynamicSolver

methods (Static)    
    
    function [save_dir] = open(basedef, counterdef, partag)
        
        if ~exist('partag', 'var'), partag = ''; end
        [save_dir] = dynamicSolver.transition(true, basedef, counterdef, partag);
        
    end
    
    
    function [save_dir] = closed(basedef, counterdef, partag)
        
        if ~exist('partag', 'var'), partag = ''; end
        [save_dir] = dynamicSolver.transition(false, basedef, counterdef, partag);
        
    end
    
end


methods (Static, Access = private)
    
    function [save_dir] = transition(isopen, basedef, counterdef, partag)
        
        % Unpack deep parameters from baseline definition
        beta  = basedef.beta ;
        gamma = basedef.gamma;
        sigma = basedef.sigma;
        deep_params = [beta, gamma, sigma];
        
        % Unpack policy parameters from counterfactual definition        
        plan = counterdef.plan;
        isbase = strcmp(plan, 'base');
        
        gcut = counterdef.gcut;
        
        
        % Identify working directories
        param_dir = dirFinder.param;
        if isopen
            save_dir = dirFinder.open  (deep_params, plan      );
        else
            save_dir = dirFinder.closed(deep_params, plan, gcut);
        end
        
        % Append parallelization tag to save directory name
        save_dir = [save_dir, partag];
        
        % Clear or create save directory
        if exist(save_dir, 'dir')
            rmdir(save_dir, 's')
        end
        mkdir(save_dir)
        
        
        
        % Load global parameters
        s = load(fullfile(param_dir, 'param_global.mat'));
        
        T   = s.T;
        Tss = s.Tss;
        
        kgrid = s.kgrid;
        bgrid = [0; s.bgrid(2:end)];            % (Lower bound on average earnings set to 0)
        
        Vbeq = s.phi1.*((1+kgrid./s.phi2).^(1-s.phi3));
        
        ndem = s.ndem;
        nb   = s.nb;
        nk   = s.nk;
        nz   = s.nz;
        
        mpci = s.mpci;
        rpci = s.rpci;
        
        ss_tax_cred = s.ss_tax_cred;
        
        A           = s.A;
        alp         = s.alp;
        d           = s.d;
        deduc_scale = s.deduc_scale;
        surv        = [s.surv(1:T-1), 0];       % (Survival probability at age T set to 0)
        tr_z        = s.tr_z;
        z           = s.z;
        
        MU2 = s.demdist_2015 * ( s.Mu2/sum(s.Mu2) );
        MU3 = repmat(1-surv, [ndem,1]) .* MU2;
        
        
        % Load social security parameters
        s = load(fullfile(param_dir, 'param_socsec.mat'));
        
        NRA        = s.NRA;
        ss_benefit = s.ss_benefit(:,1:Tss);
        ss_tax     = s.ss_tax(1:Tss);
        taxmax     = s.taxmax(1:Tss);
        
        
        % Load CBO parameters
        s = load(fullfile(param_dir, 'param_cbo.mat'));
        
        fedgovtnis     = s.fedgovtnis(1:Tss);
        rate_govs_cbo  = 1 + s.r_cbo(1:Tss);
        
        
        % Load income tax parameters
        s = load(fullfile(param_dir, sprintf('param_inctax_%s.mat', plan)));
        
        avg_deduc = deduc_scale * s.avg_deduc;
        clinton   = s.clinton;
        coefs     = s.coefs;
        limit     = s.limit;
        X         = s.X;
        
        
        % Load business tax parameters
        s = load(fullfile(param_dir, sprintf('param_bustax_%s.mat', plan)));
        
        cap_tax_share = s.cap_tax_share;
        exp_share     = s.exp_share;
        tau_cap       = s.tau_cap;
        tau_capgain   = s.tau_capgain;
        
        q_tobin = 1 - tau_cap * exp_share;
        
        s_base = load(fullfile(param_dir, 'param_bustax_base.mat'));
        tau_cap_base = s_base.tau_cap;
        q_tobin0     = 1 - tau_cap_base * s_base.exp_share;
        clear('s_base')
        
        
        % Identify reference steady state directory
        ss_dir = dirFinder.ss(deep_params);
        
        % Check for steady state solution and generate if not available
        ss_solution = fullfile(ss_dir, 'solution.mat');
        if ~exist(ss_solution, 'file')
            fprintf('\nSteady state solution not found.  Calling steady state solver.\n')
            
            if isopen
                uniquetag = sprintf('_open_plan=%s', plan);
            else
                uniquetag = sprintf('_closed_plan=%s_gcut=%+0.2f', plan, gcut);
            end
            solve_ss(deep_params, false, uniquetag);
            
            % Check for steady state solution again, possibly generated by a parallel run
            firstrun = false;
            if ~exist(ss_solution, 'file')
                % Attempt to copy to permanent directory and continue if error encountered
                % (Errors typically due to parallel copy attempt, which should result in at least one successful copy)
                try copyfile([ss_dir, uniquetag], ss_dir), firstrun = true; catch, end
            end
            
            fprintf('\nSteady state solution generated')
            if firstrun, fprintf('.\n'), else fprintf(' by a parallel run.\n'), end
            
            % Clean up temporary directory
            rmdir([ss_dir, uniquetag], 's')

        end
        
        % Load steady state solution
        s = hardyload(ss_solution);
        
        if isopen
        
            beqs        = s.beq       *ones(1,Tss);
            cap_shares  = s.cap_share *ones(1,Tss);
            debt_shares = s.debt_share*ones(1,Tss);
            rate_govs   = s.rate_gov  *ones(1,Tss);
            rate_tots   = s.rate_tot  *ones(1,Tss);

            rate_cap = ((1 - tau_cap_base)/(1 - tau_cap))*(s.rate_cap - 1) + 1;
            rate_caps = rate_cap*ones(1,Tss);

            rho = ((q_tobin*(rate_cap - 1) + d)/alp)^(1/(alp-1));
            rhos = rho*ones(1,Tss);

            wage = A*(1-alp)*rho^alp;
            wages = wage*ones(1,Tss);

            kpr_ss  = s.kpr;
            debt_ss = s.DEBTss;

            cap_series = (kpr_ss - debt_ss)*ones(1,Tss)/q_tobin0;
            
        else
            
            rho_ss  = s.rho;
            beq_ss  = s.beq;
            kpr_ss  = s.kpr;
            debt_ss = s.DEBTss;

            rhos = rho_ss *ones(1,Tss);
            beqs = beq_ss *ones(1,Tss);
            KK   = kpr_ss *ones(1,Tss);
            DD   = debt_ss*ones(1,Tss);

            cap_series = (KK - DD)/q_tobin;
            
            rate_govs = rate_govs_cbo;
            
        end
        
        % Load steady state distributions
        s = hardyload(fullfile(ss_dir, 'dist.mat'));
        
        dist_ss = s.dist;
        
        
        if isopen
            
            % Load government expenditure adjustment parameters
            s = load(fullfile(param_dir, 'param_gtilde.mat'));
            
            indplan = cellfun(@(s) strncmp(plan, s, length(s)), {'base'; 'trump'; 'clinton'; 'ryan'});
            
            revenue_percent = s.revenue_percent(indplan,:);
            if (length(revenue_percent) < Tss)
                revenue_percent = [revenue_percent, revenue_percent(end) * ones(1, Tss - length(revenue_percent))];
            end
            
        else
            
            % Identify reference open economy directories
            openbase_dir = dirFinder.open(deep_params, 'base');
            openplan_dir = dirFinder.open(deep_params, plan  );
            
            
            % Check for open economy baseline aggregates and generate if not available
            openbase_aggregates = fullfile(openbase_dir, 'aggregates.mat');
            if ~exist(openbase_aggregates, 'file')
                fprintf('\nOpen economy baseline aggregates not found.  Calling open economy solver.\n')
                
                uniquetag = sprintf('_closed_plan=%s_gcut=%+0.2f', plan, gcut);
                dynamicSolver.open( struct('beta', deep_params(1), 'gamma', deep_params(2), 'sigma', deep_params(3)), struct('plan', 'base', 'gcut', +0.00), uniquetag );
                
                % Check for open economy baseline aggregates again, possibly generated by a parallel run
                firstrun = false;
                if ~exist(openbase_aggregates, 'file')
                    % Attempt to copy to permanent directory and continue if error encountered
                    % (Errors typically due to parallel copy attempt, which should result in at least one successful copy)
                    try copyfile([openbase_dir, uniquetag], openbase_dir), firstrun = true; catch, end
                end
                
                fprintf('\nOpen economy baseline aggregates generated')
                if firstrun, fprintf('.'), else fprintf(' by a parallel run.'), end
                
                % Clean up temporary directory
                rmdir([openbase_dir, uniquetag], 's')
                fprintf('  Proceeding with closed economy dynamic aggregate generation.\n')
                
            end
            
            % Load open economy baseline aggregates
            s = hardyload(openbase_aggregates);
            Y_openbase = s.Y_total;
            
            
            % Check for open economy counterfactual aggregates and generate if not available
            openplan_aggregates = fullfile(openplan_dir, 'aggregates.mat');
            if ~exist(openplan_aggregates, 'file')
                fprintf('\nOpen economy counterfactual aggregates not found.  Calling open economy solver.\n')
                
                uniquetag = sprintf('_closed_plan=%s_gcut=%+0.2f', plan, gcut);
                dynamicSolver.open( struct('beta', deep_params(1), 'gamma', deep_params(2), 'sigma', deep_params(3)), struct('plan', plan, 'gcut', +0.00), uniquetag );
                
                % Check for open economy counterfactual aggregates again, possibly generated by a parallel run
                firstrun = false;
                if ~exist(openplan_aggregates, 'file')
                    % Attempt to copy to permanent directory and continue if error encountered
                    % (Errors typically due to parallel copy attempt, which should result in at least one successful copy)
                    try copyfile([openplan_dir, uniquetag], openplan_dir), firstrun = true; catch, end
                end
                
                fprintf('\nOpen economy counterfactual aggregates generated')
                if firstrun, fprintf('.'), else fprintf(' by a parallel run.'), end
                
                % Clean up temporary directory
                rmdir([openplan_dir, uniquetag], 's')
                fprintf('  Proceeding with closed economy dynamic aggregate generation.\n')
                
            end
            
            % Load open economy counterfactual aggregates
            s = hardyload(openplan_aggregates);
            Ttilde = s.Ttilde;
            Gtilde = s.Gtilde;
            
            % Generate government expenditure adjustments
            s = load(fullfile(param_dir, 'param_gtilde.mat'));
            GEXP_percent = s.GEXP_percent(1:Tss);
            GEXP_cut     = gcut*GEXP_percent.*Y_openbase;
            Gtilde       = Gtilde - GEXP_cut;
            
        end
        
        
        % Clear parameter loading structure
        clear('s')
        
        
        
        %% Dynamic aggregate generation
        
        % Define convergence tolerance and initialize error term sequence
        tol = 1e-3;
        eps = [];
        
        % Initialize iteration count and set maximum number of iterations
        iter    =  0;
        itermax = 25;
        
        while true
            
            % Increment iteration count
            iter = iter + 1;
            fprintf('\tIteration %2d  ...  ', iter)
            
            % Derive values
            exp_subsidys = [exp_share * max(diff(cap_series), 0), 0] ./ cap_series;
            
            if ~isopen
                
                wages        = A*(1-alp)*(rhos.^alp);
                
                cap_shares  = (KK - DD) ./ KK;
                debt_shares = 1 - cap_shares;
                
                rate_caps   = 1 + (A*alp*(rhos.^(alp-1)) - d)/q_tobin;
                
                rate_tots   = cap_shares.*rate_caps + debt_shares.*rate_govs;
                
            end
            
            
            % Initialize global aggregates
            kpr_total       = zeros(1,Tss);
            beq_total       = zeros(1,Tss);
            elab_total      = zeros(1,Tss);
            lab_total       = zeros(1,Tss);
            lfpr_total      = zeros(1,Tss);
            fedincome_total = zeros(1,Tss);
            fedit_total     = zeros(1,Tss);
            ssrev_total     = zeros(1,Tss);
            fcaptax_total   = zeros(1,Tss);
            ssexp_total     = zeros(1,Tss);
            
            
            % Define birth year range
            birthyears = (-T+1):(Tss-1);
            
            % Initialize storage structures for optimal decision values and distributions
            s_birthyear = struct('birthyear', num2cell(repmat(birthyears', [1,2])));
            opt  = s_birthyear;
            dist = s_birthyear;
            
            for idem = 1:ndem
                
                % Extract demographic adjustments
                mu2_idem = MU2(idem,:);
                mu3_idem = MU3(idem,:);
                
                % Extract steady state distributions to be used for initialization, reversing demographic adjustment
                % (Note that this currently assumes a fixed retirement age matching the one used for the steady state solution)
                N_w_ss = length(dist_ss(1,idem).dist_w);
                N_r_ss = length(dist_ss(1,idem).dist_r);
                dist_w_ss = bsxfun(@rdivide, dist_ss(1,idem).dist_w, shiftdim(mu2_idem(1:N_w_ss),          -2));
                dist_r_ss = bsxfun(@rdivide, dist_ss(1,idem).dist_r, shiftdim(mu2_idem(N_w_ss+(1:N_r_ss)), -1));
                
                
                parfor i = 1:length(birthyears)
                    
                    % Get birth year and retirement age
                    birthyear = birthyears(i);
                    Tr = NRA(i);
                    
                    
                    % Solve dynamic optimization
                    [V, Vss, kopt, koptss, labopt, bopt, ...
                     fedincome, fedincomess, ...
                     fitax,     fitaxss,     ...
                     fsstax,    fsstaxss,    ...
                     fcaptax,   fcaptaxss]   ...
                     ...
                       = solve_dynamic_optimization_mex(...
                           beqs, T, Tr, Tss, birthyear, ...
                           bgrid, kgrid, nb, nk, nz, idem, ...
                           mpci, rpci, cap_tax_share, ss_tax_cred, ...
                           surv, tau_cap, tau_capgain, ss_tax, tr_z, taxmax, z, ...
                           wages, cap_shares, debt_shares, rate_caps, rate_govs, rate_tots, exp_subsidys, q_tobin, q_tobin0, Vbeq, ...
                           avg_deduc, clinton, coefs, limit, X, ss_benefit, ...
                           beta, gamma, sigma);
                    
                    % (Duplicated variable assignments necessary for parfor)
                    opt(i,idem).V             = V             ; %#ok<PFOUS>
                    opt(i,idem).Vss           = Vss           ;
                    opt(i,idem).kopt          = kopt          ;
                    opt(i,idem).koptss        = koptss        ;
                    opt(i,idem).labopt        = labopt        ;
                    opt(i,idem).bopt          = bopt          ;
                    opt(i,idem).fedincome     = fedincome     ;
                    opt(i,idem).fedincomess   = fedincomess   ;
                    opt(i,idem).fitax         = fitax         ;
                    opt(i,idem).fitaxss       = fitaxss       ;
                    opt(i,idem).fsstax        = fsstax        ;
                    opt(i,idem).fsstaxss      = fsstaxss      ;
                    opt(i,idem).fcaptax       = fcaptax       ;
                    opt(i,idem).fcaptaxss     = fcaptaxss     ;
                    
                    % Check for unsolved dynamic optimization subproblems
                    if any(isinf([V(:); Vss(:)]))
                        warning('Infinite utility values found.  Some dynamic optimization subproblems unsolved.  Check that initial conditions satisfy constraints.')
                    end
                    
                    
                    % Generate distributions
                    % (Note that currently only the first column of ss_benefit is used to calculate social security benefits)
                    [dist_w, dist_r, N_w, N_r, Kalive, Kdead, ELab, Lab, Lfpr, Fedincome, Fedit, SSrev, Fcaptax, SSexp] ...
                     ...
                       = generate_distributions(...
                           kopt, koptss, labopt, bopt, ...
                           fedincome, fedincomess, fitax, fitaxss, fsstax, fsstaxss, fcaptax, fcaptaxss, ...
                           T, Tr, Tss, birthyear, ...
                           kgrid, bgrid, nk, nz, nb, idem, ...
                           z, tr_z, ss_benefit(:,1), ...
                           dist_w_ss, dist_r_ss, ...
                           mu2_idem, mu3_idem);
                    
                    % (Duplicated variable assignments necessary for parfor)
                    dist(i,idem).dist_w       = dist_w        ;
                    dist(i,idem).dist_r       = dist_r        ;
                    dist(i,idem).N_w          = N_w           ;
                    dist(i,idem).N_r          = N_r           ;
                    dist(i,idem).Kalive       = Kalive        ;
                    dist(i,idem).Kdead        = Kdead         ;
                    dist(i,idem).ELab         = ELab          ;
                    dist(i,idem).Lab          = Lab           ;
                    dist(i,idem).Lfpr         = Lfpr          ;
                    dist(i,idem).Fedincome    = Fedincome     ;
                    dist(i,idem).Fedit        = Fedit         ;
                    dist(i,idem).SSrev        = SSrev         ;
                    dist(i,idem).Fcaptax      = Fcaptax       ;
                    dist(i,idem).SSexp        = SSexp         ;
                    
                end
                
                for i = 1:length(birthyears)
                    
                    % Add aggregate values to global aggregates, aligning to projection years
                    T_shift = max(0, birthyears(i));
                    
                    N_w = dist(i,idem).N_w;
                    N_r = dist(i,idem).N_r;
                    
                    kpr_total      ( T_shift + (1:N_w+N_r) ) = kpr_total      ( T_shift + (1:N_w+N_r) ) + dist(i,idem).Kalive + dist(i,idem).Kdead;
                    beq_total      ( T_shift + (1:N_w+N_r) ) = beq_total      ( T_shift + (1:N_w+N_r) ) + dist(i,idem).Kdead;
                    elab_total     ( T_shift + (1:N_w)     ) = elab_total     ( T_shift + (1:N_w)     ) + dist(i,idem).ELab;
                    lab_total      ( T_shift + (1:N_w)     ) = lab_total      ( T_shift + (1:N_w)     ) + dist(i,idem).Lab;
                    lfpr_total     ( T_shift + (1:N_w)     ) = lfpr_total     ( T_shift + (1:N_w)     ) + dist(i,idem).Lfpr;
                    fedincome_total( T_shift + (1:N_w+N_r) ) = fedincome_total( T_shift + (1:N_w+N_r) ) + dist(i,idem).Fedincome;
                    fedit_total    ( T_shift + (1:N_w+N_r) ) = fedit_total    ( T_shift + (1:N_w+N_r) ) + dist(i,idem).Fedit;
                    ssrev_total    ( T_shift + (1:N_w+N_r) ) = ssrev_total    ( T_shift + (1:N_w+N_r) ) + dist(i,idem).SSrev;
                    fcaptax_total  ( T_shift + (1:N_w+N_r) ) = fcaptax_total  ( T_shift + (1:N_w+N_r) ) + dist(i,idem).Fcaptax;
                    ssexp_total    ( T_shift + N_w+(1:N_r) ) = ssexp_total    ( T_shift + N_w+(1:N_r) ) + dist(i,idem).SSexp;
                    
                end
                
            end
            
            % Calculate additional aggregates
            if isopen
                
                cap_total = rho * elab_total;
                
            else
                
                netrev_total = fedit_total + ssrev_total + fcaptax_total - ssexp_total;
                
                debt_total = [debt_ss, zeros(1,Tss-1)];
                for t = 1:Tss-1
                    debt_total(t+1) = Gtilde(t) - Ttilde(t) - netrev_total(t) + debt_total(t)*rate_govs(t);
                end
                
                cap_total = ([(kpr_ss - debt_ss)/q_tobin0, (kpr_total(1:end-1) - debt_total(2:end))/q_tobin]);
                
                rhoprs = (max([kpr_ss, kpr_total(1:end-1)] - debt_total, 0)/q_tobin) ./ elab_total;
                
            end
            
            
            % Check for convergence
            if isopen
                delta = beqs - beq_total;
            else
                delta = rhos - rhoprs;
            end
            eps = [eps, max(abs(delta))]; %#ok<AGROW>
            
            fprintf('Error term = %7.4f\n', eps(end))
            if (eps(end) < tol), break, end
            
            % Check for maximum iterations
            if (iter == itermax)
                warning('Maximum iterations reached.')
                break
            end
            
            
            % Update variables for next iteration
            beqs = beq_total;
            cap_series = cap_total;
            
            if ~isopen
                stepfactor = 0.3;
                rhos = (1-stepfactor) * rhos + stepfactor * rhoprs;
                
                KK = kpr_total;
                DD = debt_total;
            end
            
        end
        fprintf('\n')
        
        % Save log of iterations
        fid = fopen(fullfile(save_dir, 'iterations.txt'), 'wt');
        fprintf(fid, '%d iteration(s) of a maximum %d\n\nError term(s):', iter, itermax);
        for i = 1:iter
            fprintf(fid, '\n  %2d  --  %7.4f', i, eps(i));
        end
        fclose(fid);
        
        
        % Calculate additional aggregates
        if isopen
            
            domestic_cap_total  = cap_shares .* [kpr_ss, kpr_total(1:end-1)];
            foreign_cap_total   = q_tobin*cap_total - domestic_cap_total;
            
            domestic_debt_total    = (1 - cap_shares) .* kpr_total;
            
            domestic_fcaptax_total = fcaptax_total;
            foreign_fcaptax_total  = tau_cap.*(rate_caps - 1).*cap_tax_share.*foreign_cap_total;
            fcaptax_total          = domestic_fcaptax_total + foreign_fcaptax_total;
            
        else
            
            domestic_cap_total  = [q_tobin0 * cap_total(1), q_tobin * cap_total(2:end)]; %#ok<NASGU>
            foreign_cap_total   = zeros(1,Tss); %#ok<NASGU>
            
            domestic_debt_total = debt_total;
            foreign_debt_total  = zeros(1,Tss); %#ok<PREALL>
            
            domestic_fcaptax_total = fcaptax_total; %#ok<NASGU>
            foreign_fcaptax_total  = zeros(1,Tss);	%#ok<NASGU>
            
        end
        
        Y_total   = A*(max(cap_total, 0).^alp).*(elab_total.^(1-alp));
        
        labinc_total = elab_total .* wages;
        kinc_total   = q_tobin * (rate_caps - 1) .* cap_total; %#ok<NASGU>
        
        feditlab_total = fedit_total .* labinc_total ./ fedincome_total;
        fcaprev_total  = fcaptax_total + fedit_total - feditlab_total; %#ok<NASGU>
        
        
        % Save optimal decision values and distributions for baseline
        if isbase
            save(fullfile(save_dir, 'opt.mat' ), 'opt' )
            save(fullfile(save_dir, 'dist.mat'), 'dist')
        end
        
        % Save solution
        save(fullfile(save_dir, 'solution.mat'), ...
             'rhos', 'wages', 'cap_shares', 'debt_shares', 'rate_caps', 'rate_govs')
        
        % Save dynamic aggregates
        save(fullfile(save_dir, 'aggregates.mat'), ...
             'kpr_total', 'cap_total', ...
             'domestic_cap_total', 'foreign_cap_total', ...
             'domestic_debt_total', ...
             'domestic_fcaptax_total', 'foreign_fcaptax_total', ...
             'beq_total', 'elab_total', 'lab_total', 'lfpr_total', ...
             'fedincome_total', 'fedit_total', 'ssrev_total', 'fcaptax_total', 'ssexp_total', ...
             'Y_total', 'labinc_total', 'kinc_total', 'feditlab_total', 'fcaprev_total', 'exp_subsidys')
        
        
        
        %% Static aggregate generation
        
        % Identify baseline directory
        if ~isbase
            if isopen
                base_dir = dirFinder.open(deep_params, 'base');
            else
                base_dir = dirFinder.closed(deep_params, 'base', +0.00);
            end
        else
            base_dir = save_dir;    % Includes this_uniquetag
        end
        
        
        % Check for baseline aggregates and generate if not available
        if ~isbase
            
            base_aggregates = fullfile(base_dir, 'aggregates.mat');
            if ~exist(base_aggregates, 'file')
                fprintf('\nBaseline aggregates not found.  Calling solver.\n')
                
                if isopen
                    uniquetag = sprintf('_open_plan=%s', plan);
                    dynamicSolver.open( basedef, struct('plan', 'base', 'gcut', +0.00), uniquetag );
                else
                    uniquetag = sprintf('_closed_plan=%s_gcut=%+0.2f', plan, gcut);
                    dynamicSolver.closed( basedef, struct('plan', 'base', 'gcut', +0.00), uniquetag );
                end
                    
                % Check for baseline aggregates again, possibly generated by a parallel run
                firstrun = false;
                if ~exist(base_aggregates, 'file')
                    % Attempt to copy to permanent directory and continue if error encountered
                    % (Errors typically due to parallel copy attempt, which should result in at least one successful copy)
                    try copyfile([base_dir, uniquetag], base_dir), firstrun = true; catch, end
                end
                
                fprintf('\nBaseline aggregates generated')
                if firstrun, fprintf('.'), else fprintf(' by a parallel run.'), end
                
                % Clean up temporary directory
                rmdir([base_dir, uniquetag], 's')
                fprintf('  Proceeding with static aggregate generation.\n')
                
            end
            
        end
        
        if isopen
            
            % Load steady state solution
            s = hardyload(fullfile(ss_dir, 'solution.mat'));
            
            wages        = s.wage    *ones(1,Tss);
            rate_caps    = s.rate_cap*ones(1,Tss);
            exp_subsidys = zeros(1,Tss);
            
        else
            
            % Load baseline solution
            s = hardyload(fullfile(base_dir, 'solution.mat'));
            
            wages        = s.wages;
            cap_shares   = s.cap_shares;
            debt_shares  = s.debt_shares;
            rate_caps    = s.rate_caps;
            rate_govs    = s.rate_govs;
            
            % Load baseline aggregates
            s = hardyload(fullfile(base_dir, 'aggregates.mat'));
            exp_subsidys = s.exp_subsidys;
            
        end
        
        % Load optimal decision values and distributions from baseline
        s = hardyload(fullfile(base_dir, 'opt.mat'));
        base_opt = s.opt;
        
        s = hardyload(fullfile(base_dir, 'dist.mat'));
        base_dist = s.dist;
        
        
        % Clear parameter loading structure
        clear('s')
        
        
        
        % Initialize static aggregates
        fedincome_static = zeros(1,Tss);
        fedit_static     = zeros(1,Tss);
        ssrev_static     = zeros(1,Tss);
        fcaptax_static   = zeros(1,Tss);
        ssexp_static     = zeros(1,Tss);
        
        
        for idem = 1:ndem
            
            for birthyear = (-T+1):(Tss-1)
                
                % Get retirement age
                Tr = NRA(birthyear+T);
                
                
                % Extract optimal labor values
                labopt = base_opt(birthyear+T,idem).labopt;
                
                % Calculate tax terms
                [fedincome, fedincomess, fitax, fitaxss, fsstax, fsstaxss, fcaptax, fcaptaxss] ...
                  ...
                  = calculate_static_taxes_mex(...
                      T, Tr, Tss, birthyear, ...
                      kgrid, nb, nk, nz, idem, ...
                      mpci, rpci, cap_tax_share, ss_tax_cred, ...
                      tau_cap, tau_capgain, ss_tax, taxmax, z, ...
                      wages, cap_shares, debt_shares, rate_caps, rate_govs, exp_subsidys, q_tobin, q_tobin0, ...
                      avg_deduc, clinton, coefs, limit, X, ss_benefit, ...
                      labopt);
                
                
                % Extract distributions
                dist_w = base_dist(birthyear+T,idem).dist_w;
                dist_r = base_dist(birthyear+T,idem).dist_r;
                
                % Redefine effective numbers of working and retirement years
                % (Redefinition needed because distribution generator makes use of head truncation while dynamic optimization solver does not yet)
                N_w_dist = size(dist_w, 4);
                N_r_dist = size(dist_r, 3);
                
                % Calculate aggregates
                Fedincome = [sum(reshape(fedincome  (:,:,:,1:N_w_dist) .* dist_w, [], N_w_dist), 1), ...
                             sum(reshape(fedincomess(:,:,1:N_r_dist)   .* dist_r, [], N_r_dist), 1)];
                
                Fedit     = [sum(reshape(fitax      (:,:,:,1:N_w_dist) .* dist_w, [], N_w_dist), 1), ...
                             sum(reshape(fitaxss    (:,:,1:N_r_dist)   .* dist_r, [], N_r_dist), 1)];
                
                SSrev     = [sum(reshape(fsstax     (:,:,:,1:N_w_dist) .* dist_w, [], N_w_dist), 1), ...
                             sum(reshape(fsstaxss   (:,:,1:N_r_dist)   .* dist_r, [], N_r_dist), 1)];
                
                Fcaptax   = [sum(reshape(fcaptax    (:,:,:,1:N_w_dist) .* dist_w, [], N_w_dist), 1), ...
                             sum(reshape(fcaptaxss  (:,:,1:N_r_dist)   .* dist_r, [], N_r_dist), 1)];
                
                SSexp     = sum(reshape(repmat(ss_benefit(:,1)', [nk,1,N_r_dist]) .* dist_r, [], N_r_dist), 1);
                
                
                % Add aggregates to total static aggregates, aligning to projection years
                T_shift = max(0, birthyear);
                fedincome_static( T_shift + (1:N_w_dist+N_r_dist) ) = fedincome_static( T_shift + (1:N_w_dist+N_r_dist) ) + Fedincome;
                fedit_static    ( T_shift + (1:N_w_dist+N_r_dist) ) = fedit_static    ( T_shift + (1:N_w_dist+N_r_dist) ) + Fedit;
                ssrev_static    ( T_shift + (1:N_w_dist+N_r_dist) ) = ssrev_static    ( T_shift + (1:N_w_dist+N_r_dist) ) + SSrev;
                fcaptax_static  ( T_shift + (1:N_w_dist+N_r_dist) ) = fcaptax_static  ( T_shift + (1:N_w_dist+N_r_dist) ) + Fcaptax;
                ssexp_static    ( T_shift + N_w_dist+(1:N_r_dist) ) = ssexp_static    ( T_shift + N_w_dist+(1:N_r_dist) ) + SSexp;
                
            end
            
        end
        
        
        % Copy additional static aggregates from corresponding baseline dynamic aggregates
        s_base = hardyload(fullfile(base_dir, 'aggregates.mat'));
        
        elab_static             = s_base.elab_total             ; %#ok<NASGU>
        cap_static              = s_base.cap_total              ; %#ok<NASGU>
        lfpr_static             = s_base.lfpr_total             ; %#ok<NASGU>
        labinc_static           = s_base.labinc_total           ;
        kinc_static             = s_base.kinc_total             ; %#ok<NASGU>
        
        Y_static                = s_base.Y_total                ; %#ok<NASGU>
        
        domestic_cap_static     = s_base.domestic_cap_total     ; %#ok<NASGU>
        foreign_cap_static      = s_base.foreign_cap_total      ; %#ok<NASGU>
        
        domestic_debt_static    = s_base.domestic_debt_total    ; %#ok<NASGU>
        
        clear('s_base')
        
        
        % Calculate additional static aggregates
        feditlab_static         = fedit_static .* labinc_static ./ fedincome_static;
        
        domestic_fcaptax_static = fcaptax_static + fedit_static - feditlab_static;
        foreign_fcaptax_static  = zeros(1,Tss); %#ok<NASGU>
        
        fcaprev_static          = domestic_fcaptax_static; %#ok<NASGU>
        
        
        % Save static aggregates
        save(fullfile(save_dir, 'aggregates_static.mat'), ...
             'cap_static', ...
             'domestic_cap_static', 'foreign_cap_static', ...
             'domestic_debt_static', ...
             'domestic_fcaptax_static', 'foreign_fcaptax_static', ...
             'elab_static', 'lfpr_static', ...
             'fedincome_static', 'fedit_static', 'ssrev_static', 'fcaptax_static', 'ssexp_static', ...
             'Y_static', 'labinc_static', 'kinc_static', 'feditlab_static', 'fcaprev_static')
        
        
        
        %% Foreign debt calculation
        
        if isopen

            % Load baseline output
            if ~isbase
                s_base = hardyload(base_aggregates);
                Y_base = s_base.Y_total;
            else
                Y_base = Y_total;
            end

            % Calculate government expenditure adjustments
            Ttilde = revenue_percent.*Y_base - fedit_static - ssrev_static - fcaptax_static;

            if ~isbase
                Gtilde = s_base.Gtilde;
            else
                Gtilde = fedit_total + ssrev_total + fcaptax_total + Ttilde - ssexp_total - fedgovtnis.*Y_total;
            end

            % Calculate total debt
            netrev_total = fedit_total + ssrev_total + fcaptax_total - ssexp_total;

            debt_total = [debt_ss, zeros(1,Tss-1)];
            for t = 1:Tss-1
                debt_total(t+1) = Gtilde(t) - Ttilde(t) - netrev_total(t) + debt_total(t)*rate_govs_cbo(t);
            end

            % Calculate foreign debt
            foreign_debt_total = debt_total - domestic_debt_total;
            if ~isbase
                foreign_debt_static = s_base.foreign_debt_total; %#ok<NASGU>
            else
                foreign_debt_static = foreign_debt_total; %#ok<NASGU>
            end
            
        else
            
            foreign_debt_static = zeros(1,Tss); %#ok<NASGU>
            
        end

        
        % Save additional dynamic aggregates
        save(fullfile(save_dir, 'aggregates.mat'), '-append', ...
             'debt_total', 'foreign_debt_total', 'Gtilde', 'Ttilde')
        
        % Save additional static aggregate
        save(fullfile(save_dir, 'aggregates_static.mat'), '-append', ...
             'foreign_debt_static');
        
        
    end
    
end

end