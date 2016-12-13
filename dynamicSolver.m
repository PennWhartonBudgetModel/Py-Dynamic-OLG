%%
% Dynamic model solver.
% 
%%


classdef dynamicSolver

methods (Static)    
    
    
    % Solve steady state
    function [save_dir] = steady(basedef, callertag)
        
        if ~exist('callertag' , 'var'), callertag  = ''; end
        
        save_dir = dynamicSolver.solve('steady', basedef, [], callertag);
        
    end
    
    
    % Solve open economy transition path
    function [save_dir] = open(basedef, counterdef, callertag)
        
        if ~exist('counterdef', 'var'), counterdef = []; end
        if ~exist('callertag' , 'var'), callertag  = ''; end
        
        save_dir = dynamicSolver.solve('open', basedef, counterdef, callertag);
        
    end
    
    
    % Solve closed economy transition path
    function [save_dir] = closed(basedef, counterdef, callertag)
        
        if ~exist('counterdef', 'var'), counterdef = []; end
        if ~exist('callertag' , 'var'), callertag  = ''; end
        
        save_dir = dynamicSolver.solve('closed', basedef, counterdef, callertag);
        
    end
    
    
    % Generate tags for baseline and counterfactual definitions
    function [basedef_tag, counterdef_tag] = generate_tags(basedef, counterdef)
        
        % Define baseline and counterfactual parameter formats
        basedef_format    = struct( 'beta'      , '%0.3f'   , ...
                                    'gamma'     , '%0.3f'   , ...
                                    'sigma'     , '%05.2f'  );
        
        counterdef_format = struct( 'plan'      , '%s'      , ...
                                    'gcut'      , '%+0.2f'  );
        
        % Define function to construct tag from definition and format specifications
        function tag = construct_tag(def, format)
            strs = {};
            for field = fields(def)'
                strs = [strs, {sprintf([field{1}, '=', format.(field{1})], def.(field{1}))}]; %#ok<AGROW>
            end
            tag = strjoin(strs, '_');
        end
        
        % Generate baseline tag
        basedef_tag = construct_tag(basedef, basedef_format);
        
        % Generate counterfactual tag
        if (isempty(counterdef) || isempty(fields(counterdef)))
            counterdef_tag = 'baseline';
        else
            counterdef_tag = construct_tag(dynamicSolver.fill_default(counterdef), counterdef_format);
        end
        
    end
    
end


methods (Static, Access = private)
    
    % Generate counterfactual definition filled with default parameter values where necessary
    function [counterdef_filled] = fill_default(counterdef)
        
        % Define default counterfactual parameter values
        % (These are the values used for a baseline run)
        counterdef_filled = struct( 'plan'    , 'base'    , ...
                                    'gcut'    , +0.00     );
        
        % Override default parameter values with values from counterfactual definition
        for field = fields(counterdef)'
            counterdef_filled.(field{1}) = counterdef.(field{1});
        end
        
    end
    
    
    % Solve dynamic model
    function [save_dir] = solve(economy, basedef, counterdef, callertag)
        
        
        %% Initialization
        
        % Start a parallel pool if one does not already exist and JVM is enabled
        switch economy, case {'open', 'closed'}, if usejava('jvm'), gcp; end, end
        
        % Unpack parameters from baseline definition
        beta  = basedef.beta ;
        gamma = basedef.gamma;
        sigma = basedef.sigma;
        
        % Identify baseline run by empty counterfactual definition
        if isempty(counterdef), counterdef = struct(); end
        isbase = isempty(fields(counterdef));
        
        % Generate counterfactual definition filled with default parameter values where necessary
        counterdef_filled = dynamicSolver.fill_default(counterdef);
        
        % Unpack parameters from filled counterfactual definition
        plan = counterdef_filled.plan;
        gcut = counterdef_filled.gcut;
        
        
        % Identify working directories
        param_dir = dirFinder.param;
        [save_dir, ~, counterdef_tag] = dirFinder.save(economy, basedef, counterdef);
        
        % Append caller tag to save directory name and generate calling tag
        % (Obviates conflicts between parallel runs)
        save_dir   = [save_dir                                  , callertag];
        callingtag = [sprintf('^%s_%s', counterdef_tag, economy), callertag];
        
        % Clear or create save directory
        if exist(save_dir, 'dir'), rmdir(save_dir, 's'), end, mkdir(save_dir)        
        
        
        
        %% Parameter loading
        
        % Load global parameters
        s = load(fullfile(param_dir, 'param_global.mat'));
        
        T_life  = s.T_life;
        switch economy
            case 'steady'
                T_model    = 1;
                birthyears = 0;
            case {'open', 'closed'}
                T_model    = s.T_model;
                birthyears = (-T_life+1):(T_model-1);
        end
        
        kgrid = s.kgrid;
        bgrid = [0; s.bgrid(2:end)];
        
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
        proddist    = s.proddist(:,1)';
        surv        = [s.surv(1:T_life-1), 0];
        tr_z        = s.tr_z;
        z           = s.z;
        
        MU2 = s.demdist_2015 * (s.Mu2/sum(s.Mu2));
        MU3 = repmat(1-surv, [ndem,1]) .* MU2;
        
        
        % Load social security parameters
        s = load(fullfile(param_dir, 'param_socsec.mat'));
        
        Tr         = s.NRA(1);
        ss_benefit = s.ss_benefit(:,1:T_model);
        ss_tax     = s.ss_tax(1:T_model);
        taxmax     = s.taxmax(1:T_model);
        
        
        % Load CBO parameters
        s = load(fullfile(param_dir, 'param_cbo.mat'));

        D_Y          = s.FederalDebtHeldbythePublic(1)/100;
        fedgovtnis   = s.fedgovtnis(1:T_model);
        rate_cbos    = 1 + s.r_cbo(1:T_model);
        meanrate_cbo = 1 + mean(s.r_cbo);
        
        
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
        
        
        % Clear parameter loading structure
        clear('s')
        
        
        
        %% Static aggregate generation
        
        if ~isbase
            
            % Identify baseline generator and save directory
            base_generator = @() dynamicSolver.solve(economy, basedef, [], callingtag);
            base_dir = dirFinder.save(economy, basedef);
            
            % Load baseline solution
            s = hardyload('solution.mat', base_generator, base_dir);
            
            wages        = s.wages       ;
            cap_shares   = s.cap_shares  ;
            debt_shares  = s.debt_shares ;
            rate_caps    = s.rate_caps   ;
            rate_govs    = s.rate_govs   ;
            exp_subsidys = s.exp_subsidys;
            
            % Load optimal decision values and distributions from baseline
            s = hardyload('opt.mat' , base_generator, base_dir);
            base_opt  = s.opt ;
            
            s = hardyload('dist.mat', base_generator, base_dir);
            base_dist = s.dist;
            
            
            % Clear parameter loading structure
            clear('s')
            
            
            
            % Initialize static aggregates
            fedincome_static = zeros(1,T_model);
            fedit_static     = zeros(1,T_model);
            ssrev_static     = zeros(1,T_model);
            fcaptax_static   = zeros(1,T_model);
            ssexp_static     = zeros(1,T_model);
            
            
            for idem = 1:ndem %#ok<FXUP>
                
                for i = 1:length(birthyears) %#ok<FXUP>
                    
                    % Extract optimal labor values and distributions
                    labopt_static = base_opt (i,idem).lab;
                    dist_w_static = base_dist(i,idem).w  ;
                    dist_r_static = base_dist(i,idem).r  ;
                    
                    % Generate cohort aggregates
                    [~, ~, ~, N_w, N_r, ~, ~, ~, ~, ~, Fedincome, Fedit, SSrev, Fcaptax, SSexp] ...
                     ...
                       = generate_distributions(...
                           beta, gamma, sigma, T_life, Tr, T_model, T_model, birthyears(i), ...
                           kgrid, z, tr_z, bgrid, nk, nz, nb, idem, ...
                           mpci, rpci, cap_tax_share, ss_tax_cred, surv, tau_cap, tau_capgain, ss_tax, taxmax, ...
                           zeros(1,T_model), wages, cap_shares, debt_shares, rate_caps, rate_govs, zeros(1,T_model), exp_subsidys, q_tobin, q_tobin0, Vbeq, ...
                           avg_deduc, clinton, coefs, limit, X, ss_benefit, ...
                           [], [], ones(1,T_model), ones(1,T_model), ...
                           labopt_static, dist_w_static, dist_r_static);
                    
                    % Add cohort aggregates to static aggregates, aligning to model years
                    T_shift = max(0, birthyears(i));
                    fedincome_static( T_shift + (1:N_w+N_r) ) = fedincome_static( T_shift + (1:N_w+N_r) ) + Fedincome   ;
                    fedit_static    ( T_shift + (1:N_w+N_r) ) = fedit_static    ( T_shift + (1:N_w+N_r) ) + Fedit       ;
                    ssrev_static    ( T_shift + (1:N_w+N_r) ) = ssrev_static    ( T_shift + (1:N_w+N_r) ) + SSrev       ;
                    fcaptax_static  ( T_shift + (1:N_w+N_r) ) = fcaptax_static  ( T_shift + (1:N_w+N_r) ) + Fcaptax     ;
                    ssexp_static    ( T_shift + N_w+(1:N_r) ) = ssexp_static    ( T_shift + N_w+(1:N_r) ) + SSexp       ;
                    
                end
                
            end
            
            
            % Copy additional static aggregates from baseline aggregates
            s_base = hardyload('aggregates.mat', base_generator, base_dir);
            
            elab_static             = s_base.elab_total             ; %#ok<NASGU>
            cap_static              = s_base.cap_total              ; %#ok<NASGU>
            lfpr_static             = s_base.lfpr_total             ; %#ok<NASGU>
            labinc_static           = s_base.labinc_total           ;
            kinc_static             = s_base.kinc_total             ; %#ok<NASGU>
            
            Y_static                = s_base.Y_total                ; %#ok<NASGU>
            
            domestic_cap_static     = s_base.domestic_cap_total     ; %#ok<NASGU>
            foreign_cap_static      = s_base.foreign_cap_total      ; %#ok<NASGU>
            
            domestic_debt_static    = s_base.domestic_debt_total    ; %#ok<NASGU>
            foreign_debt_static     = s_base.foreign_debt_total     ; %#ok<NASGU>
            
            clear('s_base')
            
            
            % Calculate additional static aggregates
            feditlab_static         = fedit_static .* labinc_static ./ fedincome_static;
            
            domestic_fcaptax_static = fcaptax_static + fedit_static - feditlab_static;
            foreign_fcaptax_static  = zeros(1,T_model); %#ok<NASGU>
            
            fcaprev_static          = domestic_fcaptax_static; %#ok<NASGU>
            
            
            % Save static aggregates
            save(fullfile(save_dir, 'aggregates_static.mat'), ...
                 'domestic_debt_static', 'foreign_debt_static', ...
                 'cap_static', 'domestic_cap_static', 'foreign_cap_static', ...
                 'Y_static', ...
                 'elab_static', 'lfpr_static', ...
                 'fedincome_static', 'fedit_static', 'ssrev_static', ...
                 'fcaptax_static', 'domestic_fcaptax_static', 'foreign_fcaptax_static', ...
                 'ssexp_static', 'labinc_static', 'kinc_static', 'feditlab_static', 'fcaprev_static')
            
        end
        
        
        
        %% Dynamic aggregate generation
        
        switch economy
            
            case 'steady'
                
                % Load neutral solution as starting solution
                s0 = load(fullfile(param_dir, 'solution0.mat'));
                
            case {'open', 'closed'}
                
                % Identify steady state generator and save directory
                steady_generator = @() dynamicSolver.steady(basedef, callingtag);
                steady_dir = dirFinder.save('steady', basedef);
                
                % Load steady state solution as starting solution
                s0 = hardyload('solution.mat', steady_generator, steady_dir);
                
                % Load steady state distributions
                s = hardyload('dist.mat', steady_generator, steady_dir);
                
                dist_steady = s.dist;
                
        end
        
        % Unpack starting solution
        rho0        = s0.rhos        ;
        beq0        = s0.beqs        ;
        kpr0        = s0.kprs        ;
        debt0       = s0.debts       ;
        cap_share0  = s0.cap_shares  ;
        debt_share0 = s0.debt_shares ;
        rate_cap0   = s0.rate_caps   ;
        rate_gov0   = s0.rate_govs   ;
        rate_tot0   = s0.rate_tots   ;
        
        clear('s0')
        
        
        switch economy
            
            case 'open'
                
                % Load government expenditure adjustment parameters
                s = load(fullfile(param_dir, 'param_gtilde.mat'));
                
                indplan = cellfun(@(str) strncmp(plan, str, length(str)), {'base'; 'trump'; 'clinton'; 'ryan'});
                revenue_percent = s.revenue_percent(indplan,:);
                revenue_percent = [revenue_percent, revenue_percent(end)*ones(1, max(T_model-length(revenue_percent), 0))];
                
                if ~isbase
                    
                    % Calculate government expenditure adjustments
                    s_base = hardyload('aggregates.mat', base_generator, base_dir);
                    
                    Y_base = s_base.Y_total;
                    Gtilde = s_base.Gtilde - gcut * s.GEXP_percent(1:T_model) .* Y_base;
                    Ttilde = revenue_percent .* Y_base - fedit_static - ssrev_static - fcaptax_static;
                    
                    clear('s_base')
                    
                end
                
            case 'closed'
                
                % Identify open economy generator and save directory
                open_generator = @() dynamicSolver.open(basedef, counterdef, callingtag);
                open_dir = dirFinder.save('open', basedef, counterdef);
                
                % Load government expenditure adjustments
                s = hardyload('aggregates.mat', open_generator, open_dir);
                
                Gtilde = s.Gtilde;
                Ttilde = s.Ttilde;
                
        end
        
        
        % Clear parameter loading structure
        clear('s')
        
        
        
        % Define convergence tolerance and initialize error term
        tol = 1e-3;
        eps = Inf;
        
        % Initialize iteration count and set maximum number of iterations
        iter    =  0;
        itermax = 25;
        
        % Create file for logging iterations
        % (Note that this should be done after any parallel pool is started to avoid file access issues)
        iterlog = fopen(fullfile(save_dir, 'iterations.csv'), 'w');
        
        % Display header
        fprintf('\n[')
        switch economy
    
            case 'steady'
                fprintf('Steady state')
            
            case {'open', 'closed'}
                str1 = [upper(economy(1)), economy(2:end)];
                if isbase, str2 = 'baseline'; else, str2 = 'counterfactual'; end
                fprintf('%s economy %s', str1, str2)
                
        end
        fprintf(']\n')
        
        
        % Define dynamic aggregate generation function
        % (Modularization necessary for steady state labor elasticity calculation)
        function [kpr_total, beq_total, elab_total, lab_total, lfpr_total, ...
                  fedincome_total, fedit_total, ssrev_total, fcaptax_total, ssexp_total, ...
                  opt, dist] ...
                  ... 
                    = generate_aggregates(rhos, beqs, kprs, debts, caps, wages, ...
                                          cap_shares, debt_shares, rate_caps, rate_govs, rate_tots, exp_subsidys) %#ok<INUSL>
            
            % Initialize dynamic aggregates
            kpr_total       = zeros(1,T_model);
            beq_total       = zeros(1,T_model);
            elab_total      = zeros(1,T_model);
            lab_total       = zeros(1,T_model);
            lfpr_total      = zeros(1,T_model);
            fedincome_total = zeros(1,T_model);
            fedit_total     = zeros(1,T_model);
            ssrev_total     = zeros(1,T_model);
            fcaptax_total   = zeros(1,T_model);
            ssexp_total     = zeros(1,T_model);
            
            % Initialize storage structures for optimal decision values and distributions
            s_birthyear = struct('birthyear', num2cell(repmat(birthyears', [1,ndem])));
            opt     = s_birthyear;
            dist    = s_birthyear;
            cohorts = s_birthyear;
            
            % Define dynamic optimization and distribution generation time periods
            switch economy
        
                case 'steady'
                    T_opt  = 1;
                    T_dist = T_life;
                
                case {'open', 'closed'}
                    T_opt  = T_model;
                    T_dist = T_model;
                    
            end
            
            for idem = 1:ndem %#ok<FXUP>
                
                % Extract demographic adjustments
                mu2_idem = MU2(idem,:);
                mu3_idem = MU3(idem,:);
                
                % Define initial distributions
                switch economy
            
                    case 'steady'
                        dist_w0 = padarray(proddist, [nk-1, 0, nb-1], 0, 'post');
                        dist_r0 = [];
                    
                    case {'open', 'closed'}
                        % if isdynamic
                        N_w_steady = length(dist_steady(1,idem).w);
                        N_r_steady = length(dist_steady(1,idem).r);
                        dist_w0 = bsxfun(@rdivide, dist_steady(1,idem).w, shiftdim(mu2_idem(            1:N_w_steady ), -2));
                        dist_r0 = bsxfun(@rdivide, dist_steady(1,idem).r, shiftdim(mu2_idem(N_w_steady+(1:N_r_steady)), -1));
                        % else
                end
                
                parfor i = 1:length(birthyears) %#ok<FXUP>
                    
                    % Generate optimal decision values, distributions, and cohort aggregates
                    [labopt, dist_w, dist_r, N_w, N_r, Kalive, Kdead, ELab, Lab, Lfpr, Fedincome, Fedit, SSrev, Fcaptax, SSexp] ...
                     ...
                       = generate_distributions(...
                           beta, gamma, sigma, T_life, Tr, T_opt, T_dist, birthyears(i), ...
                           kgrid, z, tr_z, bgrid, nk, nz, nb, idem, ...
                           mpci, rpci, cap_tax_share, ss_tax_cred, surv, tau_cap, tau_capgain, ss_tax, taxmax, ...
                           beqs, wages, cap_shares, debt_shares, rate_caps, rate_govs, rate_tots, exp_subsidys, q_tobin, q_tobin0, Vbeq, ...
                           avg_deduc, clinton, coefs, limit, X, ss_benefit, ...
                           dist_w0, dist_r0, mu2_idem, mu3_idem, ...
                           [], [], []);
                    
                    % Store values
                    opt(i,idem).lab = labopt;
                    
                    dist(i,idem).w = dist_w;
                    dist(i,idem).r = dist_r;
                    
                    cohorts(i,idem).N_w         = N_w           ;
                    cohorts(i,idem).N_r         = N_r           ;
                    cohorts(i,idem).Kalive      = Kalive        ;
                    cohorts(i,idem).Kdead       = Kdead         ;
                    cohorts(i,idem).ELab        = ELab          ;
                    cohorts(i,idem).Lab         = Lab           ;
                    cohorts(i,idem).Lfpr        = Lfpr          ;
                    cohorts(i,idem).Fedincome   = Fedincome     ;
                    cohorts(i,idem).Fedit       = Fedit         ;
                    cohorts(i,idem).SSrev       = SSrev         ;
                    cohorts(i,idem).Fcaptax     = Fcaptax       ;
                    cohorts(i,idem).SSexp       = SSexp         ;
                    
                end
                
                % Add cohort aggregates to dynamic aggregates
                for i = 1:length(birthyears) %#ok<FXUP>
                    
                    switch economy
                        
                        case 'steady'
                            
                            kpr_total  = kpr_total  + sum(cohorts(i,idem).Kalive + cohorts(i,idem).Kdead);
                            beq_total  = beq_total  + sum(cohorts(i,idem).Kdead);
                            elab_total = elab_total + sum(cohorts(i,idem).ELab);
                            
                        case {'open', 'closed'}
                            
                            % Align aggregates to model years
                            T_shift = max(0, birthyears(i));
                            
                            N_w = cohorts(i,idem).N_w;
                            N_r = cohorts(i,idem).N_r;
                            
                            kpr_total      ( T_shift + (1:N_w+N_r) ) = kpr_total      ( T_shift + (1:N_w+N_r) ) + cohorts(i,idem).Kalive + cohorts(i,idem).Kdead;
                            beq_total      ( T_shift + (1:N_w+N_r) ) = beq_total      ( T_shift + (1:N_w+N_r) ) + cohorts(i,idem).Kdead     ;
                            elab_total     ( T_shift + (1:N_w)     ) = elab_total     ( T_shift + (1:N_w)     ) + cohorts(i,idem).ELab      ;
                            lab_total      ( T_shift + (1:N_w)     ) = lab_total      ( T_shift + (1:N_w)     ) + cohorts(i,idem).Lab       ;
                            lfpr_total     ( T_shift + (1:N_w)     ) = lfpr_total     ( T_shift + (1:N_w)     ) + cohorts(i,idem).Lfpr      ;
                            fedincome_total( T_shift + (1:N_w+N_r) ) = fedincome_total( T_shift + (1:N_w+N_r) ) + cohorts(i,idem).Fedincome ;
                            fedit_total    ( T_shift + (1:N_w+N_r) ) = fedit_total    ( T_shift + (1:N_w+N_r) ) + cohorts(i,idem).Fedit     ;
                            ssrev_total    ( T_shift + (1:N_w+N_r) ) = ssrev_total    ( T_shift + (1:N_w+N_r) ) + cohorts(i,idem).SSrev     ;
                            fcaptax_total  ( T_shift + (1:N_w+N_r) ) = fcaptax_total  ( T_shift + (1:N_w+N_r) ) + cohorts(i,idem).Fcaptax   ;
                            ssexp_total    ( T_shift + N_w+(1:N_r) ) = ssexp_total    ( T_shift + N_w+(1:N_r) ) + cohorts(i,idem).SSexp     ;
                            
                    end
                    
                end
                
            end
            
        end
        
        
        while ( eps > tol && iter < itermax )
            
            % Increment iteration count
            iter = iter + 1;
            fprintf('\tIteration %2d  ...  ', iter)
            
            
            % Define prices
            switch economy
                
                case {'steady', 'closed'}
                    
                    if (iter == 1)
                        rhos  = rho0 *ones(1,T_model);
                        beqs  = beq0 *ones(1,T_model);
                        kprs  = kpr0 *ones(1,T_model);
                        debts = debt0*ones(1,T_model);
                        caps  = (kprs - debts)/q_tobin;
                    else
                        switch economy
                            case 'steady'
                                rhos  = 0.5*rhos + 0.5*rhoprs;
                                debts = D_Y*Y_total;
                            case 'closed'
                                rhos  = 0.3*rhos + 0.7*rhoprs;
                                debts = debt_total;
                        end
                        beqs  = beq_total;
                        kprs  = kpr_total;
                        caps  = cap_total;
                    end
                    
                    cap_shares  = (kprs - debts) ./ kprs;
                    debt_shares = 1 - cap_shares;
                    rate_caps   = 1 + (A*alp*(rhos.^(alp-1)) - d)/q_tobin;
                    switch economy
                        case 'steady', rate_govs = meanrate_cbo;
                        case 'closed', rate_govs = rate_cbos;
                    end
                    rate_tots   = cap_shares.*rate_caps + debt_shares.*rate_govs;
                    
                    
                case 'open'
                    
                    if (iter == 1)
                        
                        kprs  = kpr0 *ones(1,T_model);
                        debts = debt0*ones(1,T_model);
                        caps  = (kprs - debts)/q_tobin0;
                        
                        cap_shares  = cap_share0 *ones(1,T_model);
                        debt_shares = debt_share0*ones(1,T_model);
                        rate_caps   = (((1 - tau_cap_base)/(1 - tau_cap))*(rate_cap0 - 1) + 1)*ones(1,T_model);
                        rate_govs   = rate_gov0  *ones(1,T_model);
                        rate_tots   = rate_tot0  *ones(1,T_model);
                        
                        rhos  = ((q_tobin*(rate_caps - 1) + d)/alp).^(1/(alp-1));
                        beqs  = beq0*ones(1,T_model);
                        
                    else
                        beqs = beq_total;
                        caps = cap_total;
                    end
                    
            end
            
            wages        = A*(1-alp)*(rhos.^alp);
            exp_subsidys = [exp_share * max(diff(caps), 0), 0] ./ caps;
            
            
            % Generate dynamic aggregates
            [kpr_total, beq_total, elab_total, lab_total, lfpr_total, ...
             fedincome_total, fedit_total, ssrev_total, fcaptax_total, ssexp_total, ...
             opt, dist] ...
             ...
               = generate_aggregates(rhos, beqs, kprs, debts, caps, wages, ...
                                     cap_shares, debt_shares, rate_caps, rate_govs, rate_tots, exp_subsidys); %#ok<ASGLU>
            
            
            % Calculate additional dynamic aggregates
            % (Note that open economy requires capital calculation before debt calculation while closed economy requires the reverse)
            switch economy
                
                case 'steady'
                    
                    % Calculate debt
                    debt_total = debts;
                    
                    % Calculate capital and output
                    cap_total = (kpr_total - debt_total)/q_tobin;
                    Y_total   = A*(max(cap_total, 0).^alp).*(elab_total.^(1-alp));
                    
                    % Calculate convergence value
                    rhoprs = (max(kpr_total - debt_total, 0)/q_tobin) ./ elab_total;
                    delta = rhos - rhoprs;
                    
                case 'open'
                    
                    % Calculate capital and output
                    cap_total = rhos .* elab_total;
                    Y_total   = A*(max(cap_total, 0).^alp).*(elab_total.^(1-alp));
                    
                    domestic_cap_total = cap_shares .* [kpr0, kpr_total(1:end-1)];
                    foreign_cap_total  = q_tobin*cap_total - domestic_cap_total;
                    
                    % Calculate debt
                    domestic_fcaptax_total = fcaptax_total;
                    foreign_fcaptax_total  = tau_cap.*(rate_caps - 1).*cap_tax_share.*foreign_cap_total;
                    fcaptax_total          = domestic_fcaptax_total + foreign_fcaptax_total;
                    
                    if isbase
                        Gtilde = (revenue_percent - fedgovtnis).*Y_total - ssexp_total;
                        Ttilde = revenue_percent.*Y_total - fedit_total - ssrev_total - fcaptax_total;
                    end
                    
                    netrev_total = fedit_total + ssrev_total + fcaptax_total - ssexp_total;
                    debt_total = [debt0, zeros(1,T_model-1)];
                    for t = 1:T_model-1
                        debt_total(t+1) = Gtilde(t) - Ttilde(t) - netrev_total(t) + debt_total(t)*rate_cbos(t);
                    end
                    
                    domestic_debt_total = (1 - cap_shares) .* kpr_total;
                    foreign_debt_total  = debt_total - domestic_debt_total; %#ok<NASGU>
                    
                    % Calculate income
                    labinc_total = elab_total .* wages;
                    kinc_total   = q_tobin * (rate_caps - 1) .* cap_total; %#ok<NASGU>
                    
                    feditlab_total = fedit_total .* labinc_total ./ fedincome_total;
                    fcaprev_total  = fcaptax_total + fedit_total - feditlab_total; %#ok<NASGU>
                    
                    % Calculate convergence value
                    delta = beqs - beq_total;
                    
                case 'closed'
                    
                    % Calculate debt
                    domestic_fcaptax_total = fcaptax_total;
                    foreign_fcaptax_total  = zeros(1,T_model);
                    fcaptax_total          = domestic_fcaptax_total + foreign_fcaptax_total;
                    
                    netrev_total = fedit_total + ssrev_total + fcaptax_total - ssexp_total;
                    debt_total = [debt0, zeros(1,T_model-1)];
                    for t = 1:T_model-1
                        debt_total(t+1) = Gtilde(t) - Ttilde(t) - netrev_total(t) + debt_total(t)*rate_cbos(t);
                    end
                    
                    domestic_debt_total = debt_total;   %#ok<NASGU>
                    foreign_debt_total  = zeros(1,T_model); %#ok<NASGU>
                    
                    % Calculate capital and output
                    cap_total = ([(kpr0 - debt0)/q_tobin0, (kpr_total(1:end-1) - debt_total(2:end))/q_tobin]);
                    Y_total   = A*(max(cap_total, 0).^alp).*(elab_total.^(1-alp));
                    
                    domestic_cap_total = [q_tobin0 * cap_total(1), q_tobin * cap_total(2:end)]; %#ok<NASGU>
                    foreign_cap_total  = zeros(1,T_model); %#ok<NASGU>
                    
                    % Calculate income
                    labinc_total = elab_total .* wages;
                    kinc_total   = q_tobin * (rate_caps - 1) .* cap_total; %#ok<NASGU>
                    
                    feditlab_total = fedit_total .* labinc_total ./ fedincome_total;
                    fcaprev_total  = fcaptax_total + fedit_total - feditlab_total; %#ok<NASGU>
                    
                    % Calculate convergence value
                    rhoprs = (max([kpr0, kpr_total(1:end-1)] - debt_total, 0)/q_tobin) ./ elab_total;
                    delta = rhos - rhoprs;
                    
            end
            
            
            % Calculate convergence error
            eps = max(abs(delta));
            
            fprintf('Error term = %7.4f\n', eps)
            fprintf(iterlog, '%u,%0.4f\n', iter, eps);
            
            
        end
        fprintf('\n')
        fclose(iterlog);
        
        % Issue warning if maximum iterations reached
        if (iter == itermax), warning('Maximum iterations reached.'), end
        
        
        % Save optimal decision values and distributions for baseline
        if isbase
            save(fullfile(save_dir, 'opt.mat' ), 'opt' )
            save(fullfile(save_dir, 'dist.mat'), 'dist')
        end
        
        % Save solution
        save(fullfile(save_dir, 'solution.mat'), ...
             'rhos', 'beqs', 'kprs', 'debts', 'caps', 'wages', ...
             'cap_shares', 'debt_shares', 'rate_caps', 'rate_govs', 'rate_tots', 'exp_subsidys')
        
        % Save dynamic aggregates
        switch economy
            case {'open', 'closed'}
                save(fullfile(save_dir, 'aggregates.mat'), ...
                     'beq_total', 'kpr_total', ...
                     'debt_total', 'domestic_debt_total', 'foreign_debt_total', ...
                     'cap_total', 'domestic_cap_total', 'foreign_cap_total', ...
                     'Y_total', ...
                     'elab_total', 'lab_total', 'lfpr_total', ...
                     'fedincome_total', 'fedit_total', 'ssrev_total', ...
                     'fcaptax_total', 'domestic_fcaptax_total', 'foreign_fcaptax_total', ...
                     'ssexp_total', 'labinc_total', 'kinc_total', 'feditlab_total', 'fcaprev_total', ...
                     'Gtilde', 'Ttilde')
        end
        
        
        
        %% Elasticity calculation
        
        switch economy
            case 'steady'
                
                % Calculate capital to output ratio
                K_Y = (kpr_total-debt_total)/Y_total;
                
                
                % Calculate labor elasticity
                working_mass = 0;
                frisch_total = 0;
                
                for idem = 1:ndem %#ok<FXUP>
                    
                    labopt = opt(1,idem).lab;
                    dist_w = dist(1,idem).w;
                    
                    working_ind = (labopt > 0.01);
                    
                    working_mass = working_mass + sum(dist_w(working_ind));
                    frisch_total = frisch_total + sum(dist_w(working_ind).*(1-labopt(working_ind))./labopt(working_ind))*(1-gamma*(1-sigma))/sigma;
                    
                end
                
                labor_elas = frisch_total / working_mass;
                
                
                % Calculate savings elasticity
                deviation = 0.005;
                
                rate_caps_dev = rate_caps * (1 + deviation);
                rate_govs_dev = rate_govs * (1 + deviation);
                rate_tots_dev = rate_tots * (1 + deviation);
                
                [kpr_dev] = generate_aggregates(rhos, beqs, kprs, debts, caps, wages, ...
                                                cap_shares, debt_shares, rate_caps_dev, rate_govs_dev, rate_tots_dev, exp_subsidys);
                
                savings_elas = ((kpr_dev - kpr_total)/kpr_total) / ((rate_tots_dev - rate_tots)/(rate_tots-1));
                
                
                % Package, save, and display elasticities
                elasticities = [K_Y, labor_elas, savings_elas];
                
                save(fullfile(save_dir, 'elasticities.mat'), 'K_Y', 'labor_elas', 'savings_elas')
                
                displaynames = { 'Capital-to-output ratio', 'Labor elasticity', 'Savings elasticity' };
                for i = 1:length(elasticities) %#ok<FXUP>
                    fprintf('\t%-25s= % 7.4f\n', displaynames{i}, elasticities(i))
                end
                fprintf('\n')
                
        end
        
        
    end
    
end

end




%%
% Check if file exists and generate before loading if necessary, handling parallel write and read conflicts.
% 
%%


function [s] = hardyload(filename, generator, save_dir)

% Check if file exists and generate if necessary
filepath = fullfile(save_dir, filename);
if ~exist(filepath, 'file')
    
    tagged_dir = generator();
    
    % Check for file again, possibly generated by a parallel run
    if ~exist(filepath, 'file')
        % Attempt to copy to save directory and continue if error encountered
        % (Errors are typically due to multiple simultaneous copy attempts, which should result in at least one successful copy)
        try copyfile(tagged_dir, save_dir); catch, end
    end
    
    % Clean up temporary directory
    rmdir(tagged_dir, 's')
    
end


% Turn on pause and store current pause state
pause0 = pause('on');

% Set maximum number of attempts and maximum pause time in seconds
maxtries = 200;
maxpause = 2.0;

for itry = 1:maxtries
    try
        % Attempt load
        s = load(filepath);
        break
    catch e
        if (itry == maxtries)
            error('Failed to load ''%s'' after %d attempts.\n%s', filepath, maxtries, e.message)
        end
        % Take a breather
        pause(rand * maxpause)
    end
end

% Reset pause state
pause(pause0)

end

