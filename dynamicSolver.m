%%
% Dynamic model solver.
% 
%%


classdef dynamicSolver

methods (Static)    
    
    
    % Solve open economy transition path
    function [save_dir] = open(basedef, counterdef, callertag)
        
        if ~exist('counterdef', 'var'), counterdef = []; end
        if ~exist('callertag' , 'var'), callertag  = ''; end
        
        save_dir = dynamicSolver.transition('open', basedef, counterdef, callertag);
        
    end
    
    
    % Solve closed economy transition path
    function [save_dir] = closed(basedef, counterdef, callertag)
        
        if ~exist('counterdef', 'var'), counterdef = []; end
        if ~exist('callertag' , 'var'), callertag  = ''; end
        
        save_dir = dynamicSolver.transition('closed', basedef, counterdef, callertag);
        
    end
    
    
end


methods (Static, Access = private)
    
    function [save_dir] = transition(economy, basedef, counterdef, callertag)
        
        
        %% Initialization
        
        % Identify working directories
        param_dir = dirFinder.param;
        [save_dir, callingtag] = dirFinder.save(economy, basedef, counterdef);
        
        % Append caller tag to save directory name and calling tag
        % (Obviates conflicts between parallel runs)
        save_dir   = [save_dir  , callertag];
        callingtag = [callingtag, callertag];
        
        % Clear or create save directory
        if exist(save_dir, 'dir')
            rmdir(save_dir, 's')
        end
        mkdir(save_dir)
        
        
        % Unpack parameters from baseline definition
        beta  = basedef.beta ;
        gamma = basedef.gamma;
        sigma = basedef.sigma;
        
        % Identify baseline run by absence of counterfactual definition
        isbase = isempty(counterdef) || isempty(fields(counterdef));
        if isbase, counterdef = struct(); end
        
        % Unpack parameters from counterfactual definition, setting baseline values for unspecified parameters
        if isfield(counterdef, 'plan'), plan = counterdef.plan; else plan = 'base'; end
        if isfield(counterdef, 'gcut'), gcut = counterdef.gcut; else gcut = +0.00 ; end
        
        
        
        %% Parameter loading
        
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
        
        
        % Clear parameter loading structure
        clear('s')
        
        
        
        %% Static aggregate generation
        
        if ~isbase
            
            % Identify baseline generator and save directory
            base_generator = @() dynamicSolver.transition(economy, basedef, [], callingtag);
            base_dir = dirFinder.save(economy, basedef);
            
            % Load baseline solution
            s = hardyload('solution.mat', base_generator, base_dir);
            
            wages        = s.wages      ;
            cap_shares   = s.cap_shares ;
            debt_shares  = s.debt_shares;
            rate_caps    = s.rate_caps  ;
            rate_govs    = s.rate_govs  ;
            
            % Load baseline aggregates
            s = hardyload('aggregates.mat', base_generator, base_dir);
            exp_subsidys = s.exp_subsidys;
            
            % Load optimal decision values and distributions from baseline
            s = hardyload('opt.mat' , base_generator, base_dir);
            base_opt  = s.opt ;
            
            s = hardyload('dist.mat', base_generator, base_dir);
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
                    labopt = base_opt(birthyear+T,idem).lab;
                    
                    % Calculate tax terms
                    [fedincome, fedincomess, fitax, fitaxss, fsstax, fsstaxss, fcaptax, fcaptaxss] ...
                      ...
                      = calculate_static_taxes(...
                          T, Tr, Tss, birthyear, ...
                          kgrid, nb, nk, nz, idem, ...
                          mpci, rpci, cap_tax_share, ss_tax_cred, ...
                          tau_cap, tau_capgain, ss_tax, taxmax, z, ...
                          wages, cap_shares, debt_shares, rate_caps, rate_govs, exp_subsidys, q_tobin, q_tobin0, ...
                          avg_deduc, clinton, coefs, limit, X, ss_benefit, ...
                          labopt);
                    
                    
                    % Extract distributions
                    dist_w = base_dist(birthyear+T,idem).w;
                    dist_r = base_dist(birthyear+T,idem).r;
                    
                    % Redefine effective numbers of working and retirement years
                    % (Redefinition needed because distribution generator makes use of head truncation while dynamic optimization solver does not)
                    N_w_dist = size(dist_w, 4);
                    N_r_dist = size(dist_r, 3);
                    
                    % Calculate cohort aggregates
                    Fedincome = [sum(reshape(fedincome  (:,:,:,1:N_w_dist) .* dist_w, [], N_w_dist), 1), ...
                                 sum(reshape(fedincomess(:,:,1:N_r_dist)   .* dist_r, [], N_r_dist), 1)];
                    
                    Fedit     = [sum(reshape(fitax      (:,:,:,1:N_w_dist) .* dist_w, [], N_w_dist), 1), ...
                                 sum(reshape(fitaxss    (:,:,1:N_r_dist)   .* dist_r, [], N_r_dist), 1)];
                    
                    SSrev     = [sum(reshape(fsstax     (:,:,:,1:N_w_dist) .* dist_w, [], N_w_dist), 1), ...
                                 sum(reshape(fsstaxss   (:,:,1:N_r_dist)   .* dist_r, [], N_r_dist), 1)];
                    
                    Fcaptax   = [sum(reshape(fcaptax    (:,:,:,1:N_w_dist) .* dist_w, [], N_w_dist), 1), ...
                                 sum(reshape(fcaptaxss  (:,:,1:N_r_dist)   .* dist_r, [], N_r_dist), 1)];
                    
                    SSexp     = sum(reshape(repmat(ss_benefit(:,1)', [nk,1,N_r_dist]) .* dist_r, [], N_r_dist), 1);
                    
                    
                    % Add cohort aggregates to static aggregates, aligning to projection years
                    T_shift = max(0, birthyear);
                    fedincome_static( T_shift + (1:N_w_dist+N_r_dist) ) = fedincome_static( T_shift + (1:N_w_dist+N_r_dist) ) + Fedincome;
                    fedit_static    ( T_shift + (1:N_w_dist+N_r_dist) ) = fedit_static    ( T_shift + (1:N_w_dist+N_r_dist) ) + Fedit;
                    ssrev_static    ( T_shift + (1:N_w_dist+N_r_dist) ) = ssrev_static    ( T_shift + (1:N_w_dist+N_r_dist) ) + SSrev;
                    fcaptax_static  ( T_shift + (1:N_w_dist+N_r_dist) ) = fcaptax_static  ( T_shift + (1:N_w_dist+N_r_dist) ) + Fcaptax;
                    ssexp_static    ( T_shift + N_w_dist+(1:N_r_dist) ) = ssexp_static    ( T_shift + N_w_dist+(1:N_r_dist) ) + SSexp;
                    
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
            foreign_fcaptax_static  = zeros(1,Tss); %#ok<NASGU>
            
            fcaprev_static          = domestic_fcaptax_static; %#ok<NASGU>
            
            
            % Save static aggregates
            save(fullfile(save_dir, 'aggregates_static.mat'), ...
                 'cap_static', ...
                 'domestic_cap_static', 'foreign_cap_static', ...
                 'domestic_debt_static', 'foreign_debt_static', ...
                 'domestic_fcaptax_static', 'foreign_fcaptax_static', ...
                 'elab_static', 'lfpr_static', ...
                 'fedincome_static', 'fedit_static', 'ssrev_static', 'fcaptax_static', 'ssexp_static', ...
                 'Y_static', 'labinc_static', 'kinc_static', 'feditlab_static', 'fcaprev_static')
            
        end
        
        
        
        %% Dynamic aggregate generation
        
        % Identify steady state generator and save directory
        steady_generator = @() solve_ss([beta, gamma, sigma], false, callingtag);
        steady_dir = dirFinder.ss([beta, gamma, sigma]);
        
        % Load steady state solution
        s = hardyload('solution.mat', steady_generator, steady_dir);
        
        rho_ss          = s.rho         ;
        beq_ss          = s.beq         ;
        kpr_ss          = s.kpr         ;
        debt_ss         = s.DEBTss      ;
        cap_share_ss    = s.cap_share   ;
        debt_share_ss   = s.debt_share  ;
        rate_cap_ss     = s.rate_cap    ;
        rate_gov_ss     = s.rate_gov    ;
        rate_tot_ss     = s.rate_tot    ;
        
        % Load steady state distributions
        s = hardyload('dist.mat', steady_generator, steady_dir);
        
        dist_ss = s.dist;
        
        
        switch economy
            
            case 'open'
            
                % Load government expenditure adjustment parameters
                s = load(fullfile(param_dir, 'param_gtilde.mat'));
            
                indplan = cellfun(@(str) strncmp(plan, str, length(str)), {'base'; 'trump'; 'clinton'; 'ryan'});
                revenue_percent = s.revenue_percent(indplan,:);
                revenue_percent = [revenue_percent, revenue_percent(end)*ones(1, max(Tss-length(revenue_percent), 0))];
                
                if ~isbase
                    
                    % Calculate government expenditure adjustments
                    s_base = hardyload('aggregates.mat', base_generator, base_dir);
                    
                    Y_base = s_base.Y_total;
                    Gtilde = s_base.Gtilde - gcut * s.GEXP_percent(1:Tss) .* Y_base;
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
        
        
        
        % Start a parallel pool if one does not already exist and JVM is enabled
        if usejava('jvm'), gcp; end
        
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
        str1 = [upper(economy(1)), economy(2:end)];
        if isbase, str2 = 'baseline'; else str2 = 'counterfactual'; end
        fprintf('\n[%s economy %s]\n', str1, str2)
        
        
        while ( eps > tol && iter < itermax )
            
            % Increment iteration count
            iter = iter + 1;
            fprintf('\tIteration %2d  ...  ', iter)
            
            
            % Define prices
            switch economy
                
                case 'open'
                    
                    if (iter == 1)
                        cap_total = (kpr_ss - debt_ss)*ones(1,Tss)/q_tobin0;
                        
                        cap_shares  = cap_share_ss *ones(1,Tss);
                        debt_shares = debt_share_ss*ones(1,Tss);
                        rate_caps   = (((1 - tau_cap_base)/(1 - tau_cap))*(rate_cap_ss - 1) + 1)*ones(1,Tss);
                        rate_govs   = rate_gov_ss  *ones(1,Tss);
                        rate_tots   = rate_tot_ss  *ones(1,Tss);
                        
                        rhos = ((q_tobin*(rate_caps - 1) + d)/alp).^(1/(alp-1));
                        beqs = beq_ss*ones(1,Tss);
                    else
                        beqs = beq_total;
                    end
                    
                case 'closed'
                    
                    if (iter == 1)
                        kpr_total  = kpr_ss *ones(1,Tss);
                        debt_total = debt_ss*ones(1,Tss);
                        cap_total  = (kpr_total - debt_total)/q_tobin;
                        
                        rhos = rho_ss*ones(1,Tss);
                        beqs = beq_ss*ones(1,Tss);
                    else
                        rhos = 0.3*rhos + 0.7*rhoprs;
                        beqs = beq_total;
                    end
                    
                    cap_shares  = (kpr_total - debt_total) ./ kpr_total;
                    debt_shares = 1 - cap_shares;
                    rate_caps   = 1 + (A*alp*(rhos.^(alp-1)) - d)/q_tobin;
                    rate_govs   = rate_govs_cbo;
                    rate_tots   = cap_shares.*rate_caps + debt_shares.*rate_govs;
                    
            end
            
            wages        = A*(1-alp)*(rhos.^alp);
            exp_subsidys = [exp_share * max(diff(cap_total), 0), 0] ./ cap_total;
            
            
            % Initialize dynamic aggregates
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
            s_birthyear = struct('birthyear', num2cell(repmat(birthyears', [1,ndem])));
            opt     = s_birthyear;
            dist    = s_birthyear;
            cohorts = s_birthyear;
            
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
                    
                    % Generate optimal decision values, distributions, and cohort aggregates
                    [labopt, dist_w, dist_r, N_w, N_r, Kalive, Kdead, ELab, Lab, Lfpr, Fedincome, Fedit, SSrev, Fcaptax, SSexp] ...
                     ...
                       = generate_distributions(...
                           beta, gamma, sigma, T, Tr, Tss, birthyear, ...
                           kgrid, z, tr_z, bgrid, nk, nz, nb, idem, ...
                           mpci, rpci, cap_tax_share, ss_tax_cred, surv, tau_cap, tau_capgain, ss_tax, taxmax, ...
                           beqs, wages, cap_shares, debt_shares, rate_caps, rate_govs, rate_tots, exp_subsidys, q_tobin, q_tobin0, Vbeq, ...
                           avg_deduc, clinton, coefs, limit, X, ss_benefit, ...
                           dist_w_ss, dist_r_ss, mu2_idem, mu3_idem);
                    
                    % Store values
                    opt(i,idem).lab = labopt; %#ok<PFOUS>
                    
                    dist(i,idem).w = dist_w; %#ok<PFOUS>
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
                
                for i = 1:length(birthyears)
                    
                    % Add cohort aggregates to dynamic aggregates, aligning to projection years
                    T_shift = max(0, birthyears(i));
                    
                    N_w = cohorts(i,idem).N_w;
                    N_r = cohorts(i,idem).N_r;
                    
                    kpr_total      ( T_shift + (1:N_w+N_r) ) = kpr_total      ( T_shift + (1:N_w+N_r) ) + cohorts(i,idem).Kalive + cohorts(i,idem).Kdead;
                    beq_total      ( T_shift + (1:N_w+N_r) ) = beq_total      ( T_shift + (1:N_w+N_r) ) + cohorts(i,idem).Kdead;
                    elab_total     ( T_shift + (1:N_w)     ) = elab_total     ( T_shift + (1:N_w)     ) + cohorts(i,idem).ELab;
                    lab_total      ( T_shift + (1:N_w)     ) = lab_total      ( T_shift + (1:N_w)     ) + cohorts(i,idem).Lab;
                    lfpr_total     ( T_shift + (1:N_w)     ) = lfpr_total     ( T_shift + (1:N_w)     ) + cohorts(i,idem).Lfpr;
                    fedincome_total( T_shift + (1:N_w+N_r) ) = fedincome_total( T_shift + (1:N_w+N_r) ) + cohorts(i,idem).Fedincome;
                    fedit_total    ( T_shift + (1:N_w+N_r) ) = fedit_total    ( T_shift + (1:N_w+N_r) ) + cohorts(i,idem).Fedit;
                    ssrev_total    ( T_shift + (1:N_w+N_r) ) = ssrev_total    ( T_shift + (1:N_w+N_r) ) + cohorts(i,idem).SSrev;
                    fcaptax_total  ( T_shift + (1:N_w+N_r) ) = fcaptax_total  ( T_shift + (1:N_w+N_r) ) + cohorts(i,idem).Fcaptax;
                    ssexp_total    ( T_shift + N_w+(1:N_r) ) = ssexp_total    ( T_shift + N_w+(1:N_r) ) + cohorts(i,idem).SSexp;
                    
                end
                
            end
            
            
            % Calculate additional dynamic aggregates
            % (Note that open economy requires capital calculation before debt calculation while closed economy requires the reverse)
            switch economy
                
                case 'open'
                    
                    % Calculate capital and output
                    cap_total = rhos .* elab_total;
                    Y_total   = A*(max(cap_total, 0).^alp).*(elab_total.^(1-alp));
                    
                    domestic_cap_total = cap_shares .* [kpr_ss, kpr_total(1:end-1)];
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
                    debt_total = [debt_ss, zeros(1,Tss-1)];
                    for t = 1:Tss-1
                        debt_total(t+1) = Gtilde(t) - Ttilde(t) - netrev_total(t) + debt_total(t)*rate_govs_cbo(t);
                    end
                    
                    domestic_debt_total = (1 - cap_shares) .* kpr_total;
                    foreign_debt_total  = debt_total - domestic_debt_total; %#ok<NASGU>
                    
                    % Calculate convergence series delta
                    delta = beqs - beq_total;
                    
                case 'closed'
                    
                    % Calculate debt
                    domestic_fcaptax_total = fcaptax_total;
                    foreign_fcaptax_total  = zeros(1,Tss);
                    fcaptax_total          = domestic_fcaptax_total + foreign_fcaptax_total;
                    
                    netrev_total = fedit_total + ssrev_total + fcaptax_total - ssexp_total;
                    debt_total = [debt_ss, zeros(1,Tss-1)];
                    for t = 1:Tss-1
                        debt_total(t+1) = Gtilde(t) - Ttilde(t) - netrev_total(t) + debt_total(t)*rate_govs_cbo(t);
                    end
                    
                    domestic_debt_total = debt_total;   %#ok<NASGU>
                    foreign_debt_total  = zeros(1,Tss); %#ok<NASGU>
                    
                    % Calculate capital and output
                    cap_total = ([(kpr_ss - debt_ss)/q_tobin0, (kpr_total(1:end-1) - debt_total(2:end))/q_tobin]);
                    Y_total   = A*(max(cap_total, 0).^alp).*(elab_total.^(1-alp)); %#ok<NASGU>
                    
                    domestic_cap_total = [q_tobin0 * cap_total(1), q_tobin * cap_total(2:end)]; %#ok<NASGU>
                    foreign_cap_total  = zeros(1,Tss); %#ok<NASGU>
                    
                    % Calculate convergence series delta
                    rhoprs = (max([kpr_ss, kpr_total(1:end-1)] - debt_total, 0)/q_tobin) ./ elab_total;
                    delta = rhos - rhoprs;
                    
            end
            
            labinc_total = elab_total .* wages;
            kinc_total   = q_tobin * (rate_caps - 1) .* cap_total; %#ok<NASGU>
            
            feditlab_total = fedit_total .* labinc_total ./ fedincome_total;
            fcaprev_total  = fcaptax_total + fedit_total - feditlab_total; %#ok<NASGU>
            
            
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
             'rhos', 'wages', 'cap_shares', 'debt_shares', 'rate_caps', 'rate_govs')
        
        % Save dynamic aggregates
        save(fullfile(save_dir, 'aggregates.mat'), ...
             'kpr_total', 'cap_total', 'debt_total', ...
             'Gtilde', 'Ttilde', ...
             'domestic_cap_total', 'foreign_cap_total', ...
             'domestic_debt_total', 'foreign_debt_total', ...
             'domestic_fcaptax_total', 'foreign_fcaptax_total', ...
             'beq_total', 'elab_total', 'lab_total', 'lfpr_total', ...
             'fedincome_total', 'fedit_total', 'ssrev_total', 'fcaptax_total', 'ssexp_total', ...
             'Y_total', 'labinc_total', 'kinc_total', 'feditlab_total', 'fcaprev_total', 'exp_subsidys')
        
        
    end
    
end

end




%%
% Check if file exists and generate if necessary before loading, handling parallel write and read conflicts.
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