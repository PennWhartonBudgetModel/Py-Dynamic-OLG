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
        
        counterdef_format = struct( 'taxplan'   , '%s'      , ...
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
        counterdef_filled = struct( 'taxplan' , 'base'    , ...
                                    'gcut'    , +0.00     );
        
        % Override default parameter values with values from counterfactual definition
        for field = fields(counterdef)'
            counterdef_filled.(field{1}) = counterdef.(field{1});
        end
        
    end
    
    
    % Solve dynamic model
    function [save_dir] = solve(economy, basedef, counterdef, callertag)
        
        
        %% Initialization
        
        % Start parallel pool if JVM enabled and pool does not already exist
        if usejava('jvm'), gcp; end
        
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
        taxplan = counterdef_filled.taxplan;
        gcut    = counterdef_filled.gcut   ;
        
        
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
                startyears = 0;
            case {'open', 'closed'}
                T_model    = s.T_model;
                startyears = (-T_life+1):(T_model-1);
        end
        nstartyears = length(startyears);
        
        nz   = s.nz;
        nk   = s.nk;
        nb   = s.nb;
        ndem = s.ndem;
        
        zs     = s.z;
        transz = s.tr_z;
        DISTz  = s.proddist(:,1)';
        ks     = s.kgrid;
        bs     = [0; s.bgrid(2:end)];
        
        A   = s.A;
        alp = s.alp;
        d   = s.d;
        
        surv = [s.surv(1:T_life-1), 0];
        V_beq = s.phi1.*((1+ks./s.phi2).^(1-s.phi3));
        
        mu2 = s.demdist_2015 * (s.Mu2/sum(s.Mu2));
        mu3 = repmat(1-surv, [ndem,1]) .* mu2;
        
        mpci = s.mpci;
        rpci = s.rpci;
        
        deducscale  = s.deduc_scale;
        sstaxcredit = s.ss_tax_cred;
        
        
        % Load social security parameters
        s = load(fullfile(param_dir, 'param_socsec.mat'));
        
        T_work     = s.NRA(1);
        ssbenefits = s.ss_benefit(:,1:T_model);
        sstaxs     = s.ss_tax(1:T_model);
        ssincmaxs  = s.taxmax(1:T_model);
        
        
        % Load CBO parameters
        s = load(fullfile(param_dir, 'param_cbo.mat'));
        
        D_Y         = s.FederalDebtHeldbythePublic(1)/100;
        fedgovtnis  = s.fedgovtnis(1:T_model);
        cborates    = 1 + s.r_cbo(1:T_model);
        cbomeanrate = 1 + mean(s.r_cbo);
        
        
        % Load income tax parameters
        s = load(fullfile(param_dir, sprintf('param_inctax_%s.mat', taxplan)));
        
        deduc_coefs = [deducscale * s.avg_deduc, s.coefs(1:2)];
        pit_coefs   = [s.limit, s.X];
        
        
        % Load business tax parameters
        s = load(fullfile(param_dir, sprintf('param_bustax_%s.mat', taxplan)));
        
        captaxshare = s.cap_tax_share;
        expshare    = s.exp_share;
        taucap      = s.tau_cap;
        taucapgain  = s.tau_capgain;
        
        qtobin = 1 - taucap * expshare;
        
        s_base = load(fullfile(param_dir, 'param_bustax_base.mat'));
        taucap_base = s_base.tau_cap;
        qtobin0 = 1 - taucap_base * s_base.exp_share;
        clear('s_base')
        
        
        % Clear parameter loading structure
        clear('s')
        
        
        
        %% Aggregate generation function
        
        function [kpr_total, beq_total, elab_total, lab_total, lfpr_total, ...
                  fedincome_total, fedit_total, ssrev_total, fcaptax_total, ssexp_total, ...
                  LABs, DISTs] ...
                  ...
                    = generate_aggregates(beqs, wages, capshares, debtshares, caprates, govrates, totrates, expsubsidys, ...
                                          DISTs_steady, LABs_static, DISTs_static) 
            
            % Initialize aggregates
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
            
            % Initialize data storage arrays
            LABs    = cell(nstartyears, ndem);
            DISTs   = cell(nstartyears, ndem);
            cohorts = cell(nstartyears, ndem);
            
            % Set empty values for static optimal decision values and distributions if not provided
            if isempty(LABs_static)
                LABs_static  = LABs ;
                DISTs_static = DISTs;
            end
            
            % Define dynamic optimization and distribution generation model time periods
            switch economy
        
                case 'steady'
                    T_model_opt  = 1;
                    T_model_dist = T_life;
                    
                case {'open', 'closed'}
                    T_model_opt  = T_model;
                    T_model_dist = T_model;
                    
            end
            
            
            for idem = 1:ndem %#ok<FXUP>
                
                % Extract demographic adjustments
                mu2_idem = mu2(idem,:);
                mu3_idem = mu3(idem,:);
                
                % Define initial distributions
                switch economy
                    
                    case 'steady'
                        DIST0 = padarray(DISTz', [0, nk-1, nb-1], 0, 'post');
                        
                    case {'open', 'closed'}
                        if ~isempty(DISTs_steady)
                            DIST0 = DISTs_steady{1,idem} ./ repmat(shiftdim(mu2_idem(1:size(DISTs_steady{1,idem}, 4)), -2), [nz,nk,nb,1]);
                        else
                            DIST0 = [];
                        end
                        
                end
                
                parfor i = 1:nstartyears %#ok<FXUP>
                    
                    % Extract static optimal decision values and distributions if available
                    LAB_static  = LABs_static {i,idem};
                    DIST_static = DISTs_static{i,idem};
                    
                    % Generate cohort optimal decision values, distributions, and aggregates
                    [LAB, DIST, cohort] = solve_cohort(...
                        startyears(i), T_life, T_work, T_model_opt, T_model_dist, nz, nk, nb, idem, zs, transz, ks, bs, beta, gamma, sigma, surv, V_beq, mu2_idem, mu3_idem, ...
                        mpci, rpci, sstaxcredit, ssbenefits, sstaxs, ssincmaxs, deduc_coefs, pit_coefs, captaxshare, taucap, taucapgain, qtobin, qtobin0, ...
                        beqs, wages, capshares, debtshares, caprates, govrates, totrates, expsubsidys, ...
                        DIST0, LAB_static, DIST_static);
                    
                    % Store values
                    LABs   {i,idem} = LAB   ;
                    DISTs  {i,idem} = DIST  ;
                    cohorts{i,idem} = cohort;
                    
                end
                
                % Add cohort aggregates to total aggregates
                for i = 1:nstartyears %#ok<FXUP>
                    
                    switch economy
                        
                        case 'steady'
                            
                            kpr_total  = kpr_total  + sum(cohorts{i,idem}.kalive + cohorts{i,idem}.kdead);
                            beq_total  = beq_total  + sum(cohorts{i,idem}.kdead  );
                            elab_total = elab_total + sum(cohorts{i,idem}.labeff );
                            
                        case {'open', 'closed'}
                            
                            % Align aggregates to model years
                            T_shift = max(0, startyears(i));
                            T_dist  = size(DISTs{i,idem}, 4);
                            
                            kpr_total      (T_shift+(1:T_dist)) = kpr_total      (T_shift+(1:T_dist)) + cohorts{i,idem}.kalive + cohorts{i,idem}.kdead;
                            beq_total      (T_shift+(1:T_dist)) = beq_total      (T_shift+(1:T_dist)) + cohorts{i,idem}.kdead ;
                            elab_total     (T_shift+(1:T_dist)) = elab_total     (T_shift+(1:T_dist)) + cohorts{i,idem}.labeff;
                            lab_total      (T_shift+(1:T_dist)) = lab_total      (T_shift+(1:T_dist)) + cohorts{i,idem}.lab   ;
                            lfpr_total     (T_shift+(1:T_dist)) = lfpr_total     (T_shift+(1:T_dist)) + cohorts{i,idem}.lfpr  ;
                            fedincome_total(T_shift+(1:T_dist)) = fedincome_total(T_shift+(1:T_dist)) + cohorts{i,idem}.inc   ;
                            fedit_total    (T_shift+(1:T_dist)) = fedit_total    (T_shift+(1:T_dist)) + cohorts{i,idem}.pit   ;
                            ssrev_total    (T_shift+(1:T_dist)) = ssrev_total    (T_shift+(1:T_dist)) + cohorts{i,idem}.sst   ;
                            fcaptax_total  (T_shift+(1:T_dist)) = fcaptax_total  (T_shift+(1:T_dist)) + cohorts{i,idem}.cit   ;
                            ssexp_total    (T_shift+(1:T_dist)) = ssexp_total    (T_shift+(1:T_dist)) + cohorts{i,idem}.ben   ;
                            
                    end
                    
                end
                
            end
            
        end
        
        
        
        %% Static aggregate generation
        
        if ~isbase
            
            % Identify baseline generator and save directory
            base_generator = @() dynamicSolver.solve(economy, basedef, [], callingtag);
            base_dir = dirFinder.save(economy, basedef);
            
            % Load baseline solution
            s = hardyload('solution.mat'     , base_generator, base_dir);
            
            wages       = s.wages       ;
            capshares   = s.capshares   ;
            debtshares  = s.debtshares  ;
            caprates    = s.caprates    ;
            govrates    = s.govrates    ;
            expsubsidys = s.expsubsidys ;
            
            beqs        = zeros(1,T_model);
            totrates    = zeros(1,T_model);
            
            % Load optimal decision values and distributions from baseline
            s = hardyload('decisions.mat'    , base_generator, base_dir);
            LABs_static  = s.LABs ;
            
            s = hardyload('distributions.mat', base_generator, base_dir);
            DISTs_static = s.DISTs;
            
            
            % Clear parameter loading structure
            clear('s')
            
            
            % Generate static aggregates
            [~, ~, ~, ~, ~, ...
             fedincome_static, fedit_static, ssrev_static, fcaptax_static, ssexp_static, ...
             ~, ~] ...
             ...
               = generate_aggregates(beqs, wages, capshares, debtshares, caprates, govrates, totrates, expsubsidys, ...
                                     {}, LABs_static, DISTs_static); %#ok<ASGLU>
            
            
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
                
                DISTs_steady = {};
                
            case {'open', 'closed'}
                
                % Identify steady state generator and save directory
                steady_generator = @() dynamicSolver.steady(basedef, callingtag);
                steady_dir = dirFinder.save('steady', basedef);
                
                % Load steady state solution as starting solution
                s0 = hardyload('solution.mat'     , steady_generator, steady_dir);
                
                % Load steady state distributions
                s  = hardyload('distributions.mat', steady_generator, steady_dir);
                
                DISTs_steady = s.DISTs;
                
        end
        
        % Unpack starting solution
        rho0       = s0.rhos       ;
        beq0       = s0.beqs       ;
        kpr0       = s0.kprs       ;
        debt0      = s0.debts      ;
        capshare0  = s0.capshares  ;
        debtshare0 = s0.debtshares ;
        caprate0   = s0.caprates   ;
        govrate0   = s0.govrates   ;
        totrate0   = s0.totrates   ;
        
        clear('s0')
        
        
        switch economy
            
            case 'open'
                
                % Load government expenditure adjustment parameters
                s = load(fullfile(param_dir, 'param_gtilde.mat'));
                
                indtaxplan = cellfun(@(str) strncmp(taxplan, str, length(str)), {'base'; 'trump'; 'clinton'; 'ryan'});
                revenue_percent = s.revenue_percent(indtaxplan,:);
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
                        caps  = (kprs - debts)/qtobin;
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
                    
                    capshares  = (kprs - debts) ./ kprs;
                    debtshares = 1 - capshares;
                    caprates   = 1 + (A*alp*(rhos.^(alp-1)) - d)/qtobin;
                    switch economy
                        case 'steady', govrates = cbomeanrate;
                        case 'closed', govrates = cborates   ;
                    end
                    totrates   = capshares.*caprates + debtshares.*govrates;
                    
                    
                case 'open'
                    
                    if (iter == 1)
                        
                        kprs  = kpr0 *ones(1,T_model);
                        debts = debt0*ones(1,T_model);
                        caps  = (kprs - debts)/qtobin0;
                        
                        capshares  = capshare0 *ones(1,T_model);
                        debtshares = debtshare0*ones(1,T_model);
                        caprates   = (((1 - taucap_base)/(1 - taucap))*(caprate0 - 1) + 1)*ones(1,T_model);
                        govrates   = govrate0  *ones(1,T_model);
                        totrates   = totrate0  *ones(1,T_model);
                        
                        rhos  = ((qtobin*(caprates - 1) + d)/alp).^(1/(alp-1));
                        beqs  = beq0*ones(1,T_model);
                        
                    else
                        beqs = beq_total;
                        caps = cap_total;
                    end
                    
            end
            
            wages        = A*(1-alp)*(rhos.^alp);
            expsubsidys = [expshare * max(diff(caps), 0), 0] ./ caps;
            
            
            % Generate dynamic aggregates
            [kpr_total, beq_total, elab_total, lab_total, lfpr_total, ...
             fedincome_total, fedit_total, ssrev_total, fcaptax_total, ssexp_total, ...
             LABs, DISTs] ...
             ...
               = generate_aggregates(beqs, wages, capshares, debtshares, caprates, govrates, totrates, expsubsidys, ...
                                     DISTs_steady, {}, {}); %#ok<ASGLU>
            
            
            % Calculate additional dynamic aggregates
            % (Note that open economy requires capital calculation before debt calculation while closed economy requires the reverse)
            switch economy
                
                case 'steady'
                    
                    % Calculate debt
                    debt_total = debts;
                    
                    % Calculate capital and output
                    cap_total = (kpr_total - debt_total)/qtobin;
                    Y_total   = A*(max(cap_total, 0).^alp).*(elab_total.^(1-alp));
                    
                    % Calculate convergence value
                    rhoprs = (max(kpr_total - debt_total, 0)/qtobin) ./ elab_total;
                    delta = rhos - rhoprs;
                    
                case 'open'
                    
                    % Calculate capital and output
                    cap_total = rhos .* elab_total;
                    Y_total   = A*(max(cap_total, 0).^alp).*(elab_total.^(1-alp));
                    
                    domestic_cap_total = capshares .* [kpr0, kpr_total(1:end-1)];
                    foreign_cap_total  = qtobin*cap_total - domestic_cap_total;
                    
                    % Calculate debt
                    domestic_fcaptax_total = fcaptax_total;
                    foreign_fcaptax_total  = taucap.*(caprates - 1).*captaxshare.*foreign_cap_total;
                    fcaptax_total          = domestic_fcaptax_total + foreign_fcaptax_total;
                    
                    if isbase
                        Gtilde = (revenue_percent - fedgovtnis).*Y_total - ssexp_total;
                        Ttilde = revenue_percent.*Y_total - fedit_total - ssrev_total - fcaptax_total;
                    end
                    
                    netrev_total = fedit_total + ssrev_total + fcaptax_total - ssexp_total;
                    debt_total = [debt0, zeros(1,T_model-1)];
                    for t = 1:T_model-1
                        debt_total(t+1) = Gtilde(t) - Ttilde(t) - netrev_total(t) + debt_total(t)*cborates(t);
                    end
                    
                    domestic_debt_total = (1 - capshares) .* kpr_total;
                    foreign_debt_total  = debt_total - domestic_debt_total; %#ok<NASGU>
                    
                    % Calculate income
                    labinc_total = elab_total .* wages;
                    kinc_total   = qtobin * (caprates - 1) .* cap_total; %#ok<NASGU>
                    
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
                        debt_total(t+1) = Gtilde(t) - Ttilde(t) - netrev_total(t) + debt_total(t)*cborates(t);
                    end
                    
                    domestic_debt_total = debt_total;   %#ok<NASGU>
                    foreign_debt_total  = zeros(1,T_model); %#ok<NASGU>
                    
                    % Calculate capital and output
                    cap_total = ([(kpr0 - debt0)/qtobin0, (kpr_total(1:end-1) - debt_total(2:end))/qtobin]);
                    Y_total   = A*(max(cap_total, 0).^alp).*(elab_total.^(1-alp));
                    
                    domestic_cap_total = [qtobin0 * cap_total(1), qtobin * cap_total(2:end)]; %#ok<NASGU>
                    foreign_cap_total  = zeros(1,T_model); %#ok<NASGU>
                    
                    % Calculate income
                    labinc_total = elab_total .* wages;
                    kinc_total   = qtobin * (caprates - 1) .* cap_total; %#ok<NASGU>
                    
                    feditlab_total = fedit_total .* labinc_total ./ fedincome_total;
                    fcaprev_total  = fcaptax_total + fedit_total - feditlab_total; %#ok<NASGU>
                    
                    % Calculate convergence value
                    rhoprs = (max([kpr0, kpr_total(1:end-1)] - debt_total, 0)/qtobin) ./ elab_total;
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
            save(fullfile(save_dir, 'decisions.mat'    ), 'LABs' )
            save(fullfile(save_dir, 'distributions.mat'), 'DISTs')
        end
        
        % Save solution
        save(fullfile(save_dir, 'solution.mat'), ...
             'rhos', 'beqs', 'kprs', 'debts', 'caps', 'wages', ...
             'capshares', 'debtshares', 'caprates', 'govrates', 'totrates', 'expsubsidys')
        
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
                K_Y = (kpr_total - debt_total) / Y_total;
                
                
                % Calculate labor elasticity
                working_mass = 0;
                frisch_total = 0;
                
                for idem = 1:ndem %#ok<FXUP>
                    
                    LAB  = LABs {1,idem};
                    DIST = DISTs{1,idem};
                    
                    working_ind = (LAB > 0.01);
                    
                    working_mass = working_mass + sum(DIST(working_ind));
                    frisch_total = frisch_total + sum(DIST(working_ind).*(1-LAB(working_ind))./LAB(working_ind))*(1-gamma*(1-sigma))/sigma;
                    
                end
                
                labor_elas = frisch_total / working_mass;
                
                
                % Calculate savings elasticity
                deviation = 0.005;
                
                caprates_dev = caprates * (1 + deviation);
                govrates_dev = govrates * (1 + deviation);
                totrates_dev = totrates * (1 + deviation);
                
                [kpr_dev] = generate_aggregates(beqs, wages, capshares, debtshares, caprates_dev, govrates_dev, totrates_dev, expsubsidys, ...
                                                DISTs_steady, {}, {});
                
                savings_elas = ((kpr_dev - kpr_total)/kpr_total) / ((totrates_dev - totrates)/(totrates-1));
                
                
                % Save and display elasticities
                save(fullfile(save_dir, 'elasticities.mat'), 'K_Y', 'labor_elas', 'savings_elas')
                
                elasticities = { K_Y,          'Capital-to-output ratio' ;
                                 labor_elas,   'Labor elasticity'        ;
                                 savings_elas, 'Savings elasticity'      };
                for i = 1:size(elasticities, 1) %#ok<FXUP>
                    fprintf('\t%-25s= % 7.4f\n', elasticities{i, 2}, elasticities{i, 1})
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

