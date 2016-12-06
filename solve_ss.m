%%
% Steady state solver.
% 
% Arguments:
% 
%   deep_params
%       beta, gamma, and sigma preference parameters collected in a 1 x 3 vector [beta, gamma, sigma].
% 
%   callertag (optional | '' by default)
%       String used to save solution into a unique directory, used to avoid conflicts between parallel runs.
% 
% Outputs:
% 
%   save_dir
%       Directory where results are saved.
% 
%   elasticities
%       Capital-to-output ratio, labor elasticity, and savings elasticity collected in a 1 x 3 vector [K_Y, labor_elas, savings_elas].
% 
%%


function [save_dir, elasticities] = solve_ss(basedef, callertag)


%% Initialization

economy = 'steady';
if ~exist('counterdef', 'var'), counterdef = []; end
if ~exist('callertag' , 'var'), callertag  = ''; end

% Identify working directories
param_dir = dirFinder.param;
[save_dir, callingtag] = dirFinder.save('steady', basedef);

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




% Load global parameters
s = load(fullfile(param_dir, 'param_global.mat'));

T_life  = s.T_life;
T_model = s.T_model; switch economy, case 'steady', T_model = 1; end

kgrid = s.kgrid;
bgrid = [0; s.bgrid(2:end)];

Vbeq = s.phi1.*((1+kgrid./s.phi2).^(1-s.phi3));

ndem = s.ndem;
nb   = s.nb;
nk   = s.nk;
nz   = s.nz;

mpci = s.mpci;
rpci = s.rpci;

ss_tax_cred   = s.ss_tax_cred;

A           = s.A;
alp         = s.alp;
d           = s.d;
deduc_scale = s.deduc_scale;
proddist    = s.proddist(:,1)';
surv        = [s.surv(1:T_life-1), 0];
tr_z        = s.tr_z;
z           = s.z;

MU2 = s.demdist_2015 * ( s.Mu2/sum(s.Mu2) );
MU3 = repmat(1-surv, [ndem,1]) .* MU2;


% Load social security parameters
s = load(fullfile(param_dir, 'param_socsec.mat'));

NRA        = s.NRA;
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


% Clear parameter loading structure
clear('s')




%% Dynamic aggregate generation

switch economy
    
    case 'steady'
        
        % Load initial values for steady state solution
        s = load(fullfile(param_dir, 'solution0.mat'));
        
        rho0  = s.rho   ;
        beq0  = s.beq   ;
        kpr0  = s.kpr   ;
        debt0 = s.DEBTss;
        
end


% Clear parameter loading structure
clear('s')



% Start a parallel pool if one does not already exist and JVM is enabled
if usejava('jvm'), gcp; end

% Define tolerance for rho convergence and initialize collection of rho error terms
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
        
end
fprintf(']\n')


while ( eps > tol && iter < itermax )
    
    % Increment iteration count
    iter = iter + 1;
    fprintf('\tIteration %2d  ...  ', iter)
    
    
    % Define prices
    switch economy
        
        case 'steady'
            
            if (iter == 1)
                rhos  = rho0 *ones(1,T_model);
                beqs  = beq0 *ones(1,T_model);
                kprs  = kpr0 *ones(1,T_model);
                debts = debt0*ones(1,T_model);
                caps  = (kprs - debts)/q_tobin;
            else
                rhos  = 0.5*rhos + 0.5*rhoprs;
                beqs  = beq_total;
                kprs  = kpr_total;
                debts = D_Y*Y_total;
                caps  = cap_total;
            end
            
            wages = A*(1-alp)*(rhos.^alp);
            
            cap_shares  = (kprs - debts) ./ kprs;
            debt_shares = 1 - cap_shares;
            rate_caps   = 1 + (A*alp*(rhos.^(alp-1)) - d)/q_tobin;
            rate_govs   = meanrate_cbo;
            rate_tots   = cap_shares.*rate_caps + debt_shares.*rate_govs;
            
            exp_subsidys = [exp_share * max(diff(caps), 0), 0] ./ caps;
            
    end
    
    
    % Initialize global aggregates
    kpr_total  = zeros(1,T_model);
    beq_total  = zeros(1,T_model);
    elab_total = zeros(1,T_model);
    
    % Define birth year range
    switch economy
        
        case 'steady'
            birthyears = 0;
        
    end
    
    % Initialize storage structures for optimal decision values and distributions
    s_birthyear = struct('birthyear', num2cell(repmat(birthyears', [1,ndem])));
    opt     = s_birthyear;
    dist    = s_birthyear;
    cohorts = s_birthyear;
    
    % Define dynamic optimization and distribution generation time intervals
    switch economy
        
        case 'steady'
            T_opt  = 1;
            T_dist = T_life;
            
    end
    
    for idem = 1:ndem
        
        % Extract demographic adjustments
        mu2_idem = MU2(idem,:);
        mu3_idem = MU3(idem,:);
        
        % Define initial distributions
        switch economy
            
            case 'steady'
                dist_w0 = padarray(proddist, [nk-1, 0, nb-1], 0, 'post');
                dist_r0 = [];
                
        end
        
        parfor i = 1:length(birthyears)
            
            % Get birth year and retirement age
            birthyear = birthyears(i);
            Tr = NRA(i);
            
            % Generate optimal decision values, distributions, and cohort aggregates
            [labopt, dist_w, dist_r, N_w, N_r, Kalive, Kdead, ELab, Lab, Lfpr, Fedincome, Fedit, SSrev, Fcaptax, SSexp] ...
             ...
               = generate_distributions(...
                   beta, gamma, sigma, T_life, Tr, T_opt, T_dist, birthyear, ...
                   kgrid, z, tr_z, bgrid, nk, nz, nb, idem, ...
                   mpci, rpci, cap_tax_share, ss_tax_cred, surv, tau_cap, tau_capgain, ss_tax, taxmax, ...
                   beqs, wages, cap_shares, debt_shares, rate_caps, rate_govs, rate_tots, exp_subsidys, q_tobin, q_tobin, Vbeq, ...
                   avg_deduc, clinton, coefs, limit, X, ss_benefit, ...
                   dist_w0, dist_r0, mu2_idem, mu3_idem);
            
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
            
            switch economy
                
                case 'steady'
                    
                    % Add aggregate values to global aggregates
                    kpr_total  = kpr_total  + sum(cohorts(i,idem).Kalive + cohorts(i,idem).Kdead);
                    beq_total  = beq_total  + sum(cohorts(i,idem).Kdead);
                    elab_total = elab_total + sum(cohorts(i,idem).ELab);
                    
            end
            
        end
        
    end
    
    
    % Calculate additional dynamic aggregates
    switch economy
        
        case 'steady'
            
            % Calculate debt
            debt_total = debts;
            
            % Calculate capital and output
            cap_total = (kpr_total - debt_total)/q_tobin;
            Y_total   = A*(max(cap_total, 0).^alp).*(elab_total.^(1-alp));
            
            % Calculate convergence delta
            rhoprs = (max(kpr_total - debt_total, 0)/q_tobin) ./ elab_total;
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


 
%% Elasticity calculation

switch economy
    case 'steady'
        
        % Calculate capital to output ratio
        K_Y = (kpr_total-debt_total)/Y_total;
        
        
        % Calculate labor elasticity
        working_mass = 0;
        frisch_total = 0;
        
        for idem = 1:ndem
            
            labopt = opt(1,idem).lab;
            dist_w = dist(1,idem).w;
            
            working_ind = (labopt > 0.01);  % (Labor elasticity calculation sensitive to threshold value)
            
            working_mass = working_mass + sum(dist_w(working_ind));
            frisch_total = frisch_total + sum(dist_w(working_ind).*(1-labopt(working_ind))./labopt(working_ind))*(1-gamma*(1-sigma))/sigma;
            
        end
        
        labor_elas = frisch_total / working_mass;
        
        
        % Calculate savings elasticity
        
        % % Specify percentage change in net interest rate
        % rate_deviation = 0.005;
        % 
        % % Scale rates by deviation if provided
        % rate_cap = rate_cap * (1 + rate_deviation);
        % rate_gov = rate_gov * (1 + rate_deviation);
        % 
        % % Solve dynamic optimization and generate distributions
        % [~, ~, kpr_dev, ~, ~, ~, ~, ~, ~, ~, ~, rate_tot_dev] ...
        %  ...
        %    = solve_and_generate(...
        %        rho, beq, kpr, debt, ...
        %        T, Tr, bgrid, kgrid, ndem, nb, nk, nz, ...
        %        mpci, rpci, A, alp, cap_tax_share, d, q_tobin, Vbeq, proddist, ss_tax_cred, surv, ...
        %        tau_cap, tau_capgain, tau_ss, tr_z, v_ss_max, z, ...
        %        avg_deduc, clinton, coefs, limit, X, r_cbo, ben, ...
        %        beta, gamma, sigma, MU2, MU3, ...
        %        rate_deviation);
        % 
        % % Calculate savings elasticity
        % savings_elas = ((kpr_dev - kpr)/kpr) / ((rate_tot_dev - rate_tot)/(rate_tot-1));
        
        savings_elas = 0;
        
        
        % Package, save, and display elasticities
        elasticities = [K_Y, labor_elas, savings_elas];
        
        save(fullfile(save_dir, 'elasticities.mat'), 'K_Y', 'labor_elas', 'savings_elas')
        
        displaynames = { 'Capital-to-output ratio', 'Labor elasticity', 'Savings elasticity' };
        for i = 1:length(elasticities)
            fprintf('\t%-25s= % 7.4f\n', displaynames{i}, elasticities(i))
        end
        fprintf('\n')
        
end


end