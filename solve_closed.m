%%
% Closed economy transition path solver.
% 
% Arguments:
% 
%   deep_params
%       beta, gamma, and sigma preference parameters collected in a 1 x 3 vector [beta, gamma, sigma].
% 
%   plan
%       Tax plan string identifier {'base', 'clinton', 'trump', 'ryan'}.
% 
%   gcut
%       Percentage government expenditure reduction, with positive values corresponding to reductions.
% 
%   showmore (optional | true by default)
%       Set to true for more solver status updates.
% 
%%


function [] = solve_closed(deep_params, plan, gcut, showmore)


%% Initialization

% Extract deep parameters or set defaults if none provided
if ~exist('deep_params', 'var') || isempty(deep_params)
    deep_params = [1.005, 0.390, 06.00];
end
beta  = deep_params(1);
gamma = deep_params(2);
sigma = deep_params(3);

% Set base plan as default
if ~exist('plan', 'var') || isempty(plan)
    plan = 'base';
end

% Set 0 reduction in government expenditure as default
if ~exist('gcut', 'var') || isempty(gcut)
    gcut = 0.00;
end

% Turn on all status updates by default
if ~exist('showmore', 'var') || isempty(showmore)
    showmore = true;
end



% Display problem definition
fprintf('\nSolving closed economy transition path:  beta = %0.3f  gamma = %0.3f  sigma = %05.2f  plan = %s  gcut = %+0.2f\n', beta, gamma, sigma, plan, gcut)

% Identify working directories
param_dir = dirFinder.param;
save_dir  = dirFinder.closed(deep_params, plan, gcut);

% Clear or create save directory
if exist(save_dir, 'dir')
    rmdir(save_dir, 's')
end
mkdir(save_dir)

% Identify reference steady state and open economy plan directories
ss_dir       = dirFinder.ss  (deep_params);
openbase_dir = dirFinder.open(deep_params, 'base');
openplan_dir = dirFinder.open(deep_params, plan  );



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


% Load CBO parameters
s = load(fullfile(param_dir, 'param_cbo.mat'));

rate_govs  = 1 + s.r_cbo(1:Tss);


% Load social security parameters
s = load(fullfile(param_dir, 'param_socsec.mat'));

NRA        = s.NRA;
ss_benefit = s.ss_benefit(:,1:Tss);
ss_tax     = s.ss_tax(1:Tss);
taxmax     = s.taxmax(1:Tss);


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
q_tobin0 = 1 - s_base.tau_cap * s_base.exp_share;
clear('s_base')


% Check for steady state solution and generate if not available
ss_solution = fullfile(ss_dir, 'solution.mat');
if ~exist(ss_solution, 'file')
    fprintf('\nSteady state solution not found.  Calling steady state solver.\n')
    
    uniquetag = sprintf('_closed_plan=%s_gcut=%+0.2f', plan, gcut);
    solve_ss(deep_params, showmore, false, uniquetag);    
    
    % Check for steady state solution again, possibly generated by a parallel run
    firstrun = false;
    if ~exist(ss_solution, 'file')
        % Attempt to copy to permanent directory and continue if error encountered
        % (Errors typically due to parallel copy attempt, which should result in at least one successful copy)
        try copyfile([ss_dir, uniquetag], ss_dir), firstrun = true; catch, end
    end
    
    fprintf('\nSteady state solution generated')
    if firstrun, fprintf('.'), else fprintf(' by a parallel run.'), end
    
    % Clean up temporary directory
    rmdir([ss_dir, uniquetag], 's')
    fprintf('  Proceeding with closed economy transition path solution.\n')
    
end

% Load steady state solution
s = hardyload(ss_solution);

rho_ss  = s.rho;
beq_ss  = s.beq;
kpr_ss  = s.kpr;
debt_ss = s.DEBTss;

rhos = rho_ss *ones(1,Tss);
beqs = beq_ss *ones(1,Tss);
KK   = kpr_ss *ones(1,Tss);
DD   = debt_ss*ones(1,Tss);

cap_series = (KK - DD)/q_tobin;


% Load steady state distributions
s = hardyload(fullfile(ss_dir, 'dist.mat'));

dist_ss = s.dist;


% Check for baseline open economy aggregates and generate if not available
openbase_aggregates = fullfile(openbase_dir, 'aggregates.mat');
if ~exist(openbase_aggregates, 'file')
    fprintf('\nOpen economy baseline aggregates not found.  Calling open economy solver.\n')
    
    uniquetag = sprintf('_closed_plan=%s_gcut=%+0.2f', plan, gcut);
    solve_open(deep_params, 'base', showmore, uniquetag);
    
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
    fprintf('  Proceeding with closed economy transition path solution.\n')
    
end

% Load baseline open economy aggregates
openbase_aggregates = fullfile(openbase_dir, 'aggregates.mat');
s = hardyload(openbase_aggregates);
Y_openbase      = s.Y_total;


% Check for open economy plan aggregates and generate if not available
openplan_aggregates = fullfile(openplan_dir, 'aggregates.mat');
if ~exist(openplan_aggregates, 'file')
    fprintf('\nOpen economy plan aggregates not found.  Calling open economy solver.\n')
    
    uniquetag = sprintf('_closed_plan=%s_gcut=%+0.2f', plan, gcut);
    solve_open(deep_params, plan, showmore, uniquetag);
    
    % Check for open economy plan aggregates again, possibly generated by a parallel run
    firstrun = false;
    if ~exist(openplan_aggregates, 'file')
        % Attempt to copy to permanent directory and continue if error encountered
        % (Errors typically due to parallel copy attempt, which should result in at least one successful copy)
        try copyfile([openplan_dir, uniquetag], openplan_dir), firstrun = true; catch, end
    end
    
    fprintf('\nOpen economy plan aggregates generated')
    if firstrun, fprintf('.'), else fprintf(' by a parallel run.'), end
    
    % Clean up temporary directory
    rmdir([openplan_dir, uniquetag], 's')
    fprintf('  Proceeding with closed economy transition path solution.\n')
    
end

% Load open economy plan aggregates
openplan_aggregates = fullfile(openplan_dir, 'aggregates.mat');
s = hardyload(openplan_aggregates);
Ttilde          = s.Ttilde;
Gtilde          = s.Gtilde;

% Generate government expenditure adjustments
s = load(fullfile(param_dir, 'param_gtilde.mat'));
GEXP_percent = s.GEXP_percent(1:Tss);
GEXP_cut     = gcut*GEXP_percent.*Y_openbase;
Gtilde       = Gtilde - GEXP_cut;




% Clear parameter loading structure
clear('s')



%% Closed economy transition path calculation

% Define tolerance for rho convergence and initialize collection of rho error terms
rho_tol = 1e-3;
rho_eps = [];

% Initialize iteration count and set maximum number of iterations
rho_iter    =  0;
rho_itermax = 25;
if showmore, fprintf('\n'), end

while true
    
    % Increment iteration count
    rho_iter = rho_iter + 1;
    if showmore, fprintf('rho_iter = %2d\n', rho_iter), end
    
    % Derive values
    wages        = A*(1-alp)*(rhos.^alp);
    exp_subsidys = [exp_share*max(diff(cap_series),0), 0]./cap_series;
    
    cap_shares  = (KK - DD) ./ KK;
    debt_shares = 1 - cap_shares;
    
    rate_caps   = 1 + (A*alp*(rhos.^(alp-1)) - d)/q_tobin;
    
    rate_tots   = cap_shares.*rate_caps + debt_shares.*rate_govs;
    
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
            
            if showmore, fprintf('\tidem = %02d\tbirthyear = %+03d\n', idem, birthyear), end
            
            
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
    
    % Find new debt sequence
    netrev_total = fedit_total + ssrev_total + fcaptax_total - ssexp_total;
    
    debt_total = [debt_ss, zeros(1,Tss-1)];
    for t = 1:Tss-1
        debt_total(t+1) = Gtilde(t) - Ttilde(t) - netrev_total(t) + debt_total(t)*rate_govs(t);
    end
    
    cap_total = ([(kpr_ss-debt_ss)/q_tobin0 (kpr_total(1:end-1) - debt_total(2:end))/q_tobin]);
    
    
    % Check for convergence
    rhoprs = (max([kpr_ss, kpr_total(1:end-1)] - debt_total, 0)/q_tobin) ./ elab_total;
    
    rho_eps = [rho_eps, max(abs(rhos - rhoprs))]; %#ok<AGROW>
    if showmore, fprintf('rho_eps = %7.4f\n\n', rho_eps(end)), end
    if (rho_eps(end) < rho_tol), break, end
    
    % Check for maximum iterations
    if (rho_iter == rho_itermax)
        warning('Maximum iterations reached.')
        break
    end
    
    % Update variables for next iteration
    stepfactor = 0.3;
    rhos = (1-stepfactor) * rhos + stepfactor * rhoprs;
    
    beqs       = beq_total;
    KK         = kpr_total;
    DD         = debt_total;
    cap_series = cap_total;
    
end


% Calculate additional aggregates
Y_total   = A*(max(cap_total, 0).^alp).*(elab_total.^(1-alp)); %#ok<NASGU>

domestic_cap_total  = [cap_total(1)*q_tobin0 q_tobin*cap_total(2:end)]; %#ok<NASGU>
foreign_cap_total   = zeros(1,Tss);

domestic_debt_total = debt_total; 
foreign_debt_total  = debt_total - domestic_debt_total; %#ok<NASGU>

domestic_fcaptax_total = fcaptax_total;
foreign_fcaptax_total  = tau_cap.*(rate_caps - 1).*cap_tax_share.*foreign_cap_total;
fcaptax_total          = domestic_fcaptax_total + foreign_fcaptax_total; 

labinc_total = elab_total .* wages;
kinc_total   = (rate_caps - 1) .* (kpr_total - debt_total); %#ok<NASGU>

feditlab_total = fedit_total .* labinc_total ./ fedincome_total; %#ok<NASGU>
fcaprev_total  = fcaptax_total + fedit_total .* (1 - labinc_total ./ fedincome_total); %#ok<NASGU>


% Save optimal decision values and distributions
save(fullfile(save_dir, 'opt.mat'),  'opt' )
save(fullfile(save_dir, 'dist.mat'), 'dist')

% Save solution
save(fullfile(save_dir, 'solution.mat'), ...
     'rhos', 'wages', 'cap_shares', 'debt_shares', 'rate_caps', 'rate_govs')

% Save aggregates
save(fullfile(save_dir, 'aggregates.mat'), ...
     'kpr_total', 'debt_total', 'cap_total', ...
     'domestic_cap_total', 'foreign_cap_total', ...
     'domestic_debt_total', 'foreign_debt_total', ...
     'domestic_fcaptax_total', 'foreign_fcaptax_total', ...
     'beq_total', 'elab_total', 'lab_total', 'lfpr_total', ...
     'fedincome_total', 'fedit_total', 'ssrev_total', 'fcaptax_total', 'ssexp_total', ...
     'Y_total', 'labinc_total', 'kinc_total', 'feditlab_total', 'fcaprev_total', 'exp_subsidys', ...
     'Gtilde', 'Ttilde', 'GEXP_cut')

% Save log of rho iterations
fid = fopen(fullfile(save_dir, 'iterations.txt'), 'wt');
fprintf(fid, '%d rho iteration(s) of a maximum %d\n\nError term(s):', rho_iter, rho_itermax);
for i = 1:rho_iter
    fprintf(fid, '\n  %2d  --  %7.4f', i, rho_eps(i));
end
fclose(fid);


end