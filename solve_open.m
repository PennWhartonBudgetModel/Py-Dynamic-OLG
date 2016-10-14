%%
% Open economy transition path solver.
% 
% Arguments:
% 
%   deep_params
%       beta, gamma, and sigma preference parameters collected in a 1 x 3 vector [beta, gamma, sigma].
% 
%   plan
%       Tax plan string identifier {'base', 'clinton', 'trump', 'ryan'}.
% 
%   showmore (optional | true by default)
%       Set to true for more solver status updates.
% 
%   this_uniquetag (optional | '' by default)
%       String used to save solution into a unique directory, used to avoid conflicts between parallel runs.
% 
%%


function [] = solve_open(deep_params, plan, showmore, this_uniquetag)


%% Initialization

% Extract deep parameters or set defaults if none provided
if exist('deep_params', 'var')
    beta  = deep_params(1);
    gamma = deep_params(2);
    sigma = deep_params(3);
else
    beta  = 1.005;
    gamma = 0.390;
    sigma = 06.00;
    deep_params = [beta, gamma, sigma];
end

% Set base plan as default
if ~exist('plan', 'var')
    plan = 'base';
end

% Turn on all status updates by default
if ~exist('showmore', 'var')
    showmore = true;
end

% Set solution uniqueness tag to empty by default
if ~exist('this_uniquetag', 'var')
    this_uniquetag = '';
end



% Display problem definition
fprintf('\nSolving open economy transition path:  beta = %0.3f  gamma = %0.3f  sigma = %05.2f  plan = %s\n', beta, gamma, sigma, plan)

% Identify working directories
param_dir = dirFinder.param;
save_dir  = dirFinder.open(beta, gamma, sigma, plan);

% Append uniqueness tag to name of save directory
save_dir = [save_dir, this_uniquetag];

% Clear or create save directory
if exist(save_dir, 'dir')
    rmdir(save_dir, 's')
end
mkdir(save_dir)

% Identify reference steady state and baseline open economy directories
ss_dir   = dirFinder.ss  (beta, gamma, sigma);
base_dir = dirFinder.open(beta, gamma, sigma, 'base');



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


% Check for steady state solution and generate if not available
ss_solution = fullfile(ss_dir, 'solution.mat');
if ~exist(ss_solution, 'file')
    fprintf('\nSteady state solution not found.  Calling steady state solver.\n')
    
    uniquetag = sprintf('_open_plan=%s', plan);
    solve_ss([beta, gamma, sigma], showmore, false, uniquetag);
    
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
    fprintf('  Proceeding with open economy transition path solution.\n')
    
end

% Load steady state solution
s = hardyload(ss_solution);

beqs        = s.beq       *ones(1,Tss);
cap_shares  = s.cap_share *ones(1,Tss);
debt_shares = s.debt_share*ones(1,Tss);
rate_govs   = s.rate_gov  *ones(1,Tss);
rate_tots   = s.rate_tot  *ones(1,Tss);
cap_incs    = s.cap_inc   *ones(1,Tss);

rate_cap = ((1 - tau_cap_base)/(1 - tau_cap))*(s.rate_cap - 1) + 1;
rate_caps = rate_cap*ones(1,Tss);

rho   = ((q_tobin*(rate_cap - 1) + d)/alp)^(1/(alp-1));

wage = A*(1-alp)*rho^alp;
wages = wage*ones(1,Tss);

kpr_ss  = s.kpr;
debt_ss = s.DEBTss;

cap_series = (kpr_ss - debt_ss)*ones(1,Tss)/q_tobin0;



% Load steady state distributions
s = hardyload(fullfile(ss_dir, 'dist.mat'));

dist_ss = s.dist;


% Check for baseline open economy aggregates and generate if not available
if ~strcmp(plan, 'base')
    
    base_aggregates = fullfile(base_dir, 'aggregates.mat');
    if ~exist(base_aggregates, 'file')
        fprintf('\nBaseline open economy aggregates not found.  Calling open economy solver.\n')

        uniquetag = sprintf('_open_plan=%s', plan);
        solve_open([beta, gamma, sigma], 'base', showmore, uniquetag);

        % Check for baseline open economy aggregates again, possibly generated by a parallel run
        firstrun = false;
        if ~exist(base_aggregates, 'file')
            % Attempt to copy to permanent directory and continue if error encountered
            % (Errors typically due to parallel copy attempt, which should result in at least one successful copy)
            try copyfile([base_dir, uniquetag], base_dir), firstrun = true; catch, end
        end

        fprintf('\nBaseline open economy aggregates generated')
        if firstrun, fprintf('.'), else fprintf(' by a parallel run.'), end

        % Clean up temporary directory
        rmdir([base_dir, uniquetag], 's')
        fprintf('  Proceeding with open economy transition path solution.\n')

    end
    
end


% Load government expenditure adjustment parameters
s = load(fullfile(param_dir, 'param_gtilde.mat'));

indplan = cellfun(@(s) strncmp(plan, s, length(s)), {'base'; 'trump'; 'clinton'; 'ryan'});
revenue_percent = s.revenue_percent(indplan,:);

% Repeat last value if vector is too small
if length(revenue_percent)<Tss
    revenue_percent = [revenue_percent, revenue_percent(end).*ones(1,Tss-length(revenue_percent))];
end


% Clear parameter loading structure
clear('s')



%% Open economy transition path calculation

% Define tolerance for beq convergence and initialize collection of beq error terms
beq_tol = 1e-3;
beq_eps = [];

% Initialize iteration count and set maximum number of iterations
beq_iter    =  0;
beq_itermax = 25;
if showmore, fprintf('\n'), end

while true
    
    % Increment iteration count
    beq_iter = beq_iter + 1;
    if showmore, fprintf('beq_iter = %2d\n', beq_iter), end
    
    % Derive values
    exp_subsidys = [exp_share*max(diff(cap_series),0), 0]./cap_series;
    
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
        
        
        for birthyear = (-T+1):(Tss-1)
            
            if showmore, fprintf('\tidem = %02d\tbirthyear = %+03d\n', idem, birthyear), end
            
            % Get retirement age
            Tr = NRA(birthyear+T);
            
            
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
                   wages, cap_shares, debt_shares, rate_caps, rate_govs, rate_tots, cap_incs, exp_subsidys, q_tobin, q_tobin0, Vbeq, ...
                   avg_deduc, clinton, coefs, limit, X, ss_benefit, ...
                   beta, gamma, sigma);
            
            % Check for unsolved dynamic optimization subproblems
            if any(isinf([V(:); Vss(:)]))
                warning('Infinite utility values found.  Some dynamic optimization subproblems unsolved.  Check that initial conditions satisfy constraints.')
            end
            
            % Package optimal decision values
            opt(birthyear+T,idem).birthyear   = birthyear;    %#ok<AGROW>
            opt(birthyear+T,idem).V           = V;            %#ok<AGROW>
            opt(birthyear+T,idem).Vss         = Vss;          %#ok<AGROW>
            opt(birthyear+T,idem).kopt        = kopt;         %#ok<AGROW>
            opt(birthyear+T,idem).koptss      = koptss;       %#ok<AGROW>
            opt(birthyear+T,idem).labopt      = labopt;       %#ok<AGROW>
            opt(birthyear+T,idem).bopt        = bopt;         %#ok<AGROW>
            opt(birthyear+T,idem).fedincome   = fedincome;    %#ok<AGROW>
            opt(birthyear+T,idem).fedincomess = fedincomess;  %#ok<AGROW>
            opt(birthyear+T,idem).fcaptax     = fcaptax;      %#ok<AGROW>
            opt(birthyear+T,idem).fcaptaxss   = fcaptaxss;    %#ok<AGROW>
            opt(birthyear+T,idem).fitax       = fitax;        %#ok<AGROW>
            opt(birthyear+T,idem).fitaxss     = fitaxss;      %#ok<AGROW>
            opt(birthyear+T,idem).fsstax      = fsstax;       %#ok<AGROW>
            opt(birthyear+T,idem).fsstaxss    = fsstaxss;     %#ok<AGROW>
            
            
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
            
            % Package distributions
            dist(birthyear+T,idem).birthyear = birthyear;  %#ok<AGROW>
            dist(birthyear+T,idem).dist_w    = dist_w;     %#ok<AGROW>
            dist(birthyear+T,idem).dist_r    = dist_r;     %#ok<AGROW>
            dist(birthyear+T,idem).Kalive    = Kalive;     %#ok<AGROW>
            dist(birthyear+T,idem).Kdead     = Kdead;      %#ok<AGROW>
            dist(birthyear+T,idem).ELab      = ELab;       %#ok<AGROW>
            dist(birthyear+T,idem).Lab       = Lab;        %#ok<AGROW>
            dist(birthyear+T,idem).Lfpr      = Lfpr;       %#ok<AGROW>
            dist(birthyear+T,idem).Fedincome = Fedincome;  %#ok<AGROW>
            dist(birthyear+T,idem).Fedit     = Fedit;      %#ok<AGROW>
            dist(birthyear+T,idem).SSrev     = SSrev;      %#ok<AGROW>
            dist(birthyear+T,idem).Fcaptax   = Fcaptax;    %#ok<AGROW>
            dist(birthyear+T,idem).SSexp     = SSexp;      %#ok<AGROW>
            
            % Add aggregate values to global aggregates, aligning to projection years
            T_shift = max(0, birthyear);
            kpr_total      ( T_shift + (1:N_w+N_r) ) = kpr_total      ( T_shift + (1:N_w+N_r) ) + Kalive + Kdead;
            beq_total      ( T_shift + (1:N_w+N_r) ) = beq_total      ( T_shift + (1:N_w+N_r) ) + Kdead;
            elab_total     ( T_shift + (1:N_w)     ) = elab_total     ( T_shift + (1:N_w)     ) + ELab;
            lab_total      ( T_shift + (1:N_w)     ) = lab_total      ( T_shift + (1:N_w)     ) + Lab;
            lfpr_total     ( T_shift + (1:N_w)     ) = lfpr_total     ( T_shift + (1:N_w)     ) + Lfpr;
            fedincome_total( T_shift + (1:N_w+N_r) ) = fedincome_total( T_shift + (1:N_w+N_r) ) + Fedincome;
            fedit_total    ( T_shift + (1:N_w+N_r) ) = fedit_total    ( T_shift + (1:N_w+N_r) ) + Fedit;
            ssrev_total    ( T_shift + (1:N_w+N_r) ) = ssrev_total    ( T_shift + (1:N_w+N_r) ) + SSrev;
            fcaptax_total  ( T_shift + (1:N_w+N_r) ) = fcaptax_total  ( T_shift + (1:N_w+N_r) ) + Fcaptax;
            ssexp_total    ( T_shift + N_w+(1:N_r) ) = ssexp_total    ( T_shift + N_w+(1:N_r) ) + SSexp;
            
        end
        
    end
    
    cap_total = rho * elab_total;
    
    
    % Check for convergence
    beq_eps = [beq_eps, max(abs(beqs - beq_total))]; %#ok<AGROW>
    if showmore, fprintf('beq_eps = %7.4f\n\n', beq_eps(end)), end
    if (beq_eps(end) < beq_tol), break, end
    
    % Check for maximum iterations
    if (beq_iter == beq_itermax)
        warning('Maximum iterations reached.')
        break
    end
    
    % Update bequests
    stepfactor = 1;
    beqs = (1-stepfactor) * beqs + stepfactor * beq_total;
    
    cap_series = cap_total;
    
end

% Calculate additional aggregates
Y_total   = A*((cap_total).^alp).*(elab_total.^(1-alp));  % No need to convert into capital units because rho above already in capital units

domestic_cap_total  = cap_shares .* [kpr_ss kpr_total(1:end-1)];
foreign_cap_total   = q_tobin*cap_total - domestic_cap_total;

domestic_debt_total    = (1 - cap_shares) .* kpr_total;

domestic_fcaptax_total = fcaptax_total;
foreign_fcaptax_total  = tau_cap.*(rate_caps - 1).*cap_tax_share.*foreign_cap_total;
fcaptax_total          = domestic_fcaptax_total + foreign_fcaptax_total;

labinc_total = elab_total .* wages;
kinc_total   = q_tobin*(rate_caps - 1) .* cap_total; %#ok<NASGU>

feditlab_total = fedit_total .* labinc_total ./ fedincome_total; %#ok<NASGU>
fcaprev_total  = fcaptax_total + fedit_total .* (1 - labinc_total ./ fedincome_total); %#ok<NASGU>


% Save optimal decision values and distributions
save(fullfile(save_dir, 'opt.mat'),  'opt' )
save(fullfile(save_dir, 'dist.mat'), 'dist')


% Save aggregates
save(fullfile(save_dir, 'aggregates.mat'), ...
     'kpr_total', 'cap_total', ...
     'domestic_cap_total', 'foreign_cap_total', ...
     'domestic_debt_total', ...
     'domestic_fcaptax_total', 'foreign_fcaptax_total', ...
     'beq_total', 'elab_total', 'lab_total', 'lfpr_total', ...
     'fedincome_total', 'fedit_total', 'ssrev_total', 'fcaptax_total', 'ssexp_total', ...
     'Y_total', 'labinc_total', 'kinc_total', 'feditlab_total', 'fcaprev_total', 'exp_subsidys')
 
 


% Calculate static aggregates
generate_static_aggregates(deep_params, plan, [], this_uniquetag);

% Extract aggregates from folder
s_static = load(fullfile(save_dir, 'aggregates_static.mat'));

fedit_static   = s_static.fedit_static;
ssrev_static   = s_static.ssrev_static;
fcaptax_static = s_static.fcaptax_static;

if strcmp(plan, 'base')
    Y_base = Y_total;
else
    s = hardyload(base_aggregates);
    Y_base = s.Y_total;
end

Ttilde = revenue_percent.*Y_base - fedit_static - ssrev_static - fcaptax_static;

if strcmp(plan,'base')
    Gtilde = fedit_total + ssrev_total + fcaptax_total + Ttilde - ssexp_total - fedgovtnis.*Y_total;
else
    s = hardyload(base_aggregates);
    Gtilde = s.Gtilde;
end


% Calculate aggregate debt
netrev_total = fedit_total + ssrev_total + fcaptax_total - ssexp_total;

debt_total = [debt_ss, zeros(1,Tss-1)];
for t = 1:Tss-1
    debt_total(t+1) = Gtilde(t) - Ttilde(t) - netrev_total(t) + debt_total(t)*rate_govs_cbo(t);
end

% Solving foreign debt and appending to static aggregates because of dependency conflict
foreign_debt_total = debt_total - domestic_debt_total; %#ok<NASGU>
foreign_debt_static = foreign_cap_total; %#ok<NASGU>

% Append to static aggregates because of dependency issue
save(fullfile(save_dir, 'aggregates_static.mat'), 'foreign_debt_static', '-append');

% Append to results because of dependency of statics on other results
save(fullfile(save_dir, 'aggregates.mat'), 'foreign_debt_total', 'Gtilde', 'Ttilde', 'debt_total', '-append')

% Save log of beq iterations
fid = fopen(fullfile(save_dir, 'iterations.txt'), 'wt');
fprintf(fid, '%d beq iteration(s) of a maximum %d\n\nError term(s):', beq_iter, beq_itermax);
for i = 1:beq_iter
    fprintf(fid, '\n  %2d  --  %7.4f', i, beq_eps(i));
end
fclose(fid);


end