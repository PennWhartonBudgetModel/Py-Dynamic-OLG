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

ss_tax_cred   = s.ss_tax_cred;

A           = s.A;
alp         = s.alp;
d           = s.d;
deduc_scale = s.deduc_scale;
proddist    = s.proddist(:,1);
surv        = [s.surv(1:T-1), 0];       % (Survival probability at age T set to 0)
tr_z        = s.tr_z;
z           = s.z;

MU2 = s.demdist_2015 * ( s.Mu2/sum(s.Mu2) );
MU3 = repmat(1-surv, [ndem,1]) .* MU2;


% Load social security parameters
s = load(fullfile(param_dir, 'param_socsec.mat'));

Tr       = s.NRA(1);
ben      = s.ss_benefit(:,1);
tau_ss   = s.ss_tax(1);
v_ss_max = s.taxmax(1);


% Load CBO parameters
s = load(fullfile(param_dir, 'param_cbo.mat'));

r_cbo = mean(s.r_cbo);


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


% Load initial values for steady state solution
s = load(fullfile(param_dir, 'solution0.mat'));

rho0  = s.rho;
beq0  = s.beq;
kpr0  = s.kpr;
debt0 = s.DEBTss;


% Clear parameter loading structure
clear('s')




%% Steady state calculation

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
fprintf('\n[Steady state]\n')

while ( eps > tol && iter < itermax )
    
    % Increment iteration count
    iter = iter + 1;
    fprintf('\tIteration %2d  ...  ', iter)
    
    
    % Define prices
    switch economy
        
        case 'steady'
            
            if (iter == 1)
                rho  = rho0 ;
                beq  = beq0 ;
                kpr  = kpr0 ;
                debt = debt0;
            else
                rho  = 0.5*rho + 0.5*rhopr;
                beq  = beq_total;
                kpr  = kpr_total;
                debt = 0.74*Y;
            end
        
    end
    
    
    % Solve dynamic optimization and generate distributions
    [opt, dist, kpr_total, debt_total, elab_total, beq_total, wage, cap_share, debt_share, rate_cap, rate_gov, rate_tot] ...
     ...
       = solve_and_generate(...
           rho, beq, kpr, debt, ...
           T, Tr, Tss, bgrid, kgrid, ndem, nb, nk, nz, ...
           mpci, rpci, A, alp, cap_tax_share, d, q_tobin, Vbeq, proddist, ss_tax_cred, surv, ...
           tau_cap, tau_capgain, tau_ss, tr_z, v_ss_max, z, ...
           avg_deduc, clinton, coefs, limit, X, r_cbo, ben, ...
           beta, gamma, sigma, MU2, MU3); %#ok<ASGLU>
    
    
    % Calculate additional dynamic aggregates
    Y = A*(max((kpr_total-debt_total)/q_tobin, 0)^alp)*(elab_total^(1-alp));
    
    % Calculate convergence series delta
    rhopr = (max(kpr_total-debt_total, 0)/q_tobin)/elab_total;
    
    % Calculate convergence error
    eps = abs(rho - rhopr);
    
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
     'rho', 'beq', 'kpr', 'debt', ...
     'cap_share', 'debt_share', 'rate_cap', 'rate_gov', 'rate_tot')


% Calculate capital to output ratio
K_Y = (kpr_total-debt_total)/Y;




%% Labor elasticity calculation

working_mass = 0;
frisch_total = 0;

for idem = 1:ndem
    
    labopt = opt(1,idem).lab;
    dist_w = dist(1,idem).w;
    
    working_ind = (labopt > 0.01);  % (Labor elasticity calculation sensitive to threshold value)
    
    working_mass = working_mass + sum(dist_w(working_ind));
    frisch_total = frisch_total + sum(dist_w(working_ind).*(1-labopt(working_ind))./labopt(working_ind))*(1-gamma*(1-sigma))/sigma;
    
end

% Calculate labor elasticity
labor_elas = frisch_total / working_mass;




%% Savings elasticity calculation

% Specify percentage change in net interest rate
rate_deviation = 0.005;

% Solve dynamic optimization and generate distributions
[~, ~, kpr_dev, ~, ~, ~, ~, ~, ~, ~, ~, rate_tot_dev] ...
 ...
   = solve_and_generate(...
       rho, beq, kpr, debt, ...
       T, Tr, Tss, bgrid, kgrid, ndem, nb, nk, nz, ...
       mpci, rpci, A, alp, cap_tax_share, d, q_tobin, Vbeq, proddist, ss_tax_cred, surv, ...
       tau_cap, tau_capgain, tau_ss, tr_z, v_ss_max, z, ...
       avg_deduc, clinton, coefs, limit, X, r_cbo, ben, ...
       beta, gamma, sigma, MU2, MU3, ...
       rate_deviation);

% Calculate savings elasticity
savings_elas = ((kpr_dev - kpr)/kpr) / ((rate_tot_dev - rate_tot)/(rate_tot-1));




%%

% Package, save, and display elasticities
elasticities = [K_Y, labor_elas, savings_elas];

save(fullfile(save_dir, 'elasticities.mat'), 'K_Y', 'labor_elas', 'savings_elas')

displaynames = { 'Capital-to-output ratio', 'Labor elasticity', 'Savings elasticity' };
for i = 1:length(elasticities)
    fprintf('\t%-25s= % 7.4f\n', displaynames{i}, elasticities(i))
end
fprintf('\n')


end




%% Dynamic optimization solver and distribution generator applied to steady state
% (Generalized for use in savings elasticity calculation)

function [opt, dist, kpr_total, debt_total, elab_total, beq_total, wage, cap_share, debt_share, rate_cap, rate_gov, rate_tot] ...
          ...
            = solve_and_generate(...
                rho, beq, kpr, debt, ...
                T, Tr, Tss, bgrid, kgrid, ndem, nb, nk, nz, ...
                mpci, rpci, A, alp, cap_tax_share, d, q_tobin, Vbeq, proddist, ss_tax_cred, surv, ...
                tau_cap, tau_capgain, tau_ss, tr_z, v_ss_max, z, ...
                avg_deduc, clinton, coefs, limit, X, r_cbo, ben, ...
                beta, gamma, sigma, MU2, MU3, ...
                rate_deviation)

economy = 'steady';

wage = A*(1-alp)*(rho^alp);

cap_share  = (kpr-debt)/kpr;
debt_share = 1 - cap_share;

rate_cap = 1 + (A*alp*(rho^(alp-1)) - d)/q_tobin;
rate_gov = 1 + r_cbo;

% Scale rates by deviation if provided
if exist('rate_deviation', 'var')
	rate_cap = rate_cap * (1 + rate_deviation);
	rate_gov = rate_gov * (1 + rate_deviation);
end

rate_tot  = cap_share*rate_cap + debt_share*rate_gov;

% Initialize global aggregates
kpr_total  = 0;
beq_total  = 0;
elab_total = 0;

% Define birth year range
birthyears = 0;

% Initialize storage structures for optimal decision values and distributions
s_birthyear = struct('birthyear', num2cell(repmat(birthyears', [1,ndem])));
opt     = s_birthyear;
dist    = s_birthyear;
cohorts = s_birthyear;

for idem = 1:ndem
    
    % Extract demographic adjustments
    mu2_idem = MU2(idem,:);
    mu3_idem = MU3(idem,:);
    
    % Define distribution used for initialization
    dist_w0 = zeros(nk,nz,nb,1);
    dist_w0(1,:,1,1) = proddist;
    
    parfor i = 1:length(birthyears)
        
        % Get birth year
        birthyear = birthyears(i);
        
        % Generate optimal decision values, distributions, and cohort aggregates
        [labopt, dist_w, dist_r, N_w, N_r, Kalive, Kdead, ELab, Lab, Lfpr, Fedincome, Fedit, SSrev, Fcaptax, SSexp] ...
         ...
           = generate_distributions(...
               beta, gamma, sigma, T, Tr, T, birthyear, ...
               kgrid, z, tr_z, bgrid, nk, nz, nb, idem, ...
               mpci, rpci, cap_tax_share, ss_tax_cred, surv, tau_cap, tau_capgain, tau_ss*ones(1,T), v_ss_max*ones(1,T), ...
               beq*ones(1,T), wage*ones(1,T), cap_share*ones(1,T), debt_share*ones(1,T), rate_cap*ones(1,T), rate_gov*ones(1,T), rate_tot*ones(1,T), zeros(1,T), q_tobin, q_tobin, Vbeq, ...
               avg_deduc, clinton, coefs, limit, X, ben*ones(1,T), ...
               dist_w0, [], mu2_idem, mu3_idem);
        
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
        
        debt_total = debt;
        
end


end