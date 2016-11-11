%%
% Steady state solver.
% 
% Arguments:
% 
%   deep_params
%       beta, gamma, and sigma preference parameters collected in a 1 x 3 vector [beta, gamma, sigma].
% 
%   cleanup (optional | false by default)
%       Set to true to clean up data files and directories used for intermediary results after solver execution.
% 
%   this_uniquetag (optional | '' by default)
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


function [save_dir, elasticities] = solve_ss(deep_params, cleanup, this_uniquetag)


%% Initialization

% Extract deep parameters or set defaults if none provided
if ~exist('deep_params', 'var') || isempty(deep_params)
    deep_params = [1.005, 0.390, 06.00];
end
beta  = deep_params(1);
gamma = deep_params(2);
sigma = deep_params(3);

% Turn off file cleanup by default
if ~exist('cleanup', 'var') || isempty(cleanup)
    cleanup = false;
end

% Set solution uniqueness tag to empty by default
if ~exist('this_uniquetag', 'var') || isempty(this_uniquetag)
    this_uniquetag = '';
end



% Display problem definition
fprintf('\nSolving steady state:  beta = %0.3f  gamma = %0.3f  sigma = %05.2f\n', beta, gamma, sigma)

% Identify working directories
param_dir = dirFinder.param;
save_dir  = dirFinder.ss(deep_params);

% Append uniqueness tag to name of save directory
save_dir = [save_dir, this_uniquetag];

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


% Load CBO parameters
s = load(fullfile(param_dir, 'param_cbo.mat'));

r_cbo = mean(s.r_cbo);


% Load social security parameters
s = load(fullfile(param_dir, 'param_socsec.mat'));

Tr       = s.NRA(1);
ben      = s.ss_benefit(:,1);
tau_ss   = s.ss_tax(1);
v_ss_max = s.taxmax(1);


% Load income tax parameters, using base plan
s = load(fullfile(param_dir, 'param_inctax_base.mat'));

avg_deduc = deduc_scale * s.avg_deduc;
clinton   = s.clinton;
coefs     = s.coefs;
limit     = s.limit;
X         = s.X;


% Load business tax parameters, using base plan
s = load(fullfile(param_dir, 'param_bustax_base.mat'));

cap_tax_share = s.cap_tax_share;
exp_share     = s.exp_share;
tau_cap       = s.tau_cap;
tau_capgain   = s.tau_capgain;

q_tobin = 1 - tau_cap * exp_share;


% Load initial values for steady state solution
s = load(fullfile(param_dir, 'ss_solution0.mat'));

rho    = s.rho;
beq    = s.beq;
kpr    = s.kpr;
DEBTss = s.DEBTss;

D_Y = 0.74;
Y = DEBTss / D_Y;


% Clear parameter loading structure
clear('s')




%% Steady state calculation

% Define tolerance for rho convergence and initialize collection of rho error terms
rho_tol = 1e-3;
rho_eps = [];

% Initialize iteration count and set maximum number of iterations
rho_iter    =  0;
rho_itermax = 25;

while true
    
    % Increment iteration count
    rho_iter = rho_iter + 1;
    fprintf('\trho_iter = %2d  ...  ', rho_iter)
    
    % Solve dynamic optimization and generate distributions
    [opt, dist, kpr_total, DD, elab_total, beq_total, wage, cap_share, debt_share, rate_cap, rate_gov, rate_tot] ...
     ...
       = solve_and_generate(...
           rho, kpr, DEBTss, beq, ...
           T, Tr, Tss, bgrid, kgrid, ndem, nb, nk, nz, ...
           mpci, rpci, A, alp, cap_tax_share, d, q_tobin, Vbeq, proddist, ss_tax_cred, surv, ...
           tau_cap, tau_capgain, tau_ss, tr_z, v_ss_max, z, ...
           avg_deduc, clinton, coefs, limit, X, r_cbo, ben, ...
           beta, gamma, sigma, MU2, MU3); %#ok<ASGLU>
    
    % Check for convergence
    rhopr = (max(kpr_total-DD, 0)/q_tobin)/elab_total;
    
    rho_eps = [rho_eps, abs(rho - rhopr)]; %#ok<AGROW>
    fprintf('rho_eps = %7.4f\n', rho_eps(end))
    if (rho_eps(end) < rho_tol), break, end
    
    % Check for maximum iterations
    if (rho_iter == rho_itermax)
        warning('Maximum iterations reached.')
        break
    end
    
    % Update variables
    stepfactor = 0.5;
    rho = (1-stepfactor) * rho + stepfactor * rhopr;
    
    Y      = A*(max((kpr_total-DD)/q_tobin, 0)^alp)*(elab_total^(1-alp));
    DEBTss = D_Y * Y;
    
    beq = beq_total;
    kpr = kpr_total;
    
end
fprintf('\n')

% Save log of rho iterations
fid = fopen(fullfile(save_dir, 'iterations.txt'), 'wt');
fprintf(fid, '%d rho iteration(s) of a maximum %d\n\nError term(s):', rho_iter, rho_itermax);
for i = 1:rho_iter
    fprintf(fid, '\n  % 2d  --  %7.4f', i, rho_eps(i));
end
fclose(fid);


% Save optimal decision values and distributions
save(fullfile(save_dir, 'opt.mat'),  'opt' )
save(fullfile(save_dir, 'dist.mat'), 'dist')

% Save solution
save(fullfile(save_dir, 'solution.mat'), ...
     'rho', 'kpr', 'Y', 'DEBTss', 'beq', 'wage', ...
     'cap_share', 'debt_share', 'rate_cap', 'rate_gov', 'rate_tot')


% Calculate capital to output ratio
K_Y = (kpr-DD)/Y;




%% Labor elasticity calculation

working_mass = 0;
frisch_total = 0;

for idem = 1:ndem
    
    labopt = opt(1,idem).labopt;
    dist_w = dist(1,idem).dist_w;
    
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
       rho, kpr, DEBTss, beq, ...
       T, Tr, Tss, bgrid, kgrid, ndem, nb, nk, nz, ...
       mpci, rpci, A, alp, cap_tax_share, d, q_tobin, Vbeq, proddist, ss_tax_cred, surv, ...
       tau_cap, tau_capgain, tau_ss, tr_z, v_ss_max, z, ...
       avg_deduc, clinton, coefs, limit, X, r_cbo, ben, ...
       beta, gamma, sigma, MU2, MU3, ...
       rate_deviation);

% Calculate savings elasticity
savings_elas = ((kpr_dev - kpr)/kpr) / ((rate_tot_dev - rate_tot)/(rate_tot-1));




%%

% Clean up working directory
if cleanup
    rmdir(fullfile(save_dir, '..'), 's')
end

% Package elasticities and display
elasticities = [K_Y, labor_elas, savings_elas];

displaynames = { 'Capital-to-output ratio', 'Labor elasticity', 'Savings elasticity' };
for i = 1:length(elasticities)
    fprintf('\t%-25s= % 7.4f\n', displaynames{i}, elasticities(i))
end
fprintf('\n')


end




%% Dynamic optimization solver and distribution generator applied to steady state
% (Generalized for use in savings elasticity calculation)

function [opt, dist, kpr_total, DD, elab_total, beq_total, wage, cap_share, debt_share, rate_cap, rate_gov, rate_tot] ...
          ...
            = solve_and_generate(...
                rho, kpr, DD, beq, ...
                T, Tr, Tss, bgrid, kgrid, ndem, nb, nk, nz, ...
                mpci, rpci, A, alp, cap_tax_share, d, q_tobin, Vbeq, proddist, ss_tax_cred, surv, ...
                tau_cap, tau_capgain, tau_ss, tr_z, v_ss_max, z, ...
                avg_deduc, clinton, coefs, limit, X, r_cbo, ben, ...
                beta, gamma, sigma, MU2, MU3, ...
                rate_deviation)


wage = A*(1-alp)*(rho^alp);

cap_share  = (kpr-DD)/kpr;
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

for idem = 1:ndem
    
    % Solve dynamic optimization
    [V, Vss, kopt, koptss, labopt, bopt, ...
     fedincome, fedincomess, ...
     fitax,     fitaxss,     ...
     fsstax,    fsstaxss,    ...
     fcaptax,   fcaptaxss]   ...
     ...
       = solve_dynamic_optimization_mex(...
           beq*ones(1,Tss), T, Tr, Tss, 0, ...
           bgrid, kgrid, nb, nk, nz, idem, ...
           mpci, rpci, cap_tax_share, ss_tax_cred, ...
           surv, tau_cap, tau_capgain, tau_ss*ones(1,Tss), tr_z, v_ss_max*ones(1,Tss), z, ...
           wage*ones(1,Tss), cap_share*ones(1,Tss), debt_share*ones(1,Tss), rate_cap*ones(1,Tss), rate_gov*ones(1,Tss), rate_tot*ones(1,Tss), zeros(1,Tss), q_tobin, q_tobin, Vbeq, ...
           avg_deduc, clinton, coefs, limit, X, ben*ones(1,Tss), ...
           beta, gamma, sigma);
    
    % Check for unsolved dynamic optimization subproblems
    if any(isinf([V(:); Vss(:)]))
        warning('Infinite utility values found.  Some dynamic optimization subproblems unsolved.  Check that initial conditions satisfy constraints.')
    end
    
    % Package optimal decision values
    opt(1,idem).V           = V;            %#ok<AGROW>
    opt(1,idem).Vss         = Vss;          %#ok<AGROW>
    opt(1,idem).kopt        = kopt;         %#ok<AGROW>
    opt(1,idem).koptss      = koptss;       %#ok<AGROW>
    opt(1,idem).labopt      = labopt;       %#ok<AGROW>
    opt(1,idem).bopt        = bopt;         %#ok<AGROW>
    opt(1,idem).fedincome   = fedincome;    %#ok<AGROW>
    opt(1,idem).fedincomess = fedincomess;  %#ok<AGROW>
    opt(1,idem).fcaptax     = fcaptax;      %#ok<AGROW>
    opt(1,idem).fcaptaxss   = fcaptaxss;    %#ok<AGROW>
    opt(1,idem).fitax       = fitax;        %#ok<AGROW>
    opt(1,idem).fitaxss     = fitaxss;      %#ok<AGROW>
    opt(1,idem).fsstax      = fsstax;       %#ok<AGROW>
    opt(1,idem).fsstaxss    = fsstaxss;     %#ok<AGROW>
    
    
    % Extract demographic adjustments
    mu2_idem = MU2(idem,:);
    mu3_idem = MU3(idem,:);
    
    % Define distribution used for initialization
    dist_w0 = zeros(nk,nz,nb,1);
    dist_w0(1,:,1,1) = proddist;
    
    % Generate distributions
    [dist_w, dist_r, ~, ~, Kalive, Kdead, ELab] ...
     ...
       = generate_distributions(...
           kopt, koptss, labopt, bopt, ...
           fedincome, fedincomess, fitax, fitaxss, fsstax, fsstaxss, fcaptax, fcaptaxss, ...
           T, Tr, T, 0, ...
           kgrid, bgrid, nk, nz, nb, idem, ...
           z, tr_z, ben, ...
           dist_w0, [], ...
           mu2_idem, mu3_idem);
    
    % Package distributions
    dist(1,idem).dist_w = dist_w;  %#ok<AGROW>
    dist(1,idem).dist_r = dist_r;  %#ok<AGROW>
    dist(1,idem).Kalive = Kalive;  %#ok<AGROW>
    dist(1,idem).Kdead  = Kdead;   %#ok<AGROW>
    dist(1,idem).ELab   = ELab;    %#ok<AGROW>
    
    % Add aggregate values to global aggregates
    kpr_total  = kpr_total  + sum(Kalive + Kdead);
    beq_total  = beq_total  + sum(Kdead);
    elab_total = elab_total + sum(ELab);
    
end

end