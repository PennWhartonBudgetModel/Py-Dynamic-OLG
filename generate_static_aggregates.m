%%
% Calculate static aggregates for a transition path.  Supports both open and closed economy transition paths.
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
%       Leave unspecified or set to [] for an open economy transition path.
% 
%   this_uniquetag (optional | '' by default)
%       String used to save solution into a unique directory, used to avoid conflicts between parallel runs.
% 
%%


function [] = generate_static_aggregates(deep_params, plan, gcut, this_uniquetag)

% Extract deep parameters or set defaults if none provided
if ~exist('deep_params', 'var') || isempty(deep_params)
    deep_params = [1.005, 0.390, 06.00];
end

% Set base plan as default
if ~exist('plan', 'var') || isempty(plan)
    plan = 'base';
end

% Identify open economy transition path by absence of gcut
isopen = ~exist('gcut', 'var') || isempty(gcut);

% Set solution uniqueness tag to empty by default
if ~exist('this_uniquetag', 'var') || isempty(this_uniquetag)
    this_uniquetag = '';
end



% Identify working directories
param_dir = dirFinder.param;

if isopen
    save_dir = dirFinder.open  (deep_params, plan);
else
    save_dir = dirFinder.closed(deep_params, plan, gcut);
end
save_dir = [save_dir, this_uniquetag];


% Identify corresponding reference directories
ss_dir = dirFinder.ss(deep_params);

if strcmp(plan, 'base')
    base_dir = save_dir;    % Includes this_uniquetag
else
    if isopen
        base_dir = dirFinder.open  (deep_params, 'base');
    else
        base_dir = dirFinder.closed(deep_params, 'base', +0.00);
    end
end



% Load global parameters
s = load(fullfile(param_dir, 'param_global.mat'));

T   = s.T;
Tss = s.Tss;

kgrid = s.kgrid;

ndem = s.ndem;
nb   = s.nb;
nk   = s.nk;
nz   = s.nz;

mpci = s.mpci;
rpci = s.rpci;

ss_tax_cred = s.ss_tax_cred;

deduc_scale = s.deduc_scale;
z           = s.z;


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
tau_capgain       = s.tau_capgain;

q_tobin = 1 - tau_cap * exp_share;

s_base = load(fullfile(param_dir, 'param_bustax_base.mat'));
q_tobin0 = 1 - s_base.tau_cap * s_base.exp_share;
clear('s_base')


% Load baseline solution
% (Baseline solution is equivalent to steady state solution for open economy)
if isopen
    
    s = hardyload(fullfile(ss_dir, 'solution.mat'));

    wages        = s.wage      *ones(1,Tss);
    cap_shares   = s.cap_share *ones(1,Tss);
    debt_shares  = s.debt_share*ones(1,Tss);
    rate_caps    = s.rate_cap  *ones(1,Tss);
    rate_govs    = s.rate_gov  *ones(1,Tss);
    exp_subsidys = zeros(1,Tss);
    
else

    s = hardyload(fullfile(base_dir, 'solution.mat'));

    wages        = s.wages;
    cap_shares   = s.cap_shares;
    debt_shares  = s.debt_shares;
    rate_caps    = s.rate_caps;
    rate_govs    = s.rate_govs;
    
    s_agg = hardyload(fullfile(base_dir,'aggregates.mat'));
    exp_subsidys = s_agg.exp_subsidys;
    clear('s_agg');
    
end

% Load optimal decision values and distributions from baseline
s = hardyload(fullfile(base_dir, 'opt.mat'));
base_opt = s.opt;

s = hardyload(fullfile(base_dir, 'dist.mat'));
base_dist = s.dist;



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

        
        %% (Adapted from dynamic optimization solver)
        
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
        
        
        %% (Adapted from distribution generator)
        
        % Extract distributions
        dist_w = base_dist(birthyear+T,idem).dist_w;
        dist_r = base_dist(birthyear+T,idem).dist_r;

        % Redefine effective numbers of working and retirement years
        % (Redefinition needed because distribution generator makes use of head truncation whilst dynamic optimization solver does not yet)
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
        
        
        %% (Adapted from open economy solver)
        
        % Add aggregates to total static aggregates, aligning to projection years
        T_shift = max(0, birthyear);
        fedincome_static( T_shift + (1:N_w_dist+N_r_dist) ) = fedincome_static( T_shift + (1:N_w_dist+N_r_dist) ) + Fedincome;
        fedit_static    ( T_shift + (1:N_w_dist+N_r_dist) ) = fedit_static    ( T_shift + (1:N_w_dist+N_r_dist) ) + Fedit;
        ssrev_static    ( T_shift + (1:N_w_dist+N_r_dist) ) = ssrev_static    ( T_shift + (1:N_w_dist+N_r_dist) ) + SSrev;
        fcaptax_static  ( T_shift + (1:N_w_dist+N_r_dist) ) = fcaptax_static  ( T_shift + (1:N_w_dist+N_r_dist) ) + Fcaptax;
        ssexp_static    ( T_shift + N_w_dist+(1:N_r_dist) ) = ssexp_static    ( T_shift + N_w_dist+(1:N_r_dist) ) + SSexp;
        
        
    end

end

% Copy extra static aggregates from corresponding base plan aggregates
% (Duplicates data but enhances modularity of results)
s_base = hardyload(fullfile(base_dir, 'aggregates.mat'));

elab_static          = s_base.elab_total;           %#ok<NASGU>
cap_static           = s_base.cap_total;            %#ok<NASGU>
domestic_cap_static  = s_base.domestic_cap_total;   %#ok<NASGU>
foreign_cap_static   = s_base.foreign_cap_total;    %#ok<NASGU>
domestic_debt_static = s_base.domestic_debt_total;  %#ok<NASGU>

% Loaded in solve_open for open; loaded from results if closed.
if ~isopen
    foreign_debt_static = s_base.foreign_debt_total; %#ok<NASGU>
end

Y_static      = s_base.Y_total;         %#ok<NASGU>
lfpr_static   = s_base.lfpr_total;      %#ok<NASGU>

labinc_static = s_base.labinc_total;
kinc_static   = s_base.kinc_total;      %#ok<NASGU>

feditlab_static          = fedit_static .* labinc_static ./ fedincome_static; %#ok<NASGU>
domestic_fcaptax_static  = fcaptax_static + fedit_static .* (1 - labinc_static ./ fedincome_static);
foreign_fcaptax_static   = zeros(1,Tss);
fcaprev_static           = domestic_fcaptax_static + foreign_fcaptax_static; %#ok<NASGU>

clear('s_base')


% Save static aggregates
save(fullfile(save_dir, 'aggregates_static.mat'), ...
     'fedincome_static', 'fedit_static', 'ssrev_static', 'fcaptax_static', 'ssexp_static', ...
     'elab_static', 'cap_static', 'domestic_cap_static', 'foreign_cap_static', ...
     'Y_static', 'lfpr_static', 'labinc_static', 'kinc_static', ...
     'feditlab_static', 'fcaprev_static', 'domestic_fcaptax_static', 'foreign_fcaptax_static', 'domestic_debt_static')

if ~isopen
    save(fullfile(save_dir, 'aggregates_static.mat'), 'foreign_debt_static', '-append')
end

end