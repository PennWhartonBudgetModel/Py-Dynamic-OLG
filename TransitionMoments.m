%%
% Distribution moments generator tool.
%
%%

%% SCENARIOS

% Params for K/Y=3 with d=0.056
% params.beta =  0.98745000;
% params.gamma = 0.75000000;
% params.sigma = 1.24000000;
% params.modelunit_dollar = 4.359874681178362e-05;
% params.depreciation = 0.056;

% Params for K/Y=3 with d=0.08
params.beta =  1.003510000000000;
params.gamma = 0.680000000000000;
params.sigma = 1.500000000000000;
params.modelunit_dollar = 4.135682750000000e-05;
params.depreciation = 0.08;

% Solve for baseline steady state
scenario   = Scenario(struct('economy', 'closed', 'beta', params.beta, 'gamma', params.gamma, ...
                             'sigma', params.sigma, 'modelunit_dollar', params.modelunit_dollar, ...
                             'depreciation', params.depreciation, 'bequest_phi_1', 0, 'base_brackets', 'C'));

% Scenarios and directories
sc_steady  = scenario.currentPolicy.steady;
sc_open    = scenario.open();
sc_closed  = scenario.closed();
sc_obase   = scenario.currentPolicy.open;
sc_cbase   = scenario.currentPolicy.closed;
steady_dir = PathFinder.getWorkingDir(sc_steady);

sc_use  = sc_open;
use_dir = PathFinder.getWorkingDir(sc_use);
sc_base = sc_obase;
base_dir = PathFinder.getWorkingDir(sc_base);


%% PARAMETERS

% Define time constants
s = ParamGenerator.timing( sc_use );
T_life  = s.T_life;    % Total life years
T_work  = s.T_work;    % Retirement age
T_model = s.T_model;   % Transition path model years

% Define grids
s = ParamGenerator.grids( sc_use );
ndem = s.ndem;       % demographic types
ng   = s.ng;         % num groups
nz   = s.nz;         % num labor productivity shocks
zs   = s.zs;         % shocks grid (by demographic type and age)
nk   = s.nk;         % num asset points
nb   = s.nb;         % num avg. earnings points
% Useful later for a couple of functions
kv = s.kv;
karray = repmat(reshape(kv, [1,nk,1,1,1,1,1]),[nz,1,nb,T_life,ng,T_model+1,ndem]);
yearsv = [2017:2017+T_model];


%% DISTRIBUTION AND POLICY FUNCTIONS

% Import households distribution
DIST   = zeros(nz,nk,nb,T_life,ng,T_model+1,ndem);
s      = load( fullfile(steady_dir, 'distribution.mat' ) );
DISTss = s.DIST;
% s      = load( fullfile(use_dir, 'distribution.mat' ) );
s      = load( fullfile(use_dir, 'Static_distribution.mat' ) );
DIST(:,:,:,:,:,1    ,:) = DISTss;
% DIST(:,:,:,:,:,2:end,:) = s.DIST;
DIST(:,:,:,:,:,2:end,:) = s.Static_DIST;

% Import market variables
s       = load( fullfile(steady_dir, 'market.mat' ) );
wagesss = s.wages;
s       = load( fullfile(base_dir, 'market.mat' ) );
wages   = [wagesss s.wages];

% Import policy functions
labinc = zeros(nz,nk,nb,T_life,ng,T_model+1,ndem);
lab    = zeros(nz,nk,nb,T_life,ng,T_model+1,ndem);
pit    = zeros(nz,nk,nb,T_life,ng,T_model+1,ndem);
cit    = zeros(nz,nk,nb,T_life,ng,T_model+1,ndem);
inc    = zeros(nz,nk,nb,T_life,ng,T_model+1,ndem);
sst    = zeros(nz,nk,nb,T_life,ng,T_model+1,ndem);

f = @(X) repmat(reshape(X, [nz,nk,nb,T_life,1,1,ndem]), [1,1,1,1,ng,1,1]);
s = load( fullfile(steady_dir, 'all_decisions.mat' ) );
lab(:,:,:,:,:,1,:) = f(s.LAB);
pit(:,:,:,:,:,1,:) = f(s.PIT);
cit(:,:,:,:,:,1,:) = f(s.CIT);
inc(:,:,:,:,:,1,:) = f(s.INC);
sst(:,:,:,:,:,1,:) = f(s.SST);

f = @(X) repmat(reshape(X, [nz,nk,nb,T_life,1,T_model,ndem]), [1,1,1,1,ng,1,1]);
% s = load( fullfile(use_dir, 'all_decisions.mat' ) );
s = load( fullfile(use_dir, 'Static_all_decisions.mat' ) );
lab(:,:,:,:,:,2:end,:) = f(s.LAB);
pit(:,:,:,:,:,2:end,:) = f(s.PIT);
cit(:,:,:,:,:,2:end,:) = f(s.CIT);
inc(:,:,:,:,:,2:end,:) = f(s.INC);
sst(:,:,:,:,:,2:end,:) = f(s.SST);

labinc = lab .* repmat(reshape(zs, [nz,1,1,T_life,1,1,ndem]),[1,nk,nb,1,ng,T_model+1,1]) .* repmat(reshape(wages, [1,1,1,1,1,T_model+1,1]),[nz,nk,nb,T_life,ng,1,ndem]);
taxes1 = pit + cit;
taxes2 = pit + cit + sst;


%% CIT TAXES DISTRIBUTION BY TAXABLE INCOME

inc_dist   = zeros(5,1);
inc_cum    = zeros(5,1);
inc_share  = zeros(5,1);
inc_groups = zeros(5,1);
asset_cum  = zeros(5,1);
asset_groups = zeros(5,1);
cit_cum    = zeros(5,1);
cit_groups = zeros(5,1);
pit_cum    = zeros(5,1);
pit_groups = zeros(5,1);
linc_cum    = zeros(5,1);
linc_groups = zeros(5,1);

% Period variables (steady state = 1)
t = 2;
inc_1   = inc(:,:,:,:,:,t,:);
dist_1  = DIST(:,:,:,:,:,t,:);
asset_1 = karray(:,:,:,:,:,t,:);
cit_1   = cit(:,:,:,:,:,t,:);
pit_1   = pit(:,:,:,:,:,t,:);
linc_1  = labinc(:,:,:,:,:,t,:);

% Vectorize variables
inc_1   = inc_1(:);
dist_1  = dist_1(:);
asset_1 = asset_1(:);
cit_1   = cit_1(:);
pit_1   = pit_1(:);
linc_1  = linc_1(:);

% Sort variables
[inc_1, sortv] = sort(inc_1);
dist_1  = dist_1(sortv);
asset_1 = asset_1(sortv);
cit_1   = cit_1(sortv);
pit_1   = pit_1(sortv);
linc_1  = linc_1(sortv);

% Find quintiles
q = 1;
total_inc = sum(inc_1.*dist_1);
for quintile = 0.2:0.2:0.8
    
    % Taxable income distribution
    i = find(cumsum(dist_1) >= quintile,1);
    inc_dist(q,1)  = sum(dist_1(1:i));
    inc_cum(q,1)   = sum(inc_1(1:i).*dist_1(1:i));
    inc_share(q,1) = inc_cum(q,1)/total_inc;
    
    % Asset distributed according to taxable income
    asset_cum(q,1) = sum(asset_1(1:i).*dist_1(1:i));
    % Labor income distributed according to taxable income
    linc_cum(q,1) = sum(linc_1(1:i).*dist_1(1:i));
    
    % Taxes distributed according to taxable income
    cit_cum(q,1) = sum(cit_1(1:i).*dist_1(1:i));
    pit_cum(q,1) = sum(pit_1(1:i).*dist_1(1:i));
    
    % Counter
    q = q + 1;
    
end

% Top quintile
% Taxable income
inc_dist(q,1)  = sum(dist_1);
inc_cum(q,1)   = total_inc;
inc_share(q,1) = 1;
inc_groups = params.modelunit_dollar * diff([0 inc_cum'])';
% Assets
asset_cum(q,1) = sum(asset_1(1:end).*dist_1(1:end));
asset_groups = diff([0 asset_cum'])';
% Labor income
linc_cum(q,1) = sum(linc_1(1:end).*dist_1(1:end));
linc_groups = diff([0 linc_cum'])';
% Taxes
cit_cum(q,1) = sum(cit_1(1:end).*dist_1(1:end));
cit_groups = diff([0 cit_cum'])';
pit_cum(q,1) = sum(pit_1(1:end).*dist_1(1:end));
pit_groups = diff([0 pit_cum'])';


