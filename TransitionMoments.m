%%
% Distribution moments generator tool.
%
%%

%% SCENARIOS

% Params for K/Y=3 with d=0.056
params.beta =  0.98745000;
params.gamma = 0.75000000;
params.sigma = 1.24000000;
params.modelunit_dollar = 4.359874681178362e-05;
params.depreciation = 0.056;

% Solve for baseline steady state
scenario   = Scenario(struct('economy', 'closed', 'beta', params.beta, 'gamma', params.gamma, ...
                             'sigma', params.sigma, 'modelunit_dollar', params.modelunit_dollar, ...
                             'depreciation', params.depreciation, 'bequest_phi_1', 0, 'base_brackets', 'C'));

% Scenarios and directories
sc_steady  = scenario.currentPolicy.steady;
sc_open    = scenario.open();
sc_closed  = scenario.closed();
steady_dir = PathFinder.getWorkingDir(sc_steady);
open_dir   = PathFinder.getWorkingDir(sc_open);
closed_dir = PathFinder.getWorkingDir(sc_closed);

save_dir   = closed_dir;


%% PARAMETERS

% Define time constants
s = ParamGenerator.timing( sc_closed );
T_life  = s.T_life;    % Total life years
T_work  = s.T_work;    % Retirement age
T_model = s.T_model;   % Transition path model years

% Define grids
s = ParamGenerator.grids( sc_closed );
ndem = s.ndem;       % demographic types
ng   = s.ng;         % num groups
nz   = s.nz;         % num labor productivity shocks
zs   = s.zs;         % shocks grid (by demographic type and age)
nk   = s.nk;         % num asset points
nb   = s.nb;         % num avg. earnings points
% Useful later for a couple of functions
kv = s.kv;
karray = repmat(reshape(s.kv, [1,nk,1,1,1,1,1]),[nz,1,nb,T_life,ng,T_model,ndem]);


%% DISTRIBUTION AND POLICY FUNCTIONS

% Import households distribution
s    = load( fullfile(save_dir, 'distribution.mat' ) );
DIST = s.DIST;
dist_l(1:nz,1:nk,1:nb,1:T_work,1:ng,1:T_model,1:ndem) = DIST(1:nz,1:nk,1:nb,1:T_work,1:ng,1:T_model,1:ndem); % Working age population
dist_l(1:nz,1:nk,1:nb,T_work:T_life,1:ng,1:T_model,1:ndem) = 0; % Retired population

% Import market variables
s     = load( fullfile(save_dir, 'market.mat' ) );
wages = s.wages;

% Import policy functions
f = @(X) repmat(reshape(X, [nz,nk,nb,T_life,1,T_model,ndem]), [1,1,1,1,ng,1,1]);
s      = load( fullfile(save_dir, 'all_decisions.mat' ) );
% labinc = f(s.LAB) .* repmat(reshape(zs, [nz,1,1,T_life,1,1,ndem]),[1,nk,nb,1,ng,T_model,1]) * wages;
k      = f(s.K);
cons   = f(s.CON);
pit    = f(s.PIT);


%% PIT DISTRIBUTION BY ASSETS

pit_a        = zeros(nk,T_model);     % PIT at ik asset holdings level (in dollars)
pit_a_groups = zeros( 3,T_model);
pit_a_perc   = zeros(nk,T_model);     % percentage PIT at ik asset holdings level (% of total PIT)
pit_a_pcap   = zeros(nk,T_model);     % PIT per capita at ik asset holdings level (in dollars)

for t = 1:T_model
    
   pit_t     = pit (:,:,:,:,:,t,:);
   dist_t    = DIST(:,:,:,:,:,t,:);
   pitdist_t = pit_t .* dist_t;
   
   for ik = 1:nk
       
      pop_ik = dist_t(:,ik,:,:,:,:,:);
      pitdist_t_ik = pitdist_t(:,ik,:,:,:,:,:);
      pit_a(ik,t) = sum(pitdist_t_ik(:));
      pit_a_perc(ik,t) = sum(pitdist_t_ik(:))/sum(pitdist_t(:));
      pit_a_pcap(ik,t) = (sum(pitdist_t_ik(:))/sum(pop_ik(:)))/params.modelunit_dollar;
      
   end
   
   pit_a_graph(1,t) = sum(pit_a(1:3,t));
   pit_a_graph(2,t) = sum(pit_a(4:8,t));
   pit_a_graph(3,t) = sum(pit_a(9:12,t));
    
end

figure
plot(1:T_model,pit_a_graph(1,:),1:T_model,pit_a_graph(2,:),1:T_model,pit_a_graph(3,:),'LineWidth',2)
title('Total PIT by asset holdings group','FontSize',16)
xlabel('T model','FontSize',13)
ylabel('model units','FontSize',13)
legend({'poor', 'middle class', 'up middle class'},'Location','northwest','FontSize',13)

figure
plot(1:T_model,pit_a(1,:),1:T_model,pit_a(2,:),1:T_model,pit_a(3,:),...
     1:T_model,pit_a(4,:),1:T_model,pit_a(5,:),1:T_model,pit_a(6,:),...
     1:T_model,pit_a(7,:),1:T_model,pit_a(8,:),1:T_model,pit_a(9,:),...
     'LineWidth',2)
title('Total PIT by asset holdings','FontSize',16)
xlabel('T model','FontSize',13)
ylabel('2016 dollars','FontSize',13)
legend({sprintf('%0.2f', kv(1)/params.modelunit_dollar),sprintf('%0.2f', kv(2)/params.modelunit_dollar), ...
        sprintf('%0.2f', kv(3)/params.modelunit_dollar),sprintf('%0.2f', kv(4)/params.modelunit_dollar), ...
        sprintf('%0.2f', kv(5)/params.modelunit_dollar),sprintf('%0.2f', kv(6)/params.modelunit_dollar), ...
        sprintf('%0.2f', kv(7)/params.modelunit_dollar),sprintf('%0.2f', kv(8)/params.modelunit_dollar) ...
        sprintf('%0.2f', kv(9)/params.modelunit_dollar)},'FontSize',13)

figure
plot(1:T_model,pit_a_pcap(1,:),1:T_model,pit_a_pcap(2,:),1:T_model,pit_a_pcap(3,:),...
     1:T_model,pit_a_pcap(4,:),1:T_model,pit_a_pcap(5,:),1:T_model,pit_a_pcap(6,:),...
     1:T_model,pit_a_pcap(7,:),1:T_model,pit_a_pcap(8,:),1:T_model,pit_a_pcap(9,:),...
     'LineWidth',2)
title('PIT per capita by asset holdings','FontSize',16)
xlabel('T model','FontSize',13)
ylabel('2016 dollars','FontSize',13)
legend({sprintf('%0.2f', kv(1)/params.modelunit_dollar),sprintf('%0.2f', kv(2)/params.modelunit_dollar), ...
        sprintf('%0.2f', kv(3)/params.modelunit_dollar),sprintf('%0.2f', kv(4)/params.modelunit_dollar), ...
        sprintf('%0.2f', kv(5)/params.modelunit_dollar),sprintf('%0.2f', kv(6)/params.modelunit_dollar), ...
        sprintf('%0.2f', kv(7)/params.modelunit_dollar),sprintf('%0.2f', kv(8)/params.modelunit_dollar) ...
        sprintf('%0.2f', kv(9)/params.modelunit_dollar)},'FontSize',13)


%% MODEL ASSET DISTRIBUTIONS BY GENERATION
%{
kdist_age = zeros(T_life,T_model);

for t = 1:T_model
    
   assets = karray(:,:,:,:,:,t,:);
   dist   = DIST  (:,:,:,:,:,t,:);
   kdist  = dist .* assets;
   for age = 1:T_life
       pop_age_temp   = dist(:,:,:,age,:,:,:);
       kdist_age_temp = kdist(:,:,:,age,:,:,:);
       kdist_age(age,t) = (sum(kdist_age_temp(:))/sum(pop_age_temp(:)))/modelunit_dollar;
   end
    
end

kdist_generation = zeros(T_life/10,T_model);

kdist_generation = sum(kdist_age( 1:10,:));
kdist_generation = sum(kdist_age(11:20,:));
kdist_generation = sum(kdist_age(21:30,:));
kdist_generation = sum(kdist_age(31:40,:));
kdist_generation = sum(kdist_age(41:50,:));
kdist_generation = sum(kdist_age(51:60,:));
kdist_generation = sum(kdist_age(61:70,:));
kdist_generation = sum(kdist_age(71:80,:));
%}

