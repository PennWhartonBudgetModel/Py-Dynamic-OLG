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
steady_dir = PathFinder.getWorkingDir(sc_steady);
open_dir   = PathFinder.getWorkingDir(sc_open);
closed_dir = PathFinder.getWorkingDir(sc_closed);


%% PARAMETERS

% Define time constants
s = ParamGenerator.timing( sc_open );
T_life  = s.T_life;    % Total life years
T_work  = s.T_work;    % Retirement age
T_model = s.T_model;   % Transition path model years

% Define grids
s = ParamGenerator.grids( sc_open );
ndem = s.ndem;       % demographic types
ng   = s.ng;         % num groups
nz   = s.nz;         % num labor productivity shocks
zs   = s.zs;         % shocks grid (by demographic type and age)
nk   = s.nk;         % num asset points
nb   = s.nb;         % num avg. earnings points
% Useful later for a couple of functions
kv = s.kv;
karray = repmat(reshape(s.kv, [1,nk,1,1,1,1,1]),[nz,1,nb,T_life,ng,T_model,ndem]);
yearsv = [2017:2017+T_model];


%% DISTRIBUTION AND POLICY FUNCTIONS

% Import households distribution
DIST   = zeros(nz,nk,nb,T_life,ng,T_model+1,ndem);
s      = load( fullfile(steady_dir, 'distribution.mat' ) );
DISTss = s.DIST;
s      = load( fullfile(open_dir, 'distribution.mat' ) );
DIST(:,:,:,:,:,1    ,:) = DISTss;
DIST(:,:,:,:,:,2:end,:) = s.DIST;

% Import market variables
s       = load( fullfile(steady_dir, 'market.mat' ) );
wagesss = s.wages;
s       = load( fullfile(open_dir, 'market.mat' ) );
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
s = load( fullfile(open_dir, 'all_decisions.mat' ) );
lab(:,:,:,:,:,2:end,:) = f(s.LAB);
pit(:,:,:,:,:,2:end,:) = f(s.PIT);
cit(:,:,:,:,:,2:end,:) = f(s.CIT);
inc(:,:,:,:,:,2:end,:) = f(s.INC);
sst(:,:,:,:,:,2:end,:) = f(s.SST);

labinc = lab .* repmat(reshape(zs, [nz,1,1,T_life,1,1,ndem]),[1,nk,nb,1,ng,T_model+1,1]) .* repmat(reshape(wages, [1,1,1,1,1,T_model+1,1]),[nz,nk,nb,T_life,ng,1,ndem]);
taxes = pit + cit + sst;


%% TAXABLE INCOME PERCENTILES AT THE INITIAL PERIOD

inc_ss = inc(:,:,:,:,:,1,:);
dist_ss = DIST(:,:,:,:,:,1,:);
taxable_inc_distribution = get_moments(dist_ss(:), inc_ss(:));


%% TOTAL TAXES DISTRIBUTION BY TAXABLE INCOME

total_tax  = zeros(T_model+1,1);
tax_groups = zeros(T_model+1,5);
tax_groups_perc = zeros(T_model+1,5);

for t = 1:T_model+1
 for iz = 1:nz
  for ik = 1:nk
   for ib = 1:nb
    for age = 1:T_life
     for ig = 1:ng
      for idem = 1:ndem
                            
       if (inc(iz,ik,ib,age,ig,t,idem) <= taxable_inc_distribution.threshold(1))
         tax_groups(t,1) = tax_groups(t,1) + taxes(iz,ik,ib,age,ig,t,idem);
       elseif (taxable_inc_distribution.threshold(1) < inc(iz,ik,ib,age,ig,t,idem)) && (inc(iz,ik,ib,age,ig,t,idem) <= taxable_inc_distribution.threshold(2))
        tax_groups(t,2) = tax_groups(t,2) + taxes(iz,ik,ib,age,ig,t,idem);
       elseif (taxable_inc_distribution.threshold(2) < inc(iz,ik,ib,age,ig,t,idem)) && (inc(iz,ik,ib,age,ig,t,idem) <= taxable_inc_distribution.threshold(3))
        tax_groups(t,3) = tax_groups(t,3) + taxes(iz,ik,ib,age,ig,t,idem);
       elseif (taxable_inc_distribution.threshold(3) < inc(iz,ik,ib,age,ig,t,idem)) && (inc(iz,ik,ib,age,ig,t,idem) <= taxable_inc_distribution.threshold(4))
        tax_groups(t,4) = tax_groups(t,4) + taxes(iz,ik,ib,age,ig,t,idem);
       elseif (inc(iz,ik,ib,age,ig,t,idem) > taxable_inc_distribution.threshold(4))
        tax_groups(t,5) = tax_groups(t,5) + taxes(iz,ik,ib,age,ig,t,idem);
       end
                            
       total_tax(t,1) = total_tax(t,1) + taxes(iz,ik,ib,age,ig,t,idem);
                            
      end
     end
    end
   end
  end
 end
    
 tax_groups_perc(t,:) = tax_groups(t,:)/total_tax(t,1);
 tax_groups_perc(t,:) = tax_groups_perc(t,:)/sum(tax_groups_perc(t,:));
    
end

figure
hold on
yyaxis left
plot(yearsv,tax_groups_perc(:,1),yearsv,tax_groups_perc(:,2),yearsv,tax_groups_perc(:,3), ...
     yearsv,tax_groups_perc(:,4),'LineWidth',2)
yyaxis right
plot(yearsv,tax_groups_perc(:,5),'LineWidth',2)
title('Share of total taxes by taxable labor income quintile','FontSize',16)
xlabel('T model','FontSize',13)
ylabel('share','FontSize',13)
set(gca,'XTick',yearsv(1):4:yearsv(end))
xlim([yearsv(1) yearsv(end)])
legend({'bottom 20%', '2nd quintile', '3rd quintile', '4th quintile', 'top quintile'},'Location','northwest','FontSize',13)
hold off

figure
hold on
yyaxis left
plot(yearsv,tax_groups(:,1),yearsv,tax_groups(:,2),yearsv,tax_groups(:,3), ...
     yearsv,tax_groups(:,4),'LineWidth',2)
yyaxis right
plot(yearsv,tax_groups(:,5),'LineWidth',2)
title('Total taxes by taxable labor income quintile','FontSize',16)
xlabel('T model','FontSize',13)
ylabel('2016 dollars','FontSize',13)
set(gca,'XTick',yearsv(1):4:yearsv(end))
xlim([yearsv(1) yearsv(end)])
legend({'bottom 20%', '2nd quintile', '3rd quintile', '4th quintile', 'top quintile'},'Location','northwest','FontSize',13)
hold off


%% CIT DISTRIBUTION BY ASSETS

cit_a        = zeros(nk,T_model+1);     % CIT at ik asset holdings level (in dollars)
cit_a_groups = zeros( 3,T_model+1);
cit_a_perc   = zeros(nk,T_model+1);     % percentage CIT at ik asset holdings level (% of total CIT)
cit_a_pcap   = zeros(nk,T_model+1);     % CIT per capita at ik asset holdings level (in dollars)

for t = 1:T_model+1
    
   cit_t     = cit (:,:,:,:,:,t,:);
   dist_t    = DIST(:,:,:,:,:,t,:);
   citdist_t = cit_t .* dist_t;
   
   for ik = 1:nk
       
      pop_ik = dist_t(:,ik,:,:,:,:,:);
      citdist_t_ik = citdist_t(:,ik,:,:,:,:,:);
      cit_a(ik,t) = sum(citdist_t_ik(:));
      cit_a_perc(ik,t) = sum(citdist_t_ik(:))/sum(citdist_t(:));
      cit_a_pcap(ik,t) = (sum(citdist_t_ik(:))/sum(pop_ik(:)))/params.modelunit_dollar;
      
   end
   
   cit_a_groups(1,t) = sum(cit_a(1:3,t));
   cit_a_groups(2,t) = sum(cit_a(4:8,t));
   cit_a_groups(3,t) = sum(cit_a(9:12,t));
    
end

figure
plot(yearsv,cit_a_groups(1,:),yearsv,cit_a_groups(2,:),yearsv,cit_a_groups(3,:),'LineWidth',2)
title('Total CIT by asset holdings group','FontSize',16)
xlabel('T model','FontSize',13)
ylabel('model units','FontSize',13)
set(gca,'XTick',yearsv(1):4:yearsv(end))
xlim([yearsv(1) yearsv(end)])
legend({'poor', 'middle class', 'up middle class'},'Location','northwest','FontSize',13)

figure
plot(yearsv,cit_a(1,:),yearsv,cit_a(2,:),yearsv,cit_a(3,:),...
     yearsv,cit_a(4,:),yearsv,cit_a(5,:),yearsv,cit_a(6,:),...
     yearsv,cit_a(7,:),yearsv,cit_a(8,:),yearsv,cit_a(9,:),...
     'LineWidth',2)
title('Total CIT by asset holdings','FontSize',16)
xlabel('T model','FontSize',13)
ylabel('2016 dollars','FontSize',13)
set(gca,'XTick',yearsv(1):4:yearsv(end))
xlim([yearsv(1) yearsv(end)])
legend({sprintf('%0.2f', kv(1)/params.modelunit_dollar),sprintf('%0.2f', kv(2)/params.modelunit_dollar), ...
        sprintf('%0.2f', kv(3)/params.modelunit_dollar),sprintf('%0.2f', kv(4)/params.modelunit_dollar), ...
        sprintf('%0.2f', kv(5)/params.modelunit_dollar),sprintf('%0.2f', kv(6)/params.modelunit_dollar), ...
        sprintf('%0.2f', kv(7)/params.modelunit_dollar),sprintf('%0.2f', kv(8)/params.modelunit_dollar) ...
        sprintf('%0.2f', kv(9)/params.modelunit_dollar)},'FontSize',13)

figure
plot(yearsv,cit_a_pcap(1,:),yearsv,cit_a_pcap(2,:),yearsv,cit_a_pcap(3,:),...
     yearsv,cit_a_pcap(4,:),yearsv,cit_a_pcap(5,:),yearsv,cit_a_pcap(6,:),...
     yearsv,cit_a_pcap(7,:),yearsv,cit_a_pcap(8,:),yearsv,cit_a_pcap(9,:),...
     'LineWidth',2)
title('CIT per capita by asset holdings','FontSize',16)
xlabel('T model','FontSize',13)
ylabel('2016 dollars','FontSize',13)
set(gca,'XTick',yearsv(1):4:yearsv(end))
xlim([yearsv(1) yearsv(end)])
legend({sprintf('%0.2f', kv(1)/params.modelunit_dollar),sprintf('%0.2f', kv(2)/params.modelunit_dollar), ...
        sprintf('%0.2f', kv(3)/params.modelunit_dollar),sprintf('%0.2f', kv(4)/params.modelunit_dollar), ...
        sprintf('%0.2f', kv(5)/params.modelunit_dollar),sprintf('%0.2f', kv(6)/params.modelunit_dollar), ...
        sprintf('%0.2f', kv(7)/params.modelunit_dollar),sprintf('%0.2f', kv(8)/params.modelunit_dollar) ...
        sprintf('%0.2f', kv(9)/params.modelunit_dollar)},'FontSize',13)

% Export the data
header      = {'year', 'total_CIT_1', 'total_CIT_2', 'total_CIT_3', 'total_CIT_4', 'total_CIT_5', 'total_CIT_6', ...
               'total_CIT_7', 'total_CIT_8', 'total_CIT_9', 'total_CIT_10', 'total_CIT_11', 'total_CIT_12', ...
               'pcap_CIT_1', 'pcap_CIT_2', 'pcap_CIT_3', 'pcap_CIT_4', 'pcap_CIT_5', 'pcap_CIT_6', ...
               'pcap_CIT_7', 'pcap_CIT_8', 'pcap_CIT_9', 'pcap_CIT_10', 'pcap_CIT_11', 'pcap_CIT_12'};

CITDistSummary = table(yearsv', cit_a(1,:)', cit_a(2,:)', cit_a(3,:)', cit_a(4,:)', cit_a(5,:)', cit_a(6,:)', ...
                  cit_a(7,:)', cit_a(8,:)', cit_a(9,:)', cit_a(10,:)', cit_a(11,:)', cit_a(12,:)', ...
                  cit_a_pcap(1,:)', cit_a_pcap(2,:)', cit_a_pcap(3,:)', cit_a_pcap(4,:)', cit_a_pcap(5,:)', cit_a_pcap(6,:)', ...
                  cit_a_pcap(7,:)', cit_a_pcap(8,:)', cit_a_pcap(9,:)', cit_a_pcap(10,:)', cit_a_pcap(11,:)', cit_a_pcap(12,:)', ...
                  'VariableNames', header);
writetable(CITDistSummary, 'CITDistSummary.csv')           


%% PIT DISTRIBUTION BY ASSETS

%{

pit_a        = zeros(nk,T_model+1);     % PIT at ik asset holdings level (in dollars)
pit_a_groups = zeros( 3,T_model+1);
pit_a_perc   = zeros(nk,T_model+1);     % percentage PIT at ik asset holdings level (% of total PIT)
pit_a_pcap   = zeros(nk,T_model+1);     % PIT per capita at ik asset holdings level (in dollars)

for t = 1:T_model+1
    
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
   
   pit_a_groups(1,t) = sum(pit_a(1:3,t));
   pit_a_groups(2,t) = sum(pit_a(4:8,t));
   pit_a_groups(3,t) = sum(pit_a(9:12,t));
    
end

figure
plot(yearsv,pit_a_groups(1,:),yearsv,pit_a_groups(2,:),yearsv,pit_a_groups(3,:),'LineWidth',2)
title('Total PIT by asset holdings group','FontSize',16)
xlabel('T model','FontSize',13)
ylabel('model units','FontSize',13)
set(gca,'XTick',yearsv(1):4:yearsv(end))
xlim([yearsv(1) yearsv(end)])
legend({'poor', 'middle class', 'up middle class'},'Location','northwest','FontSize',13)

figure
plot(yearsv,pit_a(1,:),yearsv,pit_a(2,:),yearsv,pit_a(3,:),...
     yearsv,pit_a(4,:),yearsv,pit_a(5,:),yearsv,pit_a(6,:),...
     yearsv,pit_a(7,:),yearsv,pit_a(8,:),yearsv,pit_a(9,:),...
     'LineWidth',2)
title('Total PIT by asset holdings','FontSize',16)
xlabel('T model','FontSize',13)
ylabel('2016 dollars','FontSize',13)
set(gca,'XTick',yearsv(1):4:yearsv(end))
xlim([yearsv(1) yearsv(end)])
legend({sprintf('%0.2f', kv(1)/params.modelunit_dollar),sprintf('%0.2f', kv(2)/params.modelunit_dollar), ...
        sprintf('%0.2f', kv(3)/params.modelunit_dollar),sprintf('%0.2f', kv(4)/params.modelunit_dollar), ...
        sprintf('%0.2f', kv(5)/params.modelunit_dollar),sprintf('%0.2f', kv(6)/params.modelunit_dollar), ...
        sprintf('%0.2f', kv(7)/params.modelunit_dollar),sprintf('%0.2f', kv(8)/params.modelunit_dollar) ...
        sprintf('%0.2f', kv(9)/params.modelunit_dollar)},'FontSize',13)

figure
plot(yearsv,pit_a_pcap(1,:),yearsv,pit_a_pcap(2,:),yearsv,pit_a_pcap(3,:),...
     yearsv,pit_a_pcap(4,:),yearsv,pit_a_pcap(5,:),yearsv,pit_a_pcap(6,:),...
     yearsv,pit_a_pcap(7,:),yearsv,pit_a_pcap(8,:),yearsv,pit_a_pcap(9,:),...
     'LineWidth',2)
title('PIT per capita by asset holdings','FontSize',16)
xlabel('T model','FontSize',13)
ylabel('2016 dollars','FontSize',13)
set(gca,'XTick',yearsv(1):4:yearsv(end))
xlim([yearsv(1) yearsv(end)])
legend({sprintf('%0.2f', kv(1)/params.modelunit_dollar),sprintf('%0.2f', kv(2)/params.modelunit_dollar), ...
        sprintf('%0.2f', kv(3)/params.modelunit_dollar),sprintf('%0.2f', kv(4)/params.modelunit_dollar), ...
        sprintf('%0.2f', kv(5)/params.modelunit_dollar),sprintf('%0.2f', kv(6)/params.modelunit_dollar), ...
        sprintf('%0.2f', kv(7)/params.modelunit_dollar),sprintf('%0.2f', kv(8)/params.modelunit_dollar) ...
        sprintf('%0.2f', kv(9)/params.modelunit_dollar)},'FontSize',13)

% Export the data
header      = {'year', 'total_PIT_1', 'total_PIT_2', 'total_PIT_3', 'total_PIT_4', 'total_PIT_5', 'total_PIT_6', ...
               'total_PIT_7', 'total_PIT_8', 'total_PIT_9', 'total_PIT_10', 'total_PIT_11', 'total_PIT_12', ...
               'pcap_PIT_1', 'pcap_PIT_2', 'pcap_PIT_3', 'pcap_PIT_4', 'pcap_PIT_5', 'pcap_PIT_6', ...
               'pcap_PIT_7', 'pcap_PIT_8', 'pcap_PIT_9', 'pcap_PIT_10', 'pcap_PIT_11', 'pcap_PIT_12'};

PITDistSummary = table(yearsv', pit_a(1,:)', pit_a(2,:)', pit_a(3,:)', pit_a(4,:)', pit_a(5,:)', pit_a(6,:)', ...
                  pit_a(7,:)', pit_a(8,:)', pit_a(9,:)', pit_a(10,:)', pit_a(11,:)', pit_a(12,:)', ...
                  pit_a_pcap(1,:)', pit_a_pcap(2,:)', pit_a_pcap(3,:)', pit_a_pcap(4,:)', pit_a_pcap(5,:)', pit_a_pcap(6,:)', ...
                  pit_a_pcap(7,:)', pit_a_pcap(8,:)', pit_a_pcap(9,:)', pit_a_pcap(10,:)', pit_a_pcap(11,:)', pit_a_pcap(12,:)', ...
                  'VariableNames', header);
writetable(PITDistSummary, 'PITDistSummary.csv')           

%}

%% ASSET DISTRIBUTION BY GENERATION
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

%% FUNCTIONS

function [moments] = get_moments(dist,x)

% Computes selected percentiles of x distributed according to dist
% Inputs:  dist = distribution vector of x
%          x    = vector with variable of interest
% Outputs: percentiles (usually slightly above the actual quintiles)
%          thresholds
%          cumulative share of x with each quintile

    moments = table([], [], [], 'VariableNames', ...
                    {'percentile', 'threshold', 'cumulativeShare'});
            
    for perc = [0.2 0.4 0.6 0.8 0.95 0.99]
        moments = [moments; struct2table(get_percentile(perc,dist,x))]; %#ok<AGROW>
    end

end

function [perc_summary] = get_percentile(perc, dist, x)

% Computes percentile perc of x distributed according to dist
% Inputs:  perc = percentile between 0 and 1
%          dist = distribution vector of x
%          x    = vector with variable of interest
% Outputs: percentile (usually slightly above perc)
%          threshold
%          cummulative share of x below percentile perc

    if perc >= 1
        perc_summary = struct('percentile',      1, ...
                              'threshold',       NaN, ...
                              'cumulativeShare', 1);
        return;
    end

    [x, sortv] = sort(x);
    dist = dist(sortv);
    i = find(cumsum(dist) >= perc,1);

    perc_summary = struct('percentile',      sum(dist(1:i)), ...
                          'threshold',       x(i), ...
                          'cumulativeShare', sum(x(1:i).*dist(1:i))/sum(x.*dist));

end

