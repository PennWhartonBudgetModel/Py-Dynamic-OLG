%%
% Moments generator tool.
% 
%%

nz = 4;
nk = 10;
f = @(lb, ub, n) lb + (ub-lb)*((0:n-1)/(n-1))'.^2;
kv = f(1e-3, 120, nk);
nb = 5;
T_life = 80;
ng = 3;
nidem = 2;
T_work = 47;

%% Importing model generated data
% The data used here come from:
% 'M:\Repositories\danielav\Development\Testing\0.990_0.486_1.66_3.491791e-05\baseline\steady'

% Importing distribution of households
s = load('distribution.mat');
measure = zeros(nz,nk,nb,T_life,ng,1,nidem);
measure = s.DIST;

% Importing labor policy function
s = load('decisions.mat');
labor = zeros(nz,nk,nb,T_life,nidem);
labor(:,:,:,:,1) = s.LABs{1,1};
labor(:,:,:,:,2) = s.LABs{1,2};
consumption = s.OPTs_CON;

% Importing market variables
s = load('market.mat');
bequests = s.beqs;
wage = s.wages;

%% Creating wealth distribution array

W_measure_aux = zeros(nz*nk*nb*T_life*ng*nidem,2);

counter = 0;

for age=1:T_life
    for iz=1:nz
        for ik=1:nk
            for ib=1:nb
                for ig=1:ng
                    for iperm=1:nidem
                        counter = counter + 1;
                        W_measure_aux(counter,1) = measure(iz,ik,ib,age,ig,1,iperm);
                        W_measure_aux(counter,2) = kv(ik);
                    end
                end
            end
        end
    end
end

% Sorting in ascending order wrt assets
W_measure = sortrows(W_measure_aux,2);

%% Compute wealth quintiles

W_quintiles = zeros(5,3);
counter = 0;

for perc=0.2:0.2:1.0
    counter = counter + 1;
    W_quintiles(counter,:) = get_percentile(perc,W_measure);
    if counter > 1
        W_quintiles(counter,3) = W_quintiles(counter,3) - sum(W_quintiles(1:counter-1,3));
    end
end

%% Compute wealth distribution at the top

W_top = zeros(5,3);

W_top(1,:) = get_percentile(0.9,W_measure);
W_top(2,:) = get_percentile(0.95,W_measure);
W_top(3,:) = get_percentile(0.99,W_measure);
W_top(4,:) = get_percentile(0.999,W_measure);
W_top(5,:) = get_percentile(0.99999,W_measure);

% Figure 1
figure
subplot(1,2,1)
plot([1:5],W_quintiles(:,2),[1:5],W_top(:,2))
title('Thresholds per quintile')
subplot(1,2,2)
plot([1:5],W_quintiles(:,3),[1:5],W_top(:,3))
title('Total assets per quintile')
suptitle('Wealth distribution')

%% Detour before automatically importing shock matrix

 f = @(n, sig) n*sig*-diff(normpdf(norminv(linspace(0, 1, n+1))));
        
        % Determine permanent and transitory shocks
        nperm  = 2; zperm  = f(nperm , sqrt(0.2105));
        ntrans = 2; ztrans = f(ntrans, sqrt(0.0630)); 
        
        % Determine persistent shocks
        npers = 2;
        
        pep = 0.973;                        % Lagged productivity coefficient
        sigep = sqrt(0.018);                % Conditional standard deviation of productivity
        sigpers = sigep/(sqrt(1-pep^2));
        
        zpers = f(npers, sigpers);
        
        % Construct Markov transition matrix for persistent shocks by approximating an AR(1) process
        persv = sigpers*norminv(linspace(0, 1, npers+1));
        transpers = zeros(npers,npers);
        for ipers = 1:npers
            integrand = @(x) exp(-x^2/(2*sigpers^2)) * diff(normcdf((persv(ipers:ipers+1) - pep*x)/sigep));
            for jpers = 1:npers
                transpers(ipers,jpers) = (npers/(sqrt(2*pi*(sigpers^2)))) * quadgk(@(v) arrayfun(integrand, v), persv(jpers), persv(jpers+1));
            end
        end
        
        % Determine initial distribution over persistent shocks
        DISTpers = diff(normcdf(persv/sqrt(0.124)));
        
        % Define deterministic lifecycle productivities
        % (Estimated coefficients from Barro and Barnes Medicaid working paper)
        zage = polyval([-5.25e-7, 1.05e-4, -8.1467e-3, 0.2891379, -1.203521], 19+(1:T_life));
        
        % Calculate total productivity
        ndem = nperm; nz = ntrans*npers;
        
        zs = max(0, repmat(reshape(zage                       , [1 ,T_life,1   ]), [nz,1     ,ndem]) ...
                  + repmat(reshape(zperm                      , [1 ,1     ,ndem]), [nz,T_life,1   ]) ...
                  + repmat(reshape(kron(ones(1,npers), ztrans) ...
                                 + kron(zpers, ones(1,ntrans)), [nz,1     ,1   ]), [1 ,T_life,ndem]));

%% Creating labor earnings distribution array

L_measure_aux = zeros(nz*nk*nb*T_work*ng*nidem,2);

counter = 0;

for age=1:T_work
    for iz=1:nz
        for ik=1:nk
            for ib=1:nb
                for ig=1:ng
                    for iperm=1:nidem
                        counter = counter + 1;
                        L_measure_aux(counter,1) = measure(iz,ik,ib,age,ig,1,iperm);
                        L_measure_aux(counter,2) = wage*zs(iz,age,iperm)*labor(iz,ik,ib,age,iperm);
                    end
                end
            end
        end
    end
end

% Sorting in ascending order wrt assets
L_measure = sortrows(L_measure_aux,2);

% Since we do not have SS income, we do not have a measure 1 population
L_measure(:,1) = (L_measure(:,1)/sum(L_measure(:,1)));

%% Compute wealth quintiles

L_quintiles = zeros(5,3);
counter = 0;

for perc=0.2:0.2:1.0
    counter = counter + 1;
    L_quintiles(counter,:) = get_percentile(perc,L_measure);
    if counter > 1
        L_quintiles(counter,3) = L_quintiles(counter,3) - sum(L_quintiles(1:counter-1,3));
    end
end

%% Compute wealth distribution at the top

L_top = zeros(5,3);

L_top(1,:) = get_percentile(0.9,L_measure);
L_top(2,:) = get_percentile(0.95,L_measure);
L_top(3,:) = get_percentile(0.99,L_measure);
L_top(4,:) = get_percentile(0.999,L_measure);
L_top(5,:) = get_percentile(0.999999,L_measure);

figure
subplot(1,2,1)
plot([1:5],L_quintiles(:,2),[1:5],L_top(:,2))
title('Thresholds per quintile')
subplot(1,2,2)
plot([1:5],L_quintiles(:,3),[1:5],L_top(:,3))
title('Total assets per quintile')
suptitle('Labor income distribution')

%% Functions

function [perc_summary] = get_percentile(percentile, matrix)

temp_mass = 0;
cnt = 0;
check = 0;
epsilon = 1.0e-06;

while cnt <= size(matrix,1)
    cnt = cnt + 1;
    temp_mass = temp_mass + matrix(cnt,1);
    if temp_mass >= percentile && check == 0 && percentile ~= 1
        check = 1;
        perc_threshold = matrix(cnt,2);
        perc_total     = sum(matrix(1:cnt,2));
        break
    end
    if temp_mass >= (percentile - epsilon) && check == 0 && percentile == 1
        check = 1;
        perc_threshold = matrix(cnt,2);
        perc_total     = sum(matrix(1:cnt,2));
        break
    end
end

percentile = temp_mass;
perc_summary = [percentile perc_threshold perc_total];

end





