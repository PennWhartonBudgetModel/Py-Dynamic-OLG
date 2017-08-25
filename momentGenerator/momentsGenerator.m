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
T_model = 1;

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


%% Importing model generated data
% The data used here come from:
% 'M:\Repositories\danielav\Development\Testing\0.990_0.486_1.66_3.491791e-05\baseline\steady'

% Importing distribution of households
s = load('distribution.mat');
measure = s.DIST(:);
mass = s.DIST;

% Importing market variables
s = load('market.mat');
bequests = s.beqs;
wage     = s.wages;

% Importing labor policy function
s = load('all_decisions.mat');

f = @(X) repmat(reshape(X, [nz,nk,nb,T_life,1,T_model,ndem]), [1,1,1,1,ng,1,1]);

labor_inc_  = f(s.LAB) .* repmat(reshape(zs, [nz,1,1,T_life,1,1,ndem]), [1,nk,nb,1,ng,T_model,1]) * wage;
ssbenefits  = f(s.BEN);
consumption = f(s.CON);

capital     = repmat(reshape(kv, [1,nk,1,1,1,1,1]), [nz,1,nb,T_life,ng,T_model,ndem]);

labor_inc   = labor_inc_(:) + ssbenefits(:);
consumption = consumption(:);
capital     = capital    (:);

%% Compute wealth distribution

W_quintiles = get_quintile(measure,capital);
W_top = get_top(measure,capital);

%% Compute labor income distribution

L_quintiles = get_quintile(measure,labor_inc);
L_top = get_top(measure,labor_inc);

%% Compute consumption distribution

C_quintiles = get_quintile(measure,consumption);
C_top = get_top(measure,consumption);

%% Functions

function [quint_summary] = get_quintile(dist,x)

% Computes quintiles of x distributed according to dist
% Inputs:  dist = distribution vector of x
%          x    = vector with variable of interest
% Outputs: quintiles (usually slightly above the actual quintiles)
%          thresholds
%          percentage of x with each quintile

quint_summary = zeros(5,3);
counter = 0;

for perc=0.2:0.2:1.0
    counter = counter + 1;
    quint_summary(counter,:) = get_percentile(perc,dist,x);
end
for counter = size(quint_summary,1):-1:2
	quint_summary(counter,3) = quint_summary(counter,3) - quint_summary(counter-1,3);
end

end

function [top_summary] = get_top(dist,x)

% Computes distribution of x  at the top
% Inputs:  dist = distribution vector of x
%          x    = vector with variable of interest
% Outputs: top percentiles (usually slightly below the actual quintiles)
%          thresholds
%          percentage of x with each top percentile

top_summary = zeros(5,3);
counter = 0;
for top = [0.9 0.95 0.99 0.999 .9999]
    counter = counter + 1;
    top_summary(counter,:) = get_percentile(top,dist,x);
    top_summary(counter,1) = 1 - top_summary(counter,1);
    top_summary(counter,3) = 1 - top_summary(counter,3);
end

end

function [perc_summary] = get_percentile(perc, dist, x)

% Computes percentile perc of x distributed according to dist
% Inputs:  perc = percentile between 0 and 1
%          dist = distribution vector of x
%          x    = vector with variable of interest
% Outputs: percentile (usually slightly above perc)
%          threshold
%          cummulative percentage of x below percentile perc

[x, sortv] = sort(x);
dist = dist(sortv);

i = find(cumsum(dist) >= perc,1);

perc_summary = [sum(dist(1:i)) x(i) (sum(x(1:i).*dist(1:i))/sum(x.*dist))];

end





