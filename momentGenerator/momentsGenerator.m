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
income      = f(s.INC);

assets = repmat(reshape(kv, [1,nk,1,1,1,1,1]), [nz,1,nb,T_life,ng,T_model,ndem]);
agesv  = repmat(reshape([1:T_life], [1,1,1,T_life,1,1,1]), [nz,nk,nb,1,ng,T_model,ndem]);

labor_inc   = labor_inc_ (:) + ssbenefits(:);
consumption = consumption(:);
assets      = assets     (:);
income      = income     (:);
agesv       = agesv      (:);

%% Compute assets distribution

a_quintiles = get_quintile(measure,assets);
a_top = get_top(measure,assets);
a_zero = zero_x(measure,assets);

measure1 = filter_by_age(15,25,agesv,measure);
a_quintiles1 = get_quintile(measure1,assets);

[a_ginicoeff a_lorenz] = gini(measure,assets,true);

%% Compute labor income distribution

l_quintiles = get_quintile(measure,labor_inc);
l_top = get_top(measure,labor_inc);

[l_ginicoeff l_lorenz] = gini(measure,labor_inc);

%% Compute consumption distribution

c_quintiles = get_quintile(measure,consumption);
c_top = get_top(measure,consumption);

%% Compute total income distribution

i_quintiles = get_quintile(measure,income);
i_top = get_top(measure,income);

%% Graphs

% figure
% bar([1:5],a_quintiles(:,4))
% title('Assets distribution')
% xlabel('quintiles')
% ylabel('fraction of total assets')

% figure
% bar([1:5],a_top(:,4))
% title('Assets distribution in the top')
% Labels = {'10%', '6%', '5%', '1%', '0.1%'};
% set(gca, 'XTick', 1:5, 'XTickLabel', Labels);
% xlabel('top percentile')
% ylabel('fraction of held by the top')

% figure
% bar([1:5],l_quintiles(:,4))
% title('Labor income distribution')
% xlabel('quintiles')
% ylabel('fraction of total labor income')

% figure
% bar([1:5],i_quintiles(:,4))
% title('Total income distribution')
% xlabel('quintiles')
% ylabel('fraction of total income')

%% Functions

function [quint_summary] = get_quintile(dist,x)

% Computes quintiles of x distributed according to dist
% Inputs:  dist = distribution vector of x
%          x    = vector with variable of interest
% Outputs: quintiles (usually slightly above the actual quintiles)
%          thresholds
%          cumulative share of x with each quintile
%          share of x with each quintile

quint_summary = zeros(5,4);
counter = 0;

for perc=0.2:0.2:1.0
    counter = counter + 1;
    quint_summary(counter,1:3) = get_percentile(perc,dist,x);
end
quint_summary(1,4) = quint_summary(1,3);
for counter = size(quint_summary,1):-1:2
	quint_summary(counter,4) = quint_summary(counter,3) - quint_summary(counter-1,3);
end

end

function [top_summary] = get_top(dist,x)

% Computes distribution of x  at the top
% Inputs:  dist = distribution vector of x
%          x    = vector with variable of interest
% Outputs: top percentiles (usually slightly below the actual quintiles)
%          thresholds
%          share of x with each top percentile

top_summary = zeros(4,3);
counter = 0;
for top = [0.9 0.94 0.95 0.99]
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
%          cummulative share of x below percentile perc

[x, sortv] = sort(x);
dist = dist(sortv);

i = find(cumsum(dist) >= perc,1);

perc_summary = [sum(dist(1:i)) x(i) (sum(x(1:i).*dist(1:i))/sum(x.*dist))];

end

function [dist1] = filter_by_age(agemin,agemax,agesv,dist0)

% Inputs:  [agemin,agemax] = age interval of interest
%          agesv = vector of ages with the same dimension of dist0
%          dist0 = distribution vector of population at each type
% Outputs: dist1 = new distribution within age interval

assert(agemin <= agemax, 'Wrong age interval: agemin > agemax.')
dist0(agesv < agemin | agesv > agemax) = 0;
dist1 = dist0./sum(dist0);

end

function [perc_x0] = zero_x(dist,x)

% Inputs:  dist = distribution vector of x
%          x    = vector with variable of interest
% Outputs: perc_x0 = percentage of population in first point of K grid

i = find(x == 0.0010,1);
perc_x0 = sum(dist(i));

end

function [ginicoeff lorenz] = gini(dist,x,makeplot)

% Inputs:  dist     = vector of population sizes of the different types
%                     (note it doesn't need to be distribution vector of x)
%          x        = vector with variable of interest
%          makeplot = boolean indicating whether a figure of the Lorentz
%                     curve should be produced or not. Default is false.
% Outputs: ginicoeff = Gini coefficients
%          lorenz    = Lorentz curve: This is a two-column array, with the 
%                      left column representing cumulative population 
%                      shares of the different classes, sorted according to
%                      val, and the right column representing the 
%                      cumulative value share that belongs to the 
%                      population up to the given class. The Lorentz curve 
%                      is a scatter plot of the left vs the right column.

if nargin < 3
	makeplot = false;
end

% Pre-append a zero because the Lorenz curve contains (0,0) by definition
dist = [0;dist(:)]; x = [0;x(:)];

% Check for negative values - this is important if there's borrowing
assert(all(dist>=0) && all(x>=0), ...
	'gini expects nonnegative vectors (neg elements in pop = %d, in val = %d).', ...
	sum(dist<0),sum(x<0))
    
% Sort in ascending order wrt x and weight vectors
z = x .* dist;
[~,ord] = sort(x);
dist    = dist(ord);      z    = z(ord);
dist    = cumsum(dist);   z    = cumsum(z);
relpop  = dist/dist(end); relz = z/z(end);

% We compute the area below the Lorentz curve. We do this by computing the 
% average of the left and right Riemann-like sums. (Riemann-'like' because 
% we evaluate not on a uniform grid, but on the points given by the pop data).
% These are the two Rieman-like sums:
    %    leftsum  = sum(relz(1:end-1) .* diff(relpop));
    %    rightsum = sum(relz(2:end)   .* diff(relpop));
% The Gini coefficient is one minus twice the average of leftsum and
% rightsum. We can put all of this into one line.

ginicoeff = 1 - sum((relz(1:end-1)+relz(2:end)) .* diff(relpop));

% Lorentz curve
lorenz = [relpop,relz];

% Plot
if makeplot
	area(relpop,relz,'FaceColor',[0.5,0.5,1.0]); % the Lorentz curve
	hold on
	plot([0,1],[0,1],'--k');                     % 45 degree line
	axis tight                                   % ranges of abscissa and ordinate are by definition exactly [0,1]
	axis square                                  % both axes should be equally long
	set(gca,'XTick',get(gca,'YTick'))            % ensure equal ticking
	set(gca,'Layer','top');                      % grid above the shaded area
	grid on;
	title(['\bfGini coefficient = ',num2str(ginicoeff)]);
	xlabel('share of population');
	ylabel('share of value');
end

end