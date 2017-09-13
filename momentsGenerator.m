%%
% Distribution moments generator tool.
% 
%%

function [gini_summary] = momentsGenerator(scenario,do_plot1,do_plot2)

    if (nargin == 1)
        do_plot1 = false;
        do_plot2 = false;
    elseif (nargin == 2)
        do_plot2 = false;
    end

    if( ~strcmp(scenario.economy, 'steady' ) )
        error('Unable to generate income distribution moments for transition paths.')
    end
    
	%% PARAMETERS
    
    save_dir = dirFinder.save(scenario);

	% Define time constants
	s       = paramGenerator.timing(scenario);
	T_life  = s.T_life;    % Total life years
	T_model = s.T_model;   % Transition path model years

	% Discretized grids, including shock process
	s    = paramGenerator.grids(scenario);
	ndem = s.ndem;       % demographic types
	ng   = s.ng;         % num groups
	nz   = s.nz;         % num labor productivity shocks
	zs   = s.zs;         % shocks grid (by demographic type and age)
	nk   = s.nk;         % num asset points
 	kv   = s.kv;         % assets grid
	nb   = s.nb;         % num avg. earnings points

    %% DISTRIBUTION AND POLICY FUNCTIONS

    % Importing distribution of households
    s    = load( fullfile(save_dir, 'distribution.mat' ) );
    dist = s.DIST(:);
    
    % Importing market variables
    s    = load( fullfile(save_dir, 'market.mat' ) );
    beq  = s.beqs;
    wage = s.wages;
    
    % Importing policy functions
    f = @(X) repmat(reshape(X, [nz,nk,nb,T_life,1,T_model,ndem]), [1,1,1,1,ng,1,1]);
    s = load( fullfile(save_dir, 'all_decisions.mat' ) );
    labinc = f(s.LAB) .* repmat(reshape(zs, [nz,1,1,T_life,1,1,ndem]),[1,nk,nb,1,ng,T_model,1]) * wage;
    cons   = f(s.CON);
    k      = f(s.K);
    labinc = labinc(:);  % Labor income
    cons   = cons  (:);  % Consumption
    k      = k     (:);  % Asset holdings for tomorrow (k')
    
    %% Import wealth and labor earnings distribution moments
    
    a_distdata = readtable(fullfile(dirFinder.param(), ...
                 'SIM_NetPersonalWealth_distribution.csv'));
    a_distdata = [a_distdata; {100 NaN 1}];     % Append last point for graph
    
    a_ginidata = 0.857;                         % Number from SIM

    l_distdata = readtable(fullfile(dirFinder.param(), ...
                 'SIM_PreTaxLaborInc_distribution.csv'));
    l_distdata = [l_distdata; {100 NaN 1}];     % Append last point for graph
    
    l_ginidata = 0.547;                         % Number from SIM

    %% WEALTH AND INCOME DISTRIBUTIONS

    % Compute wealth distribution
    a_distmodel = get_moments(dist,k);
    
    % Gini and Lorenz curve
    [a_ginimodel a_lorenz] = gini(dist,k,false); % set 3rd argument to true to graph the Lorenz curve
    
    % Compare model distribution with data
    a_ginigap = 100*(a_ginidata / a_ginimodel - 1);
    
    if do_plot1
        figure
        plot(a_distdata.percentile,a_distmodel.cumulativeShare, ...
            a_distdata.percentile,a_distdata.cumulativeShare)
        title('Assets distribution')
        xlabel('percentiles')
        ylabel('cumulative share of total assets')
        legend('model','data','Location','northwest')
    end

    % Compute labor income distribution
    l_distmodel = get_moments(dist,labinc);
    
    % Gini and Lorenz curve
    [l_ginimodel l_lorenz] = gini(dist,labinc);

    % Compare model distribution with data
    l_ginigap = 100*(l_ginidata / l_ginimodel - 1);
    
    if do_plot2
        figure
        plot(l_distdata.percentile,l_distmodel.cumulativeShare,...
            l_distdata.percentile,l_distdata.cumulativeShare)
        title('Labor income distribution')
        xlabel('percentiles')
        ylabel('cumulative share of total labor income')
        legend('model','data','Location','northwest')
    end
    
    % momentsGenerator output
    gini_summary = table(categorical({'wealth';'lab_income'}),...
                         [a_ginidata; l_ginidata],...
                         [a_ginimodel; l_ginimodel],...
                         [a_ginigap; l_ginigap],...
                         'VariableNames',{'Gini' 'data' 'model' 'percent_gap'});
                     
                     
end

%%

function [moments] = get_moments(dist,x)

% Computes selected percentiles of x distributed according to dist
% Inputs:  dist = distribution vector of x
%          x    = vector with variable of interest
% Outputs: percentiles (usually slightly above the actual quintiles)
%          thresholds
%          cumulative share of x with each quintile

    moments = table([], [], [], 'VariableNames', ...
                    {'percentile', 'threshold', 'cumulativeShare'});
            
    for perc = [0.2 0.4 0.6 0.8 0.9 0.95 0.99 1]
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

    % Check for negative population/measure
    assert(all(dist>=0), ...
           'gini expects nonnegative population vector (%d negative elements).', ...
           sum(dist<0))
    
    % Take care of first point of the Lorenz curve
    if all(x>0)
        % Pre-append a zero because the Lorenz curve contains (0,0) by definition
        dist = [0;dist(:)]; x = [0;x(:)];
    else
        % Use the lowest x holdings to set the first point
        dist = [0;dist(:)]; x = [min(x)-1e-8;x(:)];
    end

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
