%%
% Distribution moments generator tool.
% 
%%

classdef momentsGenerator
    
    properties (Access = private)
        
        modelunit_dollar;
        a_distdata; a_distmodel; a_ginimodel; a_lorenz;
        l_distdata; l_distmodel; l_ginimodel; l_lorenz;
        DIST; T_life; karray;

    end
    

    methods
       
        % Constructor
        function this = momentsGenerator(scenario,DIST,Market,OPTs)
            
            if( ~strcmp(scenario.economy, 'steady' ) )
                error('Unable to generate income distribution moments for transition paths.')
            end

            %% PARAMETERS
    
            this.modelunit_dollar = scenario.modelunit_dollar;
            save_dir = Environment.getCurrent().save(scenario);
            param_dir  = Environment.sim_param();

            % Define time constants
            s       = paramGenerator.timing(scenario);
            T_life  = s.T_life;    % Total life years
            T_work  = s.T_work;    % Retirement age
            T_model = s.T_model;   % Transition path model years

            % Define grids
            s    = paramGenerator.grids(scenario);
            ndem = s.ndem;       % demographic types
            ng   = s.ng;         % num groups
            nz   = s.nz;         % num labor productivity shocks
            zs   = s.zs;         % shocks grid (by demographic type and age)
            nk   = s.nk;         % num asset points
            nb   = s.nb;         % num avg. earnings points
            % Useful later for a couple of functions
            this.karray = repmat(reshape(s.kv, [1,nk,1,1,1,1,1]),[nz,1,nb,T_life,ng,T_model,ndem])  ;
            this.T_life = T_life;

            %% DISTRIBUTION AND POLICY FUNCTIONS

            % Import households distribution
            if ~exist('DIST','var') || isempty(DIST)
                s    = load( fullfile(save_dir, 'distribution.mat' ) );
                DIST = s.DIST;
            end
            dist = DIST(:);
            dist_l(1:nz,1:nk,1:nb,1:T_work,1:ng,1:T_model,1:ndem) = DIST(1:nz,1:nk,1:nb,1:T_work,1:ng,1:T_model,1:ndem); % Working age population
            dist_l(1:nz,1:nk,1:nb,T_work:T_life,1:ng,1:T_model,1:ndem) = 0; % Retired population
            dist_l = dist_l(:)/sum(dist_l(:));
            % Useful later for a couple of functions
            this.DIST   = DIST;

            % Import market variables
            if ~exist('Market','var') || isempty(Market)
                s     = load( fullfile(save_dir, 'market.mat' ) );
                wages = s.wages;
            else
                wages = Market.wages;
            end
            
            % Import policy functions
            f = @(X) repmat(reshape(X, [nz,nk,nb,T_life,1,T_model,ndem]), [1,1,1,1,ng,1,1]);
            if ~exist('OPTs','var') || isempty(OPTs)
                s      = load( fullfile(save_dir, 'all_decisions.mat' ) );
                labinc = f(s.LAB) .* repmat(reshape(zs, [nz,1,1,T_life,1,1,ndem]),[1,nk,nb,1,ng,T_model,1]) * wages;
                k      = f(s.K);
            else
                labinc = f(OPTs.LAB) .* repmat(reshape(zs, [nz,1,1,T_life,1,1,ndem]),[1,nk,nb,1,ng,T_model,1]) * wages;
                k      = f(OPTs.K);
            end
            labinc = labinc(:);  % Labor income
            k      = k     (:);  % Asset holdings for tomorrow (k')
    
            %% DATA WEALTH AND INCOME DISTRIBUTIONS

            this.a_distdata = readtable(fullfile(param_dir, 'SIM_NetPersonalWealth_distribution.csv'));
            this.a_distdata = [this.a_distdata; {100 NaN 1}];       % Append last point for graph

            this.l_distdata = readtable(fullfile(param_dir, 'SIM_PreTaxLaborInc_distribution.csv'));
            this.l_distdata = [this.l_distdata; {100 NaN 1}];       % Append last point for graph

            %% MODEL WEALTH AND INCOME DISTRIBUTIONS

            % Compute wealth distribution
            this.a_distmodel = get_moments(dist,k);
            % Gini and Lorenz curve
            [this.a_ginimodel this.a_lorenz] = gini(dist,k);

            % Compute labor income distribution
            this.l_distmodel = get_moments(dist_l,labinc);
            % Gini and Lorenz curve
            [this.l_ginimodel this.l_lorenz] = gini(dist_l,labinc);


        end
        
        % Table - data and model Gini and the gap between them
        function [gini_table] = giniTable(this)
            
            a_ginidata = 0.857;                         % Number from SIM
            l_ginidata = 0.4858;                        % Number from SIM

            % Compare model distribution with data
            a_ginigap = 100*(a_ginidata / this.a_ginimodel - 1);

            % Compare model distribution with data
            l_ginigap = 100*(l_ginidata / this.l_ginimodel - 1);

            % momentsGenerator output
            gini_table = table(categorical({'wealth';'lab_income'}),...
                                 [a_ginidata; l_ginidata],...
                                 [this.a_ginimodel; this.l_ginimodel],...
                                 [a_ginigap; l_ginigap],...
                                 'VariableNames',{'Gini' 'data' 'model' 'percent_gap'});
                     
        end
        
        % Table - assets and labor income cumulative shares for selected percentiles
        function [cumShare_table] = cumShareTable(this)
            
            cumShare_table = table(this.a_distdata.percentile,...
                              this.a_distdata.cumulativeShare, this.a_distmodel.cumulativeShare,...
                              this.l_distdata.cumulativeShare, this.l_distmodel.cumulativeShare,...
                              'VariableNames',{'Percentile' 'a_data' 'a_model' ...
                                                            'l_data' 'l_model'});
                     
        end
        
        % Graph - Assets distribution: model vs. data
        function [] = plot_a_lorenz(this)
            
            datapercentile = [0; this.a_distdata.percentile];
            datacumShare   = [min(this.a_distdata.cumulativeShare)-1e-8; this.a_distdata.cumulativeShare];
                         
            figure
            plot(this.a_lorenz(:,1).*100, this.a_lorenz(:,2), ...
                 datapercentile         , datacumShare, 'LineWidth',2)
            hold on
            plot([0,100],[0,1],'--k','LineWidth',1.5);     	% 45 degree line
            title('Assets distribution','FontSize',16)
            xlabel('percentiles','FontSize',13)
            ylabel('cumulative share of total assets','FontSize',13)
            legend({sprintf('model (gini = %0.3f)', this.a_ginimodel),
                    'data    (gini = 0.857)' },'Location','northwest','FontSize',13)
            
        end

        % Graph - Labor income distribution: model vs. data
        function [] = plot_l_lorenz(this)
                         
            datapercentile = [0; this.l_distdata.percentile];
            datacumShare   = [0; this.l_distdata.cumulativeShare];

            figure
            plot(this.l_lorenz(:,1).*100, this.l_lorenz(:,2), ...
                 datapercentile         , datacumShare, 'LineWidth',2)
            hold on
            plot([0,100],[0,1],'--k','LineWidth',1.5);     	% 45 degree line
            title('Labor income distribution','FontSize',16)
            xlabel('percentiles','FontSize',13)
            ylabel('cumulative share of total labor income','FontSize',13)
            legend({sprintf('model (gini = %0.4f)', this.l_ginimodel),
                    'data    (gini = 0.4858)' },'Location','northwest','FontSize',13)
            
        end
        
        % Graph - Assets threshold in dollars: model vs. data
        function [] = plot_a_threshold(this)
                         
            figure
            plot(this.a_distdata.percentile,(this.a_distmodel.threshold/this.modelunit_dollar)/1000,...
                 this.a_distdata.percentile,this.a_distdata.threshold2016dollars/1000, 'LineWidth',2)
            title('Threshold by wealth percentile','FontSize',16)
            xlabel('percentiles','FontSize',13)
            ylabel('thousands of 2016 dollars','FontSize',13)
            legend({sprintf('model (gini = %0.3f)', this.a_ginimodel),
                    'data    (gini = 0.857)' },'Location','northwest','FontSize',13)
                        
        end
        
        % Graph - Labor income threshold in dollars: model vs. data
        function [] = plot_l_threshold(this)
                         
            figure
            plot(this.l_distdata.percentile,(this.l_distmodel.threshold/this.modelunit_dollar)/1000,...
                 this.l_distdata.percentile,this.l_distdata.threshold2016dollars/1000, 'LineWidth',2)
            title('Threshold by labor income percentile','FontSize',16)
            xlabel('percentiles','FontSize',13)
            ylabel('thousands of 2016 dollars','FontSize',13)
            legend({sprintf('model (gini = %0.4f)', this.l_ginimodel),
                    'data    (gini = 0.4858)' },'Location','northwest','FontSize',13)
                        
        end
        
        % Table - Age distribution at the bottom of the capital grid
        function [topBottom_table] = topBottomTable(this)
            
            bot_g1 = this.DIST(:,1,:, 1:19,:,:,:);
            bot_g2 = this.DIST(:,1,:,20:39,:,:,:);
            bot_g3 = this.DIST(:,1,:,40:59,:,:,:);
            bot_g4 = this.DIST(:,1,:,60:80,:,:,:);
            bottom = this.DIST(:,1,:,:,:,:,:);
            bottom = sum(bottom(:));
            bot_share = [sum(bot_g1(:)) sum(bot_g2(:)) sum(bot_g3(:)) sum(bot_g4(:))];
            
            top_g1 = this.DIST(:,15,:, 1:19,:,:,:);
            top_g2 = this.DIST(:,15,:,20:39,:,:,:);
            top_g3 = this.DIST(:,15,:,40:59,:,:,:);
            top_g4 = this.DIST(:,15,:,60:80,:,:,:);
            top    = this.DIST(:,15,:,:,:,:,:);            
            top    = sum(top(:));
            top_share = [sum(top_g1(:)) sum(top_g2(:)) sum(top_g3(:)) sum(top_g4(:))];
                                                
            topBottom_table = table(categorical({'bottom'; 'top'}),[bottom; top],...
                              [bot_share(1); top_share(1)],[bot_share(2); top_share(2)],...
                              [bot_share(3); top_share(3)],[bot_share(4); top_share(4)],...
                              'VariableNames',{'grid' 'total' 'age21to40' 'age41to60' 'age61to80' 'age81to101'});

        end
        
        % Graph - Asset holdings by age
        function [] = plot_a_age(this)
            
            kdist_age = zeros(1,this.T_life);
            kdist = this.DIST .* this.karray;
            for age = 1:this.T_life
                pop_age_temp   = this.DIST(:,:,:,age,:,:,:);
                kdist_age_temp = kdist(:,:,:,age,:,:,:);
                kdist_age(age) = (sum(kdist_age_temp(:))/sum(pop_age_temp(:)))/this.modelunit_dollar;
            end
            
            figure
            plot(21:100,kdist_age, 'LineWidth',2)
            title('Average asset holdings by age','FontSize',16)
            xlabel('age','FontSize',13)
            ylabel('2016 dollars','FontSize',13)
            
        end

        % Graph - Distribution of individuals by asset holdings grid points
        function [] = plot_a_dist(this)
            
            nk     = size(this.DIST,2);
            kdist  = zeros(1,nk);
            for ik = 1:nk
                kdist_temp = this.DIST(:,ik,:,:,:,:,:);
                kdist(ik)  = sum(kdist_temp(:));
            end
            
            figure
            plot(1:nk,kdist, 'LineWidth',2)
            title('Distribution of individuals by asset holdings grid points','FontSize',16)
            xlabel('grid point number','FontSize',13)
            ylabel('share of population','FontSize',13)
            
        end
                
    end % methods
    
end % momentsGenerator

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
%                     (note it doesn't need to be distribution vector of x, but
%                      it cannot be a cumulative distribution)
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
