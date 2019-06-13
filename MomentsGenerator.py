##
# Distribution moments generator tool.

import numpy as np
import scipy.io as sio
import os
import matplotlib.pyplot as plt
import pandas as pd
import math
import PathFinder
import Scenario
import ParamGenerator

class MomentsGenerator:
    
    scenario = None
    a_distdata = None
    a_distmodel = None
    a_ginimodel = None
    a_lorenz = None
    l_distdata = None
    l_distmodel = None
    l_ginimodel = None
    l_lorenz = None
    DIST = None
    T_work = None
    T_life = None
    kv = None
    karray = None
    ben = None
    con = None
    lab = None
    
    # Constructor
    def __init__(self,scenario,DIST,Market,OPTs):
        
        if not scenario.isSteady():
            raise Exception('Unable to generate income distribution moments for transition paths.')
            
            # PARAMETERS
            pathFinder = PathFinder(scenario)
            
            self.scenario   = scenario
            save_dir        = PathFinder.getCacheDir(scenario)
            
            # Define time constants and grids
            timing      = ParamGenerator.timing(scenario)
            grids       = ParamGenerator.grids(scenario)
            T_life      = timing['T_life']        # Total life years
            T_model     = timing['T_model']       # Transition path model years
            Tmax_work   = timing['Tmax_work']     # Largest retirement age
            ng          = grids['ng']             # num groups
            nz          = grids['nz']             # num labor productivity shocks
            zs          = grids['zs']             # shocks grid (by demographic type and age)
            nk          = grids['nk']             # num asset points
            nb          = grids['nb']             # num avg. earnings points
            
            # Useful later for a couple of functions
            self.kv     = grids.kv
            self.karray = np.tile(np.reshape(grids['kv'], [1,nk,1,1,1,1]),[nz,1,nb,T_life,ng,T_model])
            self.T_work = Tmax_work
            self.T_life = T_life
            
            ## DISTRIBUTION AND POLICY FUNCTIONS
            
            # Import households distribution
            if (not 'DIST' in globals()) or DIST == None:
                s = sio.loadmat(os.path.join(save_dir, 'distribution.mat' ))
                DIST = s['DIST']
            dist = DIST[:]
            dist_l = [[[[[[]]]]]]
            dist_l[0:nz,0:nk,0:nb,0:Tmax_work,0:ng,0:T_model] = DIST[0:nz,0:nk,0:nb,0:Tmax_work,0:ng,0:T_model] # Working age population
            dist_l[0:nz,0:nk,0:nb,Tmax_work-1:T_life,0:ng,0:T_model] = 0 # Retired population
            dist_l = dist_l[:]/np.sum(dist_l[:])
            # Useful later for a couple of functions
            self.DIST = DIST
            
            # Import market variables
            if (not 'Market' in globals()) or Market == None:
                s = sio.loadmat(os.path.join(save_dir, 'market.mat'))
                wages = s['wages']
            else:
                wages = Market['wages']
            
            # Import policy functions
            f = lambda X: np.tile(np.reshape(X, [nz,nk,nb,T_life,1,T_model]), [1,1,1,1,ng,1])
            if (not 'OPTs' in globals()) or (OPTs == None):
                s = sio.loadmat(os.path.join(save_dir, 'decisions.mat'))
                s = s['OPTs']
                labinc   = f(s['LABOR']) * np.tile(np.reshape(np.transpose(zs, [3,2,1]), [nz,1,1,T_life,1,T_model]),[1,nk,nb,1,ng,1]) * wages
                k        = f(s['SAVINGS'])
                self.ben = f(s['OASI_BENEFITS'])
                self.con = f(s['CONSUMPTION'])
            else:
                labinc   = f(OPTs['LABOR']) * np.tile(np.reshape(np.transpose(zs, [3,2,1]), [nz,1,1,T_life,1,T_model]),[1,nk,nb,1,ng,1]) * wages
                k        = f(OPTs['SAVINGS'])
                self.ben = f(OPTs['OASI_BENEFITS'])
                self.lab = f(OPTs['LABOR'])
                self.con = f(OPTs['CONSUMPTION'])
            labinc = labinc[:]  # Labor income
            k      = k     [:]  # Asset holdings for tomorrow (k')
            
            # DATA WEALTH AND INCOME DISTRIBUTIONS
            file = pathFinder.getMicrosimInputPath('SIM_NetPersonalWealth_distribution')
            
            self.a_distdata = pd.from_csv(file)
            self.a_distdata.append([100, None, 1])      # Append last point for graph
            
            file = pathFinder.getMicrosimInputPath('SIM_PreTaxLaborInc_distribution')
            self.l_distdata = pd.from_csv(file)
            self.l_distdata.append([100, None, 1])       # Append last point for graph
            
            # MODEL WEALTH AND INCOME DISTRIBUTIONS
            
            # Compute wealth distribution
            self.a_distmodel = get_moments(dist,k)
            # Gini and Lorenz curve
            (self.a_ginimodel, self.a_lorenz) = gini(dist,k)
            
            # Compute labor income distribution
            self.l_distmodel = get_moments(dist_l,labinc)
            # Gini and Lorenz curve
            (self.l_ginimodel, self.l_lorenz) = gini(dist_l,labinc)
        
        # Table - data and model Gini and the gap between them
        def giniTable(self):
            
            a_ginidata = 0.857                         # Number from SIM
            l_ginidata = 0.4858                        # Number from SIM

            # Compare model distribution with data
            a_ginigap = 100*(a_ginidata / self.a_ginimodel - 1)

            # Compare model distribution with data
            l_ginigap = 100*(l_ginidata / self.l_ginimodel - 1);

            # MomentsGenerator output
            gini_table = pd.DataFrame({'Gini': ['wealth', 'lab_income'],
                                       'data': [a_ginidata, l_ginidata],
                                       'mode': [self.a_ginimodel, self.l_ginimodel],
                                       'percent_gap': [a_ginigap, l_ginigap]})
            
            return gini_table
        
        # Table - assets and labor income cumulative shares for selected percentiles
        def cumShareTable(self):
            
            cumShare_table = pd.DataFrame({'Percentile': self.a_distdata['percentile'],
                                           'a_data': self.a_distdata['cumulativeShare'],
                                           'a_model': self.a_distmodel['cumulativeShare'],
                                           'l_data': self.l_distdata['cumulativeShare'],
                                           'l_model': self.l_distmodel['cumulativeShare']})
                     
            return cumShare_table
        
        # Graph - Assets distribution: model vs. data
        def plot_a_lorenz(self):
            
            datapercentile = np.hstack((0, self.a_distdata['percentile']))
            datacumShare   = np.hstack((min(self.a_distdata['cumulativeShare'])-1e-8, self.a_distdata['cumulativeShare']))
                         
            plt.plot(self.a_lorenz[:,0] * 100, self.a_lorenz[:,1], linewidth = 2, label = 'model (gini = %0.3f)' % self.a_ginimodel)
            plt.plot(datapercentile, datacumShare, linewidth = 2, label = 'data    (gini = 0.857)')
            plt.plot([0,100],[0,1], linestyle = '--', color = 'black', linewidth = 1.5)         # 45 degree line
            plt.title('Assets distribution', fontsize = 16)
            plt.xlabel('percentiles', fontsize = 13)
            plt.ylabel('cumulative share of total assets', fontsize = 13)
            plt.legend(loc = 'upper right', fontsize = 13)
            
            
        # Graph - Labor income distribution: model vs. data
        def plot_l_lorenz(self):
                         
            datapercentile = np.hstack((0, self.l_distdata['percentile']))
            datacumShare   = np.hstack((0, self.l_distdata['cumulativeShare']))
            
            plt.plot(self.l_lorenz[:,0] * 100, self.l_lorenz[:,1], linewidth = 2, label = 'model (gini = %0.4f)' % self.l_ginimodel)
            plt.plot(datapercentile, datacumShare, linewidth = 2, label = 'data    (gini = 0.4858)')
            plt.plot([0,100],[0,1], linestyle = '--', color = 'black', linewidth = 1.5)         # 45 degree line
            plt.title('Labor income distribution', fontsize = 16)
            plt.xlabel('percentiles', fontsize = 13)
            plt.ylabel('cumulative share of total labor income', fontsize = 13)
            plt.legend(loc = 'upper right', fontsize = 13)
            
        
        # Graph - Assets threshold in dollars: model vs. data
        def plot_a_threshold(self):
                         
            plt.plot(self.a_distdata['percentile'],(self.a_distmodel['threshold']/self.scenario['modelunit_dollar'])/1000, linewidth = 2, label = 'model (gini = %0.3f)' % self.a_ginimodel)
            plt.plot(self.a_distdata['percentile'],self.a_distdata['threshold2016dollars']/1000, linewidth = 2, label = 'data    (gini = 0.857)')
            plt.title('Threshold by wealth percentile', fontsize = 16)
            plt.xlabel('percentiles', fontsize = 13)
            plt.ylabel('thousands of 2016 dollars', fontsize = 13)
            plt.legend(loc = 'upper right', fontsize = 13)

        
        # Graph - Labor income threshold in dollars: model vs. data
        def plot_l_threshold(self):
                         
            plt.plot(self.l_distdata['percentile'],(self.l_distmodel['threshold']/self.scenario['modelunit_dollar'])/1000, linewidth = 2, label = 'model (gini = %0.4f)' % self.l_ginimodel)
            plt.plot(self.l_distdata['percentile'],self.l_distdata['threshold2016dollars']/1000, linewidth = 2, label = 'data    (gini = 0.4858)')
            plt.title('Threshold by labor income percentile', fontsize = 16)
            plt.xlabel('percentiles', fontsize = 13)
            plt.ylabel('thousands of 2016 dollars', fontsize = 13)
            plt.legend(loc = 'upper right', fontsize = 13)
                        
        
        # Table - Age distribution at the bottom of the capital grid
        def topBottomTable(self):
            
            bot_g1 = self.DIST[:,0,:,0:19,:,:,:]
            bot_g2 = self.DIST[:,0,:,19:39,:,:,:]
            bot_g3 = self.DIST[:,0,:,39:59,:,:,:]
            bot_g4 = self.DIST[:,0,:,59:80,:,:,:]
            bottom = self.DIST[:,0,:,:,:,:,:]
            bottom = sum(bottom[:])
            bot_share = np.hstack((sum(bot_g1[:]),sum(bot_g2[:]),sum(bot_g3[:]),sum(bot_g4[:])))
            
            top_g1 = self.DIST[:,-1,:,0:19,:,:,:]
            top_g2 = self.DIST[:,-1,:,19:39,:,:,:]
            top_g3 = self.DIST[:,-1,:,39:59,:,:,:]
            top_g4 = self.DIST[:,-1,:,59:80,:,:,:]
            top    = self.DIST[:,-1,:,:,:,:,:]            
            top    = sum(top[:])
            top_share = np.hstack((sum(top_g1[:]), sum(top_g2[:]), sum(top_g3[:]), sum(top_g4[:])))
                                                
            topBottom_table = pd.DataFrame(['bottom', 'top'], [bottom, top],
                              [bot_share[0], top_share[0]], [bot_share[1], top_share[1]],
                              [bot_share[2], top_share[2]], [bot_share[3], top_share[3]],
                              columns = ['grid', 'total', 'age21to40', 'age41to60', 'age61to80', 'age81to101'])
            
            return topBottom_table
        
        # Graph - Asset holdings by age
        def plot_a_age(self):
            
            kdist_age = np.zeros((1,self.T_life))
            kdist = self.DIST * self.karray
            for age in range(self.T_life):
                pop_age_temp   = self.DIST[:,:,:,age,:,:,:]
                kdist_age_temp = kdist[:,:,:,age,:,:,:]
                kdist_age[age] = (sum(kdist_age_temp[:])/sum(pop_age_temp[:]))/self.scenario['modelunit_dollar']
            
            plt.plot(range(21,101), kdist_age, linewidth = 2)
            plt.title('Average asset holdings by age', fontsize = 16)
            plt.xlabel('age', fontsize = 13)
            plt.ylabel('2016 dollars', fontsize = 13)

        # Graph - Consumption by age
        def plot_c_age(self):
            
            cdist_age = np.zeros((1,self.T_life))
            cdist = self.DIST * self.con
            for age in range(self.T_life):
                pop_age_temp   = self.DIST[:,:,:,age,:,:,:]
                cdist_age_temp = cdist[:,:,:,age,:,:,:]
                cdist_age[age] = (sum(cdist_age_temp[:])/sum(pop_age_temp[:]))/self.scenario['modelunit_dollar']
            
            plt.plot(range(21,101), cdist_age, linewidth = 2)
            plt.title('Average annual consumption by age', fontsize = 16)
            plt.xlabel('age', fontsize = 13)
            plt.ylabel('2016 dollars', fontsize = 13)
            
        
        # Graph - Labor supply by age
        def plot_l_age(self):
            
            ldist_age = np.zeros((1,self.T_life))
            ldist = self.DIST * self.lab
            for age in range(self.T_life):
                pop_age_temp   = self.DIST[:,:,:,age,:,:,:]
                ldist_age_temp = ldist[:,:,:,age,:,:,:]
                ldist_age[age] = (sum(ldist_age_temp[:])) / sum(pop_age_temp[:])
            
            plt.plot(range(21,101),ldist_age, linewidth = 2)
            plt.title('Average annual labor supply by age', fontsize = 16)
            plt.xlabel('age', fontsize = 13)
            
        
        # Graph - Distribution of individuals by asset holdings grid points
        def plot_a_dist(self):
            
            kv     = self.kv
            nk     = np.size(self.DIST,axis = 1)
            kdist  = np.zeros((1,nk))
            for ik in range(nk):
                kdist_temp = self.DIST[:,ik,:,:,:,:,:]
                kdist[ik]  = sum(kdist_temp[:])
            
            plt.plot(range(nk),kdist, linewidth = 2)
            plt.title('Distribution of individuals by asset holdings grid points', fontsize = 16)
            plt.xlabel('grid point in dollars', fontsize = 13)
            ticks = range(1,nk+1)
            labels = [str(math.round(kv[0]/self.scenario['modelunit_dollar'])),
                str(math.round(kv[1]/self.scenario['modelunit_dollar'])), str(math.round(kv[2]/self.scenario['modelunit_dollar'])),
                str(math.round(kv[3]/self.scenario['modelunit_dollar'])), str(math.round(kv[4]/self.scenario['modelunit_dollar'])),
                str(math.round(kv[5]/self.scenario['modelunit_dollar'])), str(math.round(kv[6]/self.scenario['modelunit_dollar'])),
                str(math.round(kv[7]/self.scenario['modelunit_dollar'])), str(math.round(kv[8]/self.scenario['modelunit_dollar'])),
                str(math.round(kv[9]/self.scenario['modelunit_dollar'])), str(math.round(kv[10]/self.scenario['modelunit_dollar'])),
                str(math.round(kv[11]/self.scenario['modelunit_dollar']))]
            plt.xticks(ticks, labels)
            plt.ylabel('share of population', fontsize = 13)
            

        # Social Security benefits moments
        def SS_distribution(self):
            
            # Import variables common to all elements of s           
            dist_retired  = self.DIST[:,:,:,(self.T_work+1):self.T_life,:,:,:]
            ben_retired   = self.ben [:,:,:,(self.T_work+1):self.T_life,:,:,:]
            dist_retired  = dist_retired[:]
            ben_retired   = ben_retired[:]
            
            # Calculate SS outlays as a percentage of GDP
            steady_dir    = PathFinder.getCacheDir(self.scenario)
            s_dynamics    = sio.loadmat(os.join.path(steady_dir, 'dynamics.mat'))
            s['SSbentoout']  = np.sum(ben_retired * dist_retired, axis = 1)/s_dynamics['outs']
            s['SStaxtoout']  = s_dynamics['ssts']/s_dynamics['outs']
            
            # Table with distribution of Social Security benefits among retired households
            dist_retired0 = self.DIST[:,:,0,(self.T_work+1):self.T_life,:,:,:]
            dist_retired0 = dist_retired0[:]/np.sum(dist_retired)
            dist_retired  = dist_retired/np.sum(dist_retired[:])
            ben_distmodel = get_moments(dist_retired,ben_retired)
            ben0          = {'percentile': sum(dist_retired0), 'threshold': 0, 'cumulativeShare': 0}
            s['ben_dist'] = pd.DataFrame(ben0)
            s['ben_dist'].append(ben_distmodel)
            
            # Average asset holdings of retiree earning no SS benefits
            k_retired0    = self.karray[:,:,0,self.T_work:self.T_life,:,:,:]
            k_retired0    = k_retired0[:]
            k_retired0    = k_retired0 * dist_retired0
            s['k_retired0']  = sum(k_retired0)/sum(dist_retired0)/self.scenario['modelunit_dollar']
            
            # Average consumption of retiree earning no SS benefits
            c_retired0    = self.con[:,:,0,self.T_work:self.T_life,:,:,:]
            c_retired0    = c_retired0[:]
            c_retired0    = c_retired0 * dist_retired0
            s.c_retired0  = sum(c_retired0)/sum(dist_retired0)/self.scenario['modelunit_dollar']
            
            return s
        

def get_moments(dist,x):

# Computes selected percentiles of x distributed according to dist
# Inputs:  dist = distribution vector of x
#          x    = vector with variable of interest
# Outputs: percentiles (usually slightly above the actual quintiles)
#          thresholds
#          cumulative share of x with each quintile

    moments = pd.DataFrame({'percentile': [], 'threshold': [], 'cumulativeShare': []})
            
    for perc in [0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99, 1]:
        moments.append(get_percentile(perc,dist,x)) #ok<AGROW>

    return moments

def get_percentile(perc, dist, x):

# Computes percentile perc of x distributed according to dist
# Inputs:  perc = percentile between 0 and 1
#          dist = distribution vector of x
#          x    = vector with variable of interest
# Outputs: percentile (usually slightly above perc)
#          threshold
#          cummulative share of x below percentile perc

    if perc >= 1:
        perc_summary = pd.DataFrame({'percentile': 1,
                              'threshold':       None,
                              'cumulativeShare': 1})
        return perc_summary

    sortv = np.argsort(x)
    x = np.sort(x)
    dist = dist[sortv]
    i = np.where(np.cumsum(dist) >= perc)

    perc_summary = {'percentile': sum(dist[0:i]),
                          'threshold':       x[i],
                          'cumulativeShare': sum(x[0:i] * dist[0:i])/sum(x * dist)}

    return perc_summary

def gini(dist,x,makeplot):

# Inputs:  dist     = vector of population sizes of the different types
#                     (note it doesn't need to be distribution vector of x, but
#                      it cannot be a cumulative distribution)
#          x        = vector with variable of interest
#          makeplot = boolean indicating whether a figure of the Lorentz
#                     curve should be produced or not. TBD Default is false.
# Outputs: ginicoeff = Gini coefficients
#          lorenz    = Lorentz curve: This is a two-column array, with the 
#                      left column representing cumulative population 
#                      shares of the different classes, sorted according to
#                      val, and the right column representing the 
#                      cumulative value share that belongs to the 
#                      population up to the given class. The Lorentz curve 
#                      is a scatter plot of the left vs the right column.


    # Check for negative population/measure
    assert all(dist>=0), 'gini expects nonnegative population vector (%d negative elements).' % sum(dist<0)
    
    # Take care of first point of the Lorenz curve
    if all(x>0):
        # Pre-append a zero because the Lorenz curve contains (0,0) by definition
        dist = np.hstack((0,dist))
        x = np.hstack((0,x))
    else:
        # Use the lowest x holdings to set the first point
        dist = np.hstack((0,dist))
        x = np.hstack((min(x)-1e-8, x))

    # Sort in ascending order wrt x and weight vectors
    z = x * dist
    ord = np.argsort(x)
    dist = dist[ord]
    z = z[ord]
    dist = np.cumsum(dist)
    z    = np.cumsum(z)
    relpop  = dist/dist[-1]
    relz = z/z[-1]

# We compute the area below the Lorentz curve. We do this by computing the 
# average of the left and right Riemann-like sums. (Riemann-'like' because 
# we evaluate not on a uniform grid, but on the points given by the pop data).
# These are the two Rieman-like sums:
    #    leftsum  = sum(relz(1:end-1) .* diff(relpop));
    #    rightsum = sum(relz(2:end)   .* diff(relpop));
# The Gini coefficient is one minus twice the average of leftsum and
# rightsum. We can put all of this into one line.

    ginicoeff = 1 - sum((relz[0:-2]+relz[1:-1]) * np.diff(relpop))

# Lorentz curve
    lorenz = np.vstack((relpop,relz))

# Plot
    if makeplot:
        plt.fill_between(relpop, relz, color = (0.5,0.5,1.0)) # the Lorentz curve
        plt.plot([0,1],[0,1],linestyle = '--', color = 'k')       # 45 degree line
        plt.axis('tight')                                   # ranges of abscissa and ordinate are by definition exactly [0,1]
        plt.axis('equal')                                  # both axes should be equally long
        ticks, _ = plt.yticks()
        plt.xticks(ticks)            # ensure equal ticking
        plt.grid()
        plt.title('\bfGini coefficient = ' + str(ginicoeff))
        plt.xlabel('share of population')
        plt.ylabel('share of value')
    
    return (ginicoeff, lorenz)
