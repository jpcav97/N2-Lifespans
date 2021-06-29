#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 11:19:02 2021

@author: josephcavataio
"""
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from N2_functions import mean_confidence_interval
    
all_mls_data = pd.read_csv('all_mls_data.csv',index_col=0)
all_rml = pd.read_csv('saved_data/total_dataset_rmls.csv')

#%% Calculate Mean, CI, and Quartiles
all_rml = list(all_rml['All Reported mean lifespan data'])
m = np.mean(all_rml)
s = np.std(all_rml)
L = len(all_rml)
sem = s/math.sqrt(L) # Figure out the z value CI = mean ± z * mean/√n
h = sem * 1.96
print('\n############ Entire Dataset (n = {}) ############'.format(L))
print('mean = {}, lower CI = {}, upper CI = {}'.format(round(m,2),round(m-h,2),round(m+h,2)))

stat = pd.Series(all_rml).describe()
print('median = {}, 25% = {}, 75% = {}'.format(round(stat['50%'],2),round(stat['25%'],2),round(stat['75%'],2)))

labels_temp = ['15˚C','20˚C','25˚C']
rml_temp = [list(all_mls_data[x].dropna()) for x in labels_temp]
stat_temp = mean_confidence_interval(rml_temp,labels_temp,(7,5))

labels_fudr20 = ['(-)FUDR 20˚C','FUDR 20˚C']
rml_fudr20 = [list(all_mls_data[x].dropna()) for x in labels_fudr20]
stat_FUDR = mean_confidence_interval(rml_fudr20,labels_fudr20,(7,5))

labels_fudr25 = ['(-)FUDR 25˚C','FUDR 25˚C']
rml_fudr25 = [list(all_mls_data[x].dropna()) for x in labels_fudr25]
stat_FUDR25 = mean_confidence_interval(rml_fudr25,labels_fudr25,(7,5))

labels_fudr_conc = ['[FUDR] = 25µM','[FUDR] = 40µM','[FUDR] = 100µM',
                    '[FUDR] = 200µM','[FUDR] = 400µM']
rml_fudrc = [list(all_mls_data[x].dropna()) for x in labels_fudr_conc]
stat_FUDRc = mean_confidence_interval(rml_fudrc,labels_fudr_conc,(11,4))

labels_ttfp = ['No Manipulation 20˚C','Manipulation 20˚C']
rml_ttfp = [list(all_mls_data[x].dropna()) for x in labels_ttfp]
stat_ttfp = mean_confidence_interval(rml_ttfp,labels_ttfp,(7,5))

labels_cntry = ['USA','China','Germany','Republic of Korea','Japan']
rml_cntry = [list(all_mls_data[x].dropna()) for x in labels_cntry]
stat_cntry = mean_confidence_interval(rml_cntry,labels_cntry,(10,4))

labels_states = ['CA','MA','NY','WA','TX','MI']
rml_state = [list(all_mls_data[x].dropna()) for x in labels_states]
stat_states = mean_confidence_interval(rml_state,labels_states,(10,4))

#%% ANOM
# Group data by analyses
temp_crit_value = 2.39 # Based on table from paper given by Santiago
fudr20_crit_value = 2.24
fudr25_crit_value = 2.35
fudrc_crit_value = 2.65
ttfp_crit_value = 2.24
cntry_crit_value = 2.57
state_crit_value = 2.67

# Save above data as list of lists
rmls = [rml_temp,rml_fudr20,rml_fudr25,rml_fudrc,rml_ttfp,rml_cntry,rml_state]

# Save critical values into one list
crit_vals = [temp_crit_value,fudr20_crit_value,fudr25_crit_value,fudrc_crit_value,
             ttfp_crit_value,cntry_crit_value,state_crit_value]

# Save labels into one list
labels = [labels_temp,labels_fudr20,labels_fudr25,labels_fudr_conc,labels_ttfp,
          labels_cntry,labels_states]

all_CIs = [[] for _ in range(len(rmls))]
all_avgs = [[] for _ in range(len(rmls))]
all_gms = [0]*len(rmls)

for i in range(len(rmls)):
    # Calculate Grand Mean
    Ns = [len(x) for x in rmls[i]]
    N = sum(Ns)
    avg_tot = sum([np.mean(a)*len(a) for a in rmls[i]])/N
    
    # Calculate Mean Square Error
    numerator = sum([(len(x)-1)*np.var(x) for x in rmls[i]])
    denominator = N-1
    MSe = numerator/denominator
    
    # Calculate decision lines
    all_CIs[i] = [crit_vals[i]*math.sqrt(MSe)*math.sqrt((N-j)/(N*j)) for j in Ns]
    all_avgs[i] = [np.mean(m) for m in rmls[i]]
    all_gms[i] = avg_tot
    
    # Plot line representing Grand Mean
    if i == 3: # Different [FUDR] (NEEDS MORE SPACE FOR LABELS)
        figsize = (9,4)
    else:
        figsize = (6,4)
    
    m = 1
    w = 0.5
    plt.figure(figsize=figsize)
    y = np.arange(0,m*len(rmls[i])+3)
    plt.plot(y,[all_gms[i]]*len(y),'--k')
    
    xcoord=0
    for j in range(len(rmls[i])):
        # Plot line from Grand Mean to mean of group i
        if np.mean(rmls[i][j]) < all_gms[i]:
            y1 = np.arange(all_avgs[i][j],all_gms[i],0.01)
        else:
            y1 = np.arange(all_gms[i],all_avgs[i][j],0.01)
        if i == 3: # Different [FUDR]
            xcoord = xcoord + m
        else:
            xcoord = j+m
        x1 = [xcoord]*len(y1)
        plt.plot(x1,y1,'k')
        
        # Plot decision lines
        upper_dl_x = [xcoord-w,xcoord,xcoord+w]
        upper_dl_y1 = [all_gms[i]]*3
        upper_dl_y2 = [all_gms[i]+all_CIs[i][j]]*3
        lower_dl_x = [xcoord-w,xcoord,xcoord+w]
        lower_dl_y1 = [all_gms[i]-all_CIs[i][j]]*3
        lower_dl_y2 = [all_gms[i]]*3
        
        plt.fill_between(upper_dl_x,upper_dl_y1,upper_dl_y2,color='lightgrey',
                         linestyle='--',edgecolor='k')
        plt.fill_between(lower_dl_x,lower_dl_y1,lower_dl_y2,color='lightgrey',
                         linestyle='--',edgecolor='k')
        
        # Plot Dot
        plt.scatter(xcoord,np.mean(rmls[i][j]),c='k',s=50)
    
    new_labels = [labels[i][n]+'\n(n='+str(Ns[n])+')' for n in range(len(Ns))]
    if i == 5:
        plt.xticks(np.arange(0,m*len(rmls[i]),m)+m,labels=new_labels,rotation=25)
    else:
        plt.xticks(np.arange(0,m*len(rmls[i]),m)+m,labels=new_labels)
    
    plt.xlim(0,m*len(rmls[i])+m)
    
    # Find max difference between all decision lines and Grand Mean
    CI_diff = max(all_CIs[i])
    # Find max difference between all means and Grand Mean
    mean_diff = max([abs(all_gms[i]-np.mean(x)) for x in rmls[i]])
    # Set ylim to whichever is great from above differences
    ylim = max(CI_diff,mean_diff)
    if i == 0:
        plt.ylim(13,all_gms[i]+ylim*1.15)
    elif ylim < 1:
        plt.ylim(all_gms[i]-ylim*2,all_gms[i]+ylim*2)
    else:
        plt.ylim(all_gms[i]-ylim*1.15,all_gms[i]+ylim*1.15)
        
    plt.ylabel('Reported Mean Lifespan')
    plt.tight_layout()
    plt.savefig('figures/{}.eps'.format(labels[i]), format='eps')
    plt.show()
   
temp_ANOM_data = pd.DataFrame({'Averages': all_avgs[0],'Grand Means': [all_gms[0]]*len(all_avgs[0]),
                          'Decision Lines': all_CIs[0]}).to_csv('saved_data/temp_ANOM_data.csv')

fudr20_ANOM_data = pd.DataFrame({'Averages': all_avgs[1],'Grand Means': [all_gms[0]]*len(all_avgs[1]),
                          'Decision Lines': all_CIs[1]}).to_csv('saved_data/fudr20_ANOM_data.csv')

fudr25_ANOM_data = pd.DataFrame({'Averages': all_avgs[2],'Grand Means': [all_gms[0]]*len(all_avgs[2]),
                          'Decision Lines': all_CIs[2]}).to_csv('saved_data/fudr25_ANOM_data.csv')

fudrc_ANOM_data = pd.DataFrame({'Averages': all_avgs[3],'Grand Means': [all_gms[0]]*len(all_avgs[3]),
                          'Decision Lines': all_CIs[3]}).to_csv('saved_data/fudrc_ANOM_data.csv')

ttfp_ANOM_data = pd.DataFrame({'Averages': all_avgs[4],'Grand Means': [all_gms[0]]*len(all_avgs[4]),
                          'Decision Lines': all_CIs[4]}).to_csv('saved_data/ttfp_ANOM_data.csv')

cntry_ANOM_data = pd.DataFrame({'Averages': all_avgs[5],'Grand Means': [all_gms[0]]*len(all_avgs[5]),
                          'Decision Lines': all_CIs[5]}).to_csv('saved_data/cntry_ANOM_data.csv')

state_temp_ANOM_data = pd.DataFrame({'Averages': all_avgs[6],'Grand Means': [all_gms[0]]*len(all_avgs[6]),
                          'Decision Lines': all_CIs[6]}).to_csv('saved_data/state_ANOM_data.csv')
