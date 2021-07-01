#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 12:30:57 2021

@author: josephcavataio
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import os

from N2_functions import get_ttfp, make_groups
    
    #%%# Read in N2 Lifespans Data ####
path = os.getcwd()
filename = 'N2 Lifespans FINAL.xlsx'
path = os.path.join(path, filename)
data = pd.read_excel(path)
data = data.drop(columns=data.columns[0])

#%%# Assumption: Temp cultivated = temp maintained if no data is present ####
data['Reported mean lifespan'] = [-1 if x == '' or x == ' 'or x == '  ' else x for x in data['Reported mean lifespan']]
data_new = data.fillna(-1)
L_new = len(data_new)
ind_tempm = data_new.columns.get_loc('Temperature Maintained (Through L4, C)')
ind_tempc = data_new.columns.get_loc('Temperature Cultivated (Adult, C)')

for i in range(L_new):
    if data_new.iloc[(i,ind_tempm)] == -1:
        data_new.iloc[(i,ind_tempm)] = data_new.iloc[(i,ind_tempc)]
    elif data_new.iloc[(i,ind_tempc)] == -1:
        data_new.iloc[(i,ind_tempc)] = data_new.iloc[(i,ind_tempm)]
        
        
#%% Group entries by Temp, Growth Media, and Plate Transfers (data for figure 2)
grouptypes3 = ['Growth Media','Temperature Maintained (Through L4, C)',\
              'Temperature Cultivated (Adult, C)','# Transfers to Fresh Plates']

P=3
data_all3 = make_groups(data_new,grouptypes3,P)

### Group entries by Temp, Growth Media, FUDR, and Plate Transfers (data for figure 2)
grouptypes4 = ['Growth Media','Temperature Maintained (Through L4, C)',\
              'Temperature Cultivated (Adult, C)','FUDR (Yes/No?)','# Transfers to Fresh Plates']

data_all4 = make_groups(data_new,grouptypes4,P) 

days = ['% Alive on Day3','% Alive on Day5','% Alive on Day10','% Alive on Day15',
        '% Alive on Day20','% Alive on Day25','% Alive on Day30',
        '% Alive on Day40','% Alive on Day50']

column_days = ['Day3','Day5','Day10','Day15','Day20','Day25','Day30','Day40','Day50']

#%% Plate transfers
"""
                    FIGURE 3 - Worm Manipulations
                                                                 """
# Indices and their plate transfer type
# 0 = none, 1 = Every other day, 2 = Everyday, 3 = After L4, 6 = Every 2-3 days, 7 = Every 1-2 days

# Figure 3A - 
df_totnm,df_totm,data_nm,data_m,df_totnm_all,df_totm_all = get_ttfp(data_all3,'Every other day','After L4',True,20,20,False,False)
# df_totnm,df_totm = just %Alive on day# data
# data_nm,data_m = statistical data on %Alive data
# df_totnm_all,df_totm_all = % Alive on day# data plus 

fig, ((ax1, ax2)) = plt.subplots(2,1,figsize=(11,14))
fig.suptitle('Effect of Increasing Worm Transfers', fontsize=22, y=1.02)
color_cycle = ["black", "dodgerblue", "chocolate"] 
ind = [3,5,10,15,20,25,30,40,50]
                 
ax1.plot(ind,data_nm['mean'], color='k',label='mean (None)')
ax1.plot(ind,data_m['mean'], color='r',label='mean (Plate Manipulation)')
ax1.fill_between(ind, data_nm['mean']-2*data_nm['std'], data_nm['mean']+2*data_nm['std'],  
                 facecolor=color_cycle[2], label='2x std')
ax1.fill_between(ind, data_nm['mean']-data_nm['std'], data_nm['mean']+data_nm['std'],
                 facecolor=color_cycle[1], label='1x std')

ax2.plot(ind,data_m['mean'], color='k',label='mean')
ax2.fill_between(ind, data_m['mean']-2*data_m['std'], data_m['mean']+2*data_m['std'],  
                 facecolor=color_cycle[2], label='2x std')
ax2.fill_between(ind, data_m['mean']-data_m['std'], data_m['mean']+data_m['std'],
                 facecolor=color_cycle[1], label='1x std')


ax1.set_ylim(-25,130)
ax2.set_ylim(-25,130)

ax1.set_xticks(ind)
ax2.set_xticks(ind)

ax1.set_xticklabels(column_days, fontsize=12, rotation = 70)
ax2.set_xticklabels(column_days, fontsize=12, rotation = 70)

ax1.set_yticks(np.arange(-20,140,20))
ax2.set_yticks(np.arange(-20,140,20))

ax1.set_yticklabels(np.arange(-20,140,20), fontsize=12)
ax2.set_yticklabels(np.arange(-20,140,20), fontsize=12)

ax1.set_title('No Worm Manipulation', fontsize=14,x=0.5)
ax2.set_title('Worm Manipulation', fontsize=14,x=0.5)

ax1.set_ylabel('% Survival',fontsize=20)
ax2.set_ylabel('% Survival',fontsize=20)

ax1.legend(fontsize=14)
ax2.legend(fontsize=14)

plt.tight_layout()
plt.show()

#%% Get mean lifespan info for FUdR+manipulation data
isFUDR1 = True # True if sort by non-FUdR no manipulation vs FUdR no manipulation
isFUDR2 = False # True if sort by non-FUdR any manipulation vs FUdR any manipulation
df_totnm_noFUDR,df_totnm_wFUDR,data_nm_noFUDR,data_nm_wFUDR,df_totnm_all_noFUDR,df_totnm_all_wFUDR \
    = get_ttfp(data_all4,'Every other day','After L4',False,20,20,True,False)

isFUDR1 = False
isFUDR2 = True
df_totm_noFUDR,df_totm_wFUDR,data_m_noFUDR,data_m_wFUDR,df_totm_all_noFUDR,df_totm_all_wFUDR \
    = get_ttfp(data_all4,'Every other day','After L4',False,20,20,False,True)
ind_mls = df_totm_all_noFUDR.columns.get_loc('Reported mean lifespan')
#df_ttfp_list = [df_totED_20,df_totE23D_20]
df_ttfp_list = [df_totnm_noFUDR,df_totnm_wFUDR,df_totm_noFUDR,df_totm_wFUDR]

mls_ttfp2 = [[] for _ in range(len(df_ttfp_list))] # Plotting 2 different ttfp types

for j in range(len(df_ttfp_list)):
    temp = list(df_ttfp_list[j].iloc[:,ind_mls])
    temp = [x for x in temp if x > -1.0]
    #a = np.array(temp, int)
    mls_ttfp2[j] = temp
pd.DataFrame(mls_ttfp2,index=['No Manipulation No FUDR','No Manipulation w FUDR',
                             'Manipulation No FUDR','Manipulation w FUDR']).transpose().to_csv('figures/Fig3F.csv')

print('t-test between No Manipulation No FUDR vs w FUDR gives a p-value of {}'.format(
      round(stats.ttest_ind(mls_ttfp2[0],mls_ttfp2[1])[1],3)))
print('t-test between Manipulation No FUDR vs w FUDR gives a p-value of {}'.format(
      round(stats.ttest_ind(mls_ttfp2[2],mls_ttfp2[3])[1],3)))

#%% Figure 3B & 3C - LIFESPAN PLOTS FOR DIFFERENT PLATE TRANSFERS AT 20˚C
# Doing only 20C because data does not allow for 15C and 25C 

# 0 = none, 1 = Every other day, 2 = Everyday, 3 = After L4, 6 = Every 2-3 days, 7 = Every 1-2 days
df_totED_20,df_totE23D_20,data_ED_20,data_E23D_20,df_totED_20_all,df_totE23D_20_all = get_ttfp(data_all3,'Everyday','Every 2-3 days',False,20,20,False,False)

# FIGURE 3B
fig, (ax1,ax2) = plt.subplots(2,1,figsize=(11,14))
fig.suptitle('Every Day vs Every 2-3 Days (20 C)', fontsize=22, y=1.01)
color_cycle = ["black", "dodgerblue", "chocolate"] 
ind = [3,5,10,15,20,25,30,40,50]

ax1.plot(ind,data_ED_20['mean'], color='k',label='mean (Everyday)')
ax1.plot(ind,data_E23D_20['mean'], color='r',label='mean (Every 2-3 days)')
ax1.fill_between(ind, data_ED_20['mean']-2*data_ED_20['std'], data_ED_20['mean']+2*data_ED_20['std'],  
                 facecolor=color_cycle[2], label='2x std')
ax1.fill_between(ind, data_ED_20['mean']-data_ED_20['std'], data_ED_20['mean']+data_ED_20['std'],
                 facecolor=color_cycle[1], label='1x std')

ax2.plot(ind,data_E23D_20['mean'], color='k',label='mean')
ax2.fill_between(ind, data_E23D_20['mean']-2*data_E23D_20['std'], data_E23D_20['mean']+2*data_E23D_20['std'],
                 facecolor=color_cycle[2], label='2x std')
ax2.fill_between(ind, data_E23D_20['mean']-data_E23D_20['std'], data_E23D_20['mean']+data_E23D_20['std'],
                 facecolor=color_cycle[1], label='1x std')

ax1.set_ylim(-25,130)
ax2.set_ylim(-25,130)

ax1.set_xticks(ind)
ax2.set_xticks(ind)

ax1.set_xticklabels(column_days, fontsize=12, rotation = 70)
ax2.set_xticklabels(column_days, fontsize=12, rotation = 70)

ax1.set_yticks(np.arange(-20,140,20))
ax2.set_yticks(np.arange(-20,140,20))

ax1.set_yticklabels(np.arange(-20,140,20), fontsize=12)
ax2.set_yticklabels(np.arange(-20,140,20), fontsize=12)

ax1.set_title('Everyday (n = {})'.format(int(data_ED_20['count'][0])), fontsize=14,x=0.5)
ax2.set_title('Every 2-3 Days (n = {})'.format(int(data_E23D_20['count'][0])), fontsize=14,x=0.5)

ax1.set_ylabel('% Survival',fontsize=20)
ax2.set_ylabel('% Survival',fontsize=20)

ax1.legend(fontsize=14)
ax2.legend(fontsize=14)

plt.tight_layout()
plt.show()

# data_ED_20.to_csv('figures/Fig3A_everyday.csv')
# data_E23D_20.to_csv('figures/Fig3A_every 2-3 days.csv')

#%% Figure 3E - MEAN LIFESPAN PLOT FOR 2 DIFFERENT TTFP types
###(I need to ask more questions to do Figure 3D)
ind_mls = df_totnm.columns.get_loc('Reported mean lifespan')
#df_ttfp_list = [df_totED_20,df_totE23D_20]
df_ttfp_list = [df_totnm,df_totm]

mls_ttfp = [[] for _ in range(len(df_ttfp_list))] # Plotting 2 different ttfp types

for j in range(len(df_ttfp_list)):
    temp = list(df_ttfp_list[j].iloc[:,ind_mls])
    temp = [x for x in temp if x > -1.0]
    #a = np.array(temp, int)
    mls_ttfp[j] = temp

# Create a figure instance
fig = plt.figure(figsize=(8,6))
ax = fig.add_axes([0,0,1,1])

# Create the boxplot
#bp = ax.violinplot(mls_temp)
bp = ax.boxplot(mls_ttfp)
#ax.set_ylim(0,30)
ax.set_xticks([1,2])
#ax.set_xticklabels(['Everyday \n (n = {})'.format(len(mls_ttfp[0])),
                    #'Every 2-3 Days \n (n = {})'.format(len(mls_ttfp[1]))], fontsize=12)
ax.set_xticklabels(['No Manipulation \n (n = {})'.format(len(mls_ttfp[0])),
                    'Manipulation \n (n = {})'.format(len(mls_ttfp[1]))], fontsize=12)
ax.set_yticks(np.arange(0,35,5))
ax.set_yticklabels(np.arange(0,35,5), fontsize=12)
ax.set_title('Plate Transfers Everyday vs Every 2-3 Days', fontsize=16,x=0.5,y=1.02) # Change titles here
ax.set_ylabel('Average Lifespan (Days)',fontsize=12)

for line in bp['medians']:
    # get position data for median line
    x, y = line.get_xydata()[1] 
    # overlay median value
    plt.text(x, y, '%.1f' % y, verticalalignment='center') # draw above, centered

for line in bp['boxes']:
    x, y = line.get_xydata()[0] 
    plt.text(x,y, '%.1f' % y,
         verticalalignment='top',
         horizontalalignment='right')  
    x, y = line.get_xydata()[3]
    plt.text(x,y, '%.1f' % y,
         verticalalignment='bottom',
         horizontalalignment='right')

######### Saving the Data #############
# pd.DataFrame(mls_ttfp,index=['Everyday','Every 2-3 Days']).transpose().to_csv('/Users/nickurban/Desktop/N2 Lifespan Variability/Figures/mls_fig3E.csv')
# pd.DataFrame(mls_ttfp,index=['No Manipulation','Manipulation']).transpose().to_csv('figures/Fig3E_v2.csv')
