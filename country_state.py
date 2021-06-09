#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 13:23:00 2021

@author: josephcavataio
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import os
from PIL import Image
import io
import random
from sklearn.linear_model import LinearRegression

from N2_functions import get_ttfp,randomvals_and_diff,residuals,pairplot,make_groups,\
    plotlifespanspread,categorize,plotmeanlifespan, get_FUDR
    
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
        
#%% Group entries by Temp and Growth Media (data for figure 1)
P=3
grouptypes1 = ['Growth Media','Temperature Maintained (Through L4, C)',\
              'Temperature Cultivated (Adult, C)']

data_all = make_groups(data_new,grouptypes1,P)

### Group entries by Temp, Growth Media, and FUDR (data for figure 2)
grouptypes2 = ['Growth Media','Temperature Maintained (Through L4, C)',\
              'Temperature Cultivated (Adult, C)','FUDR (Yes/No?)']

data_all2 = make_groups(data_new,grouptypes2,P) 

### Group entries by Temp, Growth Media, and Plate Transfers (data for figure 2)
grouptypes3 = ['Growth Media','Temperature Maintained (Through L4, C)',\
              'Temperature Cultivated (Adult, C)','# Transfers to Fresh Plates']

data_all3 = make_groups(data_new,grouptypes3,P)

### Group entries by Temp, Growth Media, FUDR, and Plate Transfers (data for figure 2)
grouptypes4 = ['Growth Media','Temperature Maintained (Through L4, C)',\
              'Temperature Cultivated (Adult, C)','FUDR (Yes/No?)','# Transfers to Fresh Plates']

data_all4 = make_groups(data_new,grouptypes4,P) 

days = ['% Alive on Day3','% Alive on Day5','% Alive on Day10','% Alive on Day15',
        '% Alive on Day20','% Alive on Day25','% Alive on Day30',
        '% Alive on Day40','% Alive on Day50']

column_days = ['Day3','Day5','Day10','Day15','Day20','Day25','Day30','Day40','Day50']
#%% Figure 4B-D
"""
                    FIGURE 4 - % ALIVE BY COUNTRY
                                                                """
### LIFESPAN SPREAD GRAPHS FOR FUDR = 40 day 3-50                             
N = 0 # This determines the set that is plotted from whichever data_all is placed in the fxn 4 lines below
daystoplot = ['% Alive on Day15','% Alive on Day20','% Alive on Day25']

for k in range(len(daystoplot)):
    set_data,unique_group_countries,avg_set=categorize(data_all,'Lab Location (Country)', N, daystoplot[k])
    
    Lcountries = 10 # Number of countries to plot    
    perc_alive_lists = [[] for _ in range(6)]
    
    ### Create bar graph ###         
    for i in range(Lcountries):
        count20 = 0
        count40 = 0
        count60 = 0
        count80 = 0
        count100 = 0
        
        for j in range(len(set_data[i])):
            
            if int(set_data[i].iloc[j]) <= 20:
                count20 = count20 + 1
                
            elif ((int(set_data[i].iloc[j]) <= 40) and (int(set_data[i].iloc[j]) > 20)):
                count40 = count40 + 1
                
            elif ((int(set_data[i].iloc[j]) <= 60) and (int(set_data[i].iloc[j]) > 40)):
                count60 = count60 + 1
                
            elif ((int(set_data[i].iloc[j]) <= 80) and (int(set_data[i].iloc[j]) > 60)):
                count80 = count80 + 1
                
            elif ((int(set_data[i].iloc[j]) <= 100) and (int(set_data[i].iloc[j]) > 80)):
                count100 = count100 + 1
                 
        percent20 = count20/len(set_data[i])
        percent40 = count40/len(set_data[i])
        percent60 = count60/len(set_data[i])
        percent80 = count80/len(set_data[i])
        percent100 = count100/len(set_data[i])
        
        perc_alive_lists[0].append(percent20*100)
        perc_alive_lists[1].append(percent40*100)
        perc_alive_lists[2].append(percent60*100)
        perc_alive_lists[3].append(percent80*100)
        perc_alive_lists[4].append(percent100*100)
        perc_alive_lists[5].append(len(set_data[i]))
    
    df_alive = pd.DataFrame(perc_alive_lists,index=['20%','40%','60%','80%','100%','count']).transpose()
    
    ind = np.arange(Lcountries)+1  # the x locations for the groups
    width = 0.5
    
    fig = plt.figure(figsize=(14,10)) ##### CHANGE FIG1 to FIG ##################
    p1 = plt.bar(ind, df_alive['20%'], width, color='k')
    p2 = plt.bar(ind, df_alive['40%'], width, bottom=df_alive['20%'],color='chocolate')
    p3 = plt.bar(ind, df_alive['60%'], width, color='dodgerblue',\
                 bottom=[i+j for i,j in zip(df_alive['20%'],df_alive['40%'])])
    p4 = plt.bar(ind, df_alive['80%'], width, color='crimson',\
                 bottom=[i+j+k for i,j,k in zip(df_alive['20%'],df_alive['40%'],df_alive['60%'])])
    p5 = plt.bar(ind, df_alive['100%'], width, color='seagreen',\
                 bottom=[i+j+k+z for i,j,k,z in zip(df_alive['20%'],df_alive['40%'],df_alive['60%'],df_alive['80%'])])
        
    plt.ylabel('Percent of Total Entries In Set',fontsize=16)
    
    # Find out the 'n=?' for each country in the variable 'unique_group_countries'
    plt.xticks(ind, fontsize=14, labels=unique_group_countries.iloc[0:Lcountries,0],rotation=70)
    plt.yticks(np.arange(0, 110, 10),size=14)
    plt.title('{}'.format(daystoplot[k]),\
              fontsize=24,y=1.02)
    plt.legend((p1[0], p2[0], p3[0], p4[0], p5[0]), ('Entries with <20% Alive',\
               'Entries with 20-40% Alive','Entries with 40-60% Alive',\
               'Entries with 60-80% Alive','Entries with 80-100% Alive'),\
               bbox_to_anchor=(1.0, 1),loc='upper left',\
               borderaxespad=0.5, fontsize=14)
    plt.tight_layout()
    plt.show()
    
    test = stats.ttest_ind(set_data[0], set_data[1], equal_var=False)
    first = unique_group_countries.iloc[0,0]
    second = unique_group_countries.iloc[1,0]

    print('t-test between {} and {} gives a p-value of {}'.format(first, second,round(test[1],3)))
    
    if k == 0:
        df_alive.to_csv('figures/Fig4B.csv')

    if k == 1:
        df_alive.to_csv('figures/Fig4C.csv')

    if k == 2:
        df_alive.to_csv('figures/Fig4D.csv')

#%% Figure 4E - MEAN LIFESPAN PLOTS FOR TOP 5 COUNTRIES
for j in range(len(avg_set)):
    temp = list(avg_set[j])
    temp = [x for x in temp if x > -1.0]
    #a = np.array(temp, int)
    avg_set[j] = temp

# Create a figure instance
fig = plt.figure(figsize=(10,6))
ax = fig.add_axes([0,0,1,1])

# Create the boxplot
#bp = ax.violinplot(mls_temp)
bp = ax.boxplot(avg_set[0:5]) # Change how many countries are plotted here
ax.set_ylim(0,40)
ax.set_xticks([1,2,3,4,5]) # This number should match number 2 lines above and always starts from 1
ax.set_xticklabels(['{} \n (n = {})'.format(unique_group_countries.iloc[0,0],len(avg_set[0])),
                    '{} \n (n = {})'.format(unique_group_countries.iloc[1,0],len(avg_set[1])),
                    '{} \n (n = {})'.format(unique_group_countries.iloc[2,0],len(avg_set[2])),
                    '{} \n (n = {})'.format(unique_group_countries.iloc[3,0],len(avg_set[3])),
                    '{} \n (n = {})'.format(unique_group_countries.iloc[4,0],len(avg_set[4]))],
                     fontsize=12)
ax.set_yticks(np.arange(0,45,5))
ax.set_yticklabels(np.arange(0,45,5), fontsize=12)
ax.set_title('Different Countries', fontsize=16,x=0.5) # Change titles here
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
#pd.DataFrame(avg_set,index=unique_group_countries['Names']).transpose().to_csv('/Users/nickurban/Desktop/N2 Lifespan Variability/Figures/mls_fig4E.csv')
pd.DataFrame(avg_set,index=unique_group_countries['Names']).transpose().to_csv('figures/Fig4E.csv')

#%% Figure 5B-D
"""
                FIGURE 5 - % ALIVE BY STATES/PROVINCE
                                                                    """
### LIFESPAN SPREAD GRAPHS FOR FUDR = 40 day 3 - 50 
N = 0
daystoplot = ['% Alive on Day15','% Alive on Day20','% Alive on Day25']

for k in range(len(daystoplot)):
    set_data,unique_group_states,avg_set_states=categorize(data_all,'Lab Location (State/Province)', N, daystoplot[k])
    x = unique_group_states.index[unique_group_states['Names'] == '-1'].tolist()
    
    # Remove values for where State/Province == -1
    unique_group_states = unique_group_states.drop(x[0])
    set_data.pop(x[0])
    avg_set_states.pop(x[0])
    
    Lstates = 10 # Number of countries to plot    
    perc_alive_lists = [[] for _ in range(6)]
    
    ### Create bar graph ###         
    for i in range(Lstates):
        count20 = 0
        count40 = 0
        count60 = 0
        count80 = 0
        count100 = 0
        
        for j in range(len(set_data[i])):
            
            if int(set_data[i].iloc[j]) <= 20:
                count20 = count20 + 1
                
            elif ((int(set_data[i].iloc[j]) <= 40) and (int(set_data[i].iloc[j]) > 20)):
                count40 = count40 + 1
                
            elif ((int(set_data[i].iloc[j]) <= 60) and (int(set_data[i].iloc[j]) > 40)):
                count60 = count60 + 1
                
            elif ((int(set_data[i].iloc[j]) <= 80) and (int(set_data[i].iloc[j]) > 60)):
                count80 = count80 + 1
                
            elif ((int(set_data[i].iloc[j]) <= 100) and (int(set_data[i].iloc[j]) > 80)):
                count100 = count100 + 1
                 
        percent20 = count20/len(set_data[i])
        percent40 = count40/len(set_data[i])
        percent60 = count60/len(set_data[i])
        percent80 = count80/len(set_data[i])
        percent100 = count100/len(set_data[i])
        
        perc_alive_lists[0].append(percent20*100)
        perc_alive_lists[1].append(percent40*100)
        perc_alive_lists[2].append(percent60*100)
        perc_alive_lists[3].append(percent80*100)
        perc_alive_lists[4].append(percent100*100)
        perc_alive_lists[5].append(len(set_data[i]))
    
    df_alive = pd.DataFrame(perc_alive_lists,index=['20%','40%','60%','80%','100%','count']).transpose()
    
    ind = np.arange(Lstates)+1  # the x locations for the groups
    width = 0.5
    
    fig1 = plt.figure(figsize=(14,10))
    p1 = plt.bar(ind, df_alive['20%'], width, color='k')
    p2 = plt.bar(ind, df_alive['40%'], width, bottom=df_alive['20%'],color='chocolate')
    p3 = plt.bar(ind, df_alive['60%'], width, color='dodgerblue',\
                 bottom=[i+j for i,j in zip(df_alive['20%'],df_alive['40%'])])
    p4 = plt.bar(ind, df_alive['80%'], width, color='crimson',\
                 bottom=[i+j+k for i,j,k in zip(df_alive['20%'],df_alive['40%'],df_alive['60%'])])
    p5 = plt.bar(ind, df_alive['100%'], width, color='seagreen',\
                 bottom=[i+j+k+z for i,j,k,z in zip(df_alive['20%'],df_alive['40%'],df_alive['60%'],df_alive['80%'])])
        
    plt.ylabel('Percent of Total Entries In Set',fontsize=16)
    
    # Find out the 'n=?' for each country in the variable 'unique_group_states'
    plt.xticks(ind, fontsize=14, labels=unique_group_states.iloc[0:Lstates,0])
    plt.yticks(np.arange(0, 110, 10),size=14)
    plt.title('{}'.format(daystoplot[k]),\
              fontsize=24,y=1.02)
    plt.legend((p1[0], p2[0], p3[0], p4[0], p5[0]), ('Entries with <20% Alive',\
               'Entries with 20-40% Alive','Entries with 40-60% Alive',\
               'Entries with 60-80% Alive','Entries with 80-100% Alive'),\
               bbox_to_anchor=(1.0, 1),loc='upper left',\
               borderaxespad=0.5, fontsize=14)
    plt.show()

    test = stats.ttest_ind(set_data[0], set_data[1], equal_var=False)
    first = unique_group_states.iloc[0,0]
    second = unique_group_states.iloc[1,0]

    print('t-test between {} and {} gives a p-value of {}'.format(first, second,round(test[1],3)))
    
    if k == 0:
        df_alive.to_csv('figures/Fig5B.csv')

    elif k == 1:
        df_alive.to_csv('figures/Fig5C.csv')

    elif k == 2:
        df_alive.to_csv('figures/Fig5D.csv')
    
#%% Figure 5E - MENA LIFESPAN PLOTS FOR TOP 5 STATES/PROVINCES
for j in range(len(avg_set_states)):
    temp = list(avg_set_states[j])
    temp = [x for x in temp if x > -1.0]
    #a = np.array(temp, int)
    avg_set_states[j] = temp

# Create a figure instance
fig = plt.figure(figsize=(10,6))
ax = fig.add_axes([0,0,1,1])

# Create the boxplot
#bp = ax.violinplot(mls_temp)
bp = ax.boxplot(avg_set_states[0:5])
ax.set_ylim(0,40)
ax.set_xticks([1,2,3,4,5])
ax.set_xticklabels(['{} \n (n = {})'.format(unique_group_states.iloc[0,0],len(avg_set_states[0])),
                    '{} \n (n = {})'.format(unique_group_states.iloc[1,0],len(avg_set_states[1])),
                    '{} \n (n = {})'.format(unique_group_states.iloc[2,0],len(avg_set_states[2])),
                    '{} \n (n = {})'.format(unique_group_states.iloc[3,0],len(avg_set_states[3])),
                    '{} \n (n = {})'.format(unique_group_states.iloc[4,0],len(avg_set_states[4]))],
                     fontsize=12)

ax.set_yticks(np.arange(0,45,5))
ax.set_yticklabels(np.arange(0,45,5), fontsize=12)
ax.set_title('Different States/Provinces (Without FUDR)', fontsize=16,x=0.5) # Change title here
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
# pd.DataFrame(avg_set_states,index=unique_group_states['Names']).transpose().to_csv('/Users/nickurban/Desktop/N2 Lifespan Variability/Figures/mls_fig5E.csv')
# pd.DataFrame(avg_set_states,index=unique_group_states['Names']).transpose().to_csv('figures/Fig5E.csv')

#%% Multivariate Linear Regression
data_mv = data_new[data_new['Reported mean lifespan'] != -1]

data_mv[data_mv['FUDR Concentration (microM)'] == -1] = 0
FUDR_ind = data_mv.columns.get_loc('FUDR Concentration (microM)')

x = data_mv[['Temperature Cultivated (Adult, C)','FUDR Concentration (microM)']]
y = data_mv[['Reported mean lifespan']]

regr = LinearRegression()
regr.fit(x,y)

print(regr.coef_)
