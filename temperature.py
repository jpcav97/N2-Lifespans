#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 27 13:51:39 2021

@author: josephcavataio
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import os
import random

from N2_functions import randomvals_and_diff,residuals,pairplot,make_groups

#%%# Read in N2 Lifespans Data ####
path = os.getcwd()
filename = 'N2 Lifespans FINAL.xlsx'
filename = 'Supplemental File 1.xlsx'
path = os.path.join(path, filename)
data = pd.read_excel(path)
data = data.drop(columns=data.columns[0])

#### Percent of data excluded for different categories ####
""" UNCOMMENT BELOW FOR STATS ON MISSING DATA """
# strainused =  data['N2 Strain Used']
# strainsource = data['N2 Strain Source']
# growthmedia = data['Growth Media']
# rmeanls = data['Reported mean lifespan']
# rmedls = data['Reported median lifespan']
# minanimperplate = data['Min #Animals/Plate']
# maxanimperplate = data['Max #Animals/Plate']
# plateper = data['#Plates/Experiment']
# FUDRbool = data['FUDR (Yes/No?)']
# FUDRconc = data['FUDR Concentration (microM)']
# transfertofresh = data['# Transfers to Fresh Plates']
# tempm = data['Temperature Maintained (Through L4, C)']
# tempc = data['Temperature Cultivated (Adult, C)']
# loccity = data['Lab Location (City)']
# loccountry = data['Lab Location (Country)']

# missingstrainused = strainused.isna().sum()
# print("Percent of entries missing N2 Strain Used")
# print(round(100*missingstrainused/len(data),2), '%', '{}/{}'.format(missingstrainused,numberofentries))

# missingstrainsource = strainsource.isna().sum()
# print("Percent of entries missing N2 Strain Source")
# print(round(100*missingstrainsource/len(data),2), '%', '{}/{}'.format(missingstrainsource,numberofentries))

# missinggrowth = growthmedia.isna().sum()
# print("Percent of entries missing Growth Media")
# print(round(100*missinggrowth/len(data),2), '%', '{}/{}'.format(missinggrowth,numberofentries))

# missingrmeanls = rmeanls.isna().sum()
# print("Percent of entries missing Reported Mean Lifespan")
# print(round(100*missingrmeanls/len(data),2), '%', '{}/{}'.format(missingrmeanls,numberofentries))

# missingrmedls = rmedls.isna().sum()
# print("Percent of entries missing Reported Median Lifespan")
# print(round(100*missingrmedls/len(data),2), '%', '{}/{}'.format(missingrmedls,numberofentries))

# missingminanim = minanimperplate.isna().sum()
# print("Percent of entries missing Min # Animals/Plate")
# print(round(100*missingminanim/len(data),2), '%', '{}/{}'.format(missingminanim,numberofentries))

# missingmaxanim = maxanimperplate.isna().sum()
# print("Percent of entries missing Max # Animals/Plate")
# print(round(100*missingmaxanim/len(data),2), '%', '{}/{}'.format(missingmaxanim,numberofentries))

# missingplateper = plateper.isna().sum()
# print("Percent of entries missing # Plates per Experiment")
# print(round(100*missingplateper/len(data),2), '%', '{}/{}'.format(missingplateper,numberofentries))

# missingFUDRconc = FUDRconc.isna().sum()
# print("Percent of entries missing FUDR Concentration Info")
# print(round(100*missingFUDRconc/len(data),2), '%', '{}/{}'.format(missingFUDRconc,numberofentries))

# missingtransfer = transfertofresh.isna().sum()
# print("Percent of entries missing # Transfers to Fresh Plates Info")
# print(round(100*missingtransfer/len(data),2), '%', '{}/{}'.format(missingtransfer,numberofentries))

# missingtempm = tempm.isna().sum()
# print("Percent of entries missing # Temp Maintained")
# print(round(100*missingtempm/len(data),2), '%', '{}/{}'.format(missingtempm,numberofentries))

# missingtempc = tempc.isna().sum()
# print("Percent of entries missing # Temp Cultivated")
# print(round(100*missingtempc/len(data),2), '%', '{}/{}'.format(missingtempc,numberofentries))

# missingloccity = loccity.isna().sum()
# print("Percent of entries missing # Cities")
# print(round(100*missingloccity/len(data),2), '%', '{}/{}'.format(missingloccity,numberofentries))

# missingloccountry = loccountry.isna().sum()
# print("Percent of entries missing # Countries")
# print(round(100*missingloccountry/len(data),2), '%', '{}/{}'.format(missingloccountry,numberofentries))

#%%########## Array of # Entries, # of Labs, and # of Countries ##############
labs, freqslabs = np.unique(data['Lab'],return_counts = True)
cntry, freqscntry = np.unique(data['Lab Location (Country)'],return_counts = True)
growth, freqsgrowth = np.unique(data['Growth Media'].astype(str),return_counts = True)
PMID, freqsPMID = np.unique(data['PM ID'].astype(str),return_counts = True)
year, freqsyear = np.unique(data['Date Published'].astype(str),return_counts = True)

alldatastats=pd.DataFrame([len(data),len(labs),len(cntry),len(PMID),len(year)],
                         index=['# of Experiments','# of Last Authors','# of Countries',
                                '# of Papers','Unique Years of Publications']).transpose()

print(alldatastats)
#alldatastats.to_csv('/Users/nickurban/Desktop/N2 Lifespan Variability/Figures/table_1.csv')

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

days = ['% Alive on Day3','% Alive on Day5','% Alive on Day10','% Alive on Day15',
        '% Alive on Day20','% Alive on Day25','% Alive on Day30',
        '% Alive on Day40','% Alive on Day50']

column_days = ['Day3','Day5','Day10','Day15','Day20','Day25','Day30','Day40','Day50']
#%%
"""
                    Temperatures - 15C vs 20C vs 25C
                                                                     """
# Figure 1A - LIFESPAN PLOT OF ENTIRE DATASET
lists_s = [[] for _ in range(len(days))]

# Take stats for each day and save into list of lists
for i in range(len(days)):
    s = data_new[days[i]].describe()
    lists_s[i] = [s['count'],round(s['mean'],2),s['25%'],s['50%'],s['75%'],s['std']]

#Transpose list of lists
lists_s = list(map(list, zip(*lists_s))) 

# Turn list of lists into dataframe
df_totdata = pd.DataFrame(lists_s,index=['count','mean','25%','median','75%','std'], 
                          columns=column_days).transpose()

fig = plt.figure(figsize=(12,6))
color_cycle = ["black", "dodgerblue", "chocolate"] 
ind = [3,5,10,15,20,25,30,40,50]

plt.plot(ind,df_totdata['mean'], color=color_cycle[0], label='mean')
plt.fill_between(ind, df_totdata['mean']-2*df_totdata['std'], df_totdata['mean']+2*df_totdata['std'],
                 facecolor=color_cycle[2], label='2x std')
plt.fill_between(ind, df_totdata['mean']-df_totdata['std'], df_totdata['mean']+df_totdata['std'],
                 facecolor=color_cycle[1], label='1x std')

plt.ylabel('% Survival',fontsize=16)
plt.xticks(ind,df_totdata.index, fontsize=14, rotation = 70)
plt.yticks(np.arange(-40,160,20), fontsize=14)
plt.ylim(-40,140)
plt.title('Entire Dataset (n = {})'.format(int(df_totdata['count'][0])),fontsize=20,y=1.02)
plt.legend(bbox_to_anchor=(0.99, 1), loc='upper right', borderaxespad=0.2,\
            fontsize = 10)
plt.tight_layout() 
plt.show()

#df_totdata.to_csv('/Users/nickurban/Desktop/N2 Lifespan Variability/Figures/totaldatasetstats.csv')
df_totdata.to_csv('figures/Fig1A.csv')

print('n = {}'.format(df_totdata['count'][0]))

#%% Joseph's Figure
# Retrieve all the data from the first row and place into new dataframe
datalists = [[] for _ in range(len(data_all.columns)-4)] #Subtract four for the NGM, temps, and count 
for i in range(len(datalists)):
    datalists[i] = data_all.iloc[0,i+4]

df_pairs = pd.DataFrame(datalists,index = data_all.columns[4:]).transpose()
df_pairs = df_pairs[df_pairs['Reported mean lifespan'] > -1]

ind_rml = df_pairs.columns.get_loc('Reported mean lifespan')
df_plot = randomvals_and_diff(df_pairs,ind_rml,False)

idx_std_res,num_out,rsquared,linreg = residuals(df_plot)

# Pair Plot
slope_tot = pairplot(df_plot,rsquared,linreg)

#%%Figure 1B, C, and D - LIFESPAN PLOTS OF GROUPS AT 15, 20, AND 25˚C
ind_day3 = data_all.columns.get_loc('% Alive on Day3')
ind_day50 = data_all.columns.get_loc('% Alive on Day50')
ind_tempc = data_all.columns.get_loc('Temperature Cultivated (Adult, C)')

# Figure out stats for each day
set1_datalists = [[] for _ in range(len(days))]
set2_datalists = [[] for _ in range(len(days))]
set3_datalists = [[] for _ in range(len(days))]

ind_15 = data_all[(data_all['Growth Media'] == 'NGM') & 
                  (data_all['Temperature Maintained (Through L4, C)'] == 15) &
                  (data_all['Temperature Cultivated (Adult, C)'] == 15)].index[0]
ind_20 = data_all[(data_all['Growth Media'] == 'NGM') & 
                  (data_all['Temperature Maintained (Through L4, C)'] == 20) &
                  (data_all['Temperature Cultivated (Adult, C)'] == 20)].index[0]
ind_25 = data_all[(data_all['Growth Media'] == 'NGM') & 
                  (data_all['Temperature Maintained (Through L4, C)'] == 25) &
                  (data_all['Temperature Cultivated (Adult, C)'] == 25)].index[0]
temp_ind = [ind_15,ind_20,ind_25] 
for j in temp_ind:
    # Turn list of lists into dataframe
    if j == temp_ind[0]:
        for i in range(ind_day3,ind_day50+1):
            s = pd.Series(data_all.iloc[j,i]).describe()
            set1_datalists[i-ind_day3] = [s['count'],round(s['mean'],2),s['25%'],s['50%'],s['75%'],s['std']]
        
        #Transpose list of lists
        set1_datalists = list(map(list, zip(*set1_datalists))) 
        df_set1 = pd.DataFrame(set1_datalists,index=['count','mean','25%','median','75%','std'], 
                              columns=column_days).transpose()
    if j == temp_ind[1]:
        for i in range(ind_day3,ind_day50+1):
            s = pd.Series(data_all.iloc[j,i]).describe() 
            set2_datalists[i-ind_day3] = [s['count'],round(s['mean'],2),s['25%'],s['50%'],s['75%'],s['std']]
        
        #Transpose list of lists
        set2_datalists = list(map(list, zip(*set2_datalists))) 
        df_set2 = pd.DataFrame(set2_datalists,index=['count','mean','25%','median','75%','std'], 
                              columns=column_days).transpose()
        
    if j == temp_ind[2]:
        for i in range(ind_day3,ind_day50+1):
            s = pd.Series(data_all.iloc[j,i]).describe() 
            set3_datalists[i-ind_day3] = [s['count'],round(s['mean'],2),s['25%'],s['50%'],s['75%'],s['std']]
        
        #Transpose list of lists
        set3_datalists = list(map(list, zip(*set3_datalists))) 
        df_set3 = pd.DataFrame(set3_datalists,index=['count','mean','25%','median','75%','std'], 
                              columns=column_days).transpose()
 
sets = [df_set1,df_set2,df_set3]
for j in range(len(sets)):
    fig = plt.figure(figsize=(10,6))
    ax = fig.add_axes([0.1,0.15,0.8,0.8]) # main axes
    color_cycle = ["black", "dodgerblue", "chocolate"] 
    ind = [3,5,10,15,20,25,30,40,50]
    
    ax.plot(ind,sets[j]['mean'], color=color_cycle[0],label='Mean') # Change legend here
    
    minusstd = sets[j]['mean']-sets[j]['std']
    plusstd = sets[j]['mean']+sets[j]['std']
    ax.fill_between(ind, minusstd-sets[j]['std'], plusstd+sets[j]['std'], 
                     facecolor=color_cycle[2], label='2x std')
    ax.fill_between(ind, minusstd, plusstd, 
                     facecolor=color_cycle[1], label='1x std')
    
    # ax.plot(ind,df_set1['median'], color='k',label='Median') # Change legend here
    # ax.fill_between(ind, df_set1['median'], df_set1['75%'],alpha=0.3, 
    #                  facecolor=color_cycle[0], label='75%')
    # ax.fill_between(ind, df_set1['25%'], df_set1['median'],alpha=0.3, 
    #                  facecolor=color_cycle[1], label='25%')
    
    ax.set_ylim(-40,140)
    ax.set_xticks(ind)
    ax.set_xticklabels(sets[j].index, fontsize=12, rotation = 70)
    ax.set_yticks(np.arange(-40,160,20))
    ax.set_yticklabels(np.arange(-40,160,20), fontsize=12)
    ax.set_title('{}˚C (n = {})'.format(int(data_all.iloc[temp_ind[j],ind_tempc]),int(sets[j]['count'][0])), fontsize=14,x=0.5) # Change titles here
    ax.set_ylabel('% Survival',fontsize=20)
    ax.legend(fontsize=14)
    plt.show()

    print('n = {}'.format(sets[j]['count'][0]))

    if j == 0:
        sets[j].to_csv('figures/Fig1B.csv')

    if j == 1:
        sets[j].to_csv('figures/Fig1C.csv')

    if j == 6:
        sets[j].to_csv('figures/Fig1D.csv')

#%% Combined lifespan curves of 15, 20, and 25˚C
full_datalists = [[] for _ in range(len(days))]
for i in range(ind_day3,ind_day50+1):
    s = pd.Series(data_all.iloc[ind_15,i])
    s = s.append(pd.Series(data_all.iloc[ind_20,i]))
    s = s.append(pd.Series(data_all.iloc[ind_25,i]))
    s = s.describe() 
    full_datalists[i-ind_day3] = [s['count'],round(s['mean'],2),s['25%'],s['50%'],s['75%'],s['std']]

#Transpose list of lists
full_datalists = list(map(list, zip(*full_datalists))) 
df_set_full = pd.DataFrame(full_datalists,index=['count','mean','25%','median','75%','std'], 
                      columns=column_days).transpose()
df_set_full.to_csv('saved_data/15_20_25_combined_lsp_data.csv')
#%% Figure 1E  - MEAN LIFESPAN PLOTS FOR 15, 20, and 25˚C
## combine these different collections into a list
ind_mls = data_all.columns.get_loc('Reported mean lifespan')
mls_temp = [[] for _ in range(len(temp_ind))]

for j in range(len(temp_ind)):
    temp = data_all.iloc[temp_ind[j],ind_mls]
    temp = [x for x in temp if x > -1.0]
    #a = np.array(temp, int)
    mls_temp[j] = temp

# Create a figure instance
fig = plt.figure()
ax = fig.add_axes([0,0,1,1])

# Create the boxplot
#bp = ax.violinplot(mls_temp)
bp = ax.boxplot(mls_temp)
ax.set_ylim(0,40)
ax.set_xticks([1,2,3])
ax.set_xticklabels(['15˚C \n (n = {})'.format(len(mls_temp[0])),'20˚C \n (n = {})'.format(len(mls_temp[1])),'25˚C \n (n = {})'.format(len(mls_temp[2]))], fontsize=12)
ax.set_yticks(np.arange(0,45,5))
ax.set_yticklabels(np.arange(0,45,5), fontsize=12)

# EDIT TO PUT IN N=?
ax.set_title('Average Reported Lifespan at 15˚C, 20˚C, and 25˚C', fontsize=16,x=0.5,y=1.03) # Change titles here
ax.set_ylabel('Average Lifespan (Days)',fontsize=12)

for line in bp['medians']:
    # get position data for median line
    x, y = line.get_xydata()[1] # top of median line
    # overlay median value
    plt.text(x, y, '%.1f' % y, verticalalignment='center') # draw above, centered

for line in bp['boxes']:
    x, y = line.get_xydata()[0] # bottom of left line
    plt.text(x,y, '%.1f' % y,
         verticalalignment='bottom', # centered
         horizontalalignment='right')      # below
    x, y = line.get_xydata()[3] # bottom of right line
    plt.text(x,y, '%.1f' % y,
         verticalalignment='top', # centered
         horizontalalignment='right')      # below

for line in bp['caps']:
    x, y = line.get_xydata()[0] 
    plt.text(x,y, '%.1f' % y,
         verticalalignment='bottom',
         horizontalalignment='right') 
plt.show()

print('n = {}'.format(len(mls_temp[0])))
print('n = {}'.format(len(mls_temp[1])))
print('n = {}'.format(len(mls_temp[2])))

print('t-test between 20C and 25C gives a p-value of {}'.format(round(
    stats.ttest_ind(mls_temp[1],mls_temp[2])[1],3)))

######### Saving the Data #############
#pd.DataFrame(mls_temp,index=['15C','20C','25C']).transpose().to_csv('/Users/nickurban/Desktop/N2 Lifespan Variability/Figures/mls_fig1E.csv')
pd.DataFrame(mls_temp,index=['15C','20C','25C']).transpose().to_csv('figures/Fig1E.csv')

#%% New fig

""" MLS Data for new fig """ 
ind_mls = data_all.columns.get_loc('Reported mean lifespan')
ind_count = data_all.columns.get_loc('Count')

index2 = [3,7,0,0,0,1,1,1] # multiple 0's and 1's for the three random samples of both groups
mls_temp2 = [[] for _ in range(len(index2))]

m = data_all.iloc[index2[0],ind_mls]
m = [x for x in m if x > -1.0]
max_samples = len(m)
for j in range(len(index2)):
    if index2[j] == 0 or index2[j] == 1:
        vals = data_all.iloc[index2[j],ind_mls]
        vals = [x for x in vals if x > -1.0]
        temp = random.sample(vals,max_samples)
    else:
        temp = data_all.iloc[index2[j],ind_mls]
        temp = [x for x in temp if x > -1.0]
    #a = np.array(temp, int)
    mls_temp2[j] = temp
    
cols = ['20/25 (Index 3)','25/20 (Index 7)','Rand Samples 1 of 20/20 (Index 0)',\
        'Rand Samples 2 of 20/20 (Index 0)','Rand Samples 3 of 20/20 (Index 0)',\
        'Rand Samples 1 of 25/25 (Index 1)','Rand Samples 2 of 25/25 (Index 1)',\
        'Rand Samples 3 of 25/25 (Index 1)']
df_mls_data = pd.DataFrame(mls_temp2,index=cols).transpose()
# df_mls_data.to_csv('figures/New_fig_mls_data.csv')


