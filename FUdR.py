#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 11:38:19 2021

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
    
#%% Group entries by Temp, Growth Media, and FUDR (data for figure 2)
grouptypes2 = ['Growth Media','Temperature Maintained (Through L4, C)',\
              'Temperature Cultivated (Adult, C)','FUDR (Yes/No?)']

P=3
data_all2 = make_groups(data_new,grouptypes2,P) 

days = ['% Alive on Day3','% Alive on Day5','% Alive on Day10','% Alive on Day15',
        '% Alive on Day20','% Alive on Day25','% Alive on Day30',
        '% Alive on Day40','% Alive on Day50']

column_days = ['Day3','Day5','Day10','Day15','Day20','Day25','Day30','Day40','Day50']
#%% FUDR vs No FUDR
"""
                        FIGURE 2 - FUDR vs Non-FUDR
                                                             """
temps = [15,20,25]
for i in range(len(temps)):                                                    
    x = data_all2.loc[(data_all2['FUDR (Yes/No?)'] == 'No') &
                      (data_all2['Growth Media'] == 'NGM') &
                      (data_all2['Temperature Maintained (Through L4, C)'] == temps[i]) &
                      (data_all2['Temperature Cultivated (Adult, C)'] == temps[i])]
    y = data_all2.loc[(data_all2['FUDR (Yes/No?)'] == 'Yes') & 
                      (data_all2['Growth Media'] == 'NGM') &
                      (data_all2['Temperature Maintained (Through L4, C)'] == temps[i]) &
                      (data_all2['Temperature Cultivated (Adult, C)'] == temps[i])]

    # LIFESPAN PLOT to compare FUDR to w/o FUDR
    ind_day3 = x.columns.get_loc('% Alive on Day3')
    ind_day50 = x.columns.get_loc('% Alive on Day50')
    ind_count = x.columns.get_loc('Count')
    wFUDR_list = [[] for _ in range(len(days))]
    woFUDR_list = [[] for _ in range(len(days))]

    for j in range(ind_day3,ind_day50+1):
        s = pd.Series(x.iloc[0,j]).describe()
        woFUDR_list[j-ind_day3] = [s['count'],round(s['mean'],2),s['25%'],s['50%'],s['75%'],s['std']]

    #Transpose list of lists
    woFUDR_list = list(map(list, zip(*woFUDR_list))) 
    df_woFUDR = pd.DataFrame(woFUDR_list,index=['count','mean','25%','median','75%','std'], 
                          columns=column_days).transpose()

    for j in range(ind_day3,ind_day50+1):
        s = pd.Series(y.iloc[0,j]).describe()
        wFUDR_list[j-ind_day3] = [s['count'],round(s['mean'],2),s['25%'],s['50%'],s['75%'],s['std']]
    
    #Transpose list of lists
    wFUDR_list = list(map(list, zip(*wFUDR_list))) 
    df_wFUDR = pd.DataFrame(wFUDR_list,index=['count','mean','25%','median','75%','std'], 
                          columns=column_days).transpose()

    fig = plt.figure(figsize=(10,6))
    ax = fig.add_axes([0.1,0.15,0.8,0.8]) # main axes
    color_cycle = ["black", "dodgerblue", "chocolate"] 
    ind = [3,5,10,15,20,25,30,40,50]
    ax.plot(ind,df_wFUDR['mean'], color=color_cycle[0],label='Mean (with FUDR)') # Change legend here
    ax.plot(ind,df_woFUDR['mean'], color='r',label='Mean (w/o FUDR)') # Change legend here
    ax.fill_between(ind, df_wFUDR['mean']-2*df_wFUDR['std'], df_wFUDR['mean']+2*df_wFUDR['std'], 
                     facecolor=color_cycle[2], label='2x std')
    ax.fill_between(ind, df_wFUDR['mean']-df_wFUDR['std'], df_wFUDR['mean']+df_wFUDR['std'], 
                     facecolor=color_cycle[1],label='1x std')
    
    ax.set_ylim(-60,140)
    ax.set_xticks(ind)
    ax.set_xticklabels(df_wFUDR.index, fontsize=12, rotation = 70)
    ax.set_yticks(np.arange(-60,160,20))
    ax.set_yticklabels(np.arange(-60,160,20), fontsize=12)
    ax.set_title('With FUDR (n = {}) and Without FUDR (n = {})'.format(int(df_wFUDR['count'][0]),int(df_woFUDR['count'][0])), fontsize=14,x=0.5) # Change titles here
    ax.set_ylabel('% Survival',fontsize=20)
    ax.legend(fontsize=14)
    plt.show()

    print('2A (-)FUDR: n = {}'.format(df_wFUDR['count'][0]))
    print('2A FUDR: n = {}'.format(df_woFUDR['count'][0]))

    df_wFUDR.to_csv('saved_data/Temp{}_wFUdR_lsp.csv'.format(temps[i]))
    df_woFUDR.to_csv('saved_data/Temp{}_woFUdR_lsp.csv'.format(temps[i]))
    
    if i == 0:
        x15 = x
        y15 = y
    elif i == 1:
        x20 = x
        y20 = y
    elif i == 2:
        x25 = x
        y25 = y

#%% Figure 2B - MEAN LIFESPAN PLOT of FUDR vs Non-FUDR
ind_mls = x15.columns.get_loc('Reported mean lifespan')
xys = [x15,y15,x20,y20,x25,y25]
mls_FUDR = [[] for _ in range(len(xys))]

# Fill in mls blanks with -2.0
for j in range(len(xys)):
    t = xys[j].iloc[0,ind_mls]
    mls_FUDR[j] = [x for x in t if x > -1.0]

# Create a figure instance
fig = plt.figure()
ax = fig.add_axes([0,0,1,1])

# Create the boxplot
#bp = ax.violinplot(mls_temp)
bp = ax.boxplot(mls_FUDR)
ax.set_ylim(0,40)
ax.set_xticks([1,2,3,4,5,6])
ax.set_xticklabels(['(-)FUDR \n temp = 15˚C \n (n = {})'.format(len(mls_FUDR[0])),\
                    'FUDR \n temp = 15˚C \n (n = {})'.format(len(mls_FUDR[1])),\
                   '(-)FUDR \n temp = 20˚C \n (n = {})'.format(len(mls_FUDR[2])),\
                    'FUDR \n temp = 20˚C \n (n = {})'.format(len(mls_FUDR[3])),\
                   '(-)FUDR \n temp = 25˚C \n (n = {})'.format(len(mls_FUDR[4])),\
                    'FUDR \n temp = 25˚C \n (n = {})'.format(len(mls_FUDR[5]))],\
                   fontsize=10)
ax.set_yticks(np.arange(0,45,5))
ax.set_yticklabels(np.arange(0,45,5), fontsize=12)
ax.set_title('Average Reported Lifespan with and without FUDR at 15, 20, and 25˚C',
             fontsize=16,x=0.5,y=1.03)
ax.set_ylabel('Average Lifespan (Days)',fontsize=12)

for line in bp['medians']:
    # get position data for median line
    x, y = line.get_xydata()[1] 
    # overlay median value
    plt.text(x, y, '%.1f' % y, verticalalignment='center',fontsize=8) # draw above, centered

for line in bp['boxes']:
    x, y = line.get_xydata()[0] 
    plt.text(x,y, '%.1f' % y,
         verticalalignment='bottom',
         horizontalalignment='right',fontsize=8)  
    x, y = line.get_xydata()[3]
    plt.text(x,y, '%.1f' % y,
         verticalalignment='top',
         horizontalalignment='right',fontsize=8)

for line in bp['caps']:
    x, y = line.get_xydata()[0] 
    plt.text(x,y, '%.1f' % y,
         verticalalignment='bottom',
         horizontalalignment='left',fontsize=8)  

print('t-test between FUDR and no FUDR gives a p-value of {}'.format(round(
    stats.ttest_ind(mls_FUDR[0],mls_FUDR[3])[1],3)))

print('2B 15˚C (-)FUDR: n = {}'.format(len(mls_FUDR[0])))
print('2B 15˚C FUDR: n = {}'.format(len(mls_FUDR[1])))
print('2B 20˚C (-)FUDR: n = {}'.format(len(mls_FUDR[2])))
print('2B 20˚C FUDR: n = {}'.format(len(mls_FUDR[3])))
print('2B 25˚C (-)FUDR: n = {}'.format(len(mls_FUDR[4])))
print('2B 25˚C FUDR: n = {}'.format(len(mls_FUDR[5])))

######### Saving the Data #############
index = ['(-)FUDR 15˚C','FUDR 15˚C','(-)FUDR 20˚C','FUDR 20˚C','(-)FUDR 25˚C','FUDR 25˚C']
# pd.DataFrame(mls_FUDR,index=index).transpose().to_csv('/Users/nickurban/Desktop/N2 Lifespan Variability/Figures/mls_fig2B.csv')
FUdR_15_20_25_mls = pd.DataFrame(mls_FUDR,index=index).transpose()
FUdR_15_20_25_mls.to_csv('saved_data/FUdR_15_20_25_mls.csv')

#%% Figure 2C-E - LIFESPAN SPREAD GRAPHS
daystoplot = ['% Alive on Day15','% Alive on Day20','% Alive on Day25']
L = len(data_all2)
L2 = len(daystoplot)
ind_count = x15.columns.get_loc('Count')

ind_names_tempFUDR = ['(-)FUDR \n n = {}\n (15°C)'.format(x15.iloc[0,ind_count]),
                          'FUDR \n n = {}\n (15°C)'.format(y15.iloc[0,ind_count]),
                          '(-)FUDR \n n = {}\n (20°C)'.format(x20.iloc[0,ind_count]),
                          'FUDR \n n = {}\n (20°C)'.format(y20.iloc[0,ind_count]),
                          '(-)FUDR \n n = {}\n (25°C)'.format(x25.iloc[0,ind_count]),
                          'FUDR  \n n = {}\n (25°C)'.format(y25.iloc[0,ind_count]),]

for ii in range(L2):
    df_alive,lifespan_lists = plotlifespanspread(data_all2, daystoplot[ii])
    columns = ['Set1','Set2','Set3','Set4','Set5','Set6','Set7','Set8','Set9','Set10',
               'Set11','Set12','Set13','Set14','Set15','Set16','Set17','Set18','Set19',
               'Set20','Set21','Set22','Set23','Set24','Set25','Set26','Set27','Set28',
               'Set29','Set30','Set31']
    df_alive.columns = columns
    df_alive = df_alive.transpose()
    
    ind_tempFUDR = [h.index[-1] for h in xys]
    df_tempFUDR = df_alive.iloc[ind_tempFUDR,:]
    df_tempFUDR.index = ind_names_tempFUDR
    
    ind = [1,2,3.5,4.5,6,7]  # the x locations for the groups
    width = 0.7
    
    fig = plt.figure(figsize=(18,10))
    p1 = plt.bar(ind, df_tempFUDR['20%'], width, color='black')
    p2 = plt.bar(ind, df_tempFUDR['40%'], width, bottom=df_tempFUDR['20%'],color='chocolate')
    p3 = plt.bar(ind, df_tempFUDR['60%'], width, color='dodgerblue',\
                 bottom=[i+j for i,j in zip(df_tempFUDR['20%'],df_tempFUDR['40%'])])
    p4 = plt.bar(ind, df_tempFUDR['80%'], width, color='crimson',\
                 bottom=[i+j+k for i,j,k in zip(df_tempFUDR['20%'],df_tempFUDR['40%'],df_tempFUDR['60%'])])
    p5 = plt.bar(ind, df_tempFUDR['100%'], width, color='seagreen',\
                 bottom=[i+j+k+z for i,j,k,z in zip(df_tempFUDR['20%'],df_tempFUDR['40%'],df_tempFUDR['60%'],df_tempFUDR['80%'])])
    
    plt.ylabel('Percent of Total Entries In Set',fontsize=18)
    plt.ylim(0,100)
    plt.xlim(0.5,7.5)
    plt.xticks(ind, fontsize=18, labels=df_tempFUDR.index)
    plt.yticks(np.arange(0, 110, 10),size=14)
    plt.title('{}'.format(daystoplot[ii]),\
              fontsize=24,y=1.02)
    plt.legend((p1[0], p2[0], p3[0], p4[0], p5[0]), ('Entries with <20% Alive',\
               'Entries with 20-40% Alive','Entries with 40-60% Alive',\
               'Entries with 60-80% Alive','Entries with 80-100% Alive'),\
               bbox_to_anchor=(1.0, 1),loc='upper left',\
               borderaxespad=0.5, fontsize=14)
    plt.tight_layout()
    plt.show()

    [print("%'s n = {}".format(x)) for x in df_tempFUDR['count']]
    if ii == 0:
        df_tempFUDR.to_csv('saved_data/Fig2C.csv')

    elif ii == 1:
        df_tempFUDR.to_csv('saved_data/Fig2D.csv')

    elif ii == 2:
        df_tempFUDR.to_csv('saved_data/Fig2E.csv')


#%% Figure 2F
# Histogram of FUDR concentrations at temp = 20C
fig = plt.figure(figsize=(10,10))
ind_FUDRconc = data_all2.columns.get_loc('FUDR Concentration (microM)')

# Find unique values of FUDR concentrations
FUDRvalues, FUDRfreq = np.unique(data_all2.iloc[(1,ind_FUDRconc)], return_counts=True)

# Find which values have > 4 entries with same FUDR concentration
FUDRtest = np.where(FUDRfreq>4)

FUDRdata = {'Values':FUDRvalues[FUDRtest],'Frequency':FUDRfreq[FUDRtest]}
df_FUDR = pd.DataFrame(FUDRdata)

# Find and delete missing FUDR concentrations
df_FUDR = df_FUDR.drop(df_FUDR[df_FUDR['Values']==-1].index)

for i in range(len(df_FUDR)):
    plt.scatter(df_FUDR.iloc[i,0],df_FUDR.iloc[i,1],label = df_FUDR.iloc[i,0])
    
plt.ylim(0,30,5)
plt.xlim(0,450,50)
plt.title('FUDR Concentrations (20˚C)', fontsize = 16)
plt.ylabel('Frequency', fontsize = 14)
plt.xlabel('Concentration (μM)', fontsize = 14)
plt.legend(bbox_to_anchor=(0.99, 1), loc='upper right', borderaxespad=0.2,
           fontsize = 10, title = 'Unit = μM')
plt.show()

# df_FUDR.to_csv('figures/Fig2F.csv')

#%% Figure 2G
# LIFESPAN PLOTS FOR DIFFERENT FUdR CONCENTRATIONS
df_totFUDR25, df_FUDR25 = get_FUDR(data_all2,25.0)
df_totFUDR40, df_FUDR40 = get_FUDR(data_all2,40.0)
df_totFUDR100, df_FUDR100 = get_FUDR(data_all2,100.0)
df_totFUDR200, df_FUDR200 = get_FUDR(data_all2,200.0)
df_totFUDR400, df_FUDR400 = get_FUDR(data_all2,400.0)

fig, ((ax1, ax2),(ax3, ax4),(ax5, ax6)) = plt.subplots(3,2,figsize=(18,14))
fig.suptitle('FUDR Concentration Comparison', fontsize=22, y=1.01)
color_cycle = ["black", "dodgerblue", "chocolate"] 
ind = [3,5,10,15,20,25,30,40,50]
                 
ax1.plot(ind,df_FUDR25['mean'], color='k',label='mean')
ax1.fill_between(ind, df_FUDR25['mean']-2*df_FUDR25['std'], df_FUDR25['mean']+2*df_FUDR25['std'],  
                 facecolor=color_cycle[2], label='2x std')
ax1.fill_between(ind, df_FUDR25['mean']-df_FUDR25['std'], df_FUDR25['mean']+df_FUDR25['std'], 
                 facecolor=color_cycle[1], label='1x std')

ax2.plot(ind,df_FUDR40['mean'], color='k',label='mean')
ax2.fill_between(ind, df_FUDR40['mean']-2*df_FUDR40['std'], df_FUDR40['mean']+2*df_FUDR40['std'], 
                 facecolor=color_cycle[2], label='2x std')
ax2.fill_between(ind, df_FUDR40['mean']-df_FUDR40['std'], df_FUDR40['mean']+df_FUDR40['std'], 
                 facecolor=color_cycle[1], label='1x std')

ax3.plot(ind,df_FUDR100['mean'], color='k',label='mean')
ax3.fill_between(ind, df_FUDR100['mean']-2*df_FUDR100['std'], df_FUDR100['mean']+2*df_FUDR100['std'],  
                 facecolor=color_cycle[2], label='2x std')
ax3.fill_between(ind, df_FUDR100['mean']-df_FUDR100['std'], df_FUDR100['mean']+df_FUDR100['std'], 
                 facecolor=color_cycle[1], label='1x std')

ax4.plot(ind,df_FUDR200['mean'], color='k',label='mean')
ax4.fill_between(ind, df_FUDR200['mean']-2*df_FUDR200['std'], df_FUDR200['mean']+2*df_FUDR200['std'], 
                 facecolor=color_cycle[2], label='2x std')
ax4.fill_between(ind, df_FUDR200['mean']-df_FUDR200['std'], df_FUDR200['mean']+df_FUDR200['std'], 
                 facecolor=color_cycle[1], label='1x std')

ax5.plot(ind,df_FUDR400['mean'], color='k',label='mean')
ax5.fill_between(ind, df_FUDR400['mean']-2*df_FUDR400['std'], df_FUDR400['mean']+2*df_FUDR400['std'], 
                 facecolor=color_cycle[2], label='2x std')
ax5.fill_between(ind, df_FUDR400['mean']-df_FUDR400['std'], df_FUDR400['mean']+df_FUDR400['std'], 
                 facecolor=color_cycle[1], label='1x std')

ax6.plot(ind,df_FUDR25['mean'], color=color_cycle[0],label='25µM')
ax6.plot(ind,df_FUDR40['mean'], color=color_cycle[1],label='40µM')
ax6.plot(ind,df_FUDR100['mean'], color=color_cycle[2],label='100µM')
ax6.plot(ind,df_FUDR200['mean'], color='purple',label='200µM')
ax6.plot(ind,df_FUDR400['mean'], color='r',label='400µM')

axs = [ax1,ax2,ax3,ax4,ax5,ax6]
[x.set_ylim(-40,130) for x in axs]
[x.set_xticks(ind) for x in axs]
[x.set_xticklabels(df_FUDR40.index, fontsize=12, rotation = 70) for x in axs]
[x.set_yticks(np.arange(-20,140,20)) for x in axs]
[x.set_yticklabels(np.arange(-20,140,20), fontsize=12) for x in axs]
[x.legend(fontsize=14) for x in axs]

ax1.set_title('25 µM', fontsize=14,x=0.5)
ax2.set_title('40 µM', fontsize=14,x=0.5)
ax3.set_title('100 µM', fontsize=14,x=0.5)
ax4.set_title('200 µM', fontsize=14,x=0.5)
ax5.set_title('400 µM', fontsize=14,x=0.5)
ax6.set_title('All Concentrations', fontsize=14,x=0.5)

ax1.set_ylabel('% Survival',fontsize=20)
ax3.set_ylabel('% Survival',fontsize=20)
ax5.set_ylabel('% Survival',fontsize=20)

plt.tight_layout()
plt.show()

df_totFUDR25.to_csv('figures/Fig2G_25.csv')
df_totFUDR40.to_csv('figures/Fig2G_40.csv')
df_totFUDR100.to_csv('figures/Fig2G_100.csv')
df_totFUDR200.to_csv('figures/Fig2G_200.csv')
df_totFUDR400.to_csv('figures/Fig2G_400.csv')

print('FUDR25 n = {}'.format(df_FUDR25['count'][0]))
print('FUDR40 n = {}'.format(df_FUDR40['count'][0]))
print('FUDR100 n = {}'.format(df_FUDR100['count'][0]))
print('FUDR200 n = {}'.format(df_FUDR200['count'][0]))
print('FUDR400 n = {}'.format(df_FUDR400['count'][0]))
#%% Figure 2H - MEAN LIFESPAN PLOTS FOR DIFFERENT FUdR CONCENTRATIONS
ind_mls = df_totFUDR25.columns.get_loc('Reported mean lifespan')
df_FUDR_list = [df_totFUDR25,df_totFUDR40,df_totFUDR100,df_totFUDR200,df_totFUDR400]
mls_FUDR = [[] for _ in range(len(df_FUDR_list))] # Plotting 5 different concentrations

for j in range(len(df_FUDR_list)):
    temp = list(df_FUDR_list[j].iloc[:,ind_mls])
    temp = [x for x in temp if x > -1.0]
    #a = np.array(temp, int)
    mls_FUDR[j] = temp

# Create a figure instance
fig = plt.figure(figsize=(10,6))
ax = fig.add_axes([0,0,1,1])

# Create the boxplot
#bp = ax.violinplot(mls_temp)
bp = ax.boxplot(mls_FUDR)
ax.set_ylim(0,40)
ax.set_xticks([1,2,3,4,5])
ax.set_xticklabels(['[FUDR] = 25µM \n (n = {})'.format(len(mls_FUDR[0])),
                    '[FUDR] = 40µM \n (n = {})'.format(len(mls_FUDR[1])),
                    '[FUDR] = 100µM \n (n = {})'.format(len(mls_FUDR[2])),
                    '[FUDR] = 200µM \n (n = {})'.format(len(mls_FUDR[3])),
                    '[FUDR] = 400µM \n (n = {})'.format(len(mls_FUDR[4]))], fontsize=12)
ax.set_yticks(np.arange(0,45,5))
ax.set_yticklabels(np.arange(0,45,5), fontsize=12)
ax.set_title('Different FUDR Concentrations', fontsize=16,x=0.5) # Change titles here
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
#pd.DataFrame(mls_FUDR,index=['25µM','40µM','100µM','200µM','400µM']).transpose().to_csv('/Users/nickurban/Desktop/N2 Lifespan Variability/Figures/mls_fig2H.csv')
pd.DataFrame(mls_FUDR,index=['25µM','40µM','100µM','200µM','400µM']).transpose().to_csv('figures/Fig2H.csv')

for i in range(len(mls_FUDR)):
    print('mls FUDR: n = {}'.format(len(mls_FUDR[i])))
    
    