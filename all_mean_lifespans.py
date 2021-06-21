#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 14:00:05 2021

@author: josephcavataio
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import os

from N2_functions import get_ttfp,randomvals_and_diff,residuals,pairplot,make_groups,\
    plotlifespanspread,categorize,plotmeanlifespan, get_FUDR,mean_confidence_interval
    
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

#%% Histogram of Entire Dataset and run Shapiro-Wilk test for normality
rml_plot = data_new[data_new['Reported mean lifespan'] > -1]
vals = stats.shapiro(rml_plot['Reported mean lifespan'])
statistics = rml_plot['Reported mean lifespan'].describe()

# Calculate Confidence Intervals

if vals.pvalue < 0.05:
    print(f'\nSufficient data to say data comes from a normal distribution (Shapiro-Wilk statistic = {vals.statistic})')
else:
    print(f'Data does not come from a normal distribution (Shapiro-Wilk statistic = {vals.statistic})')

mean = round(statistics['mean'],1)
med = statistics['50%']
plt.hist(rml_plot['Reported mean lifespan'],label = 'mean = {}\nmed = {}'.format(mean,med))
plt.xlabel('Days Alive')
plt.ylabel('Frequency')
plt.legend()
plt.show()


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

#%% Gather all assortments of mean lifespan data

# Temperatures 15,20,25˚C
ind_mls = data_all.columns.get_loc('Reported mean lifespan')
temp_ind = [4,0,1] # 0 = 20, 1 = 25, 4 = 15 ('index' = 'temp')
mls_temp = [[] for _ in range(len(temp_ind))]

for j in range(len(temp_ind)):
    temp = data_all.iloc[temp_ind[j],ind_mls]
    temp = [x for x in temp if x > -1.0]
    mls_temp[j] = temp

df1 = mls_temp[0] # 15˚C
df2 = mls_temp[1] # 20˚C
df3 = mls_temp[2] # 25˚C

ind_mls = data_all2.columns.get_loc('Reported mean lifespan')
mls_FUDR = [[] for _ in range(2)]

#%% FUDR vs No FUDR at 20˚C
for j in range(2):
    temp = data_all2.iloc[j,ind_mls]
    temp = [x for x in temp if x > -1.0]
    #a = np.array(temp, int)
    mls_FUDR[j] = temp
    
df4 = mls_FUDR[0] # no FUdR at 20˚C
df5 = mls_FUDR[1] # wFUdR at 20˚C

#%% Any plate manipulation vs no manipulation at 15˚C
df_totnm,df_totm,data_nm,data_m,df_totnm_all,df_totm_all = get_ttfp(data_all3,'Every other day','After L4',True,15,15,False,False)
ind_mls = df_totnm.columns.get_loc('Reported mean lifespan')

temp = df_totnm[df_totnm['Reported mean lifespan'] > -1]
df6 = list(temp.iloc[:,ind_mls]) # no manipulation

temp = df_totm[df_totm['Reported mean lifespan'] > -1]
df7 = list(temp.iloc[:,ind_mls]) # manipulation

#%% Any plate manipulation vs no manipulation at 20˚C
df_totnm,df_totm,data_nm,data_m,df_totnm_all,df_totm_all = get_ttfp(data_all3,'Every other day','After L4',True,20,20,False,False)

temp = df_totnm[df_totnm['Reported mean lifespan'] > -1]
df8 = list(temp.iloc[:,ind_mls]) # no manipulation

temp = df_totm[df_totm['Reported mean lifespan'] > -1]
df9 = list(temp.iloc[:,ind_mls]) # manipulation

#%% Any plate manipulation vs no manipulation at 25˚C
df_totnm,df_totm,data_nm,data_m,df_totnm_all,df_totm_all = get_ttfp(data_all3,'Every other day','After L4',True,25,25,False,False)

temp = df_totnm[df_totnm['Reported mean lifespan'] > -1]
df10 = list(temp.iloc[:,ind_mls]) # no manipulation

temp = df_totm[df_totm['Reported mean lifespan'] > -1]
df11 = list(temp.iloc[:,ind_mls]) # manipulation

#%% Everyday vs Every 2-3 day manipulation at 20˚C
df_totED_20,df_totE23D_20,data_ED_20,data_E23D_20,df_totED_20_all,df_totE23D_20_all = get_ttfp(data_all3,'Everyday','Every 2-3 days',False,20,20,False,False)

temp = df_totED_20[df_totED_20['Reported mean lifespan'] > -1]
df12 = list(temp.iloc[:,ind_mls]) # no manipulation

temp = df_totE23D_20[df_totE23D_20['Reported mean lifespan'] > -1]
df13 = list(temp.iloc[:,ind_mls]) # manipulation

#%% Non-FUdR no manipulation vs FUdR no manipulation at 20˚C
isFUDR1 = True
isFUDR2 = False 
df_totnm_noFUDR,df_totnm_wFUDR,data_nm_noFUDR,data_nm_wFUDR,df_totnm_all_noFUDR,df_totnm_all_wFUDR \
    = get_ttfp(data_all4,'Every other day','After L4',False,20,20,True,False)
    
temp = df_totnm_noFUDR[df_totnm_noFUDR['Reported mean lifespan'] > -1]
df14 = list(temp.iloc[:,ind_mls]) # no manipulation

temp = df_totnm_wFUDR[df_totnm_wFUDR['Reported mean lifespan'] > -1]
df15 = list(temp.iloc[:,ind_mls]) # manipulation

#%% Non-FUdR any manipulation vs FUdR any manipulation at 20˚C
isFUDR1 = False
isFUDR2 = True
df_totm_noFUDR,df_totm_wFUDR,data_m_noFUDR,data_m_wFUDR,df_totm_all_noFUDR,df_totm_all_wFUDR \
    = get_ttfp(data_all4,'Every other day','After L4',False,20,20,False,True)
    
temp = df_totnm_noFUDR[df_totnm_noFUDR['Reported mean lifespan'] > -1]
df16 = list(temp.iloc[:,ind_mls]) # no manipulation

temp = df_totnm_wFUDR[df_totnm_wFUDR['Reported mean lifespan'] > -1]
df17 = list(temp.iloc[:,ind_mls]) # manipulation

#%% Countries
set_data,unique_group_countries,avg_set=categorize(data_all,'Lab Location (Country)', 0, '% Alive on Day25')
for j in range(len(avg_set)):
    temp = list(avg_set[j])
    temp = [x for x in temp if x > -1.0]
    avg_set[j] = temp

df18 = avg_set[0] # USA 20C
df19 = avg_set[1] # China 20C
df20 = avg_set[2] # Germany 20C
df21 = avg_set[3] # Republic of Korea 20C
df22 = avg_set[4] # Japan 20C

#%% States
set_data,unique_group_states,avg_set=categorize(data_all,'Lab Location (State/Province)', 0, '% Alive on Day25')

# Remove values for where State/Province == -1
x = unique_group_states.index[unique_group_states['Names'] == '-1'].tolist()
unique_group_states = unique_group_states.drop(x[0])
set_data.pop(x[0])
avg_set.pop(x[0])

for j in range(len(avg_set)):
    temp = list(avg_set[j])
    temp = [x for x in temp if x > -1.0]
    avg_set[j] = temp

df23 = avg_set[0] # CA 20C
df24 = avg_set[1] # MA 20C
df25 = avg_set[2] # NY 20C
df26 = avg_set[3] # WA 20C
df27 = avg_set[4] # TX 20C
df28 = avg_set[5] # MI 20C

#%% Create large list with all data
dfs = [df1,df2,df3,df4,df5,df6,df7,df8,df9,df10,df11,df12,df13,df14,df15,df16,df17]

index = ['15˚C','20˚C','25˚C','no FUdR at 20˚C','with FUdR at 20˚C','no manipulation at 15˚C',\
        'with any manipulation at 15˚C','no manipulation at 20˚C','with manipulation at 20˚C',\
        'no manipulation at 25˚C','with manipulation at 25˚C','Everyday manipulation at 20˚C',\
        'Every 2-3 day manipulation at 20˚C','no manipulation no FUdR at 20˚C',\
        'no manipulation with FUdR at 20˚C','with manipulation no FUdR at 20˚C',\
        'with manipulation with FUdR at 20˚C']

#pd.DataFrame(dfs,index=index).transpose().to_csv('all_mls_data.csv')

#%% All mean lifespans
ind_mls = data_all.columns.get_loc('Reported mean lifespan')
rml = [[] for _ in range(len(data_all))]

for i in range(len(data_all)):
    rml[i] = [i for i in data_all.iloc[i,ind_mls] if i > -1]
    
data_rml = pd.DataFrame(rml).transpose()

ind_mls = data.columns.get_loc('Reported mean lifespan')
all_rml = pd.DataFrame({'All Reported mean lifespan data': \
                        [i for i in data['Reported mean lifespan'] if i > -1]})
all_rml = all_rml.sort_values('All Reported mean lifespan data')
print('% of entries with reported mean lifespan = {}%'\
      .format(round(100*len(all_rml)/len(data),2)))

data_rml.to_csv('saved_data/grouped_by_temp_rmls.csv')
all_rml.to_csv('saved_data/total_dataset_rmls.csv')

# Calculate Mean, CI, and Quartiles
import math
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

labels = ['15C','20C','25C']
mean_temp, CI_temp, stat_temp = mean_confidence_interval([df1,df2,df3],labels)

labels = ['noFUDR','wFUDR']
mean_FUDR, CI_FUDR, stat_FUDR = mean_confidence_interval([df4,df5],labels)

labels = ['no manipulation 20C','manipulation 20C']
mean_ttfp, CI_ttfp, stat_ttfp = mean_confidence_interval([df8,df9],labels)

labels = ['USA','China','Germany','Republic of Korea','Japan']
mean_cntry, CI_cntry, stat_cntry = mean_confidence_interval([df18,df19,df20,df21,df22],labels)

labels = ['CA','MA','NY','WA','TX','MI']
mean_states, CI_states, stat_states = mean_confidence_interval([df23,df24,df25,df26,df27,df28],labels)