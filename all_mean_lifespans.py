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

t = [15, 20, 25]

temp_ind = [0]*len(t) # 0 = 20, 1 = 25, 4 = 15 ('index' = 'temp')
for i in range(len(t)):
    temp = data_all[(data_all['Growth Media'] == 'NGM') &
                    (data_all['Temperature Maintained (Through L4, C)'] == t[i]) & 
                    (data_all['Temperature Cultivated (Adult, C)'] == t[i])]
    temp_ind[i] = temp.index[0]

mls_temp = [[] for _ in range(len(temp_ind))]

for j in range(len(temp_ind)):
    temp = data_all.iloc[temp_ind[j],ind_mls]
    temp = [x for x in temp if x > -1.0]
    mls_temp[j] = temp

df1 = list(mls_temp[0]) # 15˚C
df2 = list(mls_temp[1]) # 20˚C
df3 = list(mls_temp[2]) # 25˚C

ind_mls = data_all2.columns.get_loc('Reported mean lifespan')
mls_FUDR = [[] for _ in range(2)]

#%% FUDR vs No FUDR at 20˚C
for j in range(2):
    temp = data_all2.iloc[j,ind_mls]
    temp = [x for x in temp if x > -1.0]
    #a = np.array(temp, int)
    mls_FUDR[j] = temp
    
df4 = list(mls_FUDR[0]) # no FUdR at 20˚C
df5 = list(mls_FUDR[1]) # wFUdR at 20˚C

#%% FUDR vs No FUDR at 25˚C
for j in range(3,5):
    temp = data_all2.iloc[j,ind_mls]
    temp = [x for x in temp if x > -1.0]
    #a = np.array(temp, int)
    mls_FUDR[j-3] = temp
    
df29 = list(mls_FUDR[0]) # no FUdR at 25˚C
df30 = list(mls_FUDR[1]) # wFUdR at 25˚C

#%% FUDR Concentrations
df_totFUDR25, df_FUDR25 = get_FUDR(data_all2,25.0)
df_totFUDR40, df_FUDR40 = get_FUDR(data_all2,40.0)
df_totFUDR100, df_FUDR100 = get_FUDR(data_all2,100.0)
df_totFUDR200, df_FUDR200 = get_FUDR(data_all2,200.0)
df_totFUDR400, df_FUDR400 = get_FUDR(data_all2,400.0)
ind_mls = df_totFUDR25.columns.get_loc('Reported mean lifespan')
df_FUDR_list = [df_totFUDR25,df_totFUDR40,df_totFUDR100,df_totFUDR200,df_totFUDR400]
mls_FUDR = [[] for _ in range(len(df_FUDR_list))] # Plotting 5 different concentrations

for j in range(len(df_FUDR_list)):
    temp = list(df_FUDR_list[j].iloc[:,ind_mls])
    temp = [x for x in temp if x > -1.0]
    #a = np.array(temp, int)
    mls_FUDR[j] = temp

df31 = list(mls_FUDR[0]) # 25 FUdR at 20˚C
df32 = list(mls_FUDR[1]) # 40 FUdR at 20˚C
df33 = list(mls_FUDR[2]) # 100 FUdR at 20˚C
df34 = list(mls_FUDR[3]) # 200 FUdR at 20˚C
df35 = list(mls_FUDR[4]) # 400 FUdR at 20˚C

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
dfs = [df1,df2,df3,df4,df5,df6,df7,df8,df9,df10,df11,df12,df13,df14,df15,df16,\
       df17,df18,df19,df20,df21,df22,df23,df24,df25,df26,df27,df28,df29,df30,\
       df31,df32,df33,df34,df35]

index = ['15˚C','20˚C','25˚C','no FUdR at 20˚C','with FUdR at 20˚C','no manipulation at 15˚C',\
        'with any manipulation at 15˚C','no manipulation at 20˚C','with manipulation at 20˚C',\
        'no manipulation at 25˚C','with manipulation at 25˚C','Everyday manipulation at 20˚C',\
        'Every 2-3 day manipulation at 20˚C','no manipulation no FUdR at 20˚C',\
        'no manipulation with FUdR at 20˚C','with manipulation no FUdR at 20˚C',\
        'with manipulation with FUdR at 20˚C', 'USA','China','Germany',\
        'Rep. of Korea','Japan','CA','MA','NY','WA','TX','MI','no FUdR at 25˚C',\
        'with FUdR at 25˚C','[25] FUDR at 20˚C','[40] FUDR at 20˚C','[100] FUDR at 20˚C',\
        '[200] FUDR at 20˚C','[400] FUDR at 20˚C']

all_mls_data = pd.DataFrame(dfs,index=index).transpose()
all_mls_data.to_csv('all_mls_data.csv')

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

labels_temp = ['15˚C','20˚C','25˚C']
stat_temp = mean_confidence_interval([df1,df2,df3],labels_temp,(7,5))

labels_fudr20 = ['(-)FUDR 20˚C','FUDR 20˚C']
stat_FUDR = mean_confidence_interval([df4,df5],labels_fudr20,(7,5))

labels_fudr25 = ['(-)FUDR 25˚C','FUDR 25˚C']
stat_FUDR25 = mean_confidence_interval([df29,df30],labels_fudr25,(7,5))

labels_fudr_conc = ['[FUDR] = 25µM','[FUDR] = 40µM','[FUDR] = 100µM',
                    '[FUDR] = 200µM','[FUDR] = 400µM']
stat_FUDRc = mean_confidence_interval([df31,df32,df33,df34,df35],labels_fudr_conc,(11,4))

labels_ttfp = ['No Manipulation 20˚C','Manipulation 20˚C']
stat_ttfp = mean_confidence_interval([df8,df9],labels_ttfp,(7,5))

labels_cntry = ['USA','China','Germany','Republic of Korea','Japan']
stat_cntry = mean_confidence_interval([df18,df19,df20,df21,df22],labels_cntry,(10,4))

labels_states = ['CA','MA','NY','WA','TX','MI']
stat_states = mean_confidence_interval([df23,df24,df25,df26,df27,df28],labels_states,(10,4))

#%% ANOM
# Group data by analyses
rml_temp = [df1,df2,df3] # each 'df' is a list
temp_crit_value = 2.39 # Based on table from paper given by Santiago

rml_fudr20 = [df4,df5]
fudr20_crit_value = 2.24

rml_fudr25 = [df29,df30]
fudr25_crit_value = 2.35

rml_fudrc = [df31,df32,df33,df34,df35]
fudrc_crit_value = 2.65

rml_ttfp = [df8,df9]
ttfp_crit_value = 2.24

rml_cntry = [df18,df19,df20,df21,df22]
cntry_crit_value = 2.57

rml_state = [df23,df24,df25,df26,df27,df28]
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

    