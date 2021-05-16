#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 12:05:28 2020

@author: josephcavataio
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import os


#### Function for creating groups of entries ####
def do_N2_stats(x):
    otpt = pd.Series(x).describe()
    return otpt

def select_outs(x,l,u):
    X = np.asarray(x)
    opts = X[(X < l) | (X > u)]
    return opts

def data_stats(X):
    L = len(X)
    
    ## Day 20 Range Calculation ##
    dataX_20 = pd.DataFrame([do_N2_stats(X['% Alive on Day20'][i]) \
                          for i in range(L)])
    low_lim_20 = dataX_20['25%'] - 1.5*(dataX_20['75%']-dataX_20['25%'])
    up_lim_20 = dataX_20['75%'] + 1.5*(dataX_20['75%']-dataX_20['25%'])
    outs_20 = pd.Series([select_outs(X['% Alive on Day20'][i],low_lim_20[i],up_lim_20[i])\
                      for i in range(L)])
    len_outs_20 = pd.Series([len(outs_20[i]) for i in range(L)])
    perc_outs_20 = 100*(len_outs_20/dataX_20['count'])
    Q_20 = np.asarray([stats.shapiro(X['% Alive on Day20'][i])\
                    for i in range(L) ])
    pvals_20 = pd.Series(Q_20[:,1])
    
    ## Day 25 Range Calculation ##
    dataX_25 = pd.DataFrame([do_N2_stats(X['% Alive on Day25'][i]) \
                          for i in range(L)])
    low_lim_25 = dataX_25['25%'] - 1.5*(dataX_25['75%']-dataX_25['25%'])
    up_lim_25 = dataX_25['75%'] + 1.5*(dataX_25['75%']-dataX_25['25%'])
    outs_25 = pd.Series([select_outs(X['% Alive on Day25'][i],low_lim_25[i],up_lim_25[i])\
                      for i in range(L)])
    len_outs_25 = pd.Series([len(outs_25[i]) for i in range(L)])
    perc_outs_25 = 100*(len_outs_25/dataX_25['count'])
    Q_25 = np.asarray([stats.shapiro(X['% Alive on Day25'][i])\
                    for i in range(L) ])
    pvals_25 = pd.Series(Q_25[:,1])

    dataX_20['LL'] = low_lim_20
    dataX_20['UL'] = up_lim_20
    dataX_20['%outliers'] = perc_outs_20
    dataX_20['p-value'] = pvals_20
    
    dataX_25['LL'] = low_lim_25
    dataX_25['UL'] = up_lim_25
    dataX_25['%outliers'] = perc_outs_25
    dataX_25['p-value'] = pvals_25
    return dataX_20, dataX_25

def make_groups(data,entries): #IsBuff,IsLogTransf:
    ## Group the entries according to criteria in the groups variable ##
    data0 = data.groupby(entries)
    groupindex = list(data0.indices)
    
    lifespan3 = []
    lifespan5 = []
    lifespan10 = []
    lifespan15 = []
    lifespan20 = []
    lifespan25 = []
    lifespan30 = []
    lifespan40 = []
    lifespan50 = []
    lifespanavg = []
    FUDRconc = []
    country = []
    PM_ID = []
    trantofresh = []
    alive50 = []
    
    L = len(groupindex)
    
    for i in range(L):
    
        example = data0.get_group(groupindex[i])
        lifespan3.append(list(example['% Alive on Day3']))
        lifespan5.append(list(example['% Alive on Day5']))
        lifespan10.append(list(example['% Alive on Day10']))
        lifespan15.append(list(example['% Alive on Day15']))
        lifespan20.append(list(example['% Alive on Day20']))
        lifespan25.append(list(example['% Alive on Day25']))
        lifespan30.append(list(example['% Alive on Day30']))
        lifespan40.append(list(example['% Alive on Day40']))
        lifespan50.append(list(example['% Alive on Day50']))
        lifespanavg.append(list(example['Reported mean lifespan']))
        FUDRconc.append(list(example['FUDR Concentration (microM)']))
        country.append(list(example['Lab Location (Country)']))
        PM_ID.append(list(example['PM ID']))
        trantofresh.append(list(example['# Transfers to Fresh Plates']))
        alive50.append(list(example['50% Alive']))
    
    grp = pd.DataFrame({'Count': data0.size()}).reset_index()
    grp['% Alive on Day3'] = pd.Series(lifespan3)
    grp['% Alive on Day5'] = pd.Series(lifespan5)
    grp['% Alive on Day10'] = pd.Series(lifespan10)
    grp['% Alive on Day15'] = pd.Series(lifespan15)
    grp['% Alive on Day20'] = pd.Series(lifespan20)
    grp['% Alive on Day25'] = pd.Series(lifespan25)
    grp['% Alive on Day30'] = pd.Series(lifespan30)
    grp['% Alive on Day40'] = pd.Series(lifespan40)
    grp['% Alive on Day50'] = pd.Series(lifespan50)
    grp['Reported mean lifespan'] = pd.Series(lifespanavg)
    grp['FUDR Concentration (microM)'] = pd.Series(FUDRconc)
    grp['Lab Location (Country)'] = pd.Series(country)
    grp['PM ID'] = pd.Series(PM_ID)
    grp['# Transfers to Fresh Plates'] = pd.Series(trantofresh)
    grp['50% Alive'] = pd.Series(alive50)
    
    grp  = grp.reindex(np.argsort(grp['Count'])[::-1])
    grp  = grp.reset_index(drop=True)
    dat = grp[grp["Count"]>=3]
    
    return dat

def plotlifespanspread(data, daytoplot):
    L = len(data)
    
    lifespan_lists = [[] for _ in range(L)]
    perc_alive_lists = [[] for _ in range(5)]
    
    ### Create bar graph ###
    ind_day = data.columns.get_loc(daytoplot)
    
    for i in range(L):
        lifespan_lists[i] = data.iloc[(i,ind_day)]
        
    for i in range(L):
        count20 = 0
        count40 = 0
        count60 = 0
        count80 = 0
        count100 = 0
        
        for j in range(len(lifespan_lists[i])):
            
            if lifespan_lists[i][j] <= 20:
                count20 = count20 + 1
                
            elif ((lifespan_lists[i][j] <= 40) and (lifespan_lists[i][j] > 20)):
                count40 = count40 + 1
                
            elif ((lifespan_lists[i][j] <= 60) and (lifespan_lists[i][j] > 40)):
                count60 = count60 + 1
                
            elif ((lifespan_lists[i][j] <= 80) and (lifespan_lists[i][j] > 60)):
                count80 = count80 + 1
                
            elif ((lifespan_lists[i][j] <= 100) and (lifespan_lists[i][j] > 80)):
                count100 = count100 + 1
                
                
        percent20 = count20/len(lifespan_lists[i])
        percent40 = count40/len(lifespan_lists[i])
        percent60 = count60/len(lifespan_lists[i])
        percent80 = count80/len(lifespan_lists[i])
        percent100 = count100/len(lifespan_lists[i])
        #print(sum([count25, count50, count75, count100]))
        #print(sum([percent25,percent50,percent75,percent100]))
        
        perc_alive_lists[0].append(percent20*100)
        perc_alive_lists[1].append(percent40*100)
        perc_alive_lists[2].append(percent60*100)
        perc_alive_lists[3].append(percent80*100)
        perc_alive_lists[4].append(percent100*100)
    
    
    df_a_lists = pd.DataFrame(perc_alive_lists,index=['20%', '40%','60%','80%','100%'])
    columns = ['Set1','Set2','Set3','Set4','Set5','Set6','Set7','Set8','Set9','Set10',
               'Set11','Set12','Set13','Set14','Set15','Set16','Set17','Set18','Set19',
               'Set20','Set21','Set22','Set23','Set24','Set25']
    df_a_lists.columns = columns
    df_a_lists = df_a_lists.transpose()
    #df_a_lists = df_a_lists.sort_values(by=['25%'],ascending=False)
    
    return df_a_lists, lifespan_lists

def plotday20(data,column_name, N):
    L = len(data)

    lists = [[] for _ in range(L)]
    ind_group = data.columns.get_loc(column_name)
    
    for i in range(L):
        lists[i] = data.iloc[(i,ind_group)]
    
    group,freq_group = np.unique(data.iloc[(N,ind_group)],return_counts=True)
    unique_group=pd.DataFrame(list(zip(group,freq_group)),columns=['Names','Frequency'])
    unique_group = unique_group.reindex(np.argsort(unique_group['Frequency'])[::-1])
    unique_group  = unique_group.reset_index(drop=True)
    
    lists_set = [[] for _ in range(len(group))]
    df_lifespan = pd.DataFrame(lifespan20_list[N])
    for i in range(len(group)):
        index = np.where(np.isin(lists[N],unique_group.iloc[(i,0)]))
        d20 = df_lifespan.iloc[(index)]
        lists_set[i] = d20
    
    return lists_set, unique_group, ind_group

def plotmeanlifespan(data,column_name, N):
    ind_group = data.columns.get_loc(column_name)
    ind_meanls = data.columns.get_loc('Reported mean lifespan')
    
    # Create temporary list to hold "country" data for set N
    templist = [data.iloc[(N,ind_group)]]
    
    # Convert list of mean lifespans for set N to dataframe
    df_meanls = pd.DataFrame(data.iloc[N,ind_meanls])
    
    # Find how many times each "country" appears in set N
    group,freq_group = np.unique(data.iloc[(N,ind_group)],return_counts=True)
    unique_group=pd.DataFrame(list(zip(group,freq_group)),columns=['Names','Frequency'])
    unique_group = unique_group.reindex(np.argsort(unique_group['Frequency'])[::-1])
    unique_group = unique_group.reset_index(drop=True)
      
    lists_set = [[] for _ in range(len(group))]
    for i in range(len(group)):
        # Find where each country appears in list of "countries"
        index = np.where(np.isin(templist,unique_group.iloc[(i,0)]))
        
        # Save mean lifespans for each "country" and drop missing values
        meanls = df_meanls.iloc[(index[1])]
        lists_set[i] = meanls
        lists_set[i] = lists_set[i].reset_index(drop=True)
        ind_minus1 = np.where(np.isin(lists_set[i],-1))
        lists_set[i] = lists_set[i].drop(ind_minus1[0])
        
        # Update frequency of "country" in unique_group
        unique_group.iloc[(i,1)] = len(lists_set[i])
        
    return lists_set, unique_group, ind_group


#%%# Read in N2 Lifespans Data ####
path = os.getcwd()
filename = 'N2 Lifespans.xlsx'
path = os.path.join(path, filename)
numberofentries = 706
rowstoskip = 998 - numberofentries
data = pd.read_excel(path, index_col=0,skipfooter=rowstoskip)

#### Percent of data excluded for different categories ####
strainused =  data['N2 Strain Used']
strainsource = data['N2 Strain Source']
growthmedia = data['Growth Media']
rmeanls = data['Reported mean lifespan']
rmedls = data['Reported median lifespan']
minanimperplate = data['Min #Animals/Plate']
maxanimperplate = data['Max #Animals/Plate']
plateper = data['#Plates/Experiment']
FUDRbool = data['FUDR (Yes/No?)']
FUDRconc = data['FUDR Concentration (microM)']
transfertofresh = data['# Transfers to Fresh Plates']
tempm = data['Temperature Maintained (Through L4, C)']
tempc = data['Temperature Cultivated (Adult, C)']
loccity = data['Lab Location (City)']
loccountry = data['Lab Location (Country)']

missingstrainused = strainused.isna().sum()
print("Percent of entries missing N2 Strain Used")
print(round(100*missingstrainused/len(data),2), '%', '{}/{}'.format(missingstrainused,numberofentries))

missingstrainsource = strainsource.isna().sum()
print("Percent of entries missing N2 Strain Source")
print(round(100*missingstrainsource/len(data),2), '%', '{}/{}'.format(missingstrainsource,numberofentries))

missinggrowth = growthmedia.isna().sum()
print("Percent of entries missing Growth Media")
print(round(100*missinggrowth/len(data),2), '%', '{}/{}'.format(missinggrowth,numberofentries))

missingrmeanls = rmeanls.isna().sum()
print("Percent of entries missing Reported Mean Lifespan")
print(round(100*missingrmeanls/len(data),2), '%', '{}/{}'.format(missingrmeanls,numberofentries))

missingrmedls = rmedls.isna().sum()
print("Percent of entries missing Reported Median Lifespan")
print(round(100*missingrmedls/len(data),2), '%', '{}/{}'.format(missingrmedls,numberofentries))

missingminanim = minanimperplate.isna().sum()
print("Percent of entries missing Min # Animals/Plate")
print(round(100*missingminanim/len(data),2), '%', '{}/{}'.format(missingminanim,numberofentries))

missingmaxanim = maxanimperplate.isna().sum()
print("Percent of entries missing Max # Animals/Plate")
print(round(100*missingmaxanim/len(data),2), '%', '{}/{}'.format(missingmaxanim,numberofentries))

missingplateper = plateper.isna().sum()
print("Percent of entries missing # Plates per Experiment")
print(round(100*missingplateper/len(data),2), '%', '{}/{}'.format(missingplateper,numberofentries))

missingFUDRconc = FUDRconc.isna().sum()
print("Percent of entries missing FUDR Concentration Info")
print(round(100*missingFUDRconc/len(data),2), '%', '{}/{}'.format(missingFUDRconc,numberofentries))

missingtransfer = transfertofresh.isna().sum()
print("Percent of entries missing # Transfers to Fresh Plates Info")
print(round(100*missingtransfer/len(data),2), '%', '{}/{}'.format(missingtransfer,numberofentries))

missingtempm = tempm.isna().sum()
print("Percent of entries missing # Temp Maintained")
print(round(100*missingtempm/len(data),2), '%', '{}/{}'.format(missingtempm,numberofentries))

missingtempc = tempc.isna().sum()
print("Percent of entries missing # Temp Cultivated")
print(round(100*missingtempc/len(data),2), '%', '{}/{}'.format(missingtempc,numberofentries))

missingloccity = loccity.isna().sum()
print("Percent of entries missing # Cities")
print(round(100*missingloccity/len(data),2), '%', '{}/{}'.format(missingloccity,numberofentries))

missingloccountry = loccountry.isna().sum()
print("Percent of entries missing # Countries")
print(round(100*missingloccountry/len(data),2), '%', '{}/{}'.format(missingloccountry,numberofentries))

#%%########## Array of # Entries, # of Labs, and # of Countries ##############
labs, freqslabs = np.unique(data['Lab'],return_counts = True)
cntry, freqscntry = np.unique(data['Lab Location (Country)'],return_counts = True)
alldatastats=pd.DataFrame([len(data),len(labs),len(cntry)],
                         index=['# of Entries','# of Labs','# of Countries']).transpose()

#%%# Clean-up: Drop columns not necessary for analysis ####
# data = data.drop(['Paper'],axis=1)
# #data = data.drop(['PM ID'],axis=1)
# data = data.drop(['Lab'],axis=1)
# data = data.drop(['Name of Person That Performed the Experiments'],axis=1)
# data = data.drop(['Date Published'],axis=1)
# data = data.drop(['Date Received'],axis=1)
# #data = data.drop(['Growth Media Supplier'],axis=1)
# data = data.drop(['N2 Strain Used'],axis=1)
# data = data.drop(['# Replicate Experiments'],axis=1)
# #data = data.drop(['#Plates/Experiment'],axis=1)
# data = data.drop(['Explosion Through Vulvua Censored or Counted?'],axis=1)
# data = data.drop(['Censoring Used? Yes/No?'],axis=1)
# data = data.drop(['Software Used for Statistics (Excel, PRISM, etc.)'],axis=1)
# data = data.drop(['Notes'],axis=1)

#%%# Assumption: Temp cultivated = temp maintained if no data is present ####
data_new = data.fillna(-1)
L_new = len(data_new)
ind_tempm = data_new.columns.get_loc('Temperature Maintained (Through L4, C)')
ind_tempc = data_new.columns.get_loc('Temperature Cultivated (Adult, C)')

for i in range(L_new):
    if data_new.iloc[(i,ind_tempm)] == -1:
        data_new.iloc[(i,ind_tempm)] = data_new.iloc[(i,ind_tempc)]
    elif data_new.iloc[(i,ind_tempc)] == -1:
        data_new.iloc[(i,ind_tempc)] = data_new.iloc[(i,ind_tempm)]

#%%###################### INTRO GRAPH ########################################
days = ['% Alive on Day3','% Alive on Day5','% Alive on Day10','% Alive on Day15',
              '% Alive on Day20','% Alive on Day25','% Alive on Day30',
              '% Alive on Day40','% Alive on Day50']
lists_s = [[] for _ in range(len(days))]

# Take stats for each day and save into list of lists
for i in range(len(days)):
    s = data_new[days[i]].describe()
    lists_s[i] = [round(s['mean'],2),s['25%'],s['50%'],s['75%'],s['std']]

#Transpose list of lists
lists_s = list(map(list, zip(*lists_s))) 

# Turn list of lists into dataframe
df_totdata = pd.DataFrame(lists_s,index=['mean','25%','50%','75%','std'], 
                          columns=['Day3','Day5','Day10','Day15','Day20',
                                   'Day25','Day30','Day40','Day50']).transpose()

fig = plt.figure(figsize=(12,6))
color_cycle = ['g', 'b', 'c']
ind = np.arange(len(days))+1

plt.plot(ind,df_totdata['mean'], color='k', label='mean')
plt.fill_between(ind, df_totdata['mean']-2*df_totdata['std'], df_totdata['mean']+2*df_totdata['std'],
                 alpha=0.3, facecolor=color_cycle[0], label='2x std')
plt.fill_between(ind, df_totdata['mean']-df_totdata['std'], df_totdata['mean']+df_totdata['std'],
                 alpha=0.3, facecolor=color_cycle[1], label='1x std')

plt.ylabel('% Alive',fontsize=16)
plt.xticks(ind,df_totdata.index, fontsize=14, rotation = 70)
plt.yticks(np.arange(0,160,20), fontsize=14)
plt.ylim(-40,150)
plt.title('The Entire Data Set',fontsize=20,y=1.02)
plt.legend(bbox_to_anchor=(0.99, 1), loc='upper right', borderaxespad=0.2,\
            fontsize = 10)
plt.show()



#%%###################### LIFESPAN TO CONDITIONs #############################

## Adding a column: last day with >50% of animals alive for each entry
daystoplot = [3,5,10,15,20,25,30,40]
data_new['50% Alive'] = np.nan
ind_50 = data_new.columns.get_loc('50% Alive')

for j in range(L_new):
    count=0
    for i in range(len(daystoplot)):
        if count==0:
            daystring = '% Alive on Day{}'.format(daystoplot[i])
            ind_day = data_new.columns.get_loc(daystring)
                    
            if data_new.iloc[(j,ind_day)] <= 50:
                data_new.iloc[(j,ind_50)] = daystoplot[i]
                count+=1

#%%###################### CONDITIONS TO LIFESPAN #############################
growth, freqsgrowth = np.unique(data['Growth Media'].astype(str),return_counts = True)
unique_growth=pd.DataFrame(list(zip(growth,freqsgrowth)),columns=['names','frequency'])

app, freqsapp = np.unique(data['Min #Animals/Plate'].astype(str),return_counts = True)
unique_app=pd.DataFrame(list(zip(app,freqsapp)),columns=['names','frequency'])

#### Group entries by the following groups ####
grouptypes = ['Growth Media','Temperature Maintained (Through L4, C)',\
              'Temperature Cultivated (Adult, C)','FUDR (Yes/No?)']

data_all = make_groups(data_new,grouptypes)

# MAKE NEW DATASET HERE #
grouptypes2 = ['Growth Media','Temperature Maintained (Through L4, C)',\
              'Temperature Cultivated (Adult, C)','FUDR (Yes/No?)','# Transfers to Fresh Plates']

data_all2 = make_groups(data_new,grouptypes2) 

### Saving the Sets to Excel ###
ind_count = data_all.columns.get_loc('Count')
count = data_all.iloc[:,ind_count]
data_display = data_all.drop(labels=['Count'], axis=1,inplace = False)
data_display.insert(0, 'Count', count)

ind_FUDR = data_all.columns.get_loc('FUDR (Yes/No?)')
ind_transf = data_all.columns.get_loc('# Transfers to Fresh Plates')

coltodrop = ['% Alive on Day3','% Alive on Day5','% Alive on Day10','% Alive on Day15',
             '% Alive on Day20','% Alive on Day25','Reported mean lifespan',
             'Lab Location (Country)','PM ID','# Transfers to Fresh Plates',]
data_display.drop(columns=coltodrop,inplace=True)

data_display.to_excel('N2 Lifespans Data Sets.xlsx', index=False)

#%%# Plotting before statistical analysis ####
import random
L = len(data)
ind_20 = data_all.columns.get_loc('Reported mean lifespan')
alive20_list_tot1 = [0]*L
alive20_list_tot2 = [0]*L

for i in range(500): # Taking x samples from entire dataset
    temp = random.sample(list(data_all.iloc[0,ind_20]),2)
    alive20_list_tot1[i] = temp[0]
    alive20_list_tot2[i] = temp[1]
    
df_pairs_tot = pd.DataFrame({'value 1': alive20_list_tot1, 'value 2': alive20_list_tot2})

plt.scatter(df_pairs_tot['value 1'],df_pairs_tot['value 2'],color='k',alpha=0.3,edgecolor=None)
plt.xticks(range(0,35,5), fontsize=12)
plt.yticks(range(0,35,5), fontsize=12)
plt.xlabel('Avg Lifespan 1',fontsize=12)
plt.ylabel('Avg Lifespan 2',fontsize=12)
plt.title('Plot of Pairs of Lifespan Values',
          fontsize = 16, y=1.02)
plt.savefig('figures/Pair_plot_3vals_all_MRCs.pdf',bbox_inches='tight')

#%%
plt.figure()
plt.scatter(df_pairs_tot['value 1'],df_pairs_tot['value 2'],color='k',alpha=0.3,edgecolor=None)
plt.plot(range(0,11,1),range(0,11,1),color='r')
plt.xticks(range(0,12,2), fontsize=12)
plt.yticks(range(0,12,2), fontsize=12)
plt.xlabel('$pK_{M1}$',fontsize=12)
plt.ylabel('$pK_{M2}$',fontsize=12)
plt.title('Plot of Pairs of Km Values Taken Randomly from Entire Dataset (only MM/{} Entries)'.format(L),
          fontsize = 16, y=1.02)
######################## Compare Set 1 with Set 2 ############################
ind_day3 = data_all2.columns.get_loc('% Alive on Day3')
ind_day50 = data_all2.columns.get_loc('% Alive on Day50')

# Figure out stats for each day
set1_datalists = [[] for _ in range(len(days))]
for i in range(ind_day3,ind_day50+1):
    s = pd.Series(data_all2.iloc[5,i]).describe() # Change index here
    set1_datalists[i-ind_day3] = [s['count'],round(s['mean'],2),s['25%'],s['50%'],s['75%'],s['std']]
    
# Figure out stats for each day
set2_datalists = [[] for _ in range(len(days))]
for i in range(ind_day3,ind_day50+1):
    s = pd.Series(data_all2.iloc[6,i]).describe() # Change index here
    set2_datalists[i-ind_day3] = [s['count'],round(s['mean'],2),s['25%'],s['50%'],s['75%'],s['std']]

#Transpose list of lists
set1_datalists = list(map(list, zip(*set1_datalists))) 
set2_datalists = list(map(list, zip(*set2_datalists))) 

# Turn list of lists into dataframe
df_set1 = pd.DataFrame(set1_datalists,index=['count','mean','25%','median','75%','std'], 
                          columns=['Day3','Day5','Day10','Day15','Day20',
                                   'Day25','Day30','Day40','Day50']).transpose()

df_set2 = pd.DataFrame(set2_datalists,index=['count','mean','25%','median','75%','std'], 
                          columns=['Day3','Day5','Day10','Day15','Day20',
                                   'Day25','Day30','Day40','Day50']).transpose()


fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2,figsize=(18,14))
fig.suptitle('Set 1 vs Set 2', fontsize=22, y=1.01)
color_cycle = ['g', 'b', 'c']
ind = np.arange(len(days))+1

ax1.plot(ind,df_set1['median'], color='k',label='median (After L4)') # Change legend here
ax1.plot(ind,df_set2['median'], color='r',label='median (Every other day)') # Change legend here
ax1.fill_between(ind, df_set1['median'], df_set1['75%'],alpha=0.3, 
                 facecolor=color_cycle[0], label='75%')
ax1.fill_between(ind, df_set1['25%'], df_set1['median'],alpha=0.3, 
                 facecolor=color_cycle[1], label='25%')
                 
ax2.plot(ind,df_set1['mean'], color='k',label='mean (After L4)') # Change legend here
ax2.plot(ind,df_set2['mean'], color='r',label='mean (Every other day)') # Change legend here
ax2.fill_between(ind, df_set1['mean']-2*df_set1['std'], df_set1['mean']+2*df_set1['std'], alpha=0.3, 
                 facecolor=color_cycle[0], label='2x std')
ax2.fill_between(ind, df_set1['mean']-df_set1['std'], df_set1['mean']+df_set1['std'],alpha=0.3, 
                 facecolor=color_cycle[1], label='1x std')

ax3.plot(ind,df_set2['median'], color='k',label='median')
ax3.fill_between(ind, df_set2['median'], df_set2['75%'],alpha=0.3, 
                 facecolor=color_cycle[0], label='75%')
ax3.fill_between(ind, df_set2['25%'], df_set2['median'],alpha=0.3, 
                 facecolor=color_cycle[1], label='25%')

ax4.plot(ind,df_set2['mean'], color='k',label='mean')
ax4.fill_between(ind, df_set2['mean']-2*df_set2['std'], df_set2['mean']+2*df_set2['std'],alpha=0.3, 
                 facecolor=color_cycle[0], label='2x std')
ax4.fill_between(ind, df_set2['mean']-df_set2['std'], df_set2['mean']+df_set2['std'],alpha=0.3, 
                 facecolor=color_cycle[1], label='1x std')

ax1.set_ylim(-5,105)
ax2.set_ylim(-45,130)
ax3.set_ylim(-5,105)
ax4.set_ylim(-40,125)

ax1.set_xticks(ind)
ax2.set_xticks(ind)
ax3.set_xticks(ind)
ax4.set_xticks(ind)

ax1.set_xticklabels(df_set1.index, fontsize=12, rotation = 70)
ax2.set_xticklabels(df_set1.index, fontsize=12, rotation = 70)
ax3.set_xticklabels(df_set1.index, fontsize=12, rotation = 70)
ax4.set_xticklabels(df_set1.index, fontsize=12, rotation = 70)

ax1.set_yticks(np.arange(0,120,20))
ax2.set_yticks(np.arange(-20,140,20))
ax3.set_yticks(np.arange(0,120,20))
ax4.set_yticks(np.arange(-20,140,20))

ax1.set_yticklabels(np.arange(0,120,20), fontsize=12)
ax2.set_yticklabels(np.arange(-20,140,20), fontsize=12)
ax3.set_yticklabels(np.arange(0,120,20), fontsize=12)
ax4.set_yticklabels(np.arange(-20,140,20), fontsize=12)

ax1.set_title('After L4 (Median)', fontsize=14,x=0.5) # Change titles here
ax2.set_title('After L4 (Mean)', fontsize=14,x=0.5) # Change titles here
ax3.set_title('Every other day (Median)', fontsize=14,x=0.5) # Change titles here
ax4.set_title('Every other day (Median)', fontsize=14,x=0.5) # Change titles here
ax1.set_ylabel('% Alive',fontsize=20)
ax3.set_ylabel('% Alive',fontsize=20)
ax1.legend(fontsize=14)
ax2.legend(fontsize=14)
ax3.legend(fontsize=14)
ax4.legend(fontsize=14)
plt.tight_layout()
plt.show()

#%%################## Plotting FUDR concentrations in Set 2 ##################
fig = plt.figure(figsize=(10,10))
ind_FUDRconc = data_all.columns.get_loc('FUDR Concentration (microM)')
x = range(len(data_all.iloc[(1,ind_FUDRconc)]))
plt.scatter(x,data_all.iloc[(1,ind_FUDRconc)])
plt.ylim(-10,220)
plt.title('FUDR Concentrations in Set 2')
plt.ylabel('Concentration (microM)')
plt.xlabel('Entry Number in Set 2')
plt.show()

# Plotting FUDR concentrations in Set 2 for entries with same concentration #
fig = plt.figure(figsize=(10,10))
ind_FUDRconc = data_all.columns.get_loc('FUDR Concentration (microM)')

# Find unique values of FUDR concentrations
FUDRvalues, FUDRfreq = np.unique(data_all.iloc[(1,ind_FUDRconc)], return_counts=True)

# Find which values have > 4 entries with same FUDR concentration
FUDRtest = np.where(FUDRfreq>4)

FUDRdata = {'Values':FUDRvalues[FUDRtest],'Frequency':FUDRfreq[FUDRtest]}
df_FUDR = pd.DataFrame(FUDRdata)

# Find and delete missing FUDR concentrations
df_FUDR = df_FUDR.drop(df_FUDR[df_FUDR['Values']==-1].index)

for i in range(len(df_FUDR)):
    plt.scatter(df_FUDR.iloc[i,0],df_FUDR.iloc[i,1],label = df_FUDR.iloc[i,0])

plt.ylim(0,30,5)
plt.title('FUDR Concentrations in Set 2')
plt.ylabel('Frequency')
plt.xlabel('Concentration (microM)')
plt.legend(bbox_to_anchor=(0.99, 1), loc='upper right', borderaxespad=0.2,
           fontsize = 10, title = 'Unit = microM')
plt.show()

#%%################## Comparing FUDR Concentrations ##########################
L_FUDR = data_all.iloc[(1,ind_count)]
listFUDR = data_all.iloc[(1,ind_FUDRconc)] #Create list of FUDR concentrations for entries
ind_FUDR40 = [i for i, x in enumerate(listFUDR) if x == 40.0] # Find where conc = 40
ind_FUDR200 = [i for i, x in enumerate(listFUDR) if x == 200.0] # Find where conc = 200

ind_day3 = data_all.columns.get_loc('% Alive on Day3')
ind_day50 = data_all.columns.get_loc('% Alive on Day50')

# Find and save % Alive values for FUDR = 40 entries for day 3 - 50
FUDR40_lists = [[] for _ in range(len(days))]
for i in range(ind_day3,ind_day50+1):
    for j in ind_FUDR40:
        FUDR40_lists[i-ind_day3].append(data_all.iloc[(1,i)][j])
             
# Find and save % Alive values for FUDR = 200 entries for day 3 - 50
FUDR200_lists = [[] for _ in range(len(days))]
for i in range(ind_day3,ind_day50+1):
    for j in ind_FUDR200:
        FUDR200_lists[i-ind_day3].append(data_all.iloc[(1,i)][j])

# Turn them into DataFrames
df_totFUDR40 = pd.DataFrame(FUDR40_lists, index = days)
df_totFUDR200 = pd.DataFrame(FUDR200_lists, index = days)

# Perform t-test % Alive for each day
ttest_FUDR = []
for i in range(len(days)):
    test = stats.ttest_ind(df_totFUDR40.iloc[i,:], df_totFUDR200.iloc[i,:], equal_var=False)
    ttest_FUDR.append([test[0],test[1]])
    
df_ttestFUDR = pd.DataFrame(ttest_FUDR,columns=['t statistic','p-value'],
                            index=['Day3','Day5','Day10','Day15','Day20',
                                     'Day25','Day30','Day40','Day50'])
    
# Figure out stats for each day
FUDR40_datalists = [[] for _ in range(len(days))]
for i in range(len(df_totFUDR40)):
    s = df_totFUDR40.iloc[i,:].describe()
    FUDR40_datalists[i] = [s['count'],round(s['mean'],2),s['25%'],s['50%'],s['75%'],s['std']]
    
# Figure out stats for each day
FUDR200_datalists = [[] for _ in range(len(days))]
for i in range(len(df_totFUDR200)):
    s = df_totFUDR200.iloc[i,:].describe()
    FUDR200_datalists[i] = [s['count'],round(s['mean'],2),s['25%'],s['50%'],s['75%'],s['std']]

#Transpose list of lists
FUDR40_datalists = list(map(list, zip(*FUDR40_datalists))) 
FUDR200_datalists = list(map(list, zip(*FUDR200_datalists))) 

# Turn list of lists into dataframe
df_FUDR40 = pd.DataFrame(FUDR40_datalists,index=['count','mean','25%','median','75%','std'], 
                          columns=['Day3','Day5','Day10','Day15','Day20',
                                   'Day25','Day30','Day40','Day50']).transpose()

df_FUDR200 = pd.DataFrame(FUDR200_datalists,index=['count','mean','25%','median','75%','std'], 
                          columns=['Day3','Day5','Day10','Day15','Day20',
                                   'Day25','Day30','Day40','Day50']).transpose()


fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2,figsize=(18,14))
fig.suptitle('FUDR Concentration Comparison', fontsize=22, y=1.01)
color_cycle = ['g', 'b', 'c']
ind = np.arange(len(days))+1

ax1.plot(ind,df_FUDR40['median'], color='k',label='median (40µM)')
ax1.plot(ind,df_FUDR200['median'], color='r',label='median (200µM)')
ax1.fill_between(ind, df_FUDR40['median'], df_FUDR40['75%'],alpha=0.3, 
                 facecolor=color_cycle[0], label='75%')
ax1.fill_between(ind, df_FUDR40['25%'], df_FUDR40['median'],alpha=0.3, 
                 facecolor=color_cycle[1], label='25%')
                 
ax2.plot(ind,df_FUDR40['mean'], color='k',label='mean (40µM)')
ax2.plot(ind,df_FUDR200['mean'], color='r',label='mean (200µM)')
ax2.fill_between(ind, df_FUDR40['mean']-2*df_FUDR40['std'], df_FUDR40['mean']+2*df_FUDR40['std'], alpha=0.3, 
                 facecolor=color_cycle[0], label='2x std')
ax2.fill_between(ind, df_FUDR40['mean']-df_FUDR40['std'], df_FUDR40['mean']+df_FUDR40['std'],alpha=0.3, 
                 facecolor=color_cycle[1], label='1x std')

ax3.plot(ind,df_FUDR200['median'], color='k',label='median')
ax3.fill_between(ind, df_FUDR200['median'], df_FUDR200['75%'],alpha=0.3, 
                 facecolor=color_cycle[0], label='75%')
ax3.fill_between(ind, df_FUDR200['25%'], df_FUDR200['median'],alpha=0.3, 
                 facecolor=color_cycle[1], label='25%')

ax4.plot(ind,df_FUDR200['mean'], color='k',label='mean')
ax4.fill_between(ind, df_FUDR200['mean']-2*df_FUDR200['std'], df_FUDR200['mean']+2*df_FUDR200['std'],alpha=0.3, 
                 facecolor=color_cycle[0], label='2x std')
ax4.fill_between(ind, df_FUDR200['mean']-df_FUDR200['std'], df_FUDR200['mean']+df_FUDR200['std'],alpha=0.3, 
                 facecolor=color_cycle[1], label='1x std')

ax1.set_ylim(-5,105)
ax2.set_ylim(-45,130)
ax3.set_ylim(-5,105)
ax4.set_ylim(-40,125)

ax1.set_xticks(ind)
ax2.set_xticks(ind)
ax3.set_xticks(ind)
ax4.set_xticks(ind)

ax1.set_xticklabels(df_FUDR40.index, fontsize=12, rotation = 70)
ax2.set_xticklabels(df_FUDR40.index, fontsize=12, rotation = 70)
ax3.set_xticklabels(df_FUDR40.index, fontsize=12, rotation = 70)
ax4.set_xticklabels(df_FUDR40.index, fontsize=12, rotation = 70)

ax1.set_yticks(np.arange(0,120,20))
ax2.set_yticks(np.arange(-20,140,20))
ax3.set_yticks(np.arange(0,120,20))
ax4.set_yticks(np.arange(-20,140,20))

ax1.set_yticklabels(np.arange(0,120,20), fontsize=12)
ax2.set_yticklabels(np.arange(-20,140,20), fontsize=12)
ax3.set_yticklabels(np.arange(0,120,20), fontsize=12)
ax4.set_yticklabels(np.arange(-20,140,20), fontsize=12)

ax1.set_title('40 µM', fontsize=14,x=0.5)
ax2.set_title('40 µM', fontsize=14,x=0.5)
ax3.set_title('200 µM', fontsize=14,x=0.5)
ax4.set_title('200 µM', fontsize=14,x=0.5)
ax1.set_ylabel('% Alive',fontsize=20)
ax3.set_ylabel('% Alive',fontsize=20)
ax1.legend(fontsize=14)
ax2.legend(fontsize=14)
ax3.legend(fontsize=14)
ax4.legend(fontsize=14)
plt.tight_layout()
plt.show()

#%%##################### Day 20 % Alive Bar Graph ############################
daystoplot = ['% Alive on Day3','% Alive on Day5','% Alive on Day10','% Alive on Day15',
              '% Alive on Day20','% Alive on Day25']
L = len(data_all)
L2 = len(daystoplot)
    
for ii in range(L2):
    df_alive,lifespan_lists = plotlifespanspread(data_all, daystoplot[ii])
    ind = np.arange(L)+1  # the x locations for the groups
    width = 0.5
    
    fig1 = plt.figure(figsize=(14,10))
    p1 = plt.bar(ind, df_alive['20%'], width, color='b')
    p2 = plt.bar(ind, df_alive['40%'], width, bottom=df_alive['20%'],color='r')
    p3 = plt.bar(ind, df_alive['60%'], width, color='g',\
                 bottom=[i+j for i,j in zip(df_alive['20%'],df_alive['40%'])])
    p4 = plt.bar(ind, df_alive['80%'], width, color='k',\
                 bottom=[i+j+k for i,j,k in zip(df_alive['20%'],df_alive['40%'],df_alive['60%'])])
    p5 = plt.bar(ind, df_alive['100%'], width, color='y',\
                 bottom=[i+j+k+z for i,j,k,z in zip(df_alive['20%'],df_alive['40%'],df_alive['60%'],df_alive['80%'])])
    
    plt.ylabel('Percent of Total Entries In Set',fontsize=16)
    plt.xlabel('Set Number',fontsize=16)
    plt.xticks(ind, fontsize=14, labels=df_alive.index)
    plt.yticks(np.arange(0, 110, 10),size=14)
    plt.title('Standardized Spread of {} for Each Set'.format(daystoplot[ii]),\
              fontsize=24,y=1.02)
    plt.legend((p1[0], p2[0], p3[0], p4[0], p5[0]), ('Entries with <20% Alive',\
               'Entries with 20-40% Alive','Entries with 40-60% Alive',\
               'Entries with 60-80% Alive','Entries with 80-100% Alive'),\
               bbox_to_anchor=(1.0, 1),loc='upper left',\
               borderaxespad=0.5, fontsize=14)
    plt.show()
    
    if daystoplot[ii] == '% Alive on Day20':
        lifespan20_list = lifespan_lists.copy()

test = stats.ttest_ind(lifespan20_list[0], lifespan20_list[1], equal_var=False)
print(test)

N = 2 # Number of Sets to Analyze
meanls_missing = [0]*N # One value per set analyzed

#################### Mean lifespan Histogram for each set ####################
for k in range(N):
    ind_meanls = data_all.columns.get_loc('Reported mean lifespan')
    
    meanls = [round(num) for num in data_all.iloc[k,ind_meanls]]
    group,freq_group = np.unique(meanls,return_counts=True)
    unique_group = pd.DataFrame(list(zip(group,freq_group)),columns=['Mean Lifespan','Frequency'])
    unique_group = unique_group.reindex(np.argsort(unique_group['Mean Lifespan'])[::-1])
    unique_group = unique_group.reset_index(drop=True)  
    meanls_missing[k] = unique_group.iloc[-1,1]
    set_stats = pd.Series(unique_group.iloc[:-1,0]).describe()
    
    fig = plt.figure(figsize=(14,10))
    plt.scatter(round(unique_group.iloc[:-1,0]),unique_group.iloc[:-1,1],
                label = '{} ± {}'.format(round(set_stats['mean'],2),round(set_stats['std'],2)))
    plt.ylabel('Frequency',fontsize=16)
    plt.xlabel('Mean Lifespan',fontsize=16)
    plt.xticks(np.arange(0,40,5), fontsize=14) # Don't want to include missing values
    plt.title('Histogram of Reported Mean Lifespan in Set {} ({} of {} missing)'
              .format(k+1,meanls_missing[k],data_all.iloc[k,ind_count]),fontsize=24,y=1.02)
    plt.legend(bbox_to_anchor=(0.99, 1), loc='upper right', borderaxespad=0.2,\
            fontsize = 10, title = 'Mean ± Std Dev')
    plt.show()

######################## % Alive on Day 20 by Country #########################
for k in range(N):
    set_data,unique_group,ind_group=plotday20(data_all,'Lab Location (Country)', k)

    ind = np.arange(len(unique_group))    # the x locations for the groups    
    fig = plt.figure(figsize=(14,8))
    stats_list = [[] for _ in range(len(ind))]
    
    for i in ind:
        x = [i]*unique_group.iloc[i,1]
        if i < 3:
            set_stats = pd.Series(set_data[i][0]).describe()
            stats_list[i] = [round(set_stats['mean'],2),round(set_stats['std'],2)]
            plt.scatter(x,set_data[i],label = '{} ± {}'.format(stats_list[i][0],stats_list[i][1]))
        else:
            plt.scatter(x,set_data[i])

    plt.ylabel('% Alive on Day 20',fontsize=16)
    plt.xlabel('Countries in Set {}'.format(k+1),fontsize=16)
    plt.xticks(ind,list(unique_group.iloc[:,0]), fontsize=14,rotation = 70)
    plt.yticks(np.arange(0, 110, 10),size=14)
    plt.title('% Alive on Day 20 for Each Country in Set {}'.format(k+1),fontsize=24,y=1.02)
    plt.legend(bbox_to_anchor=(0.99, 1), loc='upper right', borderaxespad=0.2,\
            fontsize = 14, title = 'Mean ± Std Dev')
    plt.show()
   
###################### Reported Mean Lifespan by Country ######################
for k in range(N):
    set_data,unique_group,ind_group=plotmeanlifespan(data_all,'Lab Location (Country)', k)

    ind = np.arange(len(unique_group))    # the x locations for the groups    
    fig = plt.figure(figsize=(14,8))
    stats_list = [[] for _ in range(len(ind))]
    
    for i in ind:
        x = [i]*unique_group.iloc[i,1]
        if i < 3:
            set_stats = pd.Series(set_data[i][0]).describe()
            stats_list[i] = [round(set_stats['mean'],2),round(set_stats['std'],2)]
            plt.scatter(x,set_data[i],label = '{} ± {}'.format(stats_list[i][0],stats_list[i][1]))
        else:
            plt.scatter(x,set_data[i])

    plt.ylabel('Reported Mean Lifespan (days)',fontsize=16)
    plt.xlabel('Countries in Set {}'.format(k+1),fontsize=16)
    plt.xticks(ind,list(unique_group.iloc[:,0]), fontsize=14,rotation = 70)
    plt.yticks(np.arange(0,40,5),size=14)
    plt.title('Reported Mean Lifespan for Each Country in Set {} ({} of {} missing)'
              .format(k+1,meanls_missing[k],data_all.iloc[k,ind_count]),fontsize=24,y=1.02)
    plt.legend(bbox_to_anchor=(0.99, 1), loc='upper right', borderaxespad=0.2,\
            fontsize = 14, title = 'Mean ± Std Dev')
    plt.show()

#%%########### Analysis by Country after combining set 1 and 2 ###############
ind_3 = data_all.columns.get_loc("% Alive on Day3")
ind_50 = data_all.columns.get_loc("% Alive on Day50")
df_cntry = data_all['Lab Location (Country)']
data_all.drop(labels=['Lab Location (Country)'], axis=1,inplace = True)
data_all.insert(ind_50+1, 'Lab Location (Country)', df_cntry)
ind_cntry = data_all.columns.get_loc("Lab Location (Country)")

indices = list(range(ind_3,ind_cntry+1))
data_geo = [[] for _ in range(len(indices))]
for i in indices:
    for j in range(2):
        if j==0:
            data_geo[i-ind_3] = data_all.iloc[(j,i)]
        elif j==1:
            data_geo[i-ind_3].extend(data_all.iloc[(j,i)])

df_datageo = pd.DataFrame(data_geo,index=['Day3','Day5','Day10','Day15','Day20',
                                            'Day25','Day30','Day40','Day50',
                                            'Lab Location (Country)']).transpose()

geo,freq_geo = np.unique(df_datageo['Lab Location (Country)'],return_counts=True)
unique_geo = pd.DataFrame(list(zip(geo,freq_geo)),columns=['Names','Frequency'])
unique_geo = unique_geo.reindex(np.argsort(unique_geo['Frequency'])[::-1])
unique_geo = unique_geo.reset_index(drop=True)

for i in range(len(unique_geo)):
    index = np.where(np.isin(df_datageo.iloc[:,-1],unique_geo.iloc[(i,0)]))
    geo = df_datageo.iloc[index[0],:]
    geo = geo.drop(['Lab Location (Country)'],axis=1)
    geo = geo.astype('float64')
    
    # Figure out stats for each day
    geo2 = [[] for _ in range(len(days))]
    for j in range(len(days)):
        s = geo.iloc[:,j].describe()
        geo2[j] = [s['count'],round(s['mean'],2),s['25%'],s['50%'],s['75%'],s['std']]
    
    #Transpose list of lists
    geo2 = list(map(list, zip(*geo2))) 

    # Turn list of lists into dataframe
    df_geo2 = pd.DataFrame(geo2,index=['count','mean','25%','median','75%','std'], 
                          columns=['Day3','Day5','Day10','Day15','Day20',
                                   'Day25','Day30','Day40','Day50']).transpose()
    
    if df_geo2['count'][0]>13:
        fig, (ax1, ax2) = plt.subplots(1,2,figsize=(14,8),sharey=True)
        fig.suptitle('Lifespan Data for Entries in {} (n = {})'.format(unique_geo.iloc[(i,0)],int(df_geo2['count'][0])), 
                     fontsize=22, y=1.02)
        color_cycle = ['g', 'b','c']
        ind = np.arange(len(days))+1
    
        ax1.plot(ind,df_geo2['median'], color='k',label='median')
        ax1.fill_between(ind, df_geo2['median'], df_geo2['75%'],alpha=0.3,facecolor='g',label='75%')
        ax1.fill_between(ind, df_geo2['25%'], df_geo2['median'],alpha=0.3,facecolor='b',label='25%')
                     
        ax2.plot(ind,df_geo2['mean'], color='k',label='mean (40µM)')
        ax2.fill_between(ind, df_geo2['mean']-2*df_geo2['std'], 
                         df_geo2['mean']+2*df_geo2['std'], alpha=0.3, 
                         facecolor='g',label='2x std')
        ax2.fill_between(ind, df_geo2['mean']-df_geo2['std'], 
                         df_geo2['mean']+df_geo2['std'],alpha=0.3,
                         facecolor='b',label='1x std')
    
        ax1.set_ylim(-45,130)
        ax2.set_ylim(-45,130)
        
        ax1.set_xticks(ind)
        ax2.set_xticks(ind)
        
        ax1.set_xticklabels(df_FUDR40.index, fontsize=12, rotation = 70)
        ax2.set_xticklabels(df_FUDR40.index, fontsize=12, rotation = 70)
        
        ax1.set_yticks(np.arange(-20,140,20))
        ax2.set_yticks(np.arange(-20,140,20))
        
        ax1.set_yticklabels(np.arange(-20,140,20), fontsize=12)
        ax2.set_yticklabels(np.arange(-20,140,20), fontsize=12)
        
        ax1.set_title('Median + Quartiles', fontsize=14,x=0.5)
        ax2.set_title('Mean + Standard Deviations', fontsize=14,x=0.5)
        
        ax1.set_ylabel('% Alive',fontsize=20)
        ax1.legend(fontsize=14)
        ax2.legend(fontsize=14)
        
        plt.tight_layout()
        plt.show()

#%%#################### 50% Alive Histogram for each set #####################
for k in range(N):
    set_data,unique_group,ind_group=plotday20(data_all,'50% Alive', k)
    ind = np.arange(len(unique_group))    # the x locations for the groups    
    fig = plt.figure(figsize=(14,10))
    
    plt.bar(unique_group['Names'],unique_group['Frequency'])
    plt.ylabel('Frequency',fontsize=16)
    plt.xlabel('Day',fontsize=16)
    plt.xticks(unique_group['Names'],list(unique_group.iloc[:,0]), fontsize=14)
    plt.title('Last Day with ≥5% of Animals Alive in Set {}'.format(k+1),fontsize=24,y=1.02)
    plt.show()

#%%# Descriptive statistics for MM ####

stats_data_20,stats_data_25 = data_stats(data_all)

#### Distribution Set to Set ####
N = 3 # Number of sets to analyze
color_cycle = ['g', 'b', 'c', 'm', 'y', 'k']

for i in range(N):
    fig = plt.figure(figsize=(10,6))

    mean = round(stats_data_20.iloc[(i,1)],1)
    std = round(stats_data_20.iloc[(i,2)],1)
    
    values, frequency = np.unique(lifespan20_list[i], return_counts=True)
    plt.plot(values, frequency, label = '{} ± {}'.format(mean,std),color = color_cycle[i])
    
    mean_liney = np.arange(0,11,1)
    mean_linex = [mean]*11
    plt.plot(mean_linex,mean_liney,linestyle='--',linewidth=1,color = color_cycle[i])
    
    plt.xlim(0,100)
    plt.xticks(np.arange(0,110,10))
    plt.ylim(0,10)
    plt.yticks(np.arange(0,10,1))
    plt.ylabel('Frequency',fontsize = 16)
    plt.xlabel('% Alive on Day 20',fontsize = 16)
    plt.title('Histogram of % Alive on Day 20 for Set {}'.format(i+1),fontsize = 20)
    plt.legend(bbox_to_anchor=(0.99, 1), loc='upper right', borderaxespad=0.2,\
            fontsize = 14, title = 'Mean ± Std Dev')

## % Alive on Day 20 ###
maxminusmin_20 = stats_data_20["max"] - stats_data_20["min"]
maxovermin_20 = stats_data_20["max"]/stats_data_20["min"]

proc_minus_OK_20 = np.round(100*len(maxminusmin_20[maxminusmin_20<50])/L,2)
proc_minus_nOK_20  = 100 - proc_minus_OK_20

proc_over_OK_20 = np.round(100*len(maxovermin_20[maxovermin_20<=2])/L,2)
proc_over_nOK_20  = 100 - proc_over_OK_20

## % Alive on Day 25 ###
maxminusmin_25 = stats_data_25["max"] - stats_data_25["min"]
maxovermin_25 = stats_data_25["max"]/stats_data_25["min"]

proc_minus_OK_25 = np.round(100*len(maxminusmin_25[maxminusmin_25<50])/L,2)
proc_minus_nOK_25  = 100 - proc_minus_OK_25

proc_over_OK_25 = np.round(100*len(maxovermin_25[maxovermin_25<=2])/L,2)
proc_over_nOK_25  = 100 - proc_over_OK_25

