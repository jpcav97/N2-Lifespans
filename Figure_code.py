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
from PIL import Image
import io
import random
import statsmodels.api as sm
from sklearn.linear_model import LinearRegression
def get_ttfp(data,transf1,transf2,isnm,tempm,tempc,isFUDR1,isFUDR2):
    ttfpcolumns = ['% Alive on Day3','% Alive on Day5','% Alive on Day10','% Alive on Day15',
              '% Alive on Day20','% Alive on Day25','% Alive on Day30',
              '% Alive on Day40','% Alive on Day50','Reported mean lifespan',
              'FUDR Concentration (microM)','Lab Location (Country)','Lab Location (State/Province)']

    if isnm == True:
        df_nm = data.loc[(data['# Transfers to Fresh Plates'] == -1) & 
                  (data['Temperature Maintained (Through L4, C)'] == tempm) &
                  (data['Temperature Cultivated (Adult, C)'] == tempc)]
        df_m = data.loc[(data['# Transfers to Fresh Plates'] != -1) & 
                  (data['Temperature Maintained (Through L4, C)'] == tempm) &
                  (data['Temperature Cultivated (Adult, C)'] == tempc)]
    elif isFUDR1 == True:
        df_nm = data.loc[(data['# Transfers to Fresh Plates'] == -1) &
                  (data['FUDR (Yes/No?)'] == 'No') &
                  (data['Temperature Maintained (Through L4, C)'] == tempm) &
                  (data['Temperature Cultivated (Adult, C)'] == tempc)]
        df_m = data.loc[(data['# Transfers to Fresh Plates'] == -1) & 
                  (data['FUDR (Yes/No?)'] == 'Yes') &
                  (data['Temperature Maintained (Through L4, C)'] == tempm) &
                  (data['Temperature Cultivated (Adult, C)'] == tempc)]
    elif isFUDR2 == True:
        df_nm = data.loc[(data['# Transfers to Fresh Plates'] != -1) &
                  (data['FUDR (Yes/No?)'] == 'No') &
                  (data['Temperature Maintained (Through L4, C)'] == tempm) &
                  (data['Temperature Cultivated (Adult, C)'] == tempc)]
        df_m = data.loc[(data['# Transfers to Fresh Plates'] != -1) & 
                  (data['FUDR (Yes/No?)'] == 'Yes') &
                  (data['Temperature Maintained (Through L4, C)'] == tempm) &
                  (data['Temperature Cultivated (Adult, C)'] == tempc)]
    else:
        df_nm = data.loc[(data['# Transfers to Fresh Plates'] == transf1) & 
                  (data['Temperature Maintained (Through L4, C)'] == tempm) &
                  (data['Temperature Cultivated (Adult, C)'] == tempc)]
        df_m = data.loc[(data['# Transfers to Fresh Plates'] == transf2) & 
                  (data['Temperature Maintained (Through L4, C)'] == tempm) &
                  (data['Temperature Cultivated (Adult, C)'] == tempc)]
    
    ind_start = df_nm.columns.get_loc('% Alive on Day3')
    ind_end = df_nm.columns.get_loc('Lab Location (State/Province)')
    
    # Find and save % Alive values for FUDR = 40 entries for day 3 - 50 and rml
    ttfp_lists_nm = [[] for _ in range(len(ttfpcolumns[:-3]))]
    ttfp_lists_m = [[] for _ in range(len(ttfpcolumns[:-3]))]

    for i in range(ind_start,ind_end-2):
        for j in range(len(df_nm)):
            ttfp_lists_nm[i-ind_start].extend(df_nm.iloc[(j,i)])
        
        for j in range(len(df_m)):
            ttfp_lists_m[i-ind_start].extend(df_m.iloc[(j,i)])        
        
    
    df_totnm = pd.DataFrame(ttfp_lists_nm, index = ttfpcolumns[:-3]).transpose()
    df_totm = pd.DataFrame(ttfp_lists_m, index = ttfpcolumns[:-3]).transpose()
    
    # Find and save % Alive values for FUDR = 40 entries for day 3 - 50, rml, and for country analysis
    ttfp_lists_nm_all = [[] for _ in range(len(ttfpcolumns))]
    ttfp_lists_m_all = [[] for _ in range(len(ttfpcolumns))]

    for i in range(ind_start,ind_end+1):
        for j in range(len(df_nm)):
            ttfp_lists_nm_all[i-ind_start].extend(df_nm.iloc[(j,i)])
        
        for j in range(len(df_m)):
            ttfp_lists_m_all[i-ind_start].extend(df_m.iloc[(j,i)])
    
    df_totnm_all = pd.DataFrame(ttfp_lists_nm_all, index = ttfpcolumns).transpose()
    df_totm_all = pd.DataFrame(ttfp_lists_m_all, index = ttfpcolumns).transpose()
    
    # Figure out stats for each day
    Ldays = len(ttfpcolumns)-4
    datalists_nm = [[] for _ in range(Ldays)]
    datalists_m = [[] for _ in range(Ldays)]

    for i in range(Ldays):
        s1 = df_totnm.iloc[:,i].describe()
        s2 = df_totm.iloc[:,i].describe()        
        datalists_nm[i] = [s1['count'],round(s1['mean'],2),s1['25%'],s1['50%'],s1['75%'],s1['std']]
        datalists_m[i] = [s2['count'],round(s2['mean'],2),s2['25%'],s2['50%'],s2['75%'],s2['std']]

    datalists_nm = list(map(list,zip(*datalists_nm)))
    datalists_m = list(map(list,zip(*datalists_m)))
    
    # Turn list of lists into dataframe
    data_nm = pd.DataFrame(datalists_nm,index=['count','mean','25%','median','75%','std'], 
                              columns=['Day3','Day5','Day10','Day15','Day20',
                                       'Day25','Day30','Day40','Day50']).transpose()
    data_m = pd.DataFrame(datalists_m,index=['count','mean','25%','median','75%','std'], 
                              columns=['Day3','Day5','Day10','Day15','Day20',
                                       'Day25','Day30','Day40','Day50']).transpose()
    del ttfp_lists_nm,ttfp_lists_m,ttfp_lists_nm_all,ttfp_lists_m_all
    return df_totnm,df_totm,data_nm,data_m,df_totnm_all,df_totm_all

def randomvals_and_diff(data,ind_rml,lab):
    L = int(len(data)/2)
    alive_list1 = [0]*L
    alive_list2 = [0]*L
    
    if lab == True:
        # Pick random value for Labs with more than one value 
        df_fix = data[data['Count']>1]
        inds = list(df_fix.index)
        for i in range(len(df_fix)):
            df_fix.iloc[i,ind_rml] = random.sample(data.iloc[i,ind_rml],1)
        
        data.iloc[inds] = df_fix.iloc[inds]
        
        df_fix2 = data[data['Count']==1]
        inds2 = list(df_fix2.index)
        for i in range(len(inds2)):
            data.iloc[inds2[i],ind_rml] = df_fix2.iloc[i,ind_rml][0]
    
    
    j=0
    while j < L:
        temp = random.sample(list(data.iloc[:,ind_rml]),2)
        
        if (temp[0] != -2 and temp[1] != -2):
            alive_list1[j] = temp[0]
            alive_list2[j] = temp[1]
            j = j+1
    
    df_pairs = pd.DataFrame({'value 1': alive_list1, 'value 2': alive_list2})
        
    return df_pairs

def residuals(data):
    # Calculate standardized residuals and find outliers
    val1 = np.array(data['value 1']).reshape((-1, 1))
    val2 = np.array(data['value 2'])
    
    linreg = LinearRegression(fit_intercept=False)
    linreg = linreg.fit(val1, val2)
    
    rsquared = linreg.score(val1,val2)
    
    val1 = sm.add_constant(val1)
    #val2 = sm.add_constant(val2)
    
    #fit linear regression model
    model = sm.OLS(val2, val1).fit()
    
    #create instance of influence
    influence = model.get_influence()
    
    #obtain standardized residuals
    standardized_residuals = influence.resid_studentized_internal
    
    #display standardized residuals
    idx_std_res = np.where(standardized_residuals > 3)[0]
    idx_std_res = np.append(idx_std_res,np.where(standardized_residuals < -3)[0])
    print('The number of outliers = {}'.format(len(idx_std_res)))
    
    num_out = len(idx_std_res)
    #rsquared = model.rsquared
    
    return idx_std_res,num_out,rsquared,linreg

def pairplot(data,rsquared,linreg):
    plt.figure()
    plt.scatter(data['value 1'],data['value 2'],color='k',alpha=0.3,edgecolor=None)
    
    # Fit a line to data with np.polyfit
    x = np.linspace(0,40,1000)
    m = linreg.coef_
    plt.plot(x, m*x,'g',label = 'R^2 = {}'.format(round(rsquared,2)))
    print('R^2 value = {}'.format(round(rsquared,3)))
    
    plt.xticks(range(0,45,5), fontsize=12)
    plt.yticks(range(0,45,5), fontsize=12)
    plt.ylim(0,40)
    plt.xlim(0,40)
    plt.xlabel('Mean Lifespan (days)',fontsize=12)
    plt.ylabel('Mean Lifespan (days)',fontsize=12)
    plt.legend()
    plt.show()

    return m

def make_groups(data,entries,P): #IsBuff,IsLogTransf:
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
    state = []
    PM_ID = []
    lab = []
    
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
        state.append(list(example['Lab Location (State/Province)']))
        PM_ID.append(list(example['PM ID']))
        lab.append(list(example['Lab']))
    
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
    grp['Lab Location (State/Province)'] = pd.Series(state)
    grp['PM ID'] = pd.Series(PM_ID)
    grp['Lab'] = pd.Series(lab)
    
    grp  = grp.reindex(np.argsort(grp['Count'])[::-1])
    grp  = grp.reset_index(drop=True)
    dat = grp[grp["Count"]>=P]
    
    return dat

def plotlifespanspread(data, column):
    L = len(data)
    ind_count = data.columns.get_loc('Count')
    
    lifespan_lists = [[] for _ in range(L)]
    perc_alive_lists = [[] for _ in range(6)]
    
    ### Create bar graph ###
    ind_day = data.columns.get_loc(column)
    
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
        perc_alive_lists[5].append(data.iloc[i,ind_count])
    
    
    df_a_lists = pd.DataFrame(perc_alive_lists,index=['20%','40%','60%','80%','100%','count'])
    #df_a_lists = df_a_lists.sort_values(by=['25%'],ascending=False)
    
    return df_a_lists, lifespan_lists

def categorize(data, column_name, N, daytoplot):
    L = len(data)

    ind_day = data.columns.get_loc(daytoplot)
    ind_group = data.columns.get_loc(column_name)
    ind_avg = data.columns.get_loc('Reported mean lifespan')


    lists = [[] for _ in range(L)]
    avg_lists = [[] for _ in range(L)]
    lifespan_lists = [[] for _ in range(L)]
            
    for i in range(L):
        lifespan_lists[i] = data.iloc[(i,ind_day)]
        avg_lists[i] = data.iloc[(i,ind_avg)]
        lists[i] = data.iloc[(i,ind_group)]
    
    group,freq_group = np.unique(data.iloc[(N,ind_group)],return_counts=True)
    unique_group=pd.DataFrame(list(zip(group,freq_group)),columns=['Names','Frequency'])
    unique_group = unique_group.reindex(np.argsort(unique_group['Frequency'])[::-1])
    unique_group  = unique_group.reset_index(drop=True)
    
    lists_set = [[] for _ in range(len(group))]
    avg_set = [[] for _ in range(len(group))]
    df_lifespan = pd.Series(lifespan_lists[N])
    df_avg = pd.Series(avg_lists[N])
    for i in range(len(group)):
        index = np.where(np.isin(lists[N],unique_group.iloc[(i,0)]))
        d20 = df_lifespan.iloc[(index)]
        lists_set[i] = d20
        
        a20 = df_avg.iloc[(index)]
        avg_set[i] = a20
    
    return lists_set, unique_group, avg_set

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
filename = 'N2 Lifespans FINAL.xlsx'
path = os.path.join(path, filename)
data = pd.read_excel(path)
data = data.drop(columns=data.columns[0])

#### Percent of data excluded for different categories ####
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
#%%
"""
                         FIGURE 1 
                                                             """
#Figure 1A
days = ['% Alive on Day3','% Alive on Day5','% Alive on Day10','% Alive on Day15',
              '% Alive on Day20','% Alive on Day25','% Alive on Day30',
              '% Alive on Day40','% Alive on Day50']
lists_s = [[] for _ in range(len(days))]

# Take stats for each day and save into list of lists
for i in range(len(days)):
    s = data_new[days[i]].describe()
    lists_s[i] = [s['count'],round(s['mean'],2),s['25%'],s['50%'],s['75%'],s['std']]

#Transpose list of lists
lists_s = list(map(list, zip(*lists_s))) 

# Turn list of lists into dataframe
df_totdata = pd.DataFrame(lists_s,index=['count','mean','25%','median','75%','std'], 
                          columns=['Day3','Day5','Day10','Day15','Day20',
                                   'Day25','Day30','Day40','Day50']).transpose()

fig = plt.figure(figsize=(12,6))
color_cycle = ["black", "dodgerblue", "chocolate"] 
ind = [3,5,10,15,20,25,30,40,50]

plt.plot(ind,df_totdata['mean'], color=color_cycle[0], label='mean')
plt.fill_between(ind, df_totdata['mean']-2*df_totdata['std'], df_totdata['mean']+2*df_totdata['std'],
                 facecolor=color_cycle[2], label='2x std')
plt.fill_between(ind, df_totdata['mean']-df_totdata['std'], df_totdata['mean']+df_totdata['std'],
                 facecolor=color_cycle[1], label='1x std')

# plt.plot(ind,df_totdata['median'], color='k',label='Median') # Change legend here
# plt.fill_between(ind, df_totdata['median'], df_totdata['75%'],alpha=0.3, 
#                  facecolor=color_cycle[0], label='75%')
# plt.fill_between(ind, df_totdata['25%'], df_totdata['median'],alpha=0.3, 
#                  facecolor=color_cycle[1], label='25%')

plt.ylabel('% Survival',fontsize=16)
plt.xticks(ind,df_totdata.index, fontsize=14, rotation = 70)
plt.yticks(np.arange(-40,160,20), fontsize=14)
plt.ylim(-40,140)
plt.title('Entire Dataset (n = {})'.format(int(df_totdata['count'][0])),fontsize=20,y=1.02)
plt.legend(bbox_to_anchor=(0.99, 1), loc='upper right', borderaxespad=0.2,\
            fontsize = 10)
plt.tight_layout() ######################### FIX 1 ###########################
plt.show()

#df_totdata.to_csv('/Users/nickurban/Desktop/N2 Lifespan Variability/Figures/totaldatasetstats.csv')
df_totdata.to_csv('figures/Fig1A.csv')

png1 = io.BytesIO()
fig.savefig(png1, format="png")

# # Load this image into PIL
png2 = Image.open(png1)

# # Save as TIFF
png2.save('figures/Figure1A.tif')
png1.close()
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

#%%Figure 1B, C, and D
ind_day3 = data_all.columns.get_loc('% Alive on Day3')
ind_day50 = data_all.columns.get_loc('% Alive on Day50')
ind_tempc = data_all.columns.get_loc('Temperature Cultivated (Adult, C)')

# Figure out stats for each day
set1_datalists = [[] for _ in range(len(days))]
set2_datalists = [[] for _ in range(len(days))]
set3_datalists = [[] for _ in range(len(days))]

temp_ind = [6,0,1] # 0 = 20, 1 = 25, 6 = 15 ('index' = 'temp')
for j in temp_ind:
    # Turn list of lists into dataframe
    if j == temp_ind[0]:
        for i in range(ind_day3,ind_day50+1):
            s = pd.Series(data_all.iloc[j,i]).describe()
            set1_datalists[i-ind_day3] = [s['count'],round(s['mean'],2),s['25%'],s['50%'],s['75%'],s['std']]
        
        #Transpose list of lists
        set1_datalists = list(map(list, zip(*set1_datalists))) 
        df_set1 = pd.DataFrame(set1_datalists,index=['count','mean','25%','median','75%','std'], 
                              columns=['Day3','Day5','Day10','Day15','Day20',
                                       'Day25','Day30','Day40','Day50']).transpose()
    if j == temp_ind[1]:
        for i in range(ind_day3,ind_day50+1):
            s = pd.Series(data_all.iloc[j,i]).describe() # 0 = 20, 1 = 25, 6 = 15 ('index' = 'temp')
            set2_datalists[i-ind_day3] = [s['count'],round(s['mean'],2),s['25%'],s['50%'],s['75%'],s['std']]
        
        #Transpose list of lists
        set2_datalists = list(map(list, zip(*set2_datalists))) 
        df_set2 = pd.DataFrame(set2_datalists,index=['count','mean','25%','median','75%','std'], 
                              columns=['Day3','Day5','Day10','Day15','Day20',
                                       'Day25','Day30','Day40','Day50']).transpose()
        
    if j == temp_ind[2]:
        for i in range(ind_day3,ind_day50+1):
            s = pd.Series(data_all.iloc[j,i]).describe() # 0 = 20, 1 = 25, 6 = 15 ('index' = 'temp')
            set3_datalists[i-ind_day3] = [s['count'],round(s['mean'],2),s['25%'],s['50%'],s['75%'],s['std']]
        
        #Transpose list of lists
        set3_datalists = list(map(list, zip(*set3_datalists))) 
        df_set3 = pd.DataFrame(set3_datalists,index=['count','mean','25%','median','75%','std'], 
                              columns=['Day3','Day5','Day10','Day15','Day20',
                                       'Day25','Day30','Day40','Day50']).transpose()
 
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

    # Save the image
    png1 = io.BytesIO()
    fig.savefig(png1, format="png")
    
    # # Load this image into PIL
    png2 = Image.open(png1)
    print('n = {}'.format(sets[j]['count'][0]))

    if j == 0:
        # Save as TIFF
        png2.save('figures/Figure1B.tif')
        png1.close() 
        sets[j].to_csv('figures/Fig1B.csv')

    if j == 1:
    #     # Save as TIFF
        png2.save('figures/Figure1C.tif')
        png1.close()
        sets[j].to_csv('figures/Fig1C.csv')

    if j == 6:
        # Save as TIFF
        png2.save('figures/Figure1D.tif')
        png1.close()
        sets[j].to_csv('figures/Fig1D.csv')


#%% Figure 1E 
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
pd.DataFrame(mls_temp2,index=cols).transpose().to_csv('figures/New_fig_mls_data.csv')

### Vertical Spread Graphs for new fig ###
""" 

daystoplot = ['% Alive on Day15','% Alive on Day20','% Alive on Day25']
L = len(data_all)
L2 = len(daystoplot)

max_samples2 = data_all.iloc[index2[0],ind_count]
ind_names_temp = ['M 20°C, C 25°C \n {}'.format(data_all.iloc[index2[2],ind_count]),
                  'M 25°C, C 20°C \n {}'.format(data_all.iloc[index2[3],ind_count]),
                  'Rand 1 M 20°C, C 20°C \n {}'.format(max_samples2),
                  'Rand 2 M 20°C, C 20°C \n {}'.format(max_samples2),
                  'Rand 3 M 20°C, C 20°C \n {}'.format(max_samples2),
                  'Rand 1 M 25°C, C 25°C \n {}'.format(max_samples2),
                  'Rand 2 M 25°C, C 25°C \n {}'.format(max_samples2),
                  'Rand 3 M 25°C, C 25°C \n {}'.format(max_samples2),]

for ii in range(L2):
    df_alive,lifespan_lists = plotlifespanspread(data_all, daystoplot[ii])
    columns = ['Set1','Set2','Set3','Set4','Set5','Set6','Set7','Set8','Set9','Set10',
               'Set11','Set12','Set13','Set14','Set15','Set16','Set17']
    df_alive.columns = columns
    df_alive = df_alive.transpose()
    
    df_tempFUDR = df_alive.iloc[index2,:]
    df_tempFUDR.index = ind_names_temp 
    
"""

#%% FUDR vs No FUDR
"""
                         FIGURE 2 
                                                             """
# Compare FUDR to w/o FUDR
ind_day3 = data_all2.columns.get_loc('% Alive on Day3')
ind_day50 = data_all2.columns.get_loc('% Alive on Day50')
ind_count = data_all2.columns.get_loc('Count')
wFUDR_list = [[] for _ in range(len(days))]
woFUDR_list = [[] for _ in range(len(days))]

for i in range(ind_day3,ind_day50+1):
    s = pd.Series(data_all2.iloc[0,i]).describe() # 0 = without FUDR
    woFUDR_list[i-ind_day3] = [s['count'],round(s['mean'],2),s['25%'],s['50%'],s['75%'],s['std']]

#Transpose list of lists
woFUDR_list = list(map(list, zip(*woFUDR_list))) 
df_woFUDR = pd.DataFrame(woFUDR_list,index=['count','mean','25%','median','75%','std'], 
                      columns=['Day3','Day5','Day10','Day15','Day20',
                               'Day25','Day30','Day40','Day50']).transpose()

for i in range(ind_day3,ind_day50+1):
    s = pd.Series(data_all2.iloc[1,i]).describe() # 1 = with FUDR
    wFUDR_list[i-ind_day3] = [s['count'],round(s['mean'],2),s['25%'],s['50%'],s['75%'],s['std']]

#Transpose list of lists
wFUDR_list = list(map(list, zip(*wFUDR_list))) 
df_wFUDR = pd.DataFrame(wFUDR_list,index=['count','mean','25%','median','75%','std'], 
                      columns=['Day3','Day5','Day10','Day15','Day20',
                               'Day25','Day30','Day40','Day50']).transpose()

fig = plt.figure(figsize=(10,6))
ax = fig.add_axes([0.1,0.15,0.8,0.8]) # main axes
color_cycle = ["black", "dodgerblue", "chocolate"] 
ind = [3,5,10,15,20,25,30,40,50]
ax.plot(ind,df_wFUDR['mean'], color=color_cycle[0],label='Mean (with FUDR)') # Change legend here
ax.plot(ind,df_woFUDR['mean'], color='r',label='Mean (w/o FUDR)') # Change legend here
ax.fill_between(ind, df_wFUDR['mean']-2*df_set1['std'], df_set1['mean']+2*df_set1['std'], 
                 facecolor=color_cycle[2], label='2x std')
ax.fill_between(ind, df_wFUDR['mean']-df_set1['std'], df_set1['mean']+df_set1['std'], 
                 facecolor=color_cycle[1],label='1x std')

# ax.plot(ind,df_set1['median'], color='k',label='Median (with FUDR)') # Change legend here
# ax.plot(ind,df_woFUDR['median'], color='r',label='Median (w/o FUDR)') # Change legend here
# ax.fill_between(ind, df_set1['median'], df_set1['75%'],alpha=0.3, 
#                  facecolor=color_cycle[0], label='75%')
# ax.fill_between(ind, df_set1['25%'], df_set1['median'],alpha=0.3, 
#                  facecolor=color_cycle[1], label='25%')

ax.set_ylim(-60,140)
ax.set_xticks(ind)
ax.set_xticklabels(df_set1.index, fontsize=12, rotation = 70)
ax.set_yticks(np.arange(-60,160,20))
ax.set_yticklabels(np.arange(-60,160,20), fontsize=12)
ax.set_title('With FUDR (n = {}) and Without FUDR (n = {})'.format(int(df_wFUDR['count'][0]),int(df_woFUDR['count'][0])), fontsize=14,x=0.5) # Change titles here
ax.set_ylabel('% Survival',fontsize=20)
ax.legend(fontsize=14)
plt.show()

# # Save the image
# png1 = io.BytesIO()
# fig.savefig(png1, format="png")

# # # Load this image into PIL
# png2 = Image.open(png1)

# # # Save as TIFF
# png2.save('/Users/nickurban/Desktop/N2 Lifespan Variability/Figures/Figure2A.tif')
# png1.close()

print('2A (-)FUDR: n = {}'.format(df_wFUDR['count'][0]))
print('2A FUDR: n = {}'.format(df_woFUDR['count'][0]))

df_wFUDR.to_csv('figures/Fig2A_wFUDR.csv')
df_woFUDR.to_csv('figures/Fig2A_woFUDR.csv')

#%% Figure 2B
ind_mls = data_all2.columns.get_loc('Reported mean lifespan')
mls_FUDR = [[] for _ in range(2)]

# NEED TO FILL BLANKS IN MEAN LIFESPAN COLUMN WITH -2.0
for j in range(2):
    temp = data_all2.iloc[j,ind_mls]
    temp = [x for x in temp if x > -1.0]
    #a = np.array(temp, int)
    mls_FUDR[j] = temp

# Create a figure instance
fig = plt.figure()
ax = fig.add_axes([0,0,1,1])

# Create the boxplot
#bp = ax.violinplot(mls_temp)
bp = ax.boxplot(mls_FUDR)
ax.set_ylim(0,40)
ax.set_xticks([1,2])
ax.set_xticklabels(['(-)FUDR \n (n = {})'.format(len(mls_FUDR[0])),
                    'FUDR \n (n = {})'.format(len(mls_FUDR[1]))], fontsize=12)
ax.set_yticks(np.arange(0,45,5))
ax.set_yticklabels(np.arange(0,45,5), fontsize=12)
ax.set_title('Average Reported Lifespan with and without FUDR', fontsize=16,x=0.5,y=1.03) # Change titles here
ax.set_ylabel('Average Lifespan (Days)',fontsize=12)

for line in bp['medians']:
    # get position data for median line
    x, y = line.get_xydata()[1] 
    # overlay median value
    plt.text(x, y, '%.1f' % y, verticalalignment='center') # draw above, centered

for line in bp['boxes']:
    x, y = line.get_xydata()[0] 
    plt.text(x,y, '%.1f' % y,
         verticalalignment='bottom',
         horizontalalignment='right')  
    x, y = line.get_xydata()[3]
    plt.text(x,y, '%.1f' % y,
         verticalalignment='top',
         horizontalalignment='right')

for line in bp['caps']:
    x, y = line.get_xydata()[0] 
    plt.text(x,y, '%.1f' % y,
         verticalalignment='bottom',
         horizontalalignment='right')  

print('t-test between FUDR and no FUDR gives a p-value of {}'.format(round(
    stats.ttest_ind(mls_FUDR[0],mls_FUDR[1])[1],3)))

print('2B (-)FUDR: n = {}'.format(len(mls_FUDR[0])))
print('2B FUDR: n = {}'.format(len(mls_FUDR[1])))

######### Saving the Data #############
#pd.DataFrame(mls_FUDR,index=['noFUDR','wFUDR']).transpose().to_csv('/Users/nickurban/Desktop/N2 Lifespan Variability/Figures/mls_fig2B.csv')
pd.DataFrame(mls_FUDR,index=['noFUDR','wFUDR']).transpose().to_csv('figures/Fig2B.csv')

#%% Figure 2C-E
daystoplot = ['% Alive on Day15','% Alive on Day20','% Alive on Day25']
L = len(data_all2)
L2 = len(daystoplot)
ind_count = data_all2.columns.get_loc('Count')

ind_tempFUDR = [13,0,2,10,1,3] # Indices of no FUDR (15C,20C,25C) and w/ FUDR (15C,20C,25C)
ind_names_tempFUDR = ['(-)FUDR \n {} (15°C)'.format(data_all2.iloc[13,ind_count]),
                      'FUDR \n {} (15°C)'.format(data_all2.iloc[10,ind_count]),
                      '(-)FUDR \n {} (20°C)'.format(data_all2.iloc[0,ind_count]),
                      'FUDR \n {} (20°C)'.format(data_all2.iloc[1,ind_count]),
                      '(-)FUDR \n {} (25°C)'.format(data_all2.iloc[2,ind_count]),
                      'FUDR  \n {} (25°C)'.format(data_all2.iloc[3,ind_count]),]

for ii in range(L2):
    df_alive,lifespan_lists = plotlifespanspread(data_all2, daystoplot[ii])
    columns = ['Set1','Set2','Set3','Set4','Set5','Set6','Set7','Set8','Set9','Set10',
               'Set11','Set12','Set13','Set14','Set15','Set16','Set17','Set18','Set19',
               'Set20','Set21','Set22','Set23','Set24','Set25','Set26','Set27','Set28',
               'Set29','Set30','Set31']
    df_alive.columns = columns
    df_alive = df_alive.transpose()
    
    df_tempFUDR = df_alive.iloc[ind_tempFUDR,:]
    df_tempFUDR.index = ind_names_tempFUDR
    
    ind = [1,2,3.5,4.5,6,7]  # the x locations for the groups
    width = 0.7
    
    fig = plt.figure(figsize=(16,10))
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

    # Save the image
    # png1 = io.BytesIO()
    # fig.savefig(png1, format="png")
    
    # # Load this image into PIL
    # png2 = Image.open(png1)
    #[print("%'s n = {}".format(x)) for x in df_tempFUDR['count']]
    if ii == 0:
        # Save as TIFF
        #png2.save('/Users/nickurban/Desktop/N2 Lifespan Variability/Figures/Figure2C.tif')
        png1.close()
        df_tempFUDR.to_csv('figures/Fig2C.csv')

    elif ii == 1:
        # Save as TIFF
        #png2.save('/Users/nickurban/Desktop/N2 Lifespan Variability/Figures/Figure2D.tif')
        png1.close()
        df_tempFUDR.to_csv('figures/Fig2D.csv')

    elif ii == 2:
        # Save as TIFF
        #png2.save('/Users/nickurban/Desktop/N2 Lifespan Variability/Figures/Figure2E.tif')
        png1.close()
        df_tempFUDR.to_csv('figures/Fig2E.csv')


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

# # Save the image
# png1 = io.BytesIO()
# fig.savefig(png1, format="png")

# # # Load this image into PIL
# png2 = Image.open(png1)

# png2.save('/Users/nickurban/Desktop/N2 Lifespan Variability/Figures/Figure2F.tif')
# png1.close()

df_FUDR.to_csv('figures/Fig2F.csv')

#%% Figure 2G
# Comparing FUDR Concentrations
def get_FUDR(data,conc):
    FUDRcolumns = ['% Alive on Day3','% Alive on Day5','% Alive on Day10','% Alive on Day15',
              '% Alive on Day20','% Alive on Day25','% Alive on Day30',
              '% Alive on Day40','% Alive on Day50','Reported mean lifespan']
    ind_day3 = data.columns.get_loc('% Alive on Day3')
    ind_mls = data.columns.get_loc('Reported mean lifespan')
    ind_FUDRconc = data.columns.get_loc('FUDR Concentration (microM)')
    listFUDR = data.iloc[(1,ind_FUDRconc)] #Create list of FUDR concentrations for entries
    ind_FUDR = [i for i, x in enumerate(listFUDR) if x == conc] # Find where conc = 40
    
    # Find and save % Alive values for FUDR = 40 entries for day 3 - 50 and rml
    FUDR_lists = [[] for _ in range(len(FUDRcolumns))]
    for i in range(ind_day3,ind_mls+1):
        for j in ind_FUDR:
            FUDR_lists[i-ind_day3].append(data.iloc[(1,i)][j])
    
    df_totFUDR = pd.DataFrame(FUDR_lists, index = FUDRcolumns)

    # Figure out stats for each day
    Ldays = len(FUDRcolumns)-1
    FUDR_datalists = [[] for _ in range(Ldays)]
    for i in range(Ldays):
        s = df_totFUDR.iloc[i,:].describe()
        FUDR_datalists[i] = [s['count'],round(s['mean'],2),s['25%'],s['50%'],s['75%'],s['std']]

    FUDR_datalists = list(map(list,zip(*FUDR_datalists)))
    
    # Turn list of lists into dataframe
    df_FUDR = pd.DataFrame(FUDR_datalists,index=['count','mean','25%','median','75%','std'], 
                              columns=['Day3','Day5','Day10','Day15','Day20',
                                       'Day25','Day30','Day40','Day50']).transpose()

    return df_totFUDR.transpose(),df_FUDR

df_totFUDR25, df_FUDR25 = get_FUDR(data_all2,25.0)
df_totFUDR40, df_FUDR40 = get_FUDR(data_all2,40.0)
df_totFUDR100, df_FUDR100 = get_FUDR(data_all2,100.0)
df_totFUDR200, df_FUDR200 = get_FUDR(data_all2,200.0)
df_totFUDR400, df_FUDR400 = get_FUDR(data_all2,400.0)

# Perform t-test % Alive for each day
# ttest_FUDR = []
# for i in range(len(FUDRcolumns-1)):
#     test = stats.ttest_ind(df_totFUDR40.iloc[i,:], df_totFUDR200.iloc[i,:], equal_var=False)
#     ttest_FUDR.append([test[0],test[1]])
    
# df_ttestFUDR = pd.DataFrame(ttest_FUDR,columns=['t statistic','p-value'],
#                             index=['Day3','Day5','Day10','Day15','Day20',
#                                      'Day25','Day30','Day40','Day50'])

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

# # Save the image
# png1 = io.BytesIO()
# fig.savefig(png1, format="png")

# # # Load this image into PIL
# png2 = Image.open(png1)

# png2.save('/Users/nickurban/Desktop/N2 Lifespan Variability/Figures/Figure2G.tif')
# png1.close()

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
#%% Figure 2H
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
#%% Plate transfers
"""
                         FIGURE 3 
                                                             """
# Indices and their plate transfer type
# 0 = none, 1 = Every other day, 2 = Everyday, 3 = After L4, 5 = Every 2-3 days, 7 = Every 1-2 days

# Figure 3A
df_totnm,df_totm,data_nm,data_m,df_totnm_all,df_totm_all = get_ttfp(data_all3,'Every other day','After L4',True,20,20,False,False)

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

ax1.set_xticklabels(df_FUDR40.index, fontsize=12, rotation = 70)
ax2.set_xticklabels(df_FUDR40.index, fontsize=12, rotation = 70)

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

# png1 = io.BytesIO()
# fig.savefig(png1, format="png")

# # Load this image into PIL
# png2 = Image.open(png1)

# png2.save('/Users/nickurban/Desktop/N2 Lifespan Variability/Figures/Figure3A.tif')
# png1.close()

df_totnm_noFUDR,df_totnm_wFUDR,data_nm_noFUDR,data_nm_wFUDR,df_totnm_all_noFUDR,df_totnm_all_wFUDR = get_ttfp(data_all4,'Every other day','After L4',False,20,20,True,False)
df_totm_noFUDR,df_totm_wFUDR,data_m_noFUDR,data_m_wFUDR,df_totm_all_noFUDR,df_totm_all_wFUDR = get_ttfp(data_all4,'Every other day','After L4',False,20,20,False,True)
ind_mls = df_totm_noFUDR.columns.get_loc('Reported mean lifespan')
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

#%% Figure 3B & 3C
# Doing only 20C because data does not allow for 15C and 25C 

# 0 = none, 1 = Every other day, 2 = Everyday, 3 = After L4, 5 = Every 2-3 days, 7 = Every 1-2 days
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

ax1.set_xticklabels(df_FUDR40.index, fontsize=12, rotation = 70)
ax2.set_xticklabels(df_FUDR40.index, fontsize=12, rotation = 70)

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

# # Save the image
#png1 = io.BytesIO()
#fig.savefig(png1, format="png")

# # Load this image into PIL
#png2 = Image.open(png1)

#png2.save('/Users/nickurban/Desktop/N2 Lifespan Variability/Figures/Figure3B&C.tif')
#png1.close()

data_ED_20.to_csv('figures/Fig3A_everyday.csv')
data_E23D_20.to_csv('figures/Fig3A_every 2-3 days.csv')

#%% Figure 3E (I need to ask more questions to do Figure 3D)
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
#pd.DataFrame(mls_ttfp,index=['Everyday','Every 2-3 Days']).transpose().to_csv('/Users/nickurban/Desktop/N2 Lifespan Variability/Figures/mls_fig3E.csv')
pd.DataFrame(mls_ttfp,index=['No Manipulation','Manipulation']).transpose().to_csv('figures/Fig3E_v2.csv')

#%% Figure 4B-D
"""
                        FIGURE 4
                                                        """
                                                        
######################## % Alive on by Country ################################
N = 0 # This determines the set that is plotted from whichever data_all is placed in the fxn 4 lines below
daystoplot = ['% Alive on Day15','% Alive on Day20','% Alive on Day25']
# Find and save % Alive values for FUDR = 40 entries for day 3 - 50
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
 
    # Save the image
    png1 = io.BytesIO()
    fig.savefig(png1, format="png")
    png2 = Image.open(png1)
    
    if k == 0:
    #     # Save as TIFF
        png2.save('figures/Figure4B.tif')
        png1.close()
        df_alive.to_csv('figures/Fig4B.csv')

    if k == 1:
    #     # Save as TIFF
        png2.save('figures/Figure4C.tif')
        png1.close()
        df_alive.to_csv('figures/Fig4C.csv')

    if k == 2:
    #     # Save as TIFF
        png2.save('figures/Figure4D.tif')
        png1.close()
        df_alive.to_csv('figures/Fig4D.csv')

    
#%% Figure 4E
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
                        FIGURE 5
                                                        """
                                                        
######################## % Alive on by State/Province ################################
N = 0
daystoplot = ['% Alive on Day15','% Alive on Day20','% Alive on Day25']
# Find and save % Alive values for FUDR = 40 entries for day 3 - 50
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
    #     Save as TIFF
        #png2.save('/Users/nickurban/Desktop/N2 Lifespan Variability/Figures/Figure5B.tif')
        png1.close()
        df_alive.to_csv('figures/Fig5B.csv')

    elif k == 1:
    #     Save as TIFF
        #png2.save('/Users/nickurban/Desktop/N2 Lifespan Variability/Figures/Figure5C.tif')
        png1.close()
        df_alive.to_csv('figures/Fig5C.csv')

    elif k == 2:
    #     Save as TIFF
        #png2.save('/Users/nickurban/Desktop/N2 Lifespan Variability/Figures/Figure5D.tif')
        png1.close()
        df_alive.to_csv('figures/Fig5D.csv')


#### Graph to show distribution and number of entries for each state/province
    # ind = np.arange(len(unique_group_states))    # the x locations for the groups    
    # fig = plt.figure(figsize=(14,8))
    # stats_list = [[] for _ in range(len(ind))]
    
    # for i in ind:
    #     x = [i]*unique_group_states.iloc[i,1]
    #     if i < 3:
    #         set_stats = pd.Series(set_data[i][0]).describe()
    #         stats_list[i] = [round(set_stats['mean'],2),round(set_stats['std'],2)]
    #         plt.scatter(x,set_data[i],label = '{} ± {}'.format(stats_list[i][0],stats_list[i][1]))
    #     else:
    #         plt.scatter(x,set_data[i])

    # plt.ylabel(daystoplot[k],fontsize=16)
    # plt.xlabel('Countries in Set 1',fontsize=16)
    # plt.xticks(ind,list(unique_group_states.iloc[:,0]), fontsize=14)
    # plt.yticks(np.arange(0, 110, 10),size=14)
    # plt.title('{} for Each Country in Set 1'.format(daystoplot[k]),fontsize=24,y=1.02)
    # plt.legend(bbox_to_anchor=(0.99, 1), loc='upper right', borderaxespad=0.2,\
    #         fontsize = 14, title = 'Mean ± Std Dev')
    # plt.show()
    
#%% Figure 5E
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
#pd.DataFrame(avg_set_states,index=unique_group_states['Names']).transpose().to_csv('/Users/nickurban/Desktop/N2 Lifespan Variability/Figures/mls_fig5E.csv')
pd.DataFrame(avg_set_states,index=unique_group_states['Names']).transpose().to_csv('figures/Fig5E.csv')

#%% Multivariate Linear Regression
data_mv = data_new[data_new['Reported mean lifespan'] != -1]

data_mv[data_mv['FUDR Concentration (microM)'] == -1] = 0
FUDR_ind = data_mv.columns.get_loc('FUDR Concentration (microM)')

x = data_mv[['Temperature Cultivated (Adult, C)','FUDR Concentration (microM)']]
y = data_mv[['Reported mean lifespan']]

regr = LinearRegression()
regr.fit(x,y)

print(regr.coef_)
