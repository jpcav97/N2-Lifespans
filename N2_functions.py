#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 17:20:37 2021

@author: josephcavataio
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
import statsmodels.api as sm
from sklearn.linear_model import LinearRegression

column_days = ['Day3','Day5','Day10','Day15','Day20','Day25','Day30','Day40','Day50']

def get_countrystate(data,tempm,tempc,cntrystate):
    columns = ['% Alive on Day3','% Alive on Day5','% Alive on Day10','% Alive on Day15',
              '% Alive on Day20','% Alive on Day25','% Alive on Day30',
              '% Alive on Day40','% Alive on Day50']
    
    df = data.loc[(data['Growth Media'] == 'NGM') & 
                  (data['Temperature Maintained (Through L4, C)'] == tempm) &
                  (data['Temperature Cultivated (Adult, C)'] == tempc)]
    
    if cntrystate == True:
        ind_group = df.columns.get_loc('Lab Location (Country)')
    else:
        ind_group = df.columns.get_loc('Lab Location (State/Province)')
    
    group,freq_group = np.unique(data.iloc[(0,ind_group)],return_counts=True)
    unique_group=pd.DataFrame(list(zip(group,freq_group)),columns=['Names','Frequency'])
    unique_group = unique_group.reindex(np.argsort(unique_group['Frequency'])[::-1])
    unique_group  = unique_group.reset_index(drop=True)
    
    lifespan_lists = [[] for _ in range(len(columns))]
    ind_start = df.columns.get_loc('% Alive on Day3')
    for i in range(len(columns)):
        lifespan_lists[i] = df.iloc[(0,ind_start+i)]
        
    list_of_data = [_ for _ in range(len(group))]
    group_list = df.iloc[0,ind_group]
    df_lifespan = pd.DataFrame(lifespan_lists).transpose()
    
    for i in range(len(group)):
        index = np.where(np.isin(group_list,unique_group.iloc[(i,0)]))
        d20 = df_lifespan.iloc[(index)]     
        
        # Figure out stats for each day
        Ldays = len(columns)
        datalists = [[] for _ in range(Ldays)]

        for j in range(Ldays):
            s1 = d20.iloc[:,j].describe()
            datalists[j] = [s1['count'],round(s1['mean'],2),s1['25%'],s1['50%'],s1['75%'],s1['std']]
    
        datalists = list(map(list,zip(*datalists)))
        
        # Turn list of lists into dataframe
        data = pd.DataFrame(datalists,index=['count','mean','25%','median','75%','std'], 
                              columns=['Day3','Day5','Day10','Day15','Day20',
                                       'Day25','Day30','Day40','Day50']).transpose()
        list_of_data[i] = data
        
    return list_of_data,unique_group

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
                              columns=column_days).transpose()
    data_m = pd.DataFrame(datalists_m,index=['count','mean','25%','median','75%','std'], 
                              columns=column_days).transpose()
    del ttfp_lists_nm,ttfp_lists_m,ttfp_lists_nm_all,ttfp_lists_m_all
    return df_totnm,df_totm,data_nm,data_m,df_totnm_all,df_totm_all

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
                              columns=column_days).transpose()

    return df_totFUDR.transpose(),df_FUDR

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
    
    print('The total number of measurements = {}'.format(sum(grp['Count'])))
    grp  = grp.reindex(np.argsort(grp['Count'])[::-1])
    grp  = grp.reset_index(drop=True)
    dat = grp[grp["Count"]>=P]
    print('The number of measurements after 3+ restriction = {}'.format(sum(dat['Count'])))
    
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


