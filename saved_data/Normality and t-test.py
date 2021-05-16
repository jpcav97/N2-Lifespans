#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 22 16:35:28 2021

@author: josephcavataio
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats as st
import statistics as stats
import os
from PIL import Image
import io
import random
import statsmodels.api as sm
from sklearn.linear_model import LinearRegression
import math

def measuresofquality(data):
    # Mean
    mean = stats.mean(data)
    print("\n The mean is {}".format(round(mean,3)))

    # Variance
    var = stats.variance(data)
    print("\n The variance is {}".format(round(var,3)))
    
    # Standard Deviation
    stdev = stats.stdev(data)
    print("\n The Standard Deviation is {}".format(round(stdev,3)))
    
    # Shapiro-Wilk Test
    s_test = st.shapiro(data)
    
    if s_test.pvalue < 0.05:
        print('\n Measurements come from a normal distribution')
    else:
        print('Measurements do not come from a normal distribution\n')
    
    Q1 = mean - 0.675*stdev
    Q3 = mean + 0.675*stdev
    print("Q1 = {} and Q3 = {}".format(round(Q1,2),round(Q3,2)))
    return


mls_data = pd.read_csv('Fig2B.csv',index_col=0)
noFUDR = mls_data['noFUDR']
wFUDR = mls_data['wFUDR'].dropna()

plt.figure()
d1,freq1 = np.unique(noFUDR,return_counts=True)
plt.hist(noFUDR,bins = np.arange(min(d1),max(d1),1), color= 'white',edgecolor='k')

plt.xlabel('Mean Lifespan')
plt.ylabel('Frequency')
plt.title('No FUDR')
plt.tight_layout()
plt.show()

plt.figure()
d2,freq2 = np.unique(wFUDR,return_counts=True)
plt.hist(wFUDR, bins = np.arange(min(d2),max(d2),1),color= 'white',edgecolor='k')

plt.xlabel('Mean Lifespan')
plt.ylabel('Frequency')
plt.title('With FUDR')
plt.tight_layout()
plt.show()

measuresofquality(noFUDR)
measuresofquality(wFUDR)
