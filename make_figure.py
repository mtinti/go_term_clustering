# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 15:01:40 2017

@author: mtinti-x
"""

import os
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import statsmodels.api as sm
lowess = sm.nonparametric.lowess

combined_go = pd.DataFrame.from_csv('reduced_go_terms.csv')
for col in combined_go.columns:
    combined_go[col]=combined_go[col]/combined_go[col].max()
    
combined_random_go = pd.DataFrame.from_csv('reduced_random_go_terms.csv')
for col in combined_random_go.columns:
    combined_random_go[col]=combined_go[col]/combined_random_go[col]

combined_random_go.columns = [n+'f' for n in combined_random_go.columns]
combined_pr= pd.DataFrame.from_csv('reduced_pearson.csv')
#print combined_go.head()
#print combined_pr.head()
#print combined_random_go.head()


combined_go.plot()
combined_random_go.plot()
combined_pr.plot()

combined = pd.concat([combined_go,combined_random_go,combined_pr],1)
combined = combined.fillna(0)
print combined.head()
combined.plot()
plt.xlim(0,5)




'''
#combined_norm = combined[['CC','MF','BP']]
combined_go = combined_go[['MF','BP']]
#print combined_norm.head()
for col in combined_go.columns:
    combined_go[col]=combined_go[col]/combined_go[col].max()
#print combined_norm.head()
combined_go['median']=combined_go.median(1)
combined_go['median'] = combined_go['median']/combined_go['median'].max()
combined_go['lowess'] = lowess(combined_go['median'], combined_go.index.values, frac = 0.1, return_sorted = False)
print combined_go.head()
#print np.argmax(combined_norm['median'])
plt.plot(combined_go['lowess'])
plt.xlim(0,10)
plt.ylim(0,1.1)
'''


#plt.plot(combined_norm['MF'])
#plt.plot(combined_norm['BP'])
#fig, ax = plt.subplots(nrows=1,ncols=1)