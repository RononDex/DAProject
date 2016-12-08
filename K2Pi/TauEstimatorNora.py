# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 18:50:27 2016

@author: Nora
"""

import numpy as np
import scipy.stats as stats

data = np.loadtxt('dec_lengths.dat')

n=len(data)
tau_pi=4188                                 #average decay length of pions, given
frac_pi=0.84                                #fraction of the pions in the beam

bns = 100    
n_tot, b = np.histogram(data, bins=bns)     #histogram of given data
x=b+(max(b)-min(b))/200                     #center of bins
x=x[:-1] 

N=100
tau_k_estimator=[]
for i in range(N):
    dist_pi = stats.expon.rvs(scale=tau_pi, size=int(frac_pi*n))    #exponential distribution with tau=average decay length of pions
    n_pi, b = np.histogram(dist_pi, bins=int(bns), range=[min(data),max(data)])

    n_k=n_tot-n_pi                          #values of histogram of decay length of kaons
    tau_k_est = sum(x*n_k)/sum(n_k)         #average decay length of kaons
    tau_k_estimator.append(tau_k_est)

tau_k_estimator=np.array(tau_k_estimator)   
tau_k_estimator=np.mean(tau_k_estimator)    #mean of average decay length pf kaons


print(tau_k_estimator)
