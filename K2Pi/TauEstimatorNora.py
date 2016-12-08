# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 18:50:27 2016

@author: Nora
"""

import numpy as np

import scipy.stats as stats

data = np.loadtxt('dec_lengths.dat')

n=len(data)
tau_pi=4188
frac_pi=0.84

bns = 100    
n_tot, b = np.histogram(data, bins=bns)
x=b+(max(b)-min(b))/200
x=x[:-1] 

N=100
tau_k_estimator=[]
for i in range(N):
    dist_pi = stats.expon.rvs(scale=tau_pi, size=int(frac_pi*n))
    n_pi, b = np.histogram(dist_pi, bins=int(bns), range=[min(data),max(data)])

    n_k=n_tot-n_pi
    tau_k_est = sum(x*n_k)/sum(n_k)
    tau_k_estimator.append(tau_k_est)

tau_k_estimator=np.array(tau_k_estimator)
tau_k_estimator=np.mean(tau_k_estimator)


print(tau_k_estimator)


