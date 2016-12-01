# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 21:26:03 2016

@author: Manuel
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

data = np.loadtxt('dec_lengths.dat')

eV = 1.6021766208e-19
n = 100000
P_K = 0.16
P_pi = 0.84
tau_pi = 2.6 *1e-8
tau_Pi = 4.188*1e3
p_K = 75 *1e9 *eV

'''   Reines Ausgeprobiere, kann man ignorieren
bns = 100

dist_pi = stats.expon.rvs(scale=tau_Pi, size=int(P_pi*n))

plt.subplot(3,1,1)
n_alt, b = np.histogram(data, bins=int(bns), range=[0,1000])

plt.subplot(3,1,2)
n_pi, bins_pi, p = plt.hist(dist_pi, bins=int(bns), range=[0,1000])
#print(n_pi)

n_neu = n_alt - n_pi

plt.subplot(3,1,3)
plt.plot(np.linspace(0,1000,bns),n_neu)
plt.xlim(0,1000)

tau_K = np.mean(n_neu)
print(tau_K)
'''

#tau_K estimator

bns = 100            #Anzahl Bins
N = 1000              #Anzahl Stichproben

n_alt, b = np.histogram(data, bins=int(bns), range=[0,1000])
tau_K_estimator = np.zeros(int(N))
for i in range(N):
    dist_pi = stats.expon.rvs(scale=tau_Pi, size=int(P_pi*n))
    n_pi, b = np.histogram(dist_pi, bins=int(bns), range=[0,1000])
    n_neu = n_alt-n_pi
    tau_K_estimator[i] = np.mean(n_neu)
tau_K = np.mean(tau_K_estimator)
print('tau_K ', tau_K)
plt.hist(tau_K_estimator)
    










