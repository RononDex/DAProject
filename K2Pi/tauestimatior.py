# -*- coding: utf-8 -*-
"""
Created on Sun Dec  4 15:53:55 2016

@author: nsalgo
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize


data = np.loadtxt('dec_lengths.dat')
bns = 100    
n, b = np.histogram(data, bins=bns)

y=n                                                 #number of decays
x=b+(max(b)-min(b))/200
x=x[:-1]                                            #center of bins

def f(x,A,B,C,D):       
#A = number of kaons in sample 
#B = frac of pions / frac of K 
#C = decay length of pion 
#D = the decay length of the Kaon
    return A * (B * np.exp(-1 * x / C) + np.exp(-1 * x / D) )

para, cov = scipy.optimize.curve_fit(f,x,y,bounds=([0,0,4187,0], [100000, 0.84/0.16 ,4189,4188]),p0=[16000,5.25,4188, 500])
print(para)


plt.figure()
plt.plot(x,y, label='data')                          #daten geplottet
plt.plot(x, f(x, *para),label='fited function')

plt.xlabel('decay length [m]')
plt.ylabel('number of decays')
plt.legend()
#plt.savefig('PlotTauEstimator.jpg', dpi=480)
plt.show()


