# -*- coding: utf-8 -*-
"""
Created on Sun Dec  4 15:53:55 2016

@author: nsalgo
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize



data = np.loadtxt('dec_lengths.dat')
bns = 200

xdata=np.linspace(0,40000,bns)      #instead of 40000 we tried max(data)= 43156 but the results are inconvinient
n, b = np.histogram(data, bins=int(bns))



def f1(t, tau_k, N):
    tau_pplus=4188
    a=0.84
    f=N*(a/tau_pplus*np.exp(-t/tau_pplus)+(1-a)/tau_k*np.exp(-t/tau_k))
    return f

def f2(t, tau_k,a, N):
    tau_pplus=4188
    f=N*(a/tau_pplus*np.exp(-t/tau_pplus)+(1-a)/tau_k*np.exp(-t/tau_k))
    return f
    
def f3(t, tau_k,tau_pplus,a, N):
    
    f=N*(a/tau_pplus*np.exp(-t/tau_pplus)+(1-a)/tau_k*np.exp(-t/tau_k))
    return f
    

pars1, cov1 = scipy.optimize.curve_fit(f1,xdata,n, p0=[500, 10000 ])
pars2, cov2 = scipy.optimize.curve_fit(f2,xdata,n, p0=[500, 0.84, 10000 ])
pars3, cov3 = scipy.optimize.curve_fit(f3,xdata,n, p0=[500, 4188, 0.84, 10000 ])

print('estimated tau: \n', 'f1: ', pars1[0], 'f2: ', pars2[0], 'f3: ', pars3[0])
print('estimated fraction of pions: \n', 'f2: ', pars2[1], 'f3: ', pars3[1])


plt.figure()
tau_pplus1=4188
plt.plot(xdata,n, label='data')   #daten geplottet

plt.plot(xdata, f1(xdata, *pars1),label='f1')
plt.plot(xdata, f2(xdata, *pars2),label='f2')
plt.plot(xdata, f3(xdata, *pars3),label='f3')
plt.xlabel('decay length [m]')
plt.ylabel('number of decays')
plt.legend()
#plt.hist(data,bns)
plt.savefig('plot2.jpg', dpi=480)
plt.show()
