# -*- coding: utf-8 -*-
"""
TASK 1

@author: Philipp, the MOFO
"""
"""
Als erstes haben wir die gegebenen Daten histogrammiert, um die zugrunde liegende Verteilung
abzuschätzen. Aus dem Histogramm haben wir auf eine Exponentialverteilung geschlossen, was auch intuitiv
Sinn macht (fidnet Nicola nicht).

Als zweites wollten wir die Max. Likelihood Funktion errechnen und haben dazu die ln L Funktion mit der
unteren (min(data) = 169) und oberen Grenze (596) normiert.

Anschliessend haben wir die ln L Funktion für ein Array von Tau Werten ausgerechnet und so das maximierende
Tau-Hat gefunden.

Der Erwartungswert ist nun gerade mit Tau = 128.98+169 = 298 gegeben.
"""



import numpy as np
import matplotlib.pyplot as plt
from operator import itemgetter


def read_from_file(filename):
    return np.loadtxt(filename)



def histogram(data):
    progress_data = data
    x = np.arange(169, 1000)
    y = 8.23*np.e**(-(x-297.98)/128.98)
    """183.5=200-33/2 zur normierung wegen bin)"""
    sigma_1 = 8.23*np.e**(-(x-297.98)/149.73)
    sigma_2 = 8.23*np.e**(-(x-297.98)/112.74)
    xmin, xmax = 0, 1000
    f = plt.figure(dpi = 500)
    h = plt.hist(progress_data, bins=30, range=(xmin, xmax), label = "Data")
    plt.plot(x, y, "r", linewidth = 1, label = "Maximum Likelihood Fit")
    plt.plot(x, sigma_1, "k", linewidth = 1, label = "Uncertainties")
    plt.plot(x, sigma_2, "k", linewidth = 1)
    plt.plot(297.98, 8.23, "wo", label = "tau-hat$_{shift}$")
    plt.legend()
    axes = plt.gca()
    axes.set_ylim([0, 25])
    plt.grid()
    plt.xlabel("Hours survived after first symptoms appeared")
    plt.ylabel("Number of infected individuals")
    return h, f
    
    
def dist(t, tau, data):
    normkonst = -np.e**(-max(data)/tau)+np.e**(-min(data)/tau)
    return 1/tau*np.e**(-t/tau)/normkonst

    
def maxvalue(l, i):
    key = itemgetter(1)
    return max(enumerate(sub[i] for sub in l), key = key)
    
def maximum_likelihood(data, dist, tau_array):
    LL = []
    indexlist = []
    for tau in range(len(tau_array)):
        lnL = 0.
        for i in range(len(data)):
            lnL += np.log(dist(data[i], tau_array[tau],data))
        indexlist.append([tau, lnL])
        LL.append(lnL)
    fig = plt.figure(dpi = 501)
    plt.plot(tau_array, LL, linewidth = 1, label = "Log-Likelihood")
    plt.plot(128.98, -501.34368, "wo", label = "tau-hat")
    plt.plot([112.73516, 149.72648], [-501.84368, -501.84368], "r", linestyle = "--", label = "ln L(tau-hat) - 0.5")
    plt.plot([112.73516, 112.73516], [-502.4, -501.84368], "b", linestyle = "--")
    plt.plot([149.72648, 149.72648], [-502.4, -501.84368], "b", linestyle = "--")
    plt.plot([128.98, 128.98], [-502.4, -501.34368], "b", linestyle = "--")
    plt.legend()
    plt.grid()
    plt.xlabel("tau")
    plt.ylabel("ln L(tau)")
    maxval = list(maxvalue(indexlist, 1))
    tau_hat = tau_array[maxval[0]]
    print "tau_hat: ", tau_hat
    print "new: ", max(LL)
    print "index: ", LL.index(-501.84354233514114)
    print "index: ", tau_array[LL.index(-501.84354233514114)]
    return LL, tau_hat, fig
    
"""calculate std of list"""
varianz = 0.
data = read_from_file("epidemiology_progress.dat")
temp = 0.
for i in data:
    temp += (i-np.mean(data))**2
varianz = temp/len(data)
std = varianz**(1/2.)
print(std)




data = read_from_file("epidemiology_progress.dat")
h, f = histogram(data)
plt.show()

print maximum_likelihood(data, dist, np.linspace(110, 160, 1500))
plt.show()

