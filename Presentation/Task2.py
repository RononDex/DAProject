# -*- coding: utf-8 -*-
# @Author: Marco Tresch
# @Date:   2015-09-15
# @Last Modified by:   Marco Tresch
# @Last Modified time: 2015-09-15
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

def read_from_file(filename):
    """
    input   filename, string with path and filename
    output  array with entries of file
    """
    return np.loadtxt("epidemiology_contagion.dat")

def mean(x):
	sum=0
	for i in x:
		sum+=i
	return sum/len(x)
	



data=read_from_file("")
x=data[:,1]
y=data[:,0]
ahat_0 = (mean(x**2)*mean(y)-mean(x)*mean(x*y))/float(mean(x**2)-mean(x)*mean(x))
ahat_1 = (mean(x*y)-mean(x)*mean(y))/float(mean(x**2)-mean(x)*mean(x))
print "ahat_0:", ahat_0
print "ahat_1:", ahat_1
funcx=np.linspace(0,100,1000)
funcy= ahat_0 + ahat_1*funcx

cov=1/float(len(x))*np.array([[mean(x**2), -mean(x)],[-mean(x), 1]])



print "uncertainty on ahat_0 : ", cov[0][0]**0.5
print "uncertainty on ahat_1 : ", cov[1][1]**0.5
print "covariance of ahat_0 and ahat_1 =cov(x,y): ", cov[0][1]
print "p(n=1)=", funcy[0], "%"

steepy= ahat_0-cov[0][0]**0.5+(ahat_1 + cov[1][1]**0.5)*funcx
flaty = ahat_0+cov[0][0]**0.5+(ahat_1 - cov[1][1]**0.5)*funcx



plt.scatter(data[:,1],data[:,0], label="datapoints", color="blue", marker="o")
plt.plot(funcx,funcy, label="unbinned max likelihood fit", color="blue")
plt.plot(funcx,steepy, label="steepest fit", color="red")
plt.plot(funcx,flaty, label="flattest fit", color="black")

plt.xlabel("$n$ number of people")
plt.ylabel("$p$ probability of infection in %")
plt.grid()
fig = plt.figure
axes=plt.gca()
axes.set_xlim([0,100])
axes.set_ylim([-10,50])
plt.legend(loc=4)
plt.show()

