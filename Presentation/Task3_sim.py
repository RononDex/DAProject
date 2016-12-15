print "importing..."
from sympy import Rational as frac
import matplotlib.pyplot as plt
from math import cos as cos
from math import sin as sin
from math import pi as pi
from math import *
from pylab import *
import time
import sys
import re
import numpy as np
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
p=1
oldp=0
active=[1]
t=0
deceased=[0]
deceasedold=0
deltadead=0
casualties=0
alert =0 
EW=298 #Erwartungswert
while casualties<100:
	if t>=round(EW/24.,0):
		deltadead = (active[-int(round(EW/24.,0))]-deceased[-1]+deceased[-int(round(EW/24.,0))])
		deceased.append(deltadead+deceasedold)
		deceasedold=deceased[-1]
	else:
		deceased.append(0)
	oldp=p-deltadead
	p=p-deltadead
	temp= ahat_0 + ahat_1*p 
	temp=(100-p-deceased[-1])*(temp/100.)
	p=round(temp,0)+oldp
	active.append(p)
	t+=1
	if p+deceased[-1]>=30 and alert==0:
		alert =1
		print "after ", t, " days, ", active[-1]+deceased[-1], " are infected or dead"
		print "after ", t-1, "days, ", active[-2]+deceased[-2], " are infected or dead"
	if t==50:
		break
print "..."
print deceased
print active
healthy=[]
for i in range(len(active)):
	healthy.append(100)
for h in range(len(healthy)):
	healthy[h]=100-(active[h]+deceased[h])
days=[]
for i in range(len(active)):
	days.append(i)
axes=plt.gca()
axes.set_xlim([-5,50])
axes.set_ylim([-10,150])
plt.scatter(days,healthy, color="cyan", label="healthy people")
plt.scatter(days,active, color="red", label="infected people")
plt.scatter(days,deceased, color="black", label="deceased people")
plt.xlabel("days")
plt.ylabel("number of people infected")
plt.grid()
plt.legend()
plt.show()
