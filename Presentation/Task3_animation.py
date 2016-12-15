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


xrand = np.random.randint(100,size=100)
yrand = np.random.randint(100,size=100)

p=1
oldp=0
active=[1]
t=0
deceased=[0]
deceasedold=0
deltadead=0
casualties=0
EW=298     #expected value, taken from task 1
while casualties<100:
	if t>=round(EW/24.,0):  #if the first infected people die...
		deltadead = (active[-int(round(EW/24.,0))]-deceased[-1]+deceased[-int(round(EW/24.,0))])
		deceased.append(deltadead+deceasedold)
		deceasedold=deceased[-1]
	else:
		deceased.append(0)
	oldp=p-deltadead
	p=p-deltadead
	temp= ahat_0 + ahat_1*p ##calculate probability of infection using task 2
	temp=(100-p-deceased[-1])*(temp/100.)
	p=round(temp,0)+oldp
	active.append(p)
	t+=1
	if t==50:
		break
plt.ion()
axes=plt.gca()
axes.set_xlim([-10,110])
axes.set_ylim([-10,110])
days=len(active)
z=0
print ""
print len(xrand)
repeats=1   ##enter number of repeats of the animation
while z<repeats: ##begin animation
	dead=[]
	deceasedx=[]
	healthyx=[]
	infectedx=[]
	for j in range(len(xrand)):
		dead.append(0)
		deceasedx.append(0)
		healthyx.append(1)
		infectedx.append(0)
	for t in range(days):
		p=active[t]
		d=deceased[t]
		if d==0:
			for i in range(len(xrand)):
				if i<p and deceasedx[i]==0 and infectedx[i]==0:
					infectedx[i]=1
					healthyx[i]=0	
		if d>0:
			for i in range(len(xrand)):
				if i<d and healthyx[t-i-int(round(EW/24.,0))]==0:
					deceasedx[t-i-int(round(EW/24.,0))]=1
					healthyx[t-i-int(round(EW/24.,0))]=0
					infectedx[t-i-int(round(EW/24.,0))]=0
				if i<p+d and deceasedx[i]==0 and infectedx[i]==0:
					infectedx[i]=1
					healthyx[i]=0
		for i in range(len(xrand)):
			if deceasedx[i]==1:
				plt.scatter(xrand[i],yrand[i],color="black",marker="+",s=80)
			elif infectedx[i]==1:
				plt.scatter(xrand[i],yrand[i],color="red",marker="o", s=80)
			elif healthyx[i]==1:
				plt.scatter(xrand[i],yrand[i],color="cyan",marker="o",s=60)
		timer=str("T+" + str(t)+ " days")
		plt.text(-10, 130, timer, bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
		plt.text(xrand[0], yrand[0], "patient zero")
		if p==0:
			plt.text(xrand[-1],yrand[-1], "survivor")
		plt.draw()
		plt.scatter([],[],color="cyan",marker="o",label="healthy", s=60)
		plt.scatter([],[],color="red",marker="o",label="infected",s=60)
		plt.scatter([],[],color="black",marker="+",label="deceased", s=60)
		plt.legend()
		axes=plt.gca()
		axes.set_xlim([-20,120])
		axes.set_ylim([-20,140])
		if t<(days):
			plt.pause(0.1)
		else:
			plt.pause(2)
		name='image_' + str(t) + ".png"
		#plt.savefig(name)
		plt.clf()
		percent = int(p)
		bar = "#"*int(percent/10) + " "*int(10-percent/10)
		output = "Infected: " + str(percent) + "%" + "   " + "[" + bar + "]"
		sys.stdout.write('\r')
		sys.stdout.write(str(timer)+"      ")
		sys.stdout.write(output)
		sys.stdout.flush()
	
	z+=1
	print "\n"
	print infectedx
	print deceasedx

