import numpy as np
import scipy.stats as stats


"""
theta=[0,pi]
phi=[0,2pi]
gleichverteilt in K system
"""

n=1
theta = np.array(stats.uniform.rvs(scale=np.pi, size=n))
phi = np.array(stats.uniform.rvs(scale=2*np.pi, size=n))



x0= np.sin(theta)*np.cos(phi)
y0= np.sin(theta)*np.sin(phi)
z0= np.cos(theta)



"""
x,y,z=kartesische Koordinaten
E_K=Energie der pionen in K+ system
p_K=four vector of pions in K+ system
g=gammafactor
b=betafactor
"""

"""
E_K0 =0
p_0= 0 
P_K0 =[]
for i in range(len(x0)):
    P_K0.append(np.array([E_K0,p_0*x0[i],p_0*y0[i],p_0*z0[i]]))
P_K0=np.array(P_K0)

E_K_plus = 0
p_plus= 0 
P_K_plus =[]
for i in range(len(x0)):
    P_K_plus.append(np.array([E_K_plus,p_plus*(-x0[i]),p_plus*(-y0[i]),p_plus*(-z0[i])]))
P_K_plus=np.array(P_K_plus)
"""
#E_K0,p_0 wurden von Manuel berechnet, bitte alle nachrechnen

E_K0 = 245.563588 #MeV
p_0= 205.14091 #MeV/c

E_K_plus = 248.118174 #MeV
p_plus= p_0 

P_K_0=np.array([E_K0,p_0*x0,p_0*y0,p_0*z0])
P_K_plus=np.array([E_K_plus,p_plus*(-x0),p_plus*(-y0),p_plus*(-z0)])

#Manus Werte
b = 0.99997833784995
g = 151.92756392754
boost=np.array([[g,0,0,b*g],[0,1,0,0],[0,0,1,0],[b*g,0,0,g]])

P_lab_0 = np.dot(boost,P_K_0.T)
P_lab_plus = np.dot(boost,P_K_plus.T)


"""
Zerfall von K

a = Distanz zum Downastream Detektor
tau = Zerfallsstrecke des Kaons
ez = Einheitsvektor in Z-Richtung
d = Abstand zur Z-Achse auf Detektor

"""

a = 280
tau=132.97
z_length= np.array(stats.expon.rvs(loc=0, scale=tau, size=1))


def GetAngleBetweenVectors(X, Y):
    x=np.linalg.norm(X)
    y=np.linalg.norm(Y)
    return np.arccos(np.dot(X,Y)/(x*y))

p_0 = P_lab_0[1:]
p_plus = P_lab_plus[1:]
ez= np.array([0,0,1])

d_0 =(a-z_length) * np.tan(GetAngleBetweenVectors(p_0,ez))
d_plus =(a-z_length) * np.tan(GetAngleBetweenVectors(p_plus,ez))
print(d_0)
print(d_plus)


'''Das Skript in sauberer Form:'''

#functions:

def GetAngleBetweenVectors(X, Y):
    x=np.linalg.norm(X)
    y=np.linalg.norm(Y)
    return np.arccos(np.dot(X,Y)/(x*y))

def LorentzBoost(b,g):
    return np.array([[g,0,0,b*g],[0,1,0,0],[0,0,1,0],[b*g,0,0,g]])

def SimulateDecay(E_K_0, E_K_plus, p, b, g, a, tau):
    theta = np.array(stats.uniform.rvs(scale=np.pi, size=1))
    phi = np.array(stats.uniform.rvs(scale=2*np.pi, size=1))
    x0= np.sin(theta)*np.cos(phi)
    y0= np.sin(theta)*np.sin(phi)
    z0= np.cos(theta)
    
    P_K_0 = np.array([E_K0,p*x0,p*y0,p*z0])
    P_K_plus = np.array([E_K_plus,p*(-x0),p*(-y0),p*(-z0)])

    P_lab_0 = np.dot(LorentzBoost(b,g),P_K_0.T)
    P_lab_plus = np.dot(LorentzBoost(b,g),P_K_plus.T)
    
    z_length= np.array(stats.expon.rvs(loc=0, scale=tau, size=1))
    if float(z_length) >= a:
        return [0,0]
    else:
        ez= np.array([0,0,1])
        d_0 = float((a-z_length) * np.tan(GetAngleBetweenVectors(P_lab_0[1:],ez)))
        d_plus = float((a-z_length) * np.tan(GetAngleBetweenVectors(P_lab_plus[1:],ez)))  
        return [d_0, d_plus]

def successrate(E_K_0, E_K_plus, p, b, g, a, tau, n):
    Pi_0 = np.zeros(n)
    Pi_plus = np.zeros(n)
    for i in range(n):
        decay = SimulateDecay(E_K_0, E_K_plus, p, b, g, a, tau)
        Pi_0[i] = decay[0]
        Pi_plus[i] = decay[1]
    success = 0
    for i in range(n):
        if Pi_0[i] and Pi_plus[i] <= 2:
            success += 1
    return success/n
    
#parameters:

E_K_0 = 245.563588 #MeV         #Energie der neutrale Pionen in K+ system
E_K_plus = 248.118174 #MeV      #Energie der positiven Pionen in K+ system
p = 205.14091 #MeV/c            #Impulsbetrag der Pionen (der selbe für beide)
b = 0.99997833784995            #Betafaktor
g = 151.92756392754             #Gammafaktor
a = 280                         #Totale Distanz zum Detektor
tau = 132.97                    #Mittlere Zerfallsstrecke von K+
n = 1000                        #Anzahl Teilchen

#Auswertung:
'''
plt.figure()
plt.plot(range(n),Pi_0,'r')
plt.plot(range(n),Pi_plus,'g')
plt.plot([0,n],[2,2],'k:')
plt.show()
'''

A = np.linspace(0,500,500)

'''
SR = []
for i in A:
    SR.append(successrate(E_K_0, E_K_plus, p, b, g, i, tau, n))
print(SR)
'''
    
plt.figure()
plt.plot(A,successrate(E_K_0, E_K_plus, p, b, g, A, tau, n))
plt.show()
