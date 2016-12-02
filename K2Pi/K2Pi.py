import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

#Funktionen:

def GetAngleBetweenVectors(X, Y):                                   #nicht orientierter Winkel zwischen zwei Vektoren
    x=np.linalg.norm(X)
    y=np.linalg.norm(Y)
    return np.arccos(np.dot(X,Y)/(x*y))

def LorentzBoost(b,g):                                              #Lorentz-Boost, b = Beta-Faktor, g = Gamma-Faktor
    return np.array([[g,0,0,b*g],[0,1,0,0],[0,0,1,0],[b*g,0,0,g]])

def SimulateDecay(E_K_0, E_K_plus, p, b, g, a, tau):
    theta = np.array(stats.uniform.rvs(scale=np.pi, size=1))        #gleichverteilt zwischen 0 und Pi
    phi = np.array(stats.uniform.rvs(scale=2*np.pi, size=1))        #gleichverteilt zwischen 0 und 2 Pi
    x0= np.sin(theta)*np.cos(phi)                                   
    y0= np.sin(theta)*np.sin(phi)                                   #kartesische Koordinaten
    z0= np.cos(theta)
    
    P_K_0 = np.array([E_K_0,p*x0,p*y0,p*z0])                        #4-Vektoren im K+-Frame
    P_K_plus = np.array([E_K_plus,p*(-x0),p*(-y0),p*(-z0)])

    P_lab_0 = np.dot(LorentzBoost(b,g),P_K_0.T)                     #4-Vektoren im Lab-Frame
    P_lab_plus = np.dot(LorentzBoost(b,g),P_K_plus.T)
    
    z_length= np.array(stats.expon.rvs(loc=0, scale=tau, size=1))   #Ortsvektor des K+-Zerfalls
    if float(z_length) >= a:                                        #Aussortieren der K+, die hinter Detektor zerfallen
        return [0,0]
    else:
        ez= np.array([0,0,1])                                       #Einheitsvektor in Z-Richtung
        d_0 = float((a-z_length) * np.tan(GetAngleBetweenVectors(P_lab_0[1:],ez)))
        d_plus = float((a-z_length) * np.tan(GetAngleBetweenVectors(P_lab_plus[1:],ez)))  
        return [d_0, d_plus]                                        #Output: Abstand von Mittelpunkt des Detektors von Pi_0 und Pi_plus

def successrate(E_K_0, E_K_plus, p, b, g, a, tau, n):               #Zaehlung der Erfolge im Verhaeltnis zu Anzahl K+
    Pi_0 = np.zeros(n)
    Pi_plus = np.zeros(n)
    for i in range(n):
        decay = SimulateDecay(E_K_0, E_K_plus, p, b, g, a, tau)
        Pi_0[i] = decay[0]
        Pi_plus[i] = decay[1]
    success = 0
    for i in range(n):
        if Pi_0[i] <= 2 and Pi_plus[i] <= 2:
            success += 1
    return success/n

def RunExperiment(E_K_0, E_K_plus, p, b, g, a_range, tau, n):       #Plot und optimale Position des Detektors
    A = np.linspace(a_range[0],a_range[1],a_range[2])
    SR = []
    for i in A:
        SR.append(successrate(E_K_0, E_K_plus, p, b, g, i, tau, n))  
    return plt.figure(), plt.plot(A,SR)
    
#Parameter:

E_K_0 = 245.563588 #MeV         #Energie der neutrale Pionen in K+ system
E_K_plus = 248.118174 #MeV      #Energie der positiven Pionen in K+ system
p = 205.14091 #MeV/c            #Impulsbetrag der Pionen (der selbe für beide)
b = 0.99997833784995            #Betafaktor
g = 151.92756392754             #Gammafaktor
tau = 132.97                    #Mittlere Zerfallsstrecke von K+
n = 30                          #Anzahl K+
a_range = [0,500,500]           #Anfangspunkt, Endpunkt, Anzahl Schritte der Positionsbestimmung des Detektors

#Auswertung:

RunExperiment(E_K_0, E_K_plus, p, b, g, a_range, tau, n)
plt.show()
