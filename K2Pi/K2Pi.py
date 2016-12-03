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

def SimulateOneDecay(E_K_0, E_K_plus, p, b, g, tau):                #Erzeugung eines einzelnen Pion-Paares im Lab-frame
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
    return P_lab_0, P_lab_plus, z_length

def SimulateNDecays(E_K_0, E_K_plus, p, b, g, tau, n):              #Erzeugung von n Zerfaellen:
    P_lab_0 = []
    P_lab_plus = []
    z_length = []
    for i in range(n):
        decay = SimulateOneDecay(E_K_0, E_K_plus, p, b, g, tau)
        P_lab_0.append(decay[0])
        P_lab_plus.append(decay[1])
        z_length.append(decay[2])
    return P_lab_0, P_lab_plus, z_length

def HitDistance(P_lab_0, P_lab_plus, z_length, a):                  #Abstand zu Mittelpunkt des Detektors von Pi_0 und Pi_plus
    ez= np.array([0,0,1])                                           #Einheitsvektor in Z-Richtung
    if float(z_length) >= a:                                        #Aussortieren der K+, die hinter Detektor zerfallen
        return [100,100]
    else:
        d_0 = float((a-z_length) * np.tan(GetAngleBetweenVectors(P_lab_0[1:],ez)))
        d_plus = float((a-z_length) * np.tan(GetAngleBetweenVectors(P_lab_plus[1:],ez)))  
        return [d_0, d_plus]                                       

def successrate(P_lab_0, P_lab_plus, z_length, a, n):               #Zaehlung der Erfolge einer Messung im Verhaeltnis zu Anzahl K+
    success = 0
    for i in range(n):
        if HitDistance(P_lab_0[i], P_lab_plus[i], z_length[i], a)[0] <= 2 and HitDistance(P_lab_0[i], P_lab_plus[i], z_length[i], a)[1] <= 2:
            success += 1
    return success/n

def RunExperiment(E_K_0, E_K_plus, p, b, g, a_range, tau, n):       #optimale Position des Detektors und Plot
    A = np.linspace(a_range[0],a_range[1],a_range[2])
    decay = SimulateNDecays(E_K_0, E_K_plus, p, b, g, tau, n)
    P_lab_0, P_lab_plus, z_length = decay[0], decay[1], decay[2]
    SR = []
    for i in A:
        SR.append(successrate(P_lab_0, P_lab_plus, z_length, i, n))
    SR_max = max(SR)
    a_opt = 0
    for i in range(len(A)):
        if SR[i]==SR_max:
            a_opt = A[i]
    return a_opt, SR, plt.plot(A,SR)


def SimulateKDecayPoint(sx, sy, tau):                               #Funktion gibt Zerfallspunkt eines Kaons aus (x,y,z) und Vektor in Flugrichtung
    alpha = np.array(stats.norm.rvs(loc=0, scale=sx, size=1))       #Erzeugen eines zufälligen Streuwinkels in x Richtung 
    beta = np.array(stats.norm.rvs(loc=0, scale=sy, size=1))        #Erzeugen eines zufälligen Streuwinkels in y Richtung
    theta = np.arccos(np.cos(alpha)*np.cos(beta))                   #Berechne theta und phi aus Streuwinkeln (Wäre froh würde jemand Nachrechnen) 
    phi = np.arctan(np.tan(beta)/np.sin(alpha))                         
    vlen= np.array(stats.expon.rvs(loc=0, scale=tau, size=1))       #Erzeuge Flugläng eines Kaons (Exponentialverteilt, mittlere Flugweite tau) 
    x0= np.sin(theta)*np.cos(phi)*vlen                                 
    y0= np.sin(theta)*np.sin(phi)*vlen                              #kartesische Koordinaten
    z0= np.cos(theta)*vlen
    dp= np.array([[x0],[y0],[z0]]).T                                #Zerfallspunkt (dp=decaypoint) 
    ev = dp/vlen                                                    #Normierter Vektor in Flugrichtung 
    return  dp, ev


    
#Parameter:

E_K_0 = 245.563588 #MeV         #Energie der neutrale Pionen in K+ system
E_K_plus = 248.118174 #MeV      #Energie der positiven Pionen in K+ system
p = 205.14091 #MeV/c            #Impulsbetrag der Pionen (der selbe fuer beide)
b = 0.99997833784995            #Betafaktor
g = 151.92756392754             #Gammafaktor
tau = 132.97                    #Mittlere Zerfallsstrecke von K+
n = 200                         #Anzahl K+
a_range = [0,500,500]           #Anfangspunkt, Endpunkt, Anzahl Messungen
sx = 1*10**-3			#Standardabweichung xWinkel (alpha)
sy = 1*10**-3			#Standardabweichung yWinkel (beta)

#Auswertung:

plt.figure()
a_opt, SR, f = RunExperiment(E_K_0, E_K_plus, p, b, g, a_range, tau, n)
with open("data.txt", "w") as fh:       #Ausgabe der Messdaten in Datei
	fh.write(str(SR))
plt.show()
print('Optimale Position: ', a_opt)
