
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
from threading import Thread
import multiprocessing
from multiprocessing.pool import ThreadPool
import threading

#Funktionen:

def GetAngleBetweenVectors(X, Y):                                   #nicht orientierter Winkel zwischen zwei Vektoren
    x=np.linalg.norm(X)
    y=np.linalg.norm(Y)
    return np.arccos(np.dot(X,Y)/(x*y))

def RotationAroundYAxis(alpha):                                     #Rotationsmatrix, alpha um y Achse
    return np.array([[np.cos(alpha),0,np.sin(alpha)],[0,1,0],[-np.sin(alpha),0,np.cos(alpha)]])    

def RotationAroundXAxis(beta):                                      #Rotationsmatrix, beta um x Achse
    return np.array([[1,0,0],[0,np.cos(beta),-np.sin(beta)],[0,np.sin(beta),np.cos(beta)]])

def LorentzBoost(b,g):                                              #Lorentz-Boost, b = Beta-Faktor, g = Gamma-Faktor
    return np.array([[g,0,0,b*g],[0,1,0,0],[0,0,1,0],[b*g,0,0,g]])

def SimulateKDecayPoint(sx, sy, tau):                               #Funktion gibt Zerfallspunkt eines Kaons aus (x,y,z) und Vektor in Flugrichtung
    if sx==0 and sy==0:                                             #sx, sy = 0 erzeugen einen reinen Strahl in z-Richtung
        alpha = 0
        beta = 0
    else:
        alpha = np.array(stats.norm.rvs(loc=0, scale=sx, size=1))   #Erzeugen eines zufaelligen Streuwinkels in x Richtung 
        beta = np.array(stats.norm.rvs(loc=0, scale=sy, size=1))    #Erzeugen eines zufaelligen Streuwinkels in y Richtung
    vlen= np.array(stats.expon.rvs(loc=0, scale=tau, size=1))       #Erzeuge Fluglaenge eines Kaons (Exponentialverteilt, mittlere Flugweite tau) 
    x0= vlen*np.tan(alpha)*np.cos(beta)/np.sqrt(1+np.tan(alpha)**2 *np.cos(beta)**2)                                
    y0= np.sqrt(vlen**2 -(x0**2))*np.sin(beta)                      #kartesische Koordinaten
    z0= np.sqrt(vlen**2 -(x0**2))*np.cos(beta)
    dp= np.array([x0,y0,z0])                                        #Zerfallspunkt (dp=decaypoint) 
    ev = dp/vlen                                                    #Normierter Vektor in Flugrichtung 
    return  dp, ev, alpha, beta

def SimulateK2PiDecay(E_K_0, E_K_plus, p, b, g, tau):               #Erzeugung eines einzelnen Pion-Paares im Lab-frame
    theta = np.array(stats.uniform.rvs(scale=np.pi, size=1))        #gleichverteilt zwischen 0 und Pi
    phi = np.array(stats.uniform.rvs(scale=2*np.pi, size=1))        #gleichverteilt zwischen 0 und 2 Pi
    x0= np.sin(theta)*np.cos(phi)                                   
    y0= np.sin(theta)*np.sin(phi)                                   #kartesische Koordinaten
    z0= np.cos(theta)   
    P_K_0 = np.array([E_K_0,p*x0,p*y0,p*z0])                        #4-Vektoren im K+-Frame
    P_K_plus = np.array([E_K_plus,p*(-x0),p*(-y0),p*(-z0)])
    P_lab_0 = np.dot(LorentzBoost(b,g),P_K_0.T)                     #4-Vektoren im Lab-Frame
    P_lab_plus = np.dot(LorentzBoost(b,g),P_K_plus.T)
    return P_lab_0, P_lab_plus                                

def RotateDecayVectors(sx, sy, tau, E_K_0, E_K_plus, p, b, g):                 #Erzeugt in K richtung Rotierte Zerfallsvektoren
    dp,ev,alpha,beta = SimulateKDecayPoint(sx, sy, tau)                        #Rotationswinkel    
    P_lab_0,P_lab_plus = SimulateK2PiDecay(E_K_0, E_K_plus, p, b, g, tau)      #Normale Zerfallsvektore
    P_lab_0=P_lab_0[1:]
    P_lab_plus=P_lab_plus[1:]    
    P_lab_0r = np.dot(RotationAroundXAxis(beta),np.dot(RotationAroundYAxis(alpha),P_lab_0.T))   #Gedrehte Zerfallsvektoren  
    P_lab_plusr = np.dot(RotationAroundXAxis(beta),np.dot(RotationAroundYAxis(alpha),P_lab_plus.T))
    return P_lab_0r, P_lab_plusr, dp
    
def SimulateNDecays(sx, sy, tau, E_K_0, E_K_plus, p, b, g, n):          #Erzeugt n Zerfaelle
    P_lab_0 = []
    P_lab_plus = []
    dp = []
    for i in range(n):
        decay = RotateDecayVectors(sx, sy, tau, E_K_0, E_K_plus, p, b, g)
        P_lab_0.append(decay[0])
        P_lab_plus.append(decay[1])
        dp.append(decay[2])
    return P_lab_0, P_lab_plus, dp                                      #Ausgabe: Listen von 4-Vektoren und Ortsvektoren des K+-Zerfalls

def HitDistance(P_lab_0, P_lab_plus, dp, a):                            #Abstand zu Mittelpunkt des Detektors von Pi_0 und Pi_plus                                           
    if float(dp[-1]) >= a:                                              #Aussortieren der K+, die hinter Detektor zerfallen
        return [100,100]
    else:
        n_0 = float((a-dp[-1])/float(P_lab_0[-1]))                      #Berechne wieviel mal P_lab_0r an dp angehängt werden muss damit z=a
        n_plus = float((a-dp[-1])/float(P_lab_plus[-1]))
        d_0 = np.sqrt((float(dp[0])+float(n_0*P_lab_0[0]))**2+(float(dp[1])+float(n_0*P_lab_0[1]))**2)  #Berechnet Abstand zu z Achse (r=(x^2+y^2)^(1/2))
        d_plus = np.sqrt((float(dp[0])+float(n_plus*P_lab_plus[0]))**2+(float(dp[1])+float(n_plus*P_lab_plus[1]))**2)  
        return [d_0, d_plus] 

def successrate(P_lab_0, P_lab_plus, dp, a, n):                         #Zaehlung der Erfolge einer Messung im Verhaeltnis zu Anzahl K+
    success = 0
    for i in range(n):
        if HitDistance(P_lab_0[i], P_lab_plus[i], dp[i], a)[0] <= 2 and HitDistance(P_lab_0[i], P_lab_plus[i], dp[i], a)[1] <= 2:
            success += 1
    return success/n

def GraficEvaluation(a_opt, SR_max, A, SR):			        #huebsche Darstellung der Messwerte
    plt.figure()
    plt.plot(A,SR)
    plt.xlim(xmin=a_range[0],xmax=a_range[1])
    plt.ylim(ymin=0, ymax=1)
    plt.plot([a_range[0],a_range[1]], [max(SR),max(SR)],'k:')
    plt.plot([a_opt,a_opt], [0,1],'k:')
    plt.xlabel('detector position [m]')
    plt.ylabel(r'Successrate [success/$n_{K+}$]')
    plt.show()
    
def RunExperiment(sx, sy, E_K_0, E_K_plus, p, b, g, a_range, tau, n, SR, i):   #Ausfuehrung des Experiments

    A = np.linspace(*a_range)
    decay = SimulateNDecays(sx, sy, tau, E_K_0, E_K_plus, p, b, g, n)

    P_lab_0, P_lab_plus, dp = decay[0], decay[1], decay[2]
    index = 0
    for a in A:
        if SR[i][index] == 0:
            SR[i][index] = successrate(P_lab_0, P_lab_plus, dp, a, n)
        else:
            SR[i][index] = np.mean([SR[i][index], successrate(P_lab_0, P_lab_plus, dp, a, n)])
        index = index + 1
   


def RunExperimentMultiThreaded(sx, sy, E_K_0, E_K_plus, p, b, g, a_range, tau, n, enableMultiThreading):

    if (enableMultiThreading):
        number_of_cores = multiprocessing.cpu_count()
        number_of_threads = number_of_cores
        print("Running on %s cores, using %s threads" % (number_of_cores, number_of_threads))

    decay = []
    A = np.linspace(*a_range)
    threads = []    
    SR = [ [0] * (len(A)) ] * number_of_cores

    # Use Multithreading to start as many threads as we have cores on this machine
    if (enableMultiThreading):
        # The amount of simulation runs every thread should do is n / number_of_threads
        nThread = int(n / number_of_threads)

        for i in range(number_of_cores):
            # Create mew thread
            t = Thread(target=RunExperiment, args = (sx, sy, E_K_0, E_K_plus, p, b, g, a_range, tau, nThread, SR, i))
            t.start();
            threads.append(t);

    # If Multithreading turned off, just use the "normal" main thread
    else:
        SR = [ [0] * (len(A)) ] * 1
        RunExperiment(sx, sy, E_K_0, E_K_plus, p, b, g, a_range, tau, n, SR, 0)    

    # Wait for all threads to finish
    for i in range(len(threads)):
        if (threads[i].isAlive()):
            threads[i].join()
        
    avgSR = [0] * (len(A));
    for i in range (len(A)):
        values = []
        for j in range(number_of_cores):
            values.append(SR[j][i])
        avgSR[i] = np.mean(values)

    SR_max = max(avgSR)
    a_opt = 0
    for i in range(len(A)):
        if avgSR[i]==SR_max:
            a_opt = A[i]

    print('Optimale Position: ', a_opt)
    print('Maximale Erfolgsrate: ', SR_max)

    return a_opt, SR_max, SR, GraficEvaluation(a_opt, SR_max, A, avgSR)    #Ausgabe: optimale Detektorposition, maximale Erfolgsrate, Messdaten, Plot
    

#Parameter:

E_K_0 = 245.5611565 #MeV            #Energie der neutrale Pionen in K+ system
E_K_plus = 248.115779 #MeV          #Energie der positiven Pionen in K+ system
p = 205.138 #MeV/c                  #Impulsbetrag der Pionen (der selbe fuer beide)
b = 0.999978336972366               #Betafaktor
g = 151.924486603                   #Gammafaktor
tau = 576.21338881                  #Mittlere Zerfallsstrecke von K+
uncertainty_tau=2.4
n = 100000                          #Anzahl K+
a_range = [0,500,1000]  		    #Anfangspunkt, Endpunkt, Anzahl Messungen
sx = 1e-3			                #Standardabweichung xWinkel (alpha)
sy = 1e-3   			            #Standardabweichung yWinkel (beta)

enableMultiThreading = True         # Set to true to enable multithraeding, false to disable it


#Auswertung:
a_opt, SR_max, SR, f = RunExperimentMultiThreaded(sx, sy, E_K_0, E_K_plus, p, b, g, a_range, tau, n, enableMultiThreading)
with open("data.txt", "w") as fh:               #Ausgabe der Messdaten in Datei
	fh.write(str(SR))
