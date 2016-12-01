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

E_K0 =807
p_0= 301

E_K_plus = 0
p_plus= 0 

P_K0=np.array([E_K0,p_0*x0,p_0*y0,p_0*z0])
P_K_plus=np.array([E_K_plus,p_plus*(-x0),p_plus*(-y0),p_plus*(-z0)])

b=0.8
g=5/3.
boost=np.array([[g,0,0,b*g],[0,1,0,0],[0,0,1,0],[b*g,0,0,g]])

P_lab0=np.dot(boost,P_K0.T)


"""
Zerfall von K


"""

a = 280
tau=132.97
z_length= np.array(stats.expon.rvs(loc=0, scale=tau, size=1))


def GetAngleBetweenVectors(X, Y):
    x=np.linalg.norm(X)
    y=np.linalg.norm(Y)
    return np.arccos(np.dot(X,Y)/(x*y))

P = P_lab0[1:]
ez= np.array([0,0,1])

d =(a-z_length) * np.tan(GetAngleBetweenVectors(P,ez))
print(d)
