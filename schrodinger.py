from __future__ import division
import numpy as np
from scipy.integrate import odeint
from scipy.integrate import simps
import matplotlib.pyplot as plt

def wavefn(u0, r, a, b, c, d):
    u, v = u0
    dv = [v, (a/r**2)*u - b*u - (c/r)*u + d*r*u]
    return dv

def turningpoints(r, u):
    mins = 0
    min_loc = []
    maxs = 0
    max_loc = []
    for i in range(len(u)-2):
        if u[i] < u[i+1] and u[i+2] < u[i+1]:
            maxs += 1
            maxi = i+1
            max_loc = np.append(max_loc,maxi)
        elif u[i] > u[i+1] and u[i+2] > u[i+1]:
            mins += 1
            mini = i+1
            min_loc = np.append(min_loc,mini)

    maxes = np.zeros((len(max_loc),2))
    for x in range(len(max_loc)):
        hunt = int(max_loc[x])
        maxes[x,0] = r[hunt]
        maxes[x,1] = u[hunt]

    minies = np.zeros((len(min_loc),2))
    for y in range(len(min_loc)):
        hunt = int(min_loc[y])
        minies[y,0] = r[hunt]
        minies[y,1] = u[hunt]

    print "The number of maxima is %.f" % maxs
    print "The number of minima is %.f" % mins

    return maxes, minies

u0 = [0, 1] #array of u and du/dr vals at 0

r = np.linspace(0.01, 15, 1000)
step = r[1] - r[0] #step size for normalisation

m1 = 1.34
m2 = 1.34
l = 0
invmu = 1/m1 + 1/m2
mu = 1/invmu
alpha = 0.4
be1 = 1e-02
be2 = 1.5e-02
be3 = 5e-03
E = 3.0680 - 2.68

a = l*(l+1)
b = 2*mu*E
c = (8/3)*mu*alpha
d1 = 2*mu*be1
d2 = 2*mu*be2
d3 = 2*mu*be3

sol1 = odeint(wavefn, u0, r, args=(a, b, c, d1)) #unl
#sol2 = odeint(wavefn, u0, r, args=(a, b, c, d2)) #unl
#sol3 = odeint(wavefn, u0, r, args=(a, b, c, d3)) #unl

norm = 0
for n in range(len(sol1)):
    norm = norm + abs(sol1[n,0])**2

prob = norm*step
print prob

sol1 = sol1/prob

maxima, minima = turningpoints(r,sol1[:,0])
for z in range(len(maxima[:,0])):
    print "Maxima No. %.f is found at r = %.2f, u = %.5f" % (z+1,maxima[z,0],maxima[z,1])

for t in range(len(minima[:,0])):
    print "Minima No. %.f is found at r = %.2f, u = %.5f" % (t+1,minima[t,0],minima[t,1])

#print r
#print sol1

plt.plot(r, sol1[:,0], 'b')
plt.plot(r, sol1[:,1], 'r')
#plt.plot(r, turns, 'g')
#plt.plot(r, sol3, 'g')
plt.show()


