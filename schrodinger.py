from __future__ import division
import numpy as np
from scipy.integrate import odeint
from scipy.integrate import simps
import matplotlib.pyplot as plt
from myfunctions import wavefn, turningpoints, nodes

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
E = 0.0010

a = l*(l+1)
b = 2*mu*E
c = (8/3)*mu*alpha
#c = 2*mu/137
d1 = 2*mu*be1
d2 = 2*mu*be2
d3 = 2*mu*be3

sol1 = odeint(wavefn, u0, r, args=(a, b, c, d1)) #unl
sol2 = odeint(wavefn, u0, r, args=(a, b, c, d2)) #unl
sol3 = odeint(wavefn, u0, r, args=(a, b, c, d3)) #unl

norm = 0
for n in range(len(sol1)):
    norm = norm + abs(sol1[n,0])**2

prob = norm*step
print prob

sol1 = sol1/prob

soli1 = sol1[:,0]
maxima, minima, maxs, mins = turningpoints(r,soli1)
for z in range(len(maxima[:,0])):
    print "Maxima No. %.f is found at r = %.2f, u = %.5f" % (z+1,maxima[z,0],maxima[z,1])

for t in range(len(minima[:,0])):
    print "Minima No. %.f is found at r = %.2f, u = %.5f" % (t+1,minima[t,0],minima[t,1])

noddy, numbs = nodes(maxima, minima, maxs, mins, soli1, r)
if numbs == 0:
    print "There are no nodes for this function."
else:
    numbers = len(noddy)
    for ns in range(numbers):
        print "Node No. %.f is found at r = %.2f, u = %.5f" % (ns+1,noddy[ns,0], noddy[ns,1])

plt.plot(r, sol1[:,0], 'b')
#plt.plot(r, sol1[:,1], 'r')
plt.plot((r[0],r[-1]),(0,0),'black')
#plt.plot(r, turns, 'g')
plt.plot(r, sol2[:,0], 'r')
plt.plot(r, sol3[:,0], 'g')
plt.show()
