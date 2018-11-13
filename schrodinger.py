from __future__ import division
import numpy as np
from scipy.integrate import odeint
from scipy.integrate import simps
import matplotlib.pyplot as plt
from myfunctions import wavefn, turningpoints, nodes, itera, normaliser

u0 = [0, 1] #array of u and du/dr vals at 0

r = np.linspace(0.0001, 12, 1000)
step = r[1] - r[0] #step size for normalisation

m1 = 1.34
m2 = 1.34
n = float(input("n = "))
l = float(input("l = "))
invmu = 1/m1 + 1/m2
mu = 1/invmu
alpha = 0.4*4/3
beta = 0.01
E1 = 0.01
E2 = 0.015
E3 = 0.02

sol1, sol2, sol3, du1, du2, du3, En1, En2, En3 = itera(n, l, E1, E2, E3, u0, alpha, beta, mu, r,step)

maxima, minima, maxs, mins = turningpoints(r,sol2,n,l)
for z in range(len(maxima[:,0])):
    print "Maxima No. %.f is found at r = %.2f, u = %.5f" % (z+1,maxima[z,0],maxima[z,1])

for t in range(len(minima[:,0])):
    print "Minima No. %.f is found at r = %.2f, u = %.5f" % (t+1,minima[t,0],minima[t,1])

noddy, numbs = nodes(maxima, minima, maxs, mins, sol2, r)
if numbs == 0:
    print "There are no nodes for this function."
else:
    numbers = len(noddy)
    for ns in range(numbers):
        print "Node No. %.f is found at r = %.2f, u = %.5f" % (ns+1,noddy[ns,0], noddy[ns,1])

print En1
print En2
print En3 

plt.plot(r, sol1, 'b')
plt.plot((r[0],r[-1]),(0,0),'black')
plt.plot(r, sol2, 'r')
plt.plot(r, sol3, 'g')
plt.show()
