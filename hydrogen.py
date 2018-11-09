from __future__ import division
import numpy as np
import math
from scipy.integrate import odeint
from scipy.integrate import simps
import matplotlib.pyplot as plt
from myfunctions import wavefn, turningpoints, nodes, itera, normaliser

u0 = [0, 1] #array of u and du/dr vals at 0

r = np.linspace(0.01, 2900, 1000)
step = r[1] - r[0] #step size for normalisation

n = float(input("n = "))
l = float(input("l = "))
mu = 0.511
alpha = 1/137
beta = 0
if n == 1:
    E1 = 1e-5
    E2 = 2.5e-5
    E3 = 5e-5
elif n == 2:
    E1 = 1e-6
    E2 = 2.5e-6
    E3 = 5e-6
elif n >= 3:
    E1 = 0.5e-6
    E2 = 1.25e-6
    E3 = 2e-6
else:
    E1 = 1e-7
    E2 = 5e-7
    E3 = 1e-6

sol1, sol2, sol3, du1, du2, du3, En1, En2, En3 = itera(n, l, E1, E2, E3, u0, alpha, beta, mu, r, step)

maxima, minima, maxs, mins = turningpoints(r,sol3,n,l)
for z in range(len(maxima[:,0])):
    print "Maxima No. %.f is found at r = %.2f, u = %.5e" % (z+1,maxima[z,0],maxima[z,1])

for t in range(len(minima[:,0])):
    print "Minima No. %.f is found at r = %.2f, u = %.5e" % (t+1,minima[t,0],minima[t,1])

noddy, numbs = nodes(maxima, minima, maxs, mins, sol3, r)
if numbs != 0:
    numbers = len(noddy)
    for ns in range(numbers):
        print "Node No. %.f is found at r = %.2f, u = %.5e" % (ns+1,noddy[ns,0], noddy[ns,1])

print En1
print En2
print En3 

plt.plot(r, sol1, 'b')
plt.plot(r, sol2, 'r')
plt.plot(r, sol3, 'g')
## probs get gridspec in here for separate plots
plt.plot((r[0],r[-1]),(0,0),'grey')
#plt.plot(r,du1,'orange')
#plt.plot(r,du2,'purple')
#plt.plot(r,du3,'grey')
plt.show()
