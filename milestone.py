from __future__ import division
import numpy as np
from scipy.integrate import odeint
from scipy.integrate import simps
import matplotlib.pyplot as plt
from myfunctions import wavefn, turningpoints, nodes, itera

u0 = [0, 1] #array of u and du/dr vals at 0

r = np.linspace(0.01, 15, 1000)
step = r[1] - r[0] #step size for normalisation

m1 = 1.34
m2 = 1.34
n = float(input("n = "))
l = float(input("l = "))

if l == 0:
    print "For (n, l) = (%.f, %.f), the number of nodes is %.f, and the number of turning points is %.f" % (n,l,n-1,n) 
invmu = 1/m1 + 1/m2
#mu = 0.512
mu = 1/invmu
#alpha = 1/137
alpha = 0.4
beta = 0
E1 = 1e-2
E2 = 1e-3
E3 = 1e-4

sol1, sol2, sol3 = itera(n, l, E1, E2, E3, u0, alpha, beta, mu, r)

norm1 = 0
norm2 = 0
norm3 = 0
print len(sol1)
for n in range(len(sol1)-1):
    norm1 = norm1 + abs(sol1[n])**2
    norm2 = norm2 + abs(sol2[n])**2
    norm3 = norm3 + abs(sol3[n])**2

prob1 = norm1*step
prob2 = norm2*step
prob3 = norm3*step

sol1 = sol1/prob1
sol2 = sol2/prob2
sol3 = sol3/prob3

maxima, minima, maxs, mins = turningpoints(r,sol2)
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

plt.plot(r, sol1, 'b')
plt.plot(r, sol2, 'r')
plt.plot(r, sol3, 'g')
plt.show()
