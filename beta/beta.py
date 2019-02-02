from __future__ import division
import numpy as np
import math
from scipy.integrate import odeint, simps
import matplotlib.pyplot as plt
from bfns import wavefn, turningpoints, nodes, itera, normaliser, sqr, simpson, statement

u0 = [0, 1] #array of u and du/dr vals at 0

r1 = np.linspace(0.00001, 15, 1000)

step1 = r1[1] - r1[0] #step size for normalisation

m1 = 1.34
m2 = 1.34
invmu = 1/m1 + 1/m2
mu = 1/invmu
alpha = 0.4*4/3

E = 3.068 - m1 - m2

beta1 = 0.1
beta2 = 0.5
beta3 = 1.0

n = 1
l = 0

pr1, sol1, du1, beta1 = itera(n, l, E, u0, alpha, beta1, beta2, beta3, mu, r1, step1)


norm = 1
statement(sol1,n,l,r1,beta1,norm)
print E

f1 = plt.figure(1,figsize=(8,6))
plt.plot(r1, sol1, 'b',label=r"(n,l) = (1,0), " r"$\beta =$ %.3f" % beta1)
## probs get gridspec in here for separate plots
plt.plot((r1[0],r1[-1]),(0,0),'grey')

plt.legend(loc=7, fontsize=40)
plt.title(r"Calculation of $\beta$",fontsize='xx-large')
plt.xlabel("Radial Distance, $GeV^{-1}$",fontsize='xx-large')
plt.ylabel("$u_{nl}(r)$",fontsize='xx-large')

#f2 = plt.figure(2,figsize=(8,6))
#plt.plot(r1, pr1, 'b',label=r"(n,l) = (1,0), " r"$\beta =$ %.3f" % beta1)
### probs get gridspec in here for separate plots
#plt.plot((r1[0],r1[-1]),(0,0),'grey')
#
#plt.legend(loc=1, fontsize='x-large')
#plt.title("Probability Densities of Charmonium Wavefunction",fontsize='x-large')
#plt.xlabel("Radial Distance, $GeV^{-1}$",fontsize='x-large')
#plt.ylabel("$|u_{nl}(r)|^2$",fontsize='x-large')

plt.show()
