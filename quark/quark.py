from __future__ import division
import numpy as np
import math
from scipy.integrate import odeint, simps
import matplotlib.pyplot as plt
from qfns import wavefn, turningpoints, nodes, itera, normaliser, sqr, simpson, statement
#from qfns2 import wavefn, wavey, counter, energies, itera, simpson, sqr

u0 = [0, 1] #array of u and du/dr vals at 0

r1 = np.linspace(0.000001, 15, 1000)
r2 = np.linspace(0.0000001, 15, 1000)
r3 = np.linspace(0.000001, 15, 1000)

step1 = r1[1] - r1[0] #step size for normalisation
step2 = r2[1] - r2[0]
step3 = r3[1] - r3[0]

m1 = 1.34
m2 = 1.34
invmu = 1/m1 + 1/m2
mu = 1/invmu
alpha = 0.4*4/3

beta1 = 0.195
beta2 = 0.19512

E1 = 0.1
E2 = 0.5
E3 = 1.0

E11 = 0.8
E12 = 1.0
E13 = 1.2

E01 = 1.0
E02 = 1.5
E03 = 2.0

pr1, u1, du1, E1 = itera(1, 0, E1, E2, E3, u0, alpha, beta1, mu, r1, step1)
pr11, u11, du11, E11 = itera(1, 1, E11, E12, E13, u0, alpha, beta1, mu, r2, step2)
pr20, u20, du20, E20 = itera(2, 0, E01, E02, E03, u0, alpha, beta1, mu, r3, step3)

#print u1
#print u11
#print u20

#r1 = sqr(sol1)
#pr11 = sqr(sol11)
#pr20 = sqr(sol20)

#pr1, u1, v1, norm1 = simpson(pr1, sol1, du1, 1, 0, r1)
#pr11, u11, v11, norm11 = simpson(pr11, sol11, du11, 1, 1, r2)

norm = 1
statement(u1,1,0,r1,E1,norm)
statement(u11,1,1,r2,E11,norm)
statement(u20,2,0,r3,E20,norm)

print E20

f1 = plt.figure(1,figsize=(8,6))
plt.plot(r1, u1, 'b',label="(n,l) = (1,0), $E_{nl} =$ %.3f GeV" % E1)
plt.plot(r2, u11, 'g',label="(n,l) = (1,1), $E_{nl} =$ %.3f GeV" % E11)
plt.plot(r3, u20, 'r',label="(n,l) = (2,0), $E_{nl} =$ %.3f GeV" % E20)
### probs get gridspec in here for separate plots
plt.plot((r2[0],r2[-1]),(0,0),'grey')

plt.legend(loc=1, fontsize=35)
plt.title("Solutions of the Charmonium Wavefunction",fontsize='xx-large')
plt.xlabel("Radial Distance, $GeV^{-1}$",fontsize='xx-large')
plt.ylabel("$u_{nl}(r)$",fontsize='xx-large')

f2 = plt.figure(2,figsize=(8,6))
plt.plot(r1, pr1, 'b',label="(n,l) = (1,0), $E_{nl} =$ %.3f GeV" % E1)
plt.plot(r2, pr11, 'g',label="(n,l) = (1,1), $E_{nl} =$ %.3f GeV" % E11)
plt.plot(r3, pr20, 'r',label="(n,l) = (2,0), $E_{nl} =$ %.3f GeV" % E20)
### probs get gridspec in here for separate plots
plt.plot((r2[0],r2[-1]),(0,0),'grey')

plt.legend(loc=1, fontsize=25)
plt.title("Probability Densities of Charmonium Wavefunction",fontsize='xx-large')
plt.xlabel("Radial Distance, $GeV^{-1}$",fontsize='xx-large')
plt.ylabel("$|u_{nl}(r)|^2$",fontsize='xx-large')

plt.show()
