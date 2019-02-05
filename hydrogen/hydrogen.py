from __future__ import division
import numpy as np
import math
from scipy.integrate import odeint, simps
import matplotlib.pyplot as plt
from hfns import wavefn, turningpoints, nodes, itera, normaliser, sqr, simpson, statement

u0 = [0, 1] #array of u and du/dr vals at 0

r1 = np.linspace(0.001, 5200, 1000)
r2 = np.linspace(0.001, 5200, 1000)

step1 = r1[1] - r1[0] #step size for normalisation
step2 = r2[1] - r2[0]
mu = 0.511
alpha = 1/137
beta = 0

E11 = 1e-5
E12 = 2.5e-5
E13 = 5e-5

E01 = 3e-6
E02 = 3.5e-6
E03 = 4e-6

E21 = 3e-6
E22 = 3.25e-6
E23 = 3.5e-6

sol_1, du1, En_1 = itera(1, 0, E11, E12, E13, u0, alpha, beta, mu, r1, step1)
sol_20, du20, En_20 = itera(2, 0, E01, E02, E03, u0, alpha, beta, mu, r2, step2)
sol_21, du21, En_21 = itera(2, 1, E21, E22, E23, u0, alpha, beta, mu, r2, step2)

pr1 = sqr(sol_1)
pr20 = sqr(sol_20)
pr21 = sqr(sol_21)

pr1,sol1,du1,norm1 = simpson(pr1,sol_1,du1,1,0,r1)
pr20,sol20,du20,norm20 = simpson(pr20,sol_20,du20,2,0,r2)
pr21,sol21,du21,norm21 = simpson(pr21,sol_21,du21,2,1,r2)

statement(sol_1,1,0,r1,En_1,norm1)
statement(sol_20,2,0,r2,En_20,norm20)
statement(sol_21,2,1,r2,En_21,norm21)

En_1 = En_1*-1e6
En_20 = En_20*-1e6
En_21 = En_21*-1e6

f1 = plt.figure(1,figsize=(8,6))
plt.plot(r1, sol1, 'b',label="(n,l) = (1,0), $E_{nl} =$ %.1f eV" % En_1)
plt.plot(r2, sol20, 'r',label="(n,l) = (2,0), $E_{nl} =$ %.2f eV" % En_20)
plt.plot(r2, sol21, 'g',label="(n,l) = (2,1), $E_{nl} =$ %.2f eV" % En_21)
## probs get gridspec in here for separate plots
plt.plot((r2[0],r2[-1]),(0,0),'grey')

plt.legend(loc=1, fontsize=35)
plt.title("Solutions of the Hydrogen Electron Radial Wavefunction",fontsize='xx-large')
plt.xlabel("Radial Distance, $MeV^{-1}$",fontsize='xx-large')
plt.ylabel("$u_{nl}(r)$",fontsize='xx-large')

f2 = plt.figure(2,figsize=(8,6))
plt.plot(r1, pr1, 'b',label="(n,l) = (1,0), $E_{nl} =$ %.1f eV" % En_1)
plt.plot(r2, pr20, 'r',label="(n,l) = (2,0), $E_{nl} =$ %.2f eV" % En_20)
plt.plot(r2, pr21, 'g',label="(n,l) = (2,1), $E_{nl} =$ %.2f eV" % En_21)
## probs get gridspec in here for separate plots
plt.plot((r2[0],r2[-1]),(0,0),'grey')

plt.legend(loc=1, fontsize=25)
plt.title("Probability Densities of the Hydrogen Electron Radial Wavefunction",fontsize='xx-large')
plt.xlabel("Radial Distance, $MeV^{-1}$",fontsize='xx-large')
plt.ylabel("$|u_{nl}(r)|^2$",fontsize='xx-large')

plt.show()
