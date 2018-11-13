from __future__ import division
import numpy as np
import math
from scipy.integrate import odeint
from scipy.integrate import simps
import matplotlib.pyplot as plt
from myfunctions import wavefn, turningpoints, nodes, itera, normaliser, statement

u0 = [0, 1] #array of u and du/dr vals at 0

#n = float(input("n = "))
#l = float(input("l = "))
r1 = np.linspace(0.001, 2900, 1000)
r2 = np.linspace(0.001, 4500, 1000)

step1 = r1[1] - r1[0] #step size for normalisation
step2 = r2[1] - r2[0]
mu = 0.511
alpha = 1/137
beta = 0
#if n == 1:
E11 = 1e-5
E12 = 2.5e-5
E13 = 5e-5
#elif n == 2 and l == 0:
E01 = 3e-6
E02 = 3.5e-6
E03 = 4e-6
#elif n == 2 and l == 1:
E21 = 3e-6
E22 = 3.25e-6
E23 = 3.5e-6
#elif n == 3:
#    E1 = 1e-6
#    E2 = 1.5e-6
#    E3 = 2e-6
#elif n > 3 and n < 6:
#    E1 = 0.5e-6
#    E2 = 1.25e-6
#    E3 = 2e-6
#else:
#    E1 = 1e-7
#    E2 = 5e-7
#    E3 = 1e-6

sol_1, du_1, En_1, norm_1 = itera(1, 0, E11, E12, E13, u0, alpha, beta, mu, r1, step1)
sol_20, du_20, En_20, norm_20 = itera(2, 0, E01, E02, E03, u0, alpha, beta, mu, r2, step2)
sol_21, du_21, En_21, norm_21 = itera(2, 1, E21, E22, E23, u0, alpha, beta, mu, r2, step2)

statement(sol_1,1,0,r1,En_1, norm_1)
statement(sol_20,2,0,r2,En_20, norm_20)
statement(sol_21,2,1,r2,En_21,norm_21)

sol_21 = sol_21/(norm_1*20)

En_1 = En_1*-1e6
En_20 = En_20*-1e6
En_21 = En_21*-1e6

plt.plot(r1, sol_1, 'b',label="(n,l) = (1,0), $E_{nl} =$ %.1f eV" % En_1)
plt.plot(r2, sol_20, 'r',label="(n,l) = (2,0), $E_{nl} =$ %.2f eV" % En_20)
plt.plot(r2, sol_21, 'g',label="(n,l) = (2,1), $E_{nl} =$ %.2f eV" % En_21)
## probs get gridspec in here for separate plots
plt.plot((r2[0],r2[-1]),(0,0),'grey')

plt.legend(loc=1, fontsize='x-large')
plt.title("Solutions of the Hydrogen Electron Radial Wavefunction",fontsize='x-large')
plt.xlabel("Radial Distance, $MeV^{-1}$",fontsize='x-large')
plt.ylabel("$u_{nl}(r)$",fontsize='x-large')
plt.show()



#plt.plot(r,du1,'orange')
#plt.plot(r,du2,'purple')
#plt.plot(r,du3,'grey')
