from __future__ import division
import numpy as np
import math
from scipy.integrate import odeint, simps
import matplotlib.pyplot as plt
from h2_2 import wavefn, sqr, simpson, psi, counter, itera, statement

u0 = [0, 1] #array of u and du/dr vals at 0

r1 = np.linspace(0.001, 5200, 1000)
r2 = np.linspace(0.001, 6200, 1000)
r3 = np.linspace(0.001, 9500, 1000)

step1 = r1[1] - r1[0] #step size for normalisation
step2 = r2[1] - r2[0]
step3 = r3[1] - r3[0]
mu = 0.511
alpha = 1/137
beta = 0

E11 = 1e-5
E13 = 2e-5

E01 = 3e-6
E03 = 3.5e-6

E21 = 3e-6
E23 = 3.5e-6

E31 = 1.25e-6
E33 = 1.75e-6

pr1, sol1, du1, En1 = itera(1, 0, E11, E13, u0, alpha, beta, mu, r1, step1)
pr20, sol20, du20, En20 = itera(2, 0, E01, E03, u0, alpha, beta, mu, r2, step2)
pr21, sol21, du21, En21 = itera(2, 1, E21, E23, u0, alpha, beta, mu, r2, step2)
pr30, sol30, du30, En30 = itera(3, 0, E31, E33, u0, alpha, beta, mu, r3, step3)
pr31, sol31, du31, En31 = itera(3, 1, E31, E33, u0, alpha, beta, mu, r3, step3)
pr32, sol32, du32, En32 = itera(3, 2, E31, E33, u0, alpha, beta, mu, r3, step3)

norm = 1
statement(sol1,du1,1,0,r1,En1,norm)
statement(sol20,du20,2,0,r2,En20,norm)
statement(sol21,du21,2,1,r2,En21,norm)
statement(sol30,du30,3,0,r3,En30,norm)
statement(sol31,du31,3,1,r3,En31,norm)
statement(sol32,du32,3,2,r3,En32,norm)

En1 = En1*-1e6
En20 = En20*-1e6
En21 = En21*-1e6
En30 = En30*-1e6
En31 = En31*-1e6
En32 = En32*-1e6

f1 = plt.figure(1,figsize=(8,6))
ax = f1.add_subplot(111)
ax.tick_params(axis='x', which='major', labelsize=25)
ax.tick_params(axis='y', which='major', labelsize=25)
plt.plot(r1, sol1, 'b',label="(n,l) = (1,0), $E_{nl} =$ %.1f eV" % En1)
plt.plot(r2, sol20, 'r',label="(n,l) = (2,0), $E_{nl} =$ %.2f eV" % En20)
plt.plot(r2, sol21, 'g',label="(n,l) = (2,1), $E_{nl} =$ %.2f eV" % En21)
plt.plot(r3, sol30, 'purple',label="(n,l) = (3,0), $E_{nl} =$ %.2f eV" % En30)
plt.plot(r3, sol31, 'cyan',label="(n,l) = (3,1), $E_{nl} =$ %.2f eV" % En31)
plt.plot(r3, sol32, 'm',label="(n,l) = (3,2), $E_{nl} =$ %.2f eV" % En32)
## probs get gridspec in here for separate plots
plt.plot((r3[0],r3[-1]),(0,0),'grey')

plt.legend(loc=1, fontsize=35)
plt.title("Solutions of the Hydrogen Electron Radial Wavefunction",fontsize=25)
plt.xlabel("Radial Distance, $MeV^{-1}$",fontsize=25)
plt.ylabel("$u_{nl}(r)$",fontsize=25)

#f2 = plt.figure(2,figsize=(8,6))
#ax1 = f1.add_subplot(111)
#ax1.tick_params(axis='x', which='major', labelsize=15)
#ax1.tick_params(axis='y', which='major', labelsize=15)
#plt.plot(r1, pr1, 'b',label="(n,l) = (1,0), $E_{nl} =$ %.1f eV" % En1)
#plt.plot(r2, pr20, 'r',label="(n,l) = (2,0), $E_{nl} =$ %.2f eV" % En20)
#plt.plot(r2, pr21, 'g',label="(n,l) = (2,1), $E_{nl} =$ %.2f eV" % En21)
### probs get gridspec in here for separate plots
#plt.plot((r2[0],r2[-1]),(0,0),'grey')
#
#plt.legend(loc=1, fontsize=25)
#plt.title("Probability Densities of the Hydrogen Electron Radial Wavefunction",fontsize=25)
#plt.xlabel("Radial Distance, $MeV^{-1}$",fontsize=25)
#plt.ylabel("$|u_{nl}(r)|^2$",fontsize=25)

plt.show()
