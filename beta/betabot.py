from __future__ import division
import numpy as np
import math
from scipy.integrate import odeint, simps
import matplotlib.pyplot as plt
from bfns import wavefn, counter, psi, itera, sqr, simpson, statement

u0 = [0, 1] #array of u and du/dr vals at 0

r1 = np.linspace(0.00001, 15, 1000)

step1 = r1[1] - r1[0] #step size for normalisation

m1 = 4.18
m2 = 4.18
invmu = 1/m1 + 1/m2
mu = 1/invmu
alpha = 0.28*4/3

E = 9.445 - m1 - m2

beta1 = 0.1
beta2 = 0.5
beta3 = 1.0

n = 1
l = 0

pr1, sol1, du1, beta1 = itera(n, l, E, u0, alpha, beta1, beta3, mu, r1, step1)


norm = 1
statement(sol1,du1,n,l,r1,beta1,norm)
print E

f1 = plt.figure(1,figsize=(8,6))
ax = f1.add_subplot(111)
ax.tick_params(axis='x', which='major', labelsize=15)
ax.tick_params(axis='y', which='major', labelsize=15)
plt.plot(r1, sol1, 'b',label=r"$E_{10} = 0.388, \beta =$ %.3f" % beta1)
## probs get gridspec in here for separate plots
plt.plot((r1[0],r1[-1]),(0,0),'grey')

plt.legend(loc=10, fontsize=42, handlelength=0, handletextpad=0)
plt.title(r"Calculation of $\beta$",fontsize=25)
plt.xlabel("Radial Distance, $GeV^{-1}$",fontsize=25)
plt.ylabel("$u_{nl}(r)$",fontsize=25)

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
