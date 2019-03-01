from __future__ import division
import numpy as np
import math
from scipy.integrate import odeint, simps
import scipy.optimize
import matplotlib.pyplot as plt
from cfns import wavefn, turningpoints, nodes, opti, itera, normaliser, sqr, simpson, statement

u0 = [0, 1] #array of u and du/dr vals at 0

r1 = np.linspace(0.00001, 15, 1000)

step1 = r1[1] - r1[0] #step size for normalisation

m1 = 1.34
m2 = 1.34
invmu = 1/m1 + 1/m2
mu = 1/invmu

E = 3.068 - m1 - m2

vals = [6.37,0.1]

E1 = 0.1
E2 = 0.5
E3 = 1.0

n = 1
l = 0

fit = scipy.optimize.minimize(opti, vals, args=(E,n,l,E1,E2,E3,u0,mu,r1,step1))
print fit.message

c = fit.x[0]
d = fit.x[1]

prob, psi, dpsi, Ens = itera(n, l, E1, E2, E3, u0, c, d, mu, r1, step1)

norm = 1
print E
print c
print d

f1 = plt.figure(1,figsize=(8,6))
ax = f1.add_subplot(111)
ax.tick_params(axis='x', which='major', labelsize=15)
ax.tick_params(axis='y', which='major', labelsize=15)
plt.plot(r1, psi, 'b',label=r"$E_{10} = 0.388, E =$ %.3f" % Ens)
## probs get gridspec in here for separate plots
plt.plot((r1[0],r1[-1]),(0,0),'grey')

plt.legend(loc=10, fontsize=42, handlelength=0, handletextpad=0)
plt.title(r"Calculation of $\beta$",fontsize=25)
plt.xlabel("Radial Distance, $GeV^{-1}$",fontsize=25)
plt.ylabel("$u_{nl}(r)$",fontsize=25)

plt.show()
