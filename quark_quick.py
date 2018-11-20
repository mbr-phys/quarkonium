from __future__ import division
import numpy as np
import math
from scipy.integrate import odeint, simps
import matplotlib.pyplot as plt
from hfns import wavefn, turningpoints, nodes, itera, normaliser, statement, sqr, simpson

u0 = [0, 1] #array of u and du/dr vals at 0

#n = float(input("n = "))
#l = float(input("l = "))
r1 = np.linspace(0.001, 15, 1000)

step1 = r1[1] - r1[0] #step size for normalisation
m1 = 1.34
m2 = 1.34
invmu = 1/m1 + 1/m2
mu = 1/invmu
alpha = 0.4
beta = 0.01
E = 0.08
l = 0

a = l*(l+1)
b = 2*mu*E
c = 2*mu*4*alpha/3
d = 2*mu*beta

sol1 = odeint(wavefn, u0, r1, args=(a,b,c,d))[:,0]

pr1 = sqr(sol1)

f1 = plt.figure(1,figsize=(8,6))
plt.plot(r1, sol1, 'b',label="(n,l) = (1,0), $E_{nl} =$ %.1e eV" % E)
## probs get gridspec in here for separate plots
plt.plot((r1[0],r1[-1]),(0,0),'grey')

plt.legend(loc=1, fontsize='x-large')
plt.title("First bound state of Quarkonium",fontsize='x-large')
plt.xlabel("Radial Distance, $GeV^{-1}$",fontsize='x-large')
plt.ylabel("$u_{nl}(r)$",fontsize='x-large')

f2 = plt.figure(2,figsize=(8,6))
plt.plot(r1, pr1, 'b',label="(n,l) = (1,0), $E_{nl} =$ %.1e eV" % E)
## probs get gridspec in here for separate plots
plt.plot((r1[0],r1[-1]),(0,0),'grey')

plt.legend(loc=1, fontsize='x-large')
plt.title("Probability Density",fontsize='x-large')
plt.xlabel("Radial Distance, $GeV^{-1}$",fontsize='x-large')
plt.ylabel("$|u_{nl}(r)|^2$",fontsize='x-large')

plt.show()
