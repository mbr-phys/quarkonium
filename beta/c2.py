from __future__ import division
import numpy as np
import math
from scipy.integrate import odeint, simps
import scipy.optimize
import matplotlib.pyplot as plt
#from cfns import wavefn, S_hyperfine, opti, sqr, simpson, counter, itera, statement, turningpoints
from cf2 import wavefn, sqr, simpson, psi, counter, itera, statement, opti

u0 = [0, 1] #array of u and du/dr vals at 0

r1 = np.linspace(0.00001, 15, 1000)

step1 = r1[1] - r1[0] #step size for normalisation

mc = 1.34
invmu = 1/mc + 1/mc
mu = 1/invmu

E = 3.068 - 2*mc

vals = [7.00,8.00]
d = 0.1

E1 = 0.3
E3 = 0.5

n = 1
l = 0

fit = scipy.optimize.minimize(opti, vals, args=(E,n,l,E1,E3,u0,d,mu,r1,step1))
print fit.message
print fit.nit

c = fit.x[0]
d = fit.x[1]
#c = 6.37
#d = 0.1

print c, d
#prob, psi, dpsi, Ens = itera(n, l, E1, E3, u0, c, d, mu, r1, step1)
##
#norm = 1
#statement(psi,dpsi,1,0,r1,Ens,norm)
#print Ens
##
#f1 = plt.figure(1,figsize=(8,6))
#ax = f1.add_subplot(111)
#ax.tick_params(axis='x', which='major', labelsize=15)
#ax.tick_params(axis='y', which='major', labelsize=15)
#plt.plot(r1, psi, 'b',label=r"$E_{10} = 0.388, E =$ %.3f" % Ens)
### probs get gridspec in here for separate plots
#plt.plot((r1[0],r1[-1]),(0,0),'grey')
#
#plt.legend(loc=9, fontsize=30)
#plt.title(r"Calculation of $\beta$",fontsize=25)
#plt.xlabel("Radial Distance, $GeV^{-1}$",fontsize=25)
#plt.ylabel("$u_{nl}(r)$",fontsize=25)
#
#plt.show()
#
