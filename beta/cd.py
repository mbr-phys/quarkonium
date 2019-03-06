from __future__ import division
import numpy as np
import math
from scipy.integrate import odeint, simps
import scipy.optimize
import matplotlib.pyplot as plt
#from cfns import wavefn, S_hyperfine, opti, sqr, simpson, counter, itera, statement, turningpoints
from cfns import wavefn, sqr, simpson, psi, counter, itera, statement, opti,opti_2

u0 = [0, 1] #array of u and du/dr vals at 0

r1 = np.linspace(0.00001, 15, 1000)

step1 = r1[1] - r1[0] #step size for normalisation

mc = 1.34
invmu = 1/mc + 1/mc
mu = 1/invmu

E = 3.068 - 2*mc

c1 = 7.5
c3 = 10.0
yint = 8.064
d = 0.1

E1 = 0.3
E3 = 0.5

n = 1
l = 0

#fit = scipy.optimize.minimize(opti, d, args=(c,E,n,l,E1,E3,u0,mu,r1,step1))
##fit = scipy.optimize.minimize(opti_2, d, args=(c,E,n,l,u0,mu,r1,step1))
#print fit.message
#print fit.nit
#
#d = fit.x[0]
##c = 6.37
##d = 0.1
#
#print c, d
prob, psi, dpsi, cs = itera(n, l, E, u0, c1, c3, d, yint, mu, r1, step1)
#
norm = 1
statement(psi,dpsi,1,0,r1,E,norm)
print cs
#c = 7.641
#
f1 = plt.figure(1,figsize=(8,6))
ax = f1.add_subplot(111)
ax.tick_params(axis='x', which='major', labelsize=15)
ax.tick_params(axis='y', which='major', labelsize=15)
plt.plot(r1, psi, 'b',label=r"$E_{10} = 0.388, c =$ %.3f" % cs)
## probs get gridspec in here for separate plots
plt.plot((r1[0],r1[-1]),(0,0),'grey')

plt.legend(loc=9, fontsize=30)
plt.title(r"Calculation of $\beta$",fontsize=25)
plt.xlabel("Radial Distance, $GeV^{-1}$",fontsize=25)
plt.ylabel("$u_{nl}(r)$",fontsize=25)

plt.show()

