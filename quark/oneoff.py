from __future__ import division
import numpy as np
import math
from scipy.integrate import odeint, simps
import matplotlib.pyplot as plt
from qfns import wavefn, turningpoints, nodes, itera, normaliser, sqr, simpson, statement

u0 = [0, 1] #array of u and du/dr vals at 0

r3 = np.linspace(0.000001, 12, 1000)

step3 = r3[1] - r3[0]

m1 = 1.34
m2 = 1.34
invmu = 1/m1 + 1/m2
mu = 1/invmu
alpha = 0.4*4/3

beta1 = 0.195

#E1 = 1.03670047557
E1 = 1.0
E2 = 1.5
E3 = 2.0

n = 2
l = 0

a = l*(l+1)
b1 = 2*mu*E1
b2 = 2*mu*E2
b3 = 2*mu*E3
c = 2*mu*alpha
d = 2*mu*beta1

sol1 = odeint(wavefn, u0, r3, args=(a,b1,c,d))
sol2 = odeint(wavefn, u0, r3, args=(a,b2,c,d))
sol3 = odeint(wavefn, u0, r3, args=(a,b3,c,d))

u1 = sol1[:,0]
u11 = sol2[:,0]
u20 = sol3[:,0]

du1 = sol1[:,1]
du2 = sol2[:,1]
du3 = sol3[:,1]

pr1 = sqr(u1)
pr2 = sqr(u11)
pr3 = sqr(u20)

pr1, u1, du1, norm1 = simpson(pr1,u1,du1,n,l,r3)
pr2, u11, du2, norm2 = simpson(pr2,u11,du2,n,l,r3)
pr3, u20, du3, norm3 = simpson(pr3,u20,du3,n,l,r3)

tsn1, nod1 = statement(u1,n,l,r3,E1,1)
tsn2, nod2 = statement(u11,n,l,r3,E2,1)
tsn3, nod3 = statement(u20,n,l,r3,E3,1)

f1 = plt.figure(1,figsize=(8,6))
plt.plot(r3, u1, 'b',label="$E_{nl} =$ %.3f GeV" % E1)
plt.plot(r3, u11, 'g',label="$E_{nl} =$ %.3f GeV" % E2)
plt.plot(r3, u20, 'r',label="$E_{nl} =$ %.3f GeV" % E3)
plt.plot(tsn1[:,0], tsn1[:,1], 'ob')
plt.plot(tsn2[:,0], tsn2[:,1], 'og')
plt.plot(tsn3[:,0], tsn3[:,1], 'or')
plt.plot(nod1[:,0], nod1[:,1], 'xb')
plt.plot(nod2[:,0], nod2[:,1], 'xg')
plt.plot(nod3[:,0], nod3[:,1], 'xr')
### probs get gridspec in here for separate plots
plt.plot((r3[0],r3[-1]),(0,0),'grey')

plt.legend(loc=8, fontsize=25)
plt.title("Early Iteration of Bisection Method for Charmonium (2,0)",fontsize='xx-large')
plt.xlabel("Radial Distance, $GeV^{-1}$",fontsize='xx-large')
plt.ylabel("$u_{nl}(r)$",fontsize='xx-large')

#f2 = plt.figure(2,figsize=(8,6))
#plt.plot(r3, pr1, 'b',label="$E_{20} =$ %.3f GeV" % E1)
#plt.plot(r3, pr11, 'g',label="$E_{20} =$ %.3f GeV" % E11)
#plt.plot(r3, pr20, 'r',label="$E_{20} =$ %.3f GeV" % E20)
#### probs get gridspec in here for separate plots
#plt.plot((r3[0],r3[-1]),(0,0),'grey')
#
#plt.legend(loc=1, fontsize=25)
#plt.title("Probability Densities of Charmonium Wavefunction",fontsize='xx-large')
#plt.xlabel("Radial Distance, $GeV^{-1}$",fontsize='xx-large')
#plt.ylabel("$|u_{nl}(r)|^2$",fontsize='xx-large')

plt.show()
