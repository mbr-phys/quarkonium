from __future__ import division
import numpy as np
import math
from scipy.integrate import odeint, simps
import matplotlib.pyplot as plt
from qfns import S_hyperfine
from q2 import wavefn, psi, counter, itera, sqr, simpson, statement

u0 = [0, 1] #array of u and du/dr vals at 0

r1 = np.linspace(0.000001, 15, 1000)
r2 = np.linspace(0.0000001, 15, 1000)
r3 = np.linspace(0.000001, 15, 1000)
r4 = np.linspace(0.0001, 15, 1000)
r5 = np.linspace(0.00001, 15, 1000)

step1 = r1[1] - r1[0] #step size for normalisation
step2 = r2[1] - r2[0]
step3 = r3[1] - r3[0]
step4 = r4[1] - r4[0]
step5 = r5[1] - r5[0]

mc = 1.34
invmu = 1/mc + 1/mc
mu = 1/invmu
alpha_S = 0.4
alpha = 0.4*4/3

beta1 = 0.195
beta2 = 0.19512

E1 = 0.3
E2 = 0.5
E3 = 0.5

E11 = 0.8
E12 = 1.0
E13 = 1.2

E01 = 1.0
E02 = 1.5
E03 = 2.0

Ep1 = 1.3
Ep2 = 1.55
Ep3 = 1.8

E31 = 1.5
E32 = 1.4
E33 = 1.6

pr1, u1, du1, E1 = itera(1, 0, E1, E3, u0, alpha, beta1, mu, r1, step1)
pr11, u11, du11, E11 = itera(1, 1, E11, E13, u0, alpha, beta1, mu, r2, step2)
pr20, u20, du20, E20 = itera(2, 0, E01, E03, u0, alpha, beta1, mu, r3, step3)
pr21, u21, du21, E21 = itera(2, 1, Ep1, Ep3, u0, alpha, beta1, mu, r4, step4)
pr30, u30, du30, E30 = itera(3, 0, E31, E33, u0, alpha, beta1, mu, r5, step5)

norm = 1
statement(u1,du1,1,0,r1,E1,norm)
statement(u11,du11,1,1,r2,E11,norm)
statement(u20,du20,2,0,r3,E20,norm)
statement(u21,du21,2,1,r4,E21,norm)
statement(u30,du30,3,0,r3,E30,norm)

r0 = np.average(du1[0:6])
E1_3, E1_1 = S_hyperfine(alpha_S, mc, r0, E1)
M1_1 = E1_1 + 2*mc
M1_3 = E1_3 + 2*mc
print M1_1, M1_3

r20 = np.average(du20[0:6])
E2_3, E2_1 = S_hyperfine(alpha_S, mc, r20, E20)
M2_1 = E2_1 + 2*mc
M2_3 = E2_3 + 2*mc
print M2_1, M2_3

r30 = np.average(du30[0:6])
E3_3, E3_1 = S_hyperfine(alpha_S, mc, r30, E30)
M3_1 = E3_1 + 2*mc
M3_3 = E3_3 + 2*mc
print M3_1, M3_3

m10 = 2*mc + E1
m11 = 2*mc + E11
m20 = 2*mc + E20
m21 = 2*mc + E21
m30 = 2*mc + E30

#f1 = plt.figure(1,figsize=(8,6))
#ax = f1.add_subplot(111)
#ax.tick_params(axis='x', which='major', labelsize=15)
#ax.tick_params(axis='y', which='major', labelsize=15)
#plt.plot(r1, u1, 'b',label="(n,l) = (1,0), $M_{nl} =$ %.3f GeV/$c^2$" % m10)
#plt.plot(r2, u11, 'g',label="(n,l) = (1,1), $M_{nl} =$ %.3f GeV/$c^2$" % m11)
#plt.plot(r3, u20, 'r',label="(n,l) = (2,0), $M_{nl} =$ %.3f GeV/$c^2$" % m20)
#plt.plot(r4, u21, 'purple',label="(n,l) = (2,1), $M_{nl} =$ %.3f GeV/$c^2$" % m21)
#plt.plot(r5, u30, 'cyan',label="(n,l) = (3,0), $M_{nl} =$ %.3f GeV/$c^2$" % m30)
#### probs get gridspec in here for separate plots
#plt.plot((r5[0],r5[-1]),(0,0),'grey')
#
#plt.legend(loc=1, fontsize=15)
#plt.title("Solutions of the Charmonium Wavefunction",fontsize=25)
#plt.xlabel("Radial Distance, $GeV^{-1}$",fontsize=25)
#plt.ylabel("$u_{nl}(r)$",fontsize=25)

#f2 = plt.figure(2,figsize=(8,6))
#ax1 = f1.add_subplot(111)
#ax1.tick_params(axis='x', which='major', labelsize=15)
#ax1.tick_params(axis='y', which='major', labelsize=15)
#plt.plot(r1, pr1, 'b',label="(n,l) = (1,0), $M_{nl} =$ %.3f GeV/$c^2$" % m10)
#plt.plot(r2, pr11, 'g',label="(n,l) = (1,1), $M_{nl} =$ %.3f GeV/$c^2$" % m11)
#plt.plot(r3, pr20, 'r',label="(n,l) = (2,0), $M_{nl} =$ %.3f GeV/$c^2$" % m20)
#plt.plot(r4, pr21, 'purple',label="(n,l) = (2,1), $M_{nl} =$ %.3f GeV/$c^2$" % m21)
#plt.plot(r5, pr30, 'cyan',label="(n,l) = (3,0), $M_{nl} =$ %.3f GeV/$c^2$" % m30)
## probs get gridspec in here for separate plots
#plt.plot((r2[0],r2[-1]),(0,0),'grey')

#plt.legend(loc=1, fontsize=25)
#plt.title("Probability Densities of Charmonium Wavefunction",fontsize=25)
#plt.xlabel("Radial Distance, $GeV^{-1}$",fontsize=25)
#plt.ylabel("$|u_{nl}(r)|^2$",fontsize=25)

#f3 = plt.figure(3,figsize=(8,6))
plt.plot((0,1),(3.068,3.068),'blue')
plt.plot((0,1),(M1_1,M1_1),'blue')
plt.plot((0,1),(M1_3,M1_3),'blue')
plt.plot((1,2),(m20,m20),'purple')
plt.plot((1,2),(M2_3,M2_3),'purple')
plt.plot((1,2),(M2_1,M2_1),'purple')
plt.plot((2,3),(m30,m30),'cyan')
plt.plot((2,3),(M3_1,M3_1),'cyan')
plt.plot((2,3),(M3_3,M3_3),'cyan')

plt.plot((0,1),(3.096,3.096),'red')
plt.plot((0,1),(2.98,2.98),'red')
plt.plot((1,2),(3.638,3.638),'red')
plt.plot((1,2),(3.686,3.686),'red')
plt.plot((2,3),(4.039,4.039),'red')
plt.plot((2,3),(4.024,4.024),'red')

plt.show()

