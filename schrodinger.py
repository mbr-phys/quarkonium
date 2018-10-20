from __future__ import division
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def wavefn(u, r, a, b, c, d):
    du = (a/r**2)*u - b*u - (c/r)*u + d*r*u
    return du

u0 = 0

r = np.linspace(0.0000001, 15, 100)

m1 = 1.34
m2 = 1.34
l = 1
invmu = 1/m1 + 1/m2
mu = 1/invmu
alpha = 0.4
beta = 0
E = 3.0608 - 2*m1

a = l*(l+1)
b = 2*mu*E
c = (8/3)*mu*alpha
d = 2*mu*beta

sol = odeint(wavefn, u0, r, args=(a, b, c, d)) #unl

plt.plot(r, sol, 'b')
plt.show()
