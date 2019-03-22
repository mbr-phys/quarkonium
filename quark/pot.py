from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

r = np.linspace(0.0001, 15, 1000)

v1 = -1.6/(3*r) + 0.195*r
v2 = 7.641*r**0.1 - 8.064

f1 = plt.figure(1,figsize=(8,6))
plt.plot(r,v1,'r')
plt.plot(r,v2,'b')

plt.show()

