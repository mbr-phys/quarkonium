# Project Plan

## 1. Program Design

### 1. Solve SE for specific quark mass and energy values

```
from __future__ import division
import numpy as np
from scipy.integrate import odeint
import matplotlib as plt

def wavefn(u, r, a, b, c, d):
    du = [u, (a/r**2)*u - b*u - (c/r)*u + d*r*u]
    return u

u0 = 0
a = ? #l(l+1)
b = ? #2 mu E
c = ? #8 mu alpha/3
d = ? 2 mu beta

r = np.linspace(0, 15?, 1000?)

sol = odeint(wavefn, u, r, args(a, b)) #unl
soln = sol/r #Rnl
```

What is beta??

### 2. Calculate number of nodes and turning points

### 3. Guess new energy or beta and iterate for solution

### 4. Find energy or beta, whichever not found before

### 5. Plot resulting wavefunctions

## 2. Milestone Program

### Calculate energy levels and wavefns of hydrogen

#### Check for (n,l) = (1,0), (2,0), (2,1) states

#### Compare to analytic solns

### Calculate energy of (2,1) state and plot |u(nl)|^2

## 3. Program Extension

### Apply program to Charmonium

#### Calculate beta to 3dp, using spin average mass

#### Predict energy for (1,1) and (2,0)

#### Produce plots of normalised wavefns for (1,0), (1,1), (2,0) states

## 4. Further Research


