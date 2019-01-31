'''Gonna fix qfns to work hopefully'''

from __future__ import division
import numpy as np
import math
from scipy.integrate import odeint, simps

#a = l*(l+1)
#b = 2*mu*E
#c = 2*mu*alpha
#d = 2*mu*beta

def wavefn(u0, r, a, b, c, d):
    '''
        Stuff
    '''
    u, v = u0
    dv = [v, (a/r**2)*u - b*u - (c/r)*u + d*r*u]
    return dv

def wavey(u0, r, a, b, c, d):
    '''
        Stuff
    '''
    sol = odeint(wavefn, u0, r, args=(a,b,c,d))
    ur = sol[:,0]
    dudr = sol[:,1]
    return ur, dudr

#def dun(ur, dudr, r, l, a, b, c, d):
#    '''
#        Stuff
#    '''
#    dvdr = (a/r**2)*ur - b*ur - (c/r)*ur + d*r*ur
#    return dudr, dvdr

def counter(r, u, dudr):
    '''
        Counts nodes and turning points of wvfn
    '''
    nodes = 0
    turns = 0
    for i in range(len(r)-1):
        signu = u[i]/u[i+1]
        signv = dudr[i]/dudr[i+1]
        if signu < 0:
            nodes += 1
        if signv < 0:
            turns += 1
    return nodes, turns

def energies(E1, E2, E3, nodes1, nodes2, nodes3, turns1, turns2, turns3):
    '''
        Picks new energy
    '''
    if turns1 != turns2 or nodes1 != nodes2:
        Emin = E1
        Emax = E2
        test = 1
    elif turns2 != turns3 or nodes2 != nodes3:
        Emin = E2
        Emax = E3
        test = 1
    else: 
        Emin = E1
        Emax = E3
        test = 0
    return Emin, Emax, test

def itera(n, l, E1, E3, u0, alpha, beta, mu, r):
    '''
        Iterate for solution
    '''
    a = l*(l+1)
    c = 2*mu*alpha
    d = 2*mu*beta

    i = 0
    q = 0
    while i == 0:
        if q > 1000:
            print "Iterative solution took too long to find"
            break

        E2 = (E1+E3)/2

        if E1 == E2: 
            E3 = E2
        elif E3 == E2:
            E1 = E2

        b1 = 2*mu*E1
        b2 = 2*mu*E2
        b3 = 2*mu*E3

        u1, du1 = wavey(u0, r, a, b1, c, d)
        u2, du2 = wavey(u0, r, a, b2, c, d)
        u3, du3 = wavey(u0, r, a, b3, c, d)

        nds1, tpts1 = counter(r, u1, du1)
        nds2, tpts2 = counter(r, u2, du2)
        nds3, tpts3 = counter(r, u3, du3)

        E1, E3, test = energies(E1, E2, E3, nds1, nds2, nds3, tpts1, tpts2, tpts3)
        if test == 1:
            q += 1
        else:
            if nds1 == n-1 and tpts1 == n:
                print "Found iterative solution"
                break
            elif nds2 == n-1 and tpts2 == n:
                print "Found iterative solution"
                break
            elif nds3 == n-1 and tpts2 == n:
                print "Found interactive solution"
                break
            else: 
                print "No solution found with current energies"
                delta = E3 - E1
                if E1 != 0:
                    odmag = int(math.floor(math.log10(abs(E1))))
                elif E2 != 0:
                    odmag = int(math.floor(math.log10(abs(E2))))
                elif E3 != 0:
                    odmag = int(math.floor(math.log10(abs(E3))))
                else:
                    odmag = int(input("Energies converged to 0, give new order of magnitude: "))
                    E1 = 5*10**odmag
                cond = 1e-18
                if delta < cond:
                    print "Energies converged"
                    E3 = 0.5
                    boom = E3/3
                    E1 = E3 - 2*boom
                else:
                    boom = 0.25*10**odmag
                    E3 = E1 
                    E1 = E3 - 2*boom
                q += 1

    return u2, du2, E2

def simpson(u,du,r):
    '''
        Stuff
    '''
    norm = simps(u,r,even="first")
    u = u/norm
    du = du/norm

    return u, du

def sqr(u):
    '''
        Returns square of wavefunction
    '''
    pr = np.zeros(u.shape)
    for ni in range(len(u)):
        pr[ni] = abs(u[ni])**2

    return pr

