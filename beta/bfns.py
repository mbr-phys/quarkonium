'''Module containing functions made for solving schrodinger equatino and 
    other things for quarkonium project'''

from __future__ import division
import numpy as np
import math
from scipy.integrate import odeint, simps

def wavefn(u0,r,a,b,c,d):
    '''
        Return array of the two odes required to solve quarkonium schrodinger eqn
    '''
    u, v = u0
    dv = [v, (a/r**2)*u - b*u - (c/r)*u + d*r*u]
    return dv

def sqr(psi):
    '''
        Squares wavefn
    '''
    pdf = np.zeros(psi.shape)
    for ni in range(len(psi)):
        pdf[ni] = abs(psi[ni])**2

    return pdf

def simpson(pdf,psi,dur,r):
    '''
        Normalises wavefn
    '''
    norm = simps(pdf,r,even='first')
    pdf = pdf/norm
    psi = psi/np.sqrt(norm)
    dur = dur/np.sqrt(norm)

    return pdf, psi, dur

def psi(n,l,E,u0,alpha,beta,mu,r):
    '''
        Solves wavefn
    '''
    a = l*(l+1)
    b = 2*mu*E
    c = 2*mu*alpha
    d = 2*mu*beta

    sol = odeint(wavefn,u0,r,args=(a,b,c,d))
    ur = sol[:,0]
    du = sol[:,1]

    pr = sqr(ur)

    pr, ur, du = simpson(pr,ur,du,r)

    return pr, ur, du

def counter(u,du,r):
    '''
        Counts nodes and turning points
    '''
    nodes = 0
    nod_loc_r = []
    nod_loc_u = []
    turns = 0
    t_loc_r = []
    t_loc_u = []
    for i in range(len(r)-30):
        if du[i+5]/du[i+6] < 0:
            turns += 1
            rn = (r[i+5]+r[i+6])/2
            t_loc_r = np.append(t_loc_r,rn)
            un = (u[i+5]+u[i+6])/2
            t_loc_u = np.append(t_loc_u,un)

    for ns in range(len(r)-30):
#        if len(t_loc_r) > 0:
#            if r[ns] > t_loc_r[-1]:
#                break
        if u[ns+5]/u[ns+6] < 0:
            nodes += 1
            rn = (r[ns+5]+r[ns+6])/2
            nod_loc_r = np.append(nod_loc_r,rn)
            un = (u[ns+5]+u[ns+6])/2
            nod_loc_u = np.append(nod_loc_u,un)

    print "  The number of nodes is %.f" % nodes
    print "  The number of turning points is %.f" % turns

    return nodes, nod_loc_r, nod_loc_u, turns, t_loc_r, t_loc_u

def itera(n,l,e,u0,alpha,be1,be3,mu,r,step):
    '''
        Iterates to find ur
    '''
    q = 0 
    i = 0
    while i == 0:
        if q > 1000:
            print "Iterative solution took too long to find"
            break
        be2 = (be1+be3)/2

        pr1, u1, v1 = psi(n,l,e,u0,alpha,be1,mu,r)
        pr2, u2, v2 = psi(n,l,e,u0,alpha,be2,mu,r)
        pr3, u3, v3 = psi(n,l,e,u0,alpha,be3,mu,r)

        n1, nr1, nu1, t1, tr1, tu1 = counter(u1,v1,r)
        n2, nr2, nu2, t2, tr2, tu2 = counter(u2,v2,r)
        n3, nr3, nu3, t3, tr3, tu3 = counter(u3,v3,r)

#        if t1 != t2 or n1 != n2:
#            e3 = e2
#            q += 1
#        elif t2 != t3 or n2 != n3:
#            e1 = e2
#            q += 1
#        else: 
        if n1 == (n-1) and t1 == (n):
            print "Found iterative solution in %.f iterations" % q
            break 
        elif n2 == (n-1) and t2 == (n):
            print "Found iterative solution in %.f iterations" % q
            break 
        elif n3 == (n-1) and t3 == (n):
            print "Found iterative solution in %.f iterations" % q
            break 
        if t1 != t2 or n1 != n2:
            be3 = be2
            q += 1
        elif t2 != t3 or n2 != n3:
            be1 = be2
            q += 1
        else: 
            delta = be3 - be1
            if be1 != 0:
                odmag = int(math.floor(math.log10(abs(be1))))
            elif be2 != 0:
                odmag = int(math.floor(math.log10(abs(be2))))
            elif be3 != 0:
                odmag = int(math.floor(math.log10(abs(be3))))
            else:
                odmag = int(input("Energies converged to 0, give new order of magnitude: "))
                be1 = 5*10**odmag
            if delta < 1e-18:
                print "Energies converged"
                be1 = 0.1
                be3 = 0.2
            else:
                boom = 0.5*10**odmag 
                be3 = be1 
                be1 = be3 - boom

    return pr2, u2, v2, be2


def statement(wvfn,du,n,l,r,be,norm):
    '''
        Statement of maxima, minima, nodes, and energy of functions
    '''
    print "For (n,l) = (%.f,%.f):" % (n,l)
    n, nr, nu, t, tr, tu = counter(wvfn,du,r)
    for z in range(t):
        print "    Turning Point No. %.f is found at r = %.2f, u = %.5e" % (z+1,tr[z],tu[z])

    if n != 0:
        for ns in range(n):
            print "    Node No. %.f is found at r = %.2f, u = %.5e" % (ns+1,nr[ns], nu[ns])

    print "    Beta of this state is %.4e" % be

    return None

