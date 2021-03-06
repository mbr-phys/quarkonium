'''Module containing functions made for solving schrodinger equation and
    other things for quarkonium project'''

from __future__ import division
import numpy as np
import math
from scipy.integrate import odeint,simps

def wavefn(u0, r, a, b, c, d):
    '''
        Returns array of the two odes required to solve quarkonium schrodinger eqn
    '''
    u, v = u0
    dv = [v, (a/r**2)*u - b*u - (c/r)*u + d*r*u]
    return dv

def turningpoints(r, u, n, l):
    '''
        Evaluates any maxima and minima a function may have, printing how many of each and returning coords of them
    '''
    mins = 0
    min_loc = []
    maxs = 0
    max_loc = []
    if n == 1:
        tolerance = 1e-19
    elif n == 2 and l == 0:
        tolerance = 1e-12
    else:
        tolerance = 1e-18
    for i in range(len(u)-2): #maybe look at specifying this more?
        if u[i] - u[i+1] < -1*tolerance and u[i+2] - u[i+1] < -1*tolerance:
            maxs += 1
            maxi = i+1
            max_loc = np.append(max_loc,maxi)
        elif u[i] - u[i+1] > tolerance and u[i+2] - u[i+1] > tolerance:
            mins += 1
            mini = i+1
            min_loc = np.append(min_loc,mini)

    maxes = np.zeros((len(max_loc),2))
    for x in range(len(max_loc)):
        hunt = int(max_loc[x])
        maxes[x,0] = r[hunt]
        maxes[x,1] = u[hunt]

    minies = np.zeros((len(min_loc),2))
    for y in range(len(min_loc)):
        hunt = int(min_loc[y])
        minies[y,0] = r[hunt]
        minies[y,1] = u[hunt]

    print "    The number of maxima is %.f" % maxs
    print "    The number of minima is %.f" % mins

    return maxes, minies, max_loc, min_loc

def counter(u,r,mx,mn):
    '''
    Alternate notes counter
    '''
    if mn.shape == (0,2) and mx.shape == (0,2):
        tr = r[-1]
    elif mn.shape == (0,2) and mx.shape != (0,2):
        tr = mx[-1,0]
    elif mx.shape == (0,2) and mn.shape != (0,2):
        tr = mn[-1,0]
    elif mx[-1,0] > mn[-1,0]:
        tr = mx[-1,0]
    else:
        tr = mn[-1,0]

    nodes = 0 
    nr = []
    nu = [] 
    for i in range(len(r)-6):
        if r[i] > tr:
            break
        if u[i+5]/u[i+6] < 0:
            nodes += 1
            rn = (r[i+5]+r[i+6])/2
            nr = np.append(nr,rn)
            un = (u[i+5]+u[i+6])/2
            nu = np.append(nu,un)

    print "  The number of nodes is %.f" % nodes

    return nodes, nr, nu

def nodes(maxima, minima, maxs, mins, u, r):
    '''
        Function to calculate the number of nodes of a function,
        and find their coordinates
    '''
    nod = np.zeros(2)
    if len(maxs) > 0 and len(mins) > 0:
        no_nodes = len(maxs)+len(mins)-1
        print "    The number of nodes is %.f" % no_nodes
        nodes_counter = 0
        nod = np.zeros((no_nodes,2))
        if maxima[0,0] < minima[0,0]:
            q = 0
            while nod[-1,0] == 0:
                if nodes_counter%2 == 0:
                    qis = int(q - nodes_counter/2)
                    top = mins[qis]
                    btm = maxs[qis]
                    rng = int(top-btm)
                    dr = np.zeros(rng)
                    for p in range(rng):
                        o = int(btm + p)
                        dr[p] = abs(u[o])
                    nod[q,1] = min(dr)
                    for pi in range(rng):
                        if dr[pi] == nod[q,1]:
                            oi = int(btm + pi)
                            nod[q,0] = r[oi]
                            nod[q,1] = u[oi]
                    nodes_counter += 1
                    q+= 1
                else:
                    qis = int(q - (nodes_counter-1)/2)
                    top = maxs[qis]
                    btm = mins[qis-1]
                    rng = int(top-btm)
                    dr = np.zeros(rng)
                    for p in range(rng):
                        o = int(btm + p)
                        dr[p] = abs(u[o])
                    nod[q,1] = min(dr)
                    for pi in range(rng):
                        if dr[pi] == nod[q,1]:
                            oi = int(btm + pi)
                            nod[q,0] = r[oi]
                            nod[q,1] = u[oi]
                    nodes_counter += 1
                    q += 1
            return nod, no_nodes
        else:
            q = 0
            while nod[-1,0] == 0:
                if nodes_counter%2 == 0:
                    top = maxs[q]
                    btm = mins[q]
                    rng = int(top-btm)
                    dr = np.zeros(rng)
                    for p in range(rng):
                        o = int(btm + p)
                        dr[p] = abs(u[o])
                    nod[q,1] = min(dr)
                    for pi in range(rng):
                        if dr[pi] == nod[q,1]:
                            oi = int(btm + pi)
                            nod[q,0] = r[oi]
                            nod[q,1] = u[oi]
                    nodes_counter += 1
                    q += 1
                else:
                    top = mins[q]
                    btm = maxs[q]
                    rng = int(top-btm)
                    dr = np.zeros(rng)
                    for p in range(rng):
                        o = int(btm + p)
                        dr[p] = abs(u[o])
                    nod[q,1] = min(dr)
                    for pi in range(rng):
                        if dr[pi] == nod[q,1]:
                            oi = int(btm + pi)
                            nod[q,0] = r[oi]
                            nod[q,1] = u[oi]
                    nodes_counter += 1
                return nod, no_nodes
    else:
        print "    There are no nodes for this function."
        no_nodes = 0
        return nod, no_nodes

def itera(n, l, E1, E2, E3, u0, alpha, beta, mu, r, step):
    '''
        iterate wavefunction solution over energies to find correct energy for expected nodes and turning points
    '''
    a = l*(l+1)
    c = 2*mu*alpha
    d = 2*mu*beta

    q = 0
    i = 0
    while i==0:
        if q > 2000:
            print "Iterative solution took too long to find"
            break
        if E1 == E2:
            E3 = E2
        elif E3 == E2:
            E1 = E2
        b1 = -2*mu*E1
        b2 = -2*mu*E2
        b3 = -2*mu*E3

        sol1 = odeint(wavefn, u0, r, args=(a,b1,c,d))
        sol2 = odeint(wavefn, u0, r, args=(a,b2,c,d))
        sol3 = odeint(wavefn, u0, r, args=(a,b3,c,d))

        soli1 = sol1[:,0]
        soli2 = sol2[:,0]
        soli3 = sol3[:,0]

        du1 = sol1[:,1]
        du2 = sol2[:,1]
        du3 = sol3[:,1]

        maxi1, mini1, maxs1, mins1 = turningpoints(r, soli1,n,l)
        maxi2, mini2, maxs2, mins2 = turningpoints(r, soli2,n,l)
        maxi3, mini3, maxs3, mins3 = turningpoints(r, soli3,n,l)

        tps1 = len(maxi1) + len(mini1)
        tps2 = len(maxi2) + len(mini2)
        tps3 = len(maxi3) + len(mini3)

        n1, nr1, nu1 = counter(soli1, r, maxi1, mini1)
        n2, nr2, nu2 = counter(soli2, r, maxi2, mini2)
        n3, nr3, nu3 = counter(soli3, r, maxi3, mini3)

        if n1 != n2 or tps1 != tps2:
            E3 = E2
            E2 = (E1 + E3)/2
            q += 1
        elif n2 != n3 or tps2 != tps3:
            E1 = E2
            E2 = (E1 + E3)/2
            q += 1
        else:
            if n1 == (n-l-1) and tps1 == (n-l):
                print "Found iterative solution"
                break
            elif n2 == (n-l-1) and tps2 == (n-l):
                print "Found iterative solution"
                break
            elif n3 == (n-l-1) and tps3 == (n-l):
                print "Found iterative solution"
                break
            else:
                print "No solution found with E1,2,3 = %.3e, %.3e, %.3e" % (E1, E2, E3)
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
                if delta < 1e-18:
                    print "Energies converged"
                    if n == 2:
                        E1 = 5e-6
                        boom = E1/3
                    elif n == 3:
                        E1 = 3e-6
                        boom = E1/3
                    E3 = E1
                    E2 = E3 - boom
                    E1 = E2 - boom
                else:
                    boom = 0.25*10**odmag
                    E3 = E1
                    E2 = E3 - boom
                    E1 = E2 - boom
                if E1 < 0:
                    E1 = -1*E1
                if E2 < 0:
                    E2 = -1*E2
                if E3 < 0:
                    E3 = -1*E3
                q += 1

    return soli2, du2, E2

def simpson(pr,wvfn,dwv,n,l,r):
    '''
        Stuff
    '''
    norm = simps(pr,r,even='first')
    pr = pr/norm
    wvfn = wvfn/np.sqrt(norm)
    dwv = dwv/np.sqrt(norm)

    return pr, wvfn, dwv, norm

def statement(wvfn,n,l,r,E,norm):
    '''
        Statement of maxima, minima, nodes, and energy of functions
    '''
    print "For (n,l) = (%.f,%.f):" % (n,l)
    maxima, minima, maxs, mins = turningpoints(r,wvfn,n,l)
    maxima[:,1] = maxima[:,1]/norm
    minima[:,1] = minima[:,1]/norm
    for z in range(len(maxima[:,0])):
        print "    Maxima No. %.f is found at r = %.2f, u = %.5e" % (z+1,maxima[z,0],maxima[z,1])

    for t in range(len(minima[:,0])):
        print "    Minima No. %.f is found at r = %.2f, u = %.5e" % (t+1,minima[t,0],minima[t,1])

    n, nr, nu = counter(wvfn, r, maxima, minima)
    if n != 0:
        for ns in range(n):
            print "    Node No. %.f is found at r = %.2f, u = %.5e" % (ns+1,nr[ns], nu[ns])

    print "    The energy of this state is %.4e" % E

    return None

def sqr(wvfn):
    '''
        Returns square of wavefunction
    '''

    prob_dens = np.zeros(wvfn.shape)
    for ni in range(len(wvfn)):
        prob_dens[ni] = abs(wvfn[ni])**2

    return prob_dens
