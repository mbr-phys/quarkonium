'''Module containing functions made for solving schrodinger equation and
    other things for quarkonium project'''

from __future__ import division
import numpy as np

def wavefn(u0, r, a, b, c, d):
    '''
        Returns array of the two odes required to solve quarkonium schrodinger eqn
    '''
    u, v = u0
    dv = [v, (a/r**2)*u - b*u - (c/r)*u + d*r*u]
    return dv

def turningpoints(r, u):
    '''
        Evaluates any maxima and minima a function may have, printing how many of each and returning coords of them
    '''
    mins = 0
    min_loc = []
    maxs = 0
    max_loc = []
    for i in range(len(u)-2): #maybe look at specifying this more?
        if u[i] < u[i+1] and u[i+2] < u[i+1]:
            maxs += 1
            maxi = i+1
            max_loc = np.append(max_loc,maxi)
        elif u[i] > u[i+1] and u[i+2] > u[i+1]:
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

    print "The number of maxima is %.f" % maxs
    print "The number of minima is %.f" % mins

    return maxes, minies, max_loc, min_loc

def nodes(maxima, minima, maxs, mins, u, r):
    if len(maxs) > 0 and len(mins) > 0:
        no_nodes = len(maxs)+len(mins)-1
        print "The number of nodes is %.f" % no_nodes
        nodes_counter = 0
        nod = np.zeros((no_nodes,2))
        if maxima[0,0] < minima[0,0]:
            q = 0
            while nod[-1,0] == 0: 
                if nodes_counter%2 == 0:
                    print q
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
            return nod
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
                    print nodes_counter
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
                    print nodes_counter
                return nod
    else:
        print "There are no nodes for this function."
 
'''
    to do:
        make whole thing run iterate to find beta/E based on expected turning points and nodes for known values of n and l
        switch to while loops?
'''
