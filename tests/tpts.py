def turningpoints(r,u,du,n,l):
    '''
        bla bla bla
    '''
    tpts = 0 
    tpt_loc = []
    for i in range(len(u)-1):
        if abs(du[i-1]) > abs(du[i]) and abs(du[i]) < abs(du[i+1]):
            tpts += 1
            tpti = i
            tpt_loc = np.append(tpt_loc,tpti)
        
    points = np.zeros((len(tpt_loc,2))
    for x in range(len(tpt_loc)):
        hunt = int(tpt_loc[x])
        points[x,0] = r[hunt]
        points[x,1] = u[hunt]

    print "    The number of turning points is $.f" % tpts

    return points, tpt_loc

def nodes(points, tpt_loc, u, r):
    '''
        I should probably start actually writing these properly
    '''
    nod = np.zeros(2)
    if len(tpt_loc) > 1:
        no_nodes = len(tpt_loc)-1
        print "    The number of nodes is %.f" % no_nodes
        nod = np.zeros((no_nodes,2))
        q = 0
        while nod[-1,0] == 0:
            start = tpt_loc[q]
            end = tpt_loc[q+1]
            rng = int(start-end)
            dr = np.zeros(rng)
            for p in range(rng):
                o = int(start + p)
                dr[p] = abs(u[o])
            nod[q,1] = min(dr)
            for pi in range(rng):
                if dr[pi] == nod[q,1]
                    oi = int(start + pi)
                    nod[q,0] = r[oi]
                    nod[q,1] = u[oi]
            q += 1
        return nod, no_nodes
    else:
        print "    There are no nodes for this function"
        no_nodes = 0
        return nod, no_nodes
        
