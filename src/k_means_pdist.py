#!/usr/bin/env python
# Based on the code at below and modified to use a pre-computed distance matrix
# http://stackoverflow.com/questions/5529625/is-it-possible-to-specify-your-own-distance-function-using-scikit-learn-k-means
# kmeanssample 2 pass, first sample sqrt(N)

import random
import numpy as np

def pdist(X, XA, XB):
    m = XA.shape[0]
    n = XB.shape[0]
    assert XA.shape[1] == m
    assert XA.shape[1] == m
    assert X.shape[1] == m
    d = np.ones((m,n))
    for i in xrange(m):
        for j in xrange(n):
            if XA[i,:].sum() > 0 and  XB[j,:].sum() > 0:
                d[i,j] = X[XA[i,:]==1,:][:, XB[j,:]==1].mean()
    return d

#...............................................................................
def kmeans(DIST, X, centres, delta=.001, maxiter=100, verbose=0):
    """ centres, Xtocentre, distances = kmeans( X, initial centres ... )
    in:
        X N x dim
        centres k x dim: initial centres, e.g. random.sample( X, k )
        delta: relative error, iterate until the average distance to centres
            is within delta of the previous average distance
        maxiter
        verbose: 0 silent, 2 prints running distances
    out:
        centres, k x dim
        Xtocentre: each X -> its nearest centre, ints N -> k
        distances, N
    see also: kmeanssample below, class Kmeans below.
    """
    N, dim = X.shape
    k, cdim = centres.shape
    if dim != cdim:
        raise ValueError( "kmeans: X %s and centres %s must have the same number of columns" % (
            X.shape, centres.shape ))
    if verbose:
        print "kmeans: X %s  centres %s  delta=%.2g  maxiter=%d " % (
            X.shape, centres.shape, delta, maxiter)
    allx = np.arange(N)
    prevdist = 0
    for jiter in range( 1, maxiter+1 ):
        D = pdist(DIST, X, centres)  # |X| x |centres|
        xtoc = D.argmin(axis=1)  # X -> nearest centre
        distances = D[allx,xtoc]
        avdist = distances.mean()  # median ?
        if verbose >= 2:
            print "kmeans: av |X - nearest centre| = %.4g" % avdist
        if (1 - delta) * prevdist <= avdist <= prevdist or jiter == maxiter:
            break
        prevdist = avdist
        for jc in range(k):  # (1 pass in C)
            c = np.where( xtoc == jc )[0]
            if len(c) > 0:
                centres[jc] = X[c].sum( axis=0 )
    if verbose:
        print "kmeans: %d iterations  cluster sizes:" % jiter, np.bincount(xtoc)
    if verbose >= 2:
        r50 = np.zeros(k)
        r90 = np.zeros(k)
        for j in range(k):
            dist = distances[ xtoc == j ]
            if len(dist) > 0:
                r50[j], r90[j] = np.percentile( dist, (50, 90) )
        print "kmeans: cluster 50 % radius", r50.astype(int)
        print "kmeans: cluster 90 % radius", r90.astype(int)
            # scale L1 / dim, L2 / sqrt(dim) ?
    return centres, xtoc, D

def kmeanssample(DIST, X, k, nsample=0, **kwargs ):
    """ 2-pass kmeans, fast for large N:
        1) kmeans a random sample of nsample ~ sqrt(N) from X
        2) full kmeans, starting from those centres
    """
    N, dim = X.shape
    if nsample == 0:
        nsample = max(2*np.sqrt(N), k)
    Xsample = randomsample( X, int(nsample) )
    pass1centres = randomsample( X, int(k) )
    sel = pass1centres.sum(0).reshape(-1)
    samplecentres = kmeans(DIST, Xsample, pass1centres[sel==1,:], **kwargs )[0]
    return kmeans(DIST, X, samplecentres, **kwargs )

def randomsample(X, n):
    sampleix = random.sample(xrange(X.shape[0]), int(n))
    SMP = np.zeros_like(X)
    for i in sampleix:
        SMP[i,:] = X[i,:]
    return SMP

def nearestcentres(DIST, X, centres):
    D = pdist(DIST, X, centres)  # |X| x |centres|
    return D.argmin(axis=1)

if __name__ == "__main__":
    import random
    import sys
    from time import time

    N = 1000
    dim = 3
    ncluster = 10
    kmdelta = .001
    kmiter = 100
    seed = 1

    exec( "\n".join( sys.argv[1:] ))  # run this.py N= ...
    np.set_printoptions(1, threshold=200, edgeitems=5, suppress=True)
    np.random.seed(seed)
    random.seed(seed)

    print "N %d  dim %d  ncluster %d" % ( N, dim, ncluster)
    X = np.random.exponential( size=(N,dim) )

    from scipy.spatial.distance import cdist
    DIST = cdist(X, X)
    print DIST.shape

    t0 = time()
    centres, xtoc, dist = kmeanssample(DIST, np.eye(X.shape[0]), ncluster, nsample=0, delta=kmdelta, maxiter=kmiter, verbose=2 )
    print "%.0f msec" % ((time() - t0) * 1000)
    print dist.shape

