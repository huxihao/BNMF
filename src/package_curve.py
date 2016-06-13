## Simulation the Hilbert Curve
## Author: Xihao Hu <huxihao@gmail.com>

from tools import *

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def hilbert_curve(itr=2, st=(0,0), od=1):
    ''' Generate the map of Hilbert Curve (HC) with given size
        The input number is iteration number of creating HC
    '''
    x, y = st
    if itr == 1:
        ## up: 0     down: 1  __   left: 2 __    right: 3  __
        ##      |  |         |  |            |            |
        ##      |__|         |  |          __|            |__
        ##
        if od == 0: return [(x,y), (x,y-1), (x+1,y-1), (x+1,y)] ## up
        elif od==1: return [(x,y), (x,y+1), (x+1,y+1), (x+1,y)] ## down
        elif od==2: return [(x,y), (x+1,y), (x+1,y-1), (x,y-1)] ## left
        elif od==3: return [(x,y), (x-1,y), (x-1,y-1), (x,y-1)] ## right
        else: assert False ## invalid value
    else:
        ## recursively define the four shapes
        i = 2**(itr-1) - 1 ## size of one dimention
        if od == 0:
            return hilbert_curve(itr-1, (x,y), 2) + \
                   hilbert_curve(itr-1, (x,y-i-1), 0) + \
                   hilbert_curve(itr-1, (x+i+1,y-i-1), 0) + \
                   hilbert_curve(itr-1, (x+i+i+1,y), 3)[::-1]
        elif od == 1:
            return hilbert_curve(itr-1, (x,y+i), 2)[::-1] + \
                   hilbert_curve(itr-1, (x,y+i+1), 1) + \
                   hilbert_curve(itr-1, (x+i+1,y+i+1), 1) + \
                   hilbert_curve(itr-1, (x+i+i+1,y+i), 3)
        elif od == 2:
            return hilbert_curve(itr-1, (x,y), 0) + \
                   hilbert_curve(itr-1, (x+i+1,y), 2) + \
                   hilbert_curve(itr-1, (x+i+1,y-i-1), 2) + \
                   hilbert_curve(itr-1, (x,y-i-i-1), 1)[::-1]
        elif od == 3:
            return hilbert_curve(itr-1, (x-i,y), 0)[::-1] + \
                   hilbert_curve(itr-1, (x-i-1,y), 3) + \
                   hilbert_curve(itr-1, (x-i-1,y-i-1), 3) + \
                   hilbert_curve(itr-1, (x-i,y-i-i-1), 1)
        else: assert False ## invalid value

def hilbert_show(pdf, it):
    fig = plt.figure()
    plt.suptitle("2D Hilbert Curve (n=%d)"%it)
    for i in xrange(4):
        ax = fig.add_subplot(221+i)
        verts = hilbert_curve(it, od=i)
        x, y = zip(*verts)
        plt.plot(x, y, 'bo-')
        ax.set_xlim(min(x)-1, max(x)+1)
        ax.set_ylim(min(y)-1, max(y)+1)
    pdf.savefig(); plt.clf()

def peano_curve(itr=2, st=(0,0), od=0):
    ''' Generate Peano Curve (PC) with given size '''
    x, y = st
    if itr == 1:
        ## left: 0         right: 1
        ##        ------          ------
        ##       |                      |
        ##        ------          ------
        ##              |        |
        ##        ------          ------
        if od == 0:
            return [(x,y), (x+1,y), (x+2,y), (x+2,y+1), (x+1,y+1),\
                    (x,y+1), (x,y+2), (x+1,y+2), (x+2,y+2)]
        elif od == 1:
            return [(x,y), (x+1,y), (x+2,y), (x+2,y-1), (x+1,y-1),\
                    (x,y-1), (x,y-2), (x+1,y-2), (x+2,y-2)]
        else: assert False ## invalid value
    else:
        i = 3**(itr-1) - 1
        if od == 0:
            return peano_curve(itr-1, (x,y), 0) + \
                   peano_curve(itr-1, (x+i+1,y+i), 1) + \
                   peano_curve(itr-1, (x+i+1+i+1,y), 0) + \
                   peano_curve(itr-1, (x+i+1+i+1,y+i+1+i), 1)[::-1] + \
                   peano_curve(itr-1, (x+i+1,y+i+1), 0)[::-1] + \
                   peano_curve(itr-1, (x,y+i+1+i), 1)[::-1] + \
                   peano_curve(itr-1, (x,y+i+1+i+1), 0) + \
                   peano_curve(itr-1, (x+i+1,y+i+1+i+1+i), 1) + \
                   peano_curve(itr-1, (x+i+1+i+1,y+i+1+i+1), 0)
        if od == 1:
            return peano_curve(itr-1, (x,y), 1) + \
                   peano_curve(itr-1, (x+i+1,y-i), 0) + \
                   peano_curve(itr-1, (x+i+1+i+1,y), 1) + \
                   peano_curve(itr-1, (x+i+1+i+1,y-i-1-i), 0)[::-1] + \
                   peano_curve(itr-1, (x+i+1,y-i-1), 1)[::-1] + \
                   peano_curve(itr-1, (x,y-i-1-i), 0)[::-1] + \
                   peano_curve(itr-1, (x,y-i-1-i-1), 1) + \
                   peano_curve(itr-1, (x+i+1,y-i-1-i-1-i), 0) + \
                   peano_curve(itr-1, (x+i+1+i+1,y-i-1-i-1), 1)
        else: assert False ## invalid value

def peano_show(pdf, it):
    fig = plt.figure()
    plt.suptitle("Peano Curve, %d iteration"%it)
    if True:
        ax = fig.add_subplot(111)
        verts = peano_curve(it, od=0)
        x, y = zip(*verts)
        plt.plot(x, y, 'bo-')
        ax.set_xlim(min(x)-1, max(x)+1)
        ax.set_ylim(min(y)-1, max(y)+1)
    pdf.savefig(); plt.clf()

def hilbert3_curve(itr=2, st=(0,0,0), od=10):
    ''' 3 dimentional Hilbert curve '''
    x, y, z = st
    d0 = [(0,0,0), (0,1,0), (0,1,1), (0,0,1), \
          (1,0,1), (1,1,1), (1,1,0), (1,0,0)]
    ## od = 10, x expand from (0,0,0) to (1,0,0)
    ## od = 20, y expand from (0,0,0) to (0,1,0)
    ## od = 30, z expand from (0,0,0) to (0,0,1)
    ## od = 11, clock wise rotate along x-xls
    if itr > 1:
        i = 2**(itr-1) - 1
        d0 = hilbert3_curve(itr-1, (x,y,z), 21) + \
             hilbert3_curve(itr-1, (x,y+i+1,z), 31) + \
             hilbert3_curve(itr-1, (x,y+i+1,z+i+1), 31) + \
             hilbert3_curve(itr-1, (x,y+i,z+i+1+i), 11) + \
             hilbert3_curve(itr-1, (x+i+1,y+i,z+i+1+i), 11) + \
             hilbert3_curve(itr-1, (x+i+1+i,y+i+1,z+i+1), 30)[::-1] + \
             hilbert3_curve(itr-1, (x+i+1+i,y+i+1,z), 30)[::-1] + \
             hilbert3_curve(itr-1, (x+i+1+i,y,z), 20)[::-1]
    x0, y0, z0 = d0[0]
    ## rotate
    if od == 10:
        return [(x+x1-x0,y+y1-y0,z+z1-z0) for x1,y1,z1 in d0]
    elif od == 20:
        return [(x-y1+y0,y+x1-x0,z+z1-z0) for x1,y1,z1 in d0]
    elif od == 30:
        return [(x-z1+z0,y+y1-y0,z+x1-x0) for x1,y1,z1 in d0]
    elif od == 11:
        return [(x+x1-x0,y-y1+y0,z-z1+z0) for x1,y1,z1 in d0]
    elif od == 21:
        return [(x+y1-y0,y+x1-x0,z+z1-z0) for x1,y1,z1 in d0]
    elif od == 31:
        return [(x+z1-z0,y+y1-y0,z+x1-x0) for x1,y1,z1 in d0]
    else: assert False ## invalid value

def hilbert3_show(pdf, it):
    from mpl_toolkits.mplot3d import axes3d, Axes3D
    fig = plt.figure()
    plt.suptitle("3D Hilbert Curve, %d iteration"%it)
    if True:
        ax = Axes3D(fig)
        verts = hilbert3_curve(it, od=10)
        x, y, z = zip(*verts)
        ax.plot(x,y,z,'bo-')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.set_xlim3d(min(x)-1, max(x)+1)
        ax.set_ylim3d(min(y)-1, max(y)+1)
        ax.set_zlim3d(min(z)-1, max(z)+1)
    pdf.savefig(); plt.clf()

def random_curve(n, m=1000, seed=2013):
    u = [(i-n/2, +2) for i in xrange(n)]
    v = [(i-n/2, -2) for i in xrange(n)]
    import random
    random.seed(seed)
    for iter in xrange(m):
        for i in random.choice([range(n), range(n-1,-1,-1)]):
            t = random.choice([u, v])
            if i == 0:
                x1, y1 = t[i+1]
                xm, ym = random.choice([(1,0), (-1,0), (0,1), (0,-1)])
                new_x, new_y = x1+xm, y1+ym
            elif i == n-1:
                x1, y1 = t[i-1]
                xm, ym = random.choice([(1,0), (-1,0), (0,1), (0,-1)])
                new_x, new_y = x1+xm, y1+ym
            else: ## middle ones
                x1, y1 = t[i-1]
                x2, y2 = t[i+1]
                if (x1-x2)**2 + (y1-y2)**2 != 2:
                    continue
                new_x, new_y = random.choice([(x1, y2), (x2, y1)])
            ## accept or reject
            if (new_x, new_y) in u or (new_x, new_y) in v:
                continue
            old_x, old_y = t[i]
            if (old_x**2+old_y**2) < (new_x**2+new_y**2):
                if random.random() < 0.8:
                    continue ## reject
            if random.random() < 0.8: ## accept
                t[i] = (new_x, new_y)
    return u + v

def curve_show(verts, fig=None, choose=None, color='ko-'):
    from mpl_toolkits.mplot3d import Axes3D
    v = np.array(verts) 
    if choose == None:
        choose = [i for i in xrange(v.shape[0])]
    if v.shape[1] == 2:
        if fig == None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        else:
            ax = fig
        plt.plot(v[choose,0], v[choose,1], color)
        #ax.grid()
        #ax.set_xlabel('x')
        #ax.set_ylabel('y')
        ax.set_xlim(min(v[:,0])-1, max(v[:,0])+1)
        ax.set_ylim(min(v[:,1])-1, max(v[:,1])+1)
    elif v.shape[1] == 3:
        #ax = fig.add_subplot(111, projection='3d')
        if fig == None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax = Axes3D(fig)
        else:
            ax = fig
        plt.plot(v[:,0], v[:,1], v[:,2], 'k.')
        plt.plot(v[choose,0], v[choose,1], v[choose,2], color)
        #ax.grid()
        #ax.set_xlabel('x')
        #ax.set_ylabel('y')
        #ax.set_zlabel('z')
        ax.set_xlim3d(min(v[:,0])-1, max(v[:,0])+1)
        ax.set_ylim3d(min(v[:,1])-1, max(v[:,1])+1)
        ax.set_zlim3d(min(v[:,2])-1, max(v[:,2])+1)
    return fig

def sep_map_show(pdf, verts, Q):
    v = np.array(verts)
    color = ''.join('bgrcmy'*Q.shape[1])
    fig = curve_show(v)
    ax = fig.add_subplot(111)
    for i in xrange(Q.shape[1]):
        if v.shape[1] == 3:
            pdf.savefig(); plt.clf();
        c = [Q[j,i] == np.max(Q[j,:]) and bool(Q[j,i]) for j in xrange(Q.shape[0])]
        curve_show(v, ax, [j for j in xrange(len(c)) if c[j]], color[i]+'o--')
    pdf.savefig(); plt.clf();
    fig = plt.figure()
    for i in xrange(Q.shape[1]):
        ax = fig.add_subplot(221+i)
        curve_show(v, ax, color='k-')
        c = [Q[j,i] > np.mean(Q[:,i]) and bool(Q[j,i]) for j in xrange(Q.shape[0])]
        curve_show(v, ax, [j for j in xrange(len(c)) if c[j]], color[i]+'o')
        plt.title('Cluster %s'%(i+1))
    pdf.savefig(); plt.clf();

def distance_map(A,t):
    vfunc = np.vectorize(lambda a:max(0,t-a))
    return vfunc(A)

def trace_sum(A):
    n = A.shape[0]
    s = np.zeros(n)
    for i in xrange(n):
        s[i] = np.trace(A,i)
    return s/s.sum()

def decompose_dist(pdf, curve, r=None):
    'Decompose the Euc distance matrix on curve'
    from contact_map import ContactMap, EIG, NMF_main
    map1 = ContactMap()
    curve_show(curve)
    pdf.savefig(); plt.clf()
    verts = np.array(curve)
    map1.plot_map(verts, title = "Verteces", log=False)
    pdf.savefig(); plt.clf()
    from scipy.spatial.distance import pdist, squareform
    ds = squareform(pdist(verts, 'euclidean'))

    ## transform
#    V = ds.max() - ds
#    V = ds.max() / (ds + 1)
    V = ds.max() * ((ds+1)**-2)
    map1.plot_map(V, title="Synthetic Heatmap", log=True)
    pdf.savefig(); plt.clf()

    #plt.hist(np.reshape(V,(-1)), bins=100, normed=1, facecolor='blue')
    #plt.title('Distribution of map values')
    #pdf.savefig(); plt.clf()
    plt.loglog([(i+1.0)/V.shape[0] for i in range(V.shape[0])], 
               trace_sum(V), linestyle='-.')
    plt.title('Distribution of interactions along 1D')
    plt.xlabel('Ratio of linked locations to the total length')
    plt.ylabel('Number of observed links')
    pdf.savefig(); plt.clf()

    if r == None:
        r = choose_size(pdf, V, 9)
        show('Best number of dimentions is %s\n'%r)
        r = 4
    if False: ## try PCA
        U = (V-np.mean(V.T,axis=1)).T
        Q, M = EIG(np.cov(U), r)
    else:
        Q, M = EIG(V, r)
    map1.plot_map(Q, title = 'Eig. Decomp. - Q Matrix', log=False)
    pdf.savefig(); plt.clf()
    map1.plot_map(M, title = 'Eig. Decomp. - M Matrix', log=False)
    pdf.savefig(); plt.clf()
    map1.plot_map(Q*M*Q.T, title = 'Eig. Decomp. - Recovered', log=False)
    pdf.savefig(); plt.clf()
    sep_map_show(pdf, verts, Q)

    H, S, obj = NMF_main(V, J='NMF-PoissonManifoldEqual', H=Q, S=M, r=r)
    map1.plot_map(H*S*H.T, title = 'NMF Decomp. - Recovered', log=False)
    pdf.savefig(); plt.clf()
    map1.plot_map(H, title = 'NMF Decomp. - H Matrix', log=False)
    pdf.savefig(); plt.clf()
    map1.plot_map(S, title = 'NMF Decomp. - S Matrix', log=False)
    pdf.savefig(); plt.clf()
    maxp = np.argmax(np.asarray(H),0)
    srt = np.argsort(maxp)
    sep_map_show(pdf, verts, H[:,srt])

    try:
        from sklearn.cluster import KMeans
        km = KMeans(n_clusters=r)
        H = -np.matrix(km.fit_transform(V))
        S = np.matrix(np.eye(r))
        maxp = np.argmax(np.asarray(H),0)
        srt = np.argsort(maxp)
        map1.plot_map(H, title = 'K-means Decomp. - H Matrix', log=False)
        pdf.savefig(); plt.clf()
        sep_map_show(pdf, verts, H[:,srt])
    except:
        print 'Please install SK-kit to run K-means'
        pass

def main(para):
    ''' Naive way of building space filling curves and factorize matrices '''
    pdf = PdfPages(para['ExeFile']+'plot.pdf')

    hilbert_show(pdf, 1)
    hilbert_show(pdf, 2)
#    hilbert_show(pdf, 3)
    decompose_dist(pdf, hilbert_curve(4), 4)

#    peano_show(pdf, 1)
#    peano_show(pdf, 2)
#    peano_show(pdf, 3)
#    decompose_dist(pdf, peano_curve(3), 9)
#
#    hilbert3_show(pdf, 1)
#    hilbert3_show(pdf, 2)
#    hilbert3_show(pdf, 3)
#    decompose_dist(pdf, hilbert3_curve(3), 8)
    pdf.close()

if __name__ == "__main__": main_fun(main)
