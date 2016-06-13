from tools import *
from contact_map import *

def nmf_j4a_v1(X,C,lm,H,S,minimp,maxiter,eps):
    timelist = []
    timelist.append(('Start', time.time()))
    if type(X) == type(''): X = np.load(X)
    else: X = np.array(X)
    from numpy import dot, ones
    obj = [float("infinity")]
    I = ones(X.shape)
    I[np.isnan(X)] = 0 ## missing data
    X[I == 0] = 0
    E,D = mark2network(C, lm)
    G = H.copy() ## no bias
    timelist.append(('End of init', time.time()))
    for iter in xrange(maxiter):
        R = dot(dot(G,S), G.T) + eps;
        timelist.append(('Calulate R (1)', time.time()))
        obj.append((I*(R - X*np.log(R))).sum() + \
                   lm * np.trace(dot(dot(G.T,(D-E)), G)))
        timelist.append(('Calulate obj', time.time()))
        if iter == 0: print 'Initial objective is', obj[-1]
        strobj = 'Current objective is %s; '%obj[-1]
        print strobj,
        if abs(obj[-2]-obj[-1]) <= abs(minimp*obj[-1]): break
        G *= (dot(I*X/R, dot(G,S)) + dot(E,G)) / \
             (dot(I, dot(G,S)) + dot(D,G) + eps)
        timelist.append(('Calulate G', time.time()))
        ## normalized row of H and update bias
        B = dot(G,S).mean(axis=1) ## bias vector
        B /= B.sum()/(B>0).sum() ## normalize
        B[B==0] = 1
        timelist.append(('Calulate B', time.time()))
        H = G / B[:,np.newaxis]
        ## normalize column of H and update S
        h = H.mean(axis=0)
        h /= h.mean() ## this is important!
        S *= dot(h.T, h)
        H /= h[np.newaxis,:]
        timelist.append(('Calulate H', time.time()))
        ## add the bias back to H for reconstruction
        G = H * B[:,np.newaxis]
        R = dot(dot(G,S), G.T) + eps
        timelist.append(('Calulate R (2)', time.time()))
        S *= dot(dot(G.T, I*X/R), G) / dot(dot(G.T, I), G)
        timelist.append(('Calulate S', time.time()))
        print '\r'*len(strobj),
    print ''
    #return G,S,obj
    return timelist

def nmf_j4a_v2(X,C,lm,H,S,minimp,maxiter,eps):
    timelist = []
    timelist.append(('Start', time.time()))
    if type(X) == type(''): X = np.load(X)
    else: X = np.array(X)
    from numpy import dot, ones
    obj = [float("infinity")]
    I = ones(X.shape)
    I[np.isnan(X)] = 0 ## missing data
    X[I == 0] = 0
    E,D = mark2network(C, lm)
    G = H.copy() ## no bias
    timelist.append(('End of init', time.time()))
    for iter in xrange(maxiter):
        R = dot(dot(G,S), G.T) + eps;
        timelist.append(('Calulate R (1)', time.time()))
        obj.append((I*(R - X*np.log(R))).sum() + lm * np.trace(dot(dot(G.T,(D-E)), G)))
        timelist.append(('Calulate obj', time.time()))
        if iter == 0: print 'Initial objective is', obj[-1]
        strobj = 'Current objective is %s; '%obj[-1]
        print strobj,
        if abs(obj[-2]-obj[-1]) <= abs(minimp*obj[-1]): break
        GS = dot(G,S)
        G *= (dot(X/R, GS) + dot(E,G)) / (dot(I, GS) + dot(D,G) + eps)
        timelist.append(('Calulate G', time.time()))
        ## normalized row of H and update bias
        B = dot(G,S).mean(axis=1) ## bias vector
        B /= B.sum()/(B>0).sum() ## normalize
        B[B==0] = 1
        timelist.append(('Calulate B', time.time()))
        H = G / B[:,np.newaxis]
        ## normalize column of H and update S
        h = H.mean(axis=0)
        h /= h.mean() ## this is important!
        S *= dot(h.T, h)
        H /= h[np.newaxis,:]
        timelist.append(('Calulate H', time.time()))
        ## add the bias back to H for reconstruction
        G = H * B[:,np.newaxis]
        R = dot(dot(G,S), G.T) + eps
        timelist.append(('Calulate R (2)', time.time()))
        S *= dot(dot(G.T, X/R), G) / dot(dot(G.T, I), G)
        timelist.append(('Calulate S', time.time()))
        print '\r'*len(strobj),
    print ''
    #return G,S,obj
    return timelist

def nmf_j4a_v3(X,C,lm,H,S,minimp,maxiter,eps):
    timelist = []
    timelist.append(('Start', time.time()))
    if type(X) == type(''): X = np.load(X)
    else: X = np.array(X)
    from numpy import dot, ones
    obj = [float("infinity")]
    I = ones(X.shape)
    I[np.isnan(X)] = 0 ## missing data
    X[I == 0] = 0
    E,D = mark2network(C, lm)
    G = H.copy() ## no bias
    timelist.append(('End of init', time.time()))
    for iter in xrange(maxiter):
        R = dot(dot(G,S), G.T) + eps;
        timelist.append(('Calulate R (1)', time.time()))
        obj.append((I*(R - X*np.log(R))).sum() + \
                   lm * np.trace(dot(dot(G.T,(D-E)), G)))
        timelist.append(('Calulate obj', time.time()))
        if iter == 0: print 'Initial objective is', obj[-1]
        strobj = 'Current objective is %s; '%obj[-1]
        print strobj,
        if abs(obj[-2]-obj[-1]) <= abs(minimp*obj[-1]): break
        G *= dot(I*X/R, dot(G,S)) / (dot(I, dot(G,S)) + eps)
        timelist.append(('Calulate G', time.time()))
        ## normalized row of H and update bias
        B = dot(G,S).mean(axis=1) ## bias vector
        B /= B.sum()/(B>0).sum() ## normalize
        B[B==0] = 1
        timelist.append(('Calulate B', time.time()))
        H = G / B[:,np.newaxis]
        ## normalize column of H and update S
        h = H.mean(axis=0)
        h /= h.mean() ## this is important!
        S *= dot(h.T, h)
        H /= h[np.newaxis,:]
        timelist.append(('Calulate H', time.time()))
        ## add the bias back to H for reconstruction
        G = H * B[:,np.newaxis]
        R = dot(dot(G,S), G.T) + eps
        timelist.append(('Calulate R (2)', time.time()))
        S *= dot(dot(G.T, I*X/R), G) / dot(dot(G.T, I), G)
        timelist.append(('Calulate S', time.time()))
        print '\r'*len(strobj),
    print ''
    #return G,S,obj
    return timelist

def nmf_j4a_v4(X,C,lm,H,S,minimp,maxiter,eps):
    from numexpr import evaluate
    timelist = []
    timelist.append(('Start', time.time()))
    if type(X) == type(''): X = np.load(X)
    else: X = np.array(X)
    from numpy import dot, ones
    obj = [float("infinity")]
    I = ones(X.shape)
    I[np.isnan(X)] = 0 ## missing data
    X[I == 0] = 0
    E,D = mark2network(C, lm)
    L = D-E
    G = H.copy() ## no bias
    timelist.append(('End of init', time.time()))
    for iter in xrange(maxiter):
        R = dot(dot(G,S), G.T) + eps;
        timelist.append(('Calulate R (1)', time.time()))
        temp1 = evaluate('sum(I*(R - X*log(R)))')
        obj.append(temp1 + lm * np.trace(dot(dot(G.T,L), G)))
        timelist.append(('Calulate obj', time.time()))
        if iter == 0: print 'Initial objective is', obj[-1]
        strobj = 'Current objective is %s; '%obj[-1]
        print strobj,
        if abs(obj[-2]-obj[-1]) <= abs(minimp*obj[-1]): break
        GS = dot(G,S)
        temp1 = dot(X/R, GS)
        temp2 = dot(E,G)
        temp3 = dot(I, GS)
        temp4 = G * np.diag(D)[:,np.newaxis]
        evaluate('G * (temp1 + temp2) / (temp3 + temp4 + eps)', out=G)
        timelist.append(('Calulate G', time.time()))
        ## normalized row of H and update bias
        B = dot(G,S).mean(axis=1) ## bias vector
        B /= B.sum()/(B>0).sum() ## normalize
        B[B==0] = 1
        timelist.append(('Calulate B', time.time()))
        H = G / B[:,np.newaxis]
        ## normalize column of H and update S
        h = H.mean(axis=0)
        h /= h.mean() ## this is important!
        temp1 = dot(h.T, h)
        evaluate('S * temp1', out=S)
        evaluate('H / h', out=H)
        timelist.append(('Calulate H', time.time()))
        ## add the bias back to H for reconstruction
        G = H * B[:,np.newaxis]
        R = dot(dot(G,S), G.T)
        timelist.append(('Calulate R (2)', time.time()))
        XdR = evaluate('X/(R+eps)')
        temp1 = dot(dot(G.T, XdR), G)
        temp2 = dot(dot(G.T, I), G)
        evaluate('S * temp1 / temp2', out=S)
        timelist.append(('Calulate S', time.time()))
        print '\r'*len(strobj),
    print ''
    #return G,S,obj
    return timelist

def main(para):
    map1 = ContactMap()
    map1.load('demo')
    r = 50
    X = np.array(map1.contact_map)
    C = np.array(map1.frag_chr, dtype='float')
    lm = 20
    from numpy.random import rand
    eps = 1e-20
    H = rand(X.shape[0], r) + eps
    S = rand(r,r)+np.eye(r)
    minimp = 0
    maxiter = 10
    print X.shape, C.shape, H.shape, S.shape
    timeseq = []
    timeseq.append(('v1', nmf_j4a_v1(X.copy(),C,lm,H.copy(),S.copy(),minimp,maxiter,eps)))
    timeseq.append(('v2', nmf_j4a_v2(X.copy(),C,lm,H.copy(),S.copy(),minimp,maxiter,eps)))
    timeseq.append(('v3', nmf_j4a_v3(X.copy(),C,lm,H.copy(),S.copy(),minimp,maxiter,eps)))
    timeseq.append(('v4', nmf_j4a_v4(X.copy(),C,lm,H.copy(),S.copy(),minimp,maxiter,eps)))
    records = []
    for vs, ts in timeseq:
        show(vs)
        show(ts[-1][1] - ts[0][1])
        show()
        itemseq = {}
        for i in xrange(1,len(ts)):
            seq, sec = ts[i]
            used = itemseq.get(seq, 0)
            itemseq[seq] = used + sec - ts[i-1][1]
        records.append(itemseq)
    show()
    for seq in records[0]:
        show(seq)
        for itemseq in records:
            show(itemseq[seq])
        show()

if __name__ == '__main__': main_fun(main)
