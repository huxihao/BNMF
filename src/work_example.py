from contact_map import *
from tools import *

def print_mat(name='Matrix', A=[]):
    show(name)
    show('=', True)
    if len(A.shape) == 1:
        show(A.tolist(), True)
    if len(A.shape) == 2:
        for i in xrange(A.shape[0]):
            show(A[i,:].tolist(), True)

def run1(par_lam=0):
    show('''
    Produce the example discussed in the paper
    ''', True)
    H = np.matrix([[0.5,0],[0.5,0],[0,1]], dtype='float')
    S = np.matrix([[2,0],[0,1]], dtype='float')
    B = np.diagflat(np.array([2,2,1], dtype='float'))
    W = S*H.T
    A = H*W
    X = B*A*B
    print_mat('A', A)
    print_mat('B', B)
    print_mat('H', H)
    print_mat('S', S)
    print_mat('W', W)
    initH, initS = EIG(X,2)
    initH = np.abs(np.asarray(initH))+1
    initS = np.abs(np.asarray(initS))+1
    print_mat('initH', initH)
    print_mat('initS', initS)
    newG, newS, obj = nmf_j4a(
        X = np.array(X), 
        C = np.array([0,1,1]),
        lm = par_lam,
        H = np.asarray(initH),
        S = np.asarray(initS),
        minimp = 0,
        maxiter = 100,
        eps = 1e-30)
    print len(obj), 'iterations'
    newB = np.dot(newG, newS).mean(axis=1)
    newB /= newB.mean()
    newH = newG / newB[:,np.newaxis]
    newA = np.dot(np.dot(newH,newS),newH.T)
    newW = np.dot(newS, newH.T)
    print_mat('newB', np.round(newB,2))
    print_mat('newH', np.round(newH,2))
    print_mat('newS', np.round(newS,2))
    print_mat('newW', np.round(newW,2))
    print_mat('newA', np.round(newA,2))

def run2(k=5):
    pdf = PdfPages('demo-plot.pdf')
    map1 = ContactMap()
    map1.load('demo')
    map1.create_binnedmap(60e3)
    map1.clear(['_nmf.npz'])

    n = map1.contact_map.shape[0]
    A = 0; step = n/k
    for i in xrange(k+1):
        h = np.zeros(n)
        h[min(n-1,step*i):min(n,step*i+step)] = 1
        A += np.outer(h,h)
        map1.frag_chr[h==1] = i+1
    map1.contact_map = A+0.01

    map1.plot_map()
    pdf.savefig(); plt.clf()
    map1.decompose_auto(dim_num=range(2,20,1), max_iter=100, plot=pdf)
    map1.decompose_auto(beta=1, max_iter=100, plot=pdf)
    map1.plot_submap(vmax=None, vmin=None)
    pdf.savefig(); plt.clf()
    pdf.close()

def main(para):
    run1()
    run2()

if __name__ == "__main__": main_fun(main)
