from tools import *
import numpy as np

def exp_in_grp(exp='exp-Gm12878', grp='3d-GM12878', mode='Max'):
    ## file 1
    data1 = []
    exp_file = open('loci_for_'+exp+'.txt', 'r')
    head = exp_file.readline()
    for line in exp_file:
        ele = line.split('\t')
        rpkm = float(ele[4])
#        if rpkm >= 1:
        data1.append((ele[3], int(ele[0]), rpkm)) 
    exp_file.close()
    print 'We have', len(data1), 'genes in', exp
    ## file 2
    data2 = []
    grp_file = open(grp+'_groups.txt', 'r')
    for line in grp_file:
        data2.append([float(v) for v in line.split('\t')[3:]])
    grp_file.close()
    print 'We have', len(data2), 'bins with', len(data2[0]), 'groups in', grp
    ## check
    H = np.array(data2)
    if mode == 'Max':
        E = np.zeros(H.shape[0])
        for g,i,e in data1:
            if E[i] < e:
                E[i] = e
            #E[i] += e
        E = np.argsort(np.argsort(E))
        E = np.array(E+1, dtype='float')
        E /= E.max()
        val = []
        for j in xrange(H.shape[1]):
            idx = H[:,j] > np.percentile(H[:,j], 50)
            #idx = H[:,j] > 0
            val.append(E[idx].var()/E.var())
        return min(val) #[i for i in xrange(len(val)) if min(val)==val[i]]
    elif mode == 'All':
        gene, bidx, expv = zip(*data1)
        E = np.argsort(np.argsort(expv)) ## get the rank
        E = np.array(E+1, dtype='float')
        E /= E.max()
        val = []
        for j in xrange(H.shape[1]):
            idx = H[:,j] >= np.percentile(H[:,j], 90)
    #        idx = H[:,j] > 1
            for e,i in zip(E.tolist(),bidx):
                if idx[i]:
                    val.append(e)
        return np.array(val).var()/E.var()
    else:
        return -1

def highexp_in_grp(exp='highexp-Gm12878', grp='3d-GM12878'):
    ## file 1
    data1 = []
    exp_file = open('loci_for_'+exp+'.txt', 'r')
    head = exp_file.readline()
    for line in exp_file:
        ele = line.split('\t')
        data1.append((ele[3], int(ele[0])))
    exp_file.close()
    print 'We have', len(data1), 'genes in', exp
    ## file 2
    data2 = []
    grp_file = open(grp+'_groups.txt', 'r')
    for line in grp_file:
        data2.append([float(v) for v in line.split('\t')[3:]])
    grp_file.close()
    print 'We have', len(data2), 'bins with', len(data2[0]), 'groups in', grp
    ## check
    H = np.array(data2)
    gene, bidx = zip(*data1)
    val = []
    for j in xrange(H.shape[1]):
        real = np.zeros(H.shape[0])
        for i in bidx:
            real[i] = 1
        auc, X, Y, C = performance(real.tolist(), H[:,j].tolist())
        val.append(auc)
    best = np.argsort(np.array(val))[::-1]
    show(exp)
    show(np.array(val)[best[:5]].tolist())
    show()

def shared_groups(grps, cutoff=0.3):
    data = []
    for grp in grps:
        H = []
        grp_file = open(grp+'_groups.txt', 'r')
        for line in grp_file:
            H.append([float(v) for v in line.split('\t')[3:]])
        grp_file.close()
        data.append(np.array(H))
    for j in xrange(10):
        s_and = 1; s_or = 0
        for H in data:
            h = H[:,j]
            c = h > np.percentile(h, 50)
            s_and = np.logical_and(s_and, c)
            s_or  = np.logical_or(s_or, c)
        print j, (s_and>0).sum()/float((s_or>0).sum())
    ref = np.array([H.shape[1] for H in data])
    cc = ref-1
    out = []
    while True:
        s_and = 1; s_or = 0
        for i in xrange(len(cc)):
            h = data[i][:,cc[i]]
            c = (h > np.percentile(h, 50))
            s_and = np.logical_and(s_and, c)
            s_or  = np.logical_or(s_or, c)
            if s_and.sum()/float(s_or.sum()) < cutoff:
                if i < len(cc)-1:
                    for k in xrange(i+1, len(cc)):
                        cc[k] = 0
                break
        if s_and.sum()/float(s_or.sum()) >= cutoff:
            out.append([s_and.sum()/float(s_or.sum())] + cc.tolist())
        if cc.sum() == 0:
            break
        for i in xrange(len(cc)-1, -1, -1):
            if cc[i] > 0:
                cc[i] -= 1
                break
            else:
                cc[i] = ref[i]-1
    print 'We have', len(out), 'groups'
    out.sort(reverse=True)
    unique = []
    bins = []
    for val1 in out:
        exist = False
        for val2 in unique:
            for i,j in zip(val1[1:], val2[1:]):
                if i == j:
                    exist = True
        if not exist:
            unique.append(val1)
            h = []
            for i in xrange(len(val1)-1):
                h.append(data[i][:,val1[i+1]])
            bins.append(np.array(h))
            show(val1, True)
    return bins

def coexpressed(exps, sel):
    data = {}
    for exp in exps:
        exp_file = open('loci_for_'+exp+'.txt', 'r')
        head = exp_file.readline()
        for line in exp_file:
            ele = line.split('\t')
            rpkm = float(ele[4])
            if rpkm > 0:
                gene = (ele[3], int(ele[0]))
                if gene in data:
                    rps = data[gene]
                    rps.append(rpkm)
                    data[gene] = rps
                else:
                    data[gene] = [rpkm]
        exp_file.close()
    genes = []
    shared = []
    for gene in data:
        if len(data[gene]) == len(exps):
            gname, gbin = gene
            if gbin in sel:
                genes.append(gene)
                shared.append(data[gene])
                show(gene)
                show(data[gene])
                show()
    print 'We have', len(genes), 'selected genes from', len(sel), 'bins'
    v = np.array(shared)
    c = np.corrcoef(v)
    a = np.tril(c)
    return len(sel), len(genes), a.min(), a.mean(), a.max()

def run1():
    show('''
    Variance of expression values in contact groups
    ''', True)
    grps = ['3d-GM12878', '3d-K562', '3d-hESC', '3d-IMR90'] 
    exps = ['exp-Gm12878', 'exp-K562', 'exp-H1hesc', 'exp-Imr90']
    show('')
    show(exps)
    show()
    for grp in grps:
        show(grp)
        for exp in exps:
            show(exp_in_grp(exp, grp))
        show()

def run2():
    show('''
    Highly expressed genes in matched cell lines
    ''', True)
    mtch = [('3d-GM12878', 'highexp-Gm12878'),
            ('3d-K562', 'highexp-K562'),
            ('3d-hESC', 'highexp-H1hesc'),
            ('3d-IMR90', 'highexp-Imr90')]
    for grp, exp in mtch:
        highexp_in_grp(exp, grp)

def run3():
    show('''
    Conserved structure
    ''', True)
    grps = ['3d-GM12878', '3d-K562', '3d-hESC', '3d-IMR90'] 
    exps = ['exp-Gm12878', 'exp-K562', 'exp-H1hesc', 'exp-Imr90']
    out = shared_groups(grps)
    h1 = out[0]
    h2 = out[1]
    sel = np.logical_and((h1<1).sum(0)==0, (h2<1).sum(0)==0)
#    sel = (h2<1).sum(0)==0
    idx = np.arange(len(sel))[sel]
#    for i in idx:
#        show(i)
#        show(h1[:,i].tolist())
#        show(h2[:,i].tolist())
#        show()
    show(coexpressed(exps, idx), True)

def main(para):
#    run1()
#    run2()
    run3()

if __name__ == '__main__': main_fun(main)
