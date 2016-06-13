from tools import *
import numpy as np
from contact_map import ContactMap
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def run1(name='demo'):
    show('''
    Compare decomposation results during iterations
    ''', True)
    map1 = ContactMap(name)
    assert map1.load()
    map1.decompose_auto()
    n,r = map1.contact_group.shape
    print n, r
    map1.decompose('NND', dim_num=r, A=map1.get_null_map())
    map2 = map1.duplicate()
    map2.decompose(dim_num=r)
    ref = np.array(map2.contact_group)
    show('Iter\tObj.\t# out of %s'%r)
    show('Corr. Mean\tCorr. STD\tMetric\n')
    from contact_map import gini_impurity
    for i in [1, 5, 10, 50] + range(100, 1201, 100):
        show(i)
        map2 = map1.duplicate()
        show(map2.decompose(dim_num=r, max_iter=i, stop_thrd=0))
        ## match dims
        s = []
        v = []
        for j in xrange(r):
            srt, val, pval = map2.test_enrichment(ref[:,j], 'PCC')
            s.append(srt[0])
            v.append(val[srt[0]])
        show(len(set(s)))
        show(mean_std(v))
        show(gini_impurity(np.diag(map2.group_map)))
        show()

def run2(name='demo'):
    show('''
    Show the distribution of objective values
    ''', True)
    map1 = ContactMap(name)
    assert map1.load()
    map1.decompose_auto(max_iter=50)
    n,r = map1.contact_group.shape
    map2 = map1.duplicate()
    objs = []
    for i in xrange(100):
        map2.reset_solution()
        objs.append(map2.decompose(dim_num=r, max_iter=50))
    max_obj = int(max(objs))+1
    min_obj = int(min(objs))-1
    bins = range(min_obj, max_obj, int((max_obj-min_obj)/10))
    show('')
    show(bins, True)
    show('Frequency')
    show(histogram(objs, bins), True)
    show('Mean\tSTD\n')
    m, s = mean_std(objs)
    show([m, s], True)

    map2.decompose('Null', dim_num=r)
    obj1 = map2.decompose(dim_num=r)
    show('Use Null init has obj.')
    show(obj1)
    show('and Z-score is')
    show((obj1-m)/s, True)

def run3(para, name='demo'):
    show('''
    Compare objective values in NMF and average distances in PDB.
    ''', True)
    pdf = PdfPages(para['ExeFile']+'plot.pdf')
    map1 = ContactMap('Syn3D')
    map1.genome_info(para['DataPath']+'/Tjong2012GR/yeast_chr_len-Tjong.txt')
    map1.create_binnedmap(3200)
    map2 = map1.duplicate()
    map1.contact_map = np.matrix(np.load('syn_link.npy'), dtype='float')
    map1.get_interactions()
    map1.create_binnedmap(32000)
    map1.mask_diag()
    paras = map1.decompose_auto()
    r = paras[-1][0]
    map3 = map1.duplicate()
    show(r)
    show('is the selected cluster number\n')
    print map1.contact_map.shape
    idx, names = map2.get_locations(map1.output_groups(), st=0, ch=0, po=1, nm=0, add=0)
    print len(idx)
    dist_map = np.load('syn_dist.npy')
    show(dist_map.mean())
    show('is the average of all bins\n')
    dist = dist_map[idx,:][:,idx] ## distance among bins
    show(dist.mean())
    show('is the average of selected bins\n')
    inum = []; objs = []; avgs = []; objs3 = []; avgs3 = []
    show('\tObjective function values\tAverage intra-cluster distances\tcase2\n')
    map1.reset_solution()
    map1.decompose('NND', dim_num=r)
    map3.reset_solution()
#    map3.decompose('NND', dim_num=r)
    icc = 0; step = 20
    for i in xrange(100):
        icc += step
        show(icc)
        inum.append(icc)
        obj = map1.decompose(dim_num=r, par_lam=1, max_iter=step, stop_thrd=0)
        obj3 = map3.decompose(dim_num=r, par_lam=1, max_iter=step, stop_thrd=0)
        map1.sort_groups('diagnal')
        show(obj)
        avg = []; avg3 = []
        for j in xrange(r):
            idx1 = np.asarray(map1.contact_group)[:,j] > float(map1.contact_group[:,j].mean())
            D1 = dist[idx1,:][:,idx1]
            d1 = D1[np.triu_indices(D1.shape[0], k=1)]
            avg += d1.tolist()

            idx3 = np.asarray(map3.contact_group)[:,j] > float(map3.contact_group[:,j].mean())
            D3 = dist[idx3,:][:,idx3]
            d3 = D3[np.triu_indices(D3.shape[0], k=1)]
            avg3 += d3.tolist()
        show(mean_std(avg)[0])
        show(mean_std(avg3)[0])
        show()
        objs.append(obj)
        objs3.append(obj3)
        avgs.append(mean_std(avg)[0])
        avgs3.append(mean_std(avg3)[0])
    plt.plot(inum, objs, 'r-', label='NNDSVD Initialization')
#    plt.plot(inum, objs3, 'b--', label='Random Initialization')
    plt.legend()
    plt.xlabel('Number of iterations')
    plt.ylabel('Objective function values for NMF')
    pdf.savefig(); plt.clf();
    plt.plot(objs, avgs, 'r.', label='NNDSVD Initialization')
#    plt.plot(objs3, avgs3, 'b+', label='Random Initialization')
    plt.legend()
    plt.xlabel('Objective function values for NMF')
    plt.ylabel('Average intra cluster distances (nm)')
    pdf.savefig(); plt.clf();
    show('\nCorrelation of objective with the average distances\n')
    show('Pearson Coef.')
    show(correlation(objs, avgs), True)
    show('Spearman Rank Coef.')
    show(correlation(objs, avgs, rank=True), True)
    map1.plot_submap()
    pdf.savefig(); plt.clf()
    map3.plot_submap()
    pdf.savefig(); plt.clf()
    pdf.close()

def run4(name='demo'):
    show('''
    Compare the change of clusters under different resolutions and iterations
    ''', True)
    map1 = ContactMap(name)
    assert map1.load()
    dims = [10,15,20,25,30,35,40,45,50]
    map1.decompose_auto(dim_num=dims)
    memb = np.array(map1.contact_group*map1.group_map)
    show('')
    show(dims, True)
    for ratio in xrange(1,4):
        map2 = map1.duplicate()
        map2.get_interactions()
        map2.create_binnedmap(binsize=map1.get_binsize()*ratio)
        map2.mask_diag()
        paras = map2.decompose_auto(dim_num=dims)
        bins, vals = zip(*paras)
        idx, val = map1.get_locations(map2.output_groups(), st=0, ch=0, po=1)
        newv = np.array(map2.contact_group*map2.group_map)
        show(map2.get_binsize())
        show(vals, True)

    show(map1.contact_group.shape, True)
    show('\n# of Iter.')
    show(dims, True)
    for i in xrange(10):
        it = i*100
        paras = map1.decompose_auto(max_iter=it, update=True, dim_num=dims)
        bins, vals = zip(*paras)
        show(it)
        show(vals, True)

def run5(name='demo'):
    show('''
    Mapping clusters by changing the number of total clusters 
    ''', True)
    map1 = ContactMap(name)
    map2 = ContactMap(name)
    assert map1.load()
    assert map2.load()
#    dims = [10,20,30,40,50,60,70,80]
    dims = range(5, 31, 1)
    show('Bin Size\tMetric')
    map1.decompose_auto(dim_num=30)
    full = np.arange(map1.contact_group.shape[1])
    show(full.tolist(), True)
    from contact_map import gini_impurity
    for r in dims:
        show(r)
        map2.decompose_auto(dim_num=r)
        show(gini_impurity(np.diag(map2.group_map)))
        match = map1.best_cor(map2, dims=True)
        dt = {}
        for i,j in match:
            dt[i] = j
        for i in full:
            if i in dt:
                show(dt[i])
            else:
                show('')
        show()

def main(para):
    #run1()
    #run2()
    run3(para)
    run4()
    #run5()

if __name__ == '__main__': main_fun(main)
