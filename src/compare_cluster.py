from tools import *
import numpy as np
from contact_map import ContactMap, IterateCorrect
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from sklearn.cluster import KMeans

def get_syn_map(para, bin_size=3200, with_bias=True):
    pdf = PdfPages(para['ExeFile']+'plot1.pdf')
    ## prepare
    map1 = ContactMap('Syn3D')
    map1.genome_info(para['DataPath']+'/Tjong2012GR/yeast_chr_len-Tjong.txt')
    map1.create_binnedmap(3200) ## fixed
    map2 = map1.duplicate()
    map3 = map1.duplicate()

    ## obtain links from PDB
    link_map = np.load('syn_link.npy')
    if with_bias: ## add random bias
        np.random.seed(0)
        bias = np.random.random(link_map.shape[0])
        link_map *= np.outer(bias, bias)
        print link_map.min(), link_map.max(), link_map.mean()
        link_map = np.floor(link_map) ## sampling bias
    map1.contact_map = np.matrix(link_map, dtype='float')

    output = open('syn_yeast_map_bin%s%s.txt'%(bin_size, 'bias' if with_bias else ''), 'w')
    ch = map1.frag_chr
    po = (map1.frag_sta+map1.frag_end)/2
    for i in xrange(link_map.shape[0]):
        for j in xrange(link_map.shape[1]):
            if link_map[i,j] > 0:
                output.write('%s\t%s\t%s\t%s\t0\t%s\t1e-10\t1e-10\n'%(ch[i], po[i], ch[j], po[j], link_map[i,j]))
        output.write('\n')
    output.close()

    map1.get_interactions()
    map1.create_binnedmap(bin_size)
    map1.mask_diag()
    map1.plot_map(title='Heatmap for the number of links')
    pdf.savefig(); plt.clf();
    map1.decompose('NND')
    idx, names = map2.get_locations(map1.output_groups(), st=0, ch=0, po=1, nm=0, add=0)
    dist_map = np.load('syn_dist.npy')
    dist = dist_map[idx,:][:,idx]
    map1.plot_map(dist, title='Heatmap for the average distances', log=False)
    pdf.savefig(); plt.clf();
    pdf.close()
    return map1, dist

def decompose_map(map1, method, r=40, out='inter'):
    map1.reset_solution()
    if method == 'EIG':
        map1.decompose('EIG', dim_num=r)
    elif method == 'PCA':
        map1.decompose('PCA', dim_num=r)
    elif method == 'ICE':
        map1.decompose('ICE', dim_num=r)
    elif method == 'K-means':
        from k_means_pdist import kmeanssample
        DIST = -np.array(map1.contact_map) ## simi to dist
        centres, xtoc, dist = kmeanssample(DIST, np.eye(DIST.shape[0]), r, nsample=0, delta=0.001, maxiter=20, verbose=0)
        map1.contact_group = -np.matrix(dist) ## dist to simi
    elif method == '3D-K-means':
        km = KMeans(n_clusters=r)
        dfile = 'pdb.txt'
        pb, vx = map1.get_locations(dfile, st=1, ch=0, po=1, nm=2, add=0)
        pb, vy = map1.get_locations(dfile, st=1, ch=0, po=1, nm=3, add=0)
        pb, vz = map1.get_locations(dfile, st=1, ch=0, po=1, nm=4, add=0)
        X = np.zeros((map1.contact_map.shape[0], 3))
        C = np.zeros(map1.contact_map.shape[0])
        for i,x,y,z in zip(pb,vx,vy,vz):
            X[i,0] = x
            X[i,1] = y
            X[i,2] = z
            C[i] += 1
        C[C==0] = 1
        X /= C[:,np.newaxis]
        map1.contact_group = -np.matrix(km.fit_transform(X))
    elif method == 'NMF':
        map1.decompose('NND', dim_num=r)
        map1.decompose('NMF-Gaussian', dim_num=r)
        map1.contact_group = np.dot(map1.contact_group, map1.group_map)
    elif method == 'BNMF':
        map1.decompose('NND', dim_num=r)
        map1.decompose('NMF-PoissonManifoldEqual', dim_num=r, par_lam=0)
        map1.contact_group = np.dot(map1.contact_group, map1.group_map)
    elif method == 'Random':
        n = map1.contact_map.shape[0]
        map1.contact_group = np.zeros((n,r))
        from math import ceil
        size = int(ceil(n/float(r)))
        for i in xrange(n):
            map1.contact_group[i, i/size] = 1
    elif method == 'Armatus':
        from run_armatus import Armatus
        map1.save()
        map2 = Armatus('../tools/armatus2.1/armatus', name=map1.name)
        map2.load()
        map2.decompose()
        map1.contact_group = map2.contact_group
    elif method == 'TAD':
        from run_domaincall import DomainCall
        map1.save()
        map2 = DomainCall('../tools/domaincall/', name=map1.name)
        map2.load()
        map2.decompose()
        map1.contact_group = map2.contact_group
    else:
        raise ValueError('Unknow method name '+method)

def get_distances(map1, dist, case='all'):
    dist_intra = dist.copy()
    dist_inter = dist.copy()
    map1.mask_trans(dist_intra, 0)
    map1.mask_cis(dist_inter, 0)
    select = []
    H = np.asarray(map1.contact_group)
    M = np.argmax(H, 1) ## maximum membershop to assign clusters
    for j in xrange(H.shape[1]):
        val = H[:,j]
        srt = np.argsort(val)
        idx1 = (M == j)
        if case == 'all':
            D1 = dist[idx1,:][:,idx1]
        elif case == 'intra':
            D1 = dist_intra[idx1,:][:,idx1]
        elif case == 'inter':
            D1 = dist_inter[idx1,:][:,idx1]
        else:
            raise ValueError('Unknown case: %s, eg: intra'%case)
        d1 = D1[np.triu_indices(D1.shape[0], k=1)]
        select += d1[d1>0].tolist()
    return select

def plot_dists(para, methods, dists):
    pdf = PdfPages(para['ExeFile']+'plot2.pdf')
    fig = plt.figure()
    for i in xrange(len(methods)):
        axis = fig.add_subplot(len(methods)*100+11+i)
        axis.hist(dists[i], range(100, 1300, 50), color='white')
        axis.set_xlim([100, 1300])
        axis.set_title(methods[i]+' (%.f$\pm$%.f nm)'%tuple(mean_std(dists[i])))
        axis.tick_params(labelsize=8)
    fig.tight_layout()
    fig.subplots_adjust(left=0.15, right=0.85)
    fig.savefig(pdf, format='pdf')
    pdf.close()

def main_par(para):
    m = para['Method']
    r = int(para['Cluster'])
    map1, dist = get_syn_map(para, with_bias=True)
    decompose_map(map1, method=m, r=r)
    for case in ['all', 'intra', 'inter']:
        d = get_distances(map1, dist, case)
        if len(d) == 0:
            show(0)
        else:
            show(mean_std(d)[0])
    show()

def main(para):
    if 'Method' in para and 'Cluster' in para:
        main_par(para)
        return
    methods = ['PCA', 'ICE', 'K-means', '3D-K-means', 'NMF', 'BNMF']
#    methods = ['Random', 'K-means', 'NMF', 'BNMF']
    for bias in [False, True]:
        show('With bias' if bias else 'Without bias', True)
        map1, dist = get_syn_map(para, with_bias=bias)
        if bias:
            os.system('cp output.pdb.1Dto3D.temp.biased.txt pdb.txt')
        else:
            os.system('cp output.pdb.1Dto3D.temp.txt pdb.txt')
        n = map1.contact_map.shape[0]
        show('#Bins')
        show(methods, True)
        for r in range(10,201,10):
            show(r)
            for method in methods:
                decompose_map(map1, method=method, r=r)
                for case in ['all', 'intra', 'inter']:
                    d = get_distances(map1, dist, case)
                    if len(d) == 0:
                        show(0)
                    else:
                        show(mean_std(d)[0])
            show()
        show()

if __name__ == '__main__': main_fun(main)
