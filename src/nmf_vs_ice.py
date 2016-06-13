from tools import *
from contact_map import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def plot_maps_1(para, name, pdf, dim=50):
    '''Used to check the consistant'''
    path = para['DataPath'] + '/Imakaev2012NM'
    genome = para['DataPath'] + '/hg19_chr_len.txt'
    map1 = ContactMap()
    map1.genome_info(genome)
    map2 = ContactMap()
    map2.genome_info(genome)

    map1.create_densemap([path+'/'+name+'_map-res1000k-ic.hdf5'], reso=1e6)
    N1= np.nan_to_num(map1.contact_map)

    idx, val = map1.get_locations(path+'/'+name+'_map-res1000k-ic-bias.txt', st=0, ch=0, po=1, nm=3, add=0)
    b0 = np.ones(map1.contact_map.shape[0])
    b0[idx] = np.array([float(v) for v in val])
    map1.contact_map = np.multiply(map1.contact_map, np.outer(b0, b0))

    map2.create_densemap([path+'/'+name+'_map-res1000k.hdf5'], reso=1e6)
    map2.mask_diag()
#    map2.mask_short(0.5)
#    map2.mask_low(1)
    mask_idx = np.asarray(map1.contact_map).sum(1) == 0
    map2.contact_map[mask_idx,:] = 0
    map2.contact_map[:,mask_idx] = 0
    map2.trim_high(0.05)

    ## Check two raw maps
    M1 = np.nan_to_num(map1.contact_map)
    M2 = np.nan_to_num(map2.contact_map)
    show('Differences among two maps')
    show((M1-M2).max())
    show((M1-M2).mean())
    show((M1-M2).min(), True)
#    map1.plot_map(M1-M2, log=False, title='Non Factorable Parts')
#    pdf.savefig(); plt.clf();
#    plt.hist(np.asarray((M1-M2)[np.abs(M1-M2)>10]).reshape(-1))
#    pdf.savefig(); plt.clf();

#    map1.decompose(method='Correct')
#    map2.decompose(method='Correct')

    ## Decompose
    map1.decompose(method='NND', dim_num=dim)
    map1.decompose(method='NMF-PoissonManifoldEqual', dim_num=dim, max_iter=50)

    map2.decompose(method='NND', dim_num=dim)
    map2.decompose(method='NMF-PoissonManifoldEqual', dim_num=dim, max_iter=50)

    ## Compare results
    b1 = np.array(map1.bias_vector)
    b2 = np.array(map2.bias_vector)

    r1 = correlation(b0.tolist(), b1.tolist())
    r2 = correlation(b0.tolist(), b2.tolist())
    r3 = correlation(b1.tolist(), b2.tolist())
    show(['B0 v.s. B1', r1], True)
    show(['B0 v.s. B2', r2], True)
    show(['B1 v.s. B2', r3], True)
    plt.plot(b0, b2, '+b', label='Raw Contact Map (r=%.3f)'%r2)
    plt.plot(b0, b1, '.r', label='Factorable Map (r=%.3f)'%r1)
    plt.legend(loc='upper left')
    plt.xlabel('Bias vector from ICE')
    plt.ylabel('Bias vectors from BNMF')
    pdf.savefig(); plt.clf();

def plot_maps_2(pdf, dim=100):
    path = '../data/Imakaev2012NM'
    genome = '../data/hg19_chr_len.txt'
    map1 = ContactMap()
    map1.genome_info(genome)
    map1.create_densemap([path+'/SRR027956_map-res1000k.hdf5'], reso=1e6)
    map1.mask_diag()
    map1.mask_short(0.5)
    map1.mask_low(1)
    map1.trim_high(0.05)
    map2 = map1.duplicate()

    map1.decompose(method='Correct')
    map2.decompose(method='NND', dim_num=dim, A=map2.get_null_map())
    map2.decompose(method='NMF-PoissonManifoldEqual', dim_num=dim, max_iter=10)

    ## Compare results
    b1 = np.array(map1.bias_vector)
    b2 = np.array(map2.bias_vector)

    cor = correlation(b1.tolist(), b2.tolist())
    show(['B1 v.s. B2', cor], True)
    plt.plot(b1, b2, '.b', label='PCC=%.3f'%cor)
    plt.legend(loc='upper left')
    plt.xlabel('Position bias vector from ICE')
    plt.ylabel('Position bias vector from CNMF')
    pdf.savefig(); plt.clf();

def bias_vector(name='', icefile=''):
    show('''
    Check the bias from ICE and the bias from NMF
    ''', True)
    map1 = ContactMap(name)
    assert map1.load()
    map1.decompose_auto(dim_num=map1.contact_group.shape[1])
    map1.sort_groups()
    idx, val = map1.get_locations('../data/Imakaev2012NM/'+icefile+'-bias.txt', st=0, ch=0, po=1, nm=3, add=0)
    bias = np.zeros(map1.contact_map.shape[0])
    bias[idx] = np.array([float(v) for v in val])
    show('Pearson coeff')
    show(correlation(map1.bias_vector, bias), True)
    show('Spearman rank coeff')
    show(correlation(map1.bias_vector, bias, rank=True), True)

    idx, val = map1.get_locations('../data/Imakaev2012NM/'+icefile+'-eig.txt', st=0, ch=0, po=1, nm=3, add=0)
    eig1 = np.zeros(map1.contact_map.shape[0])
    eig1[idx] = np.array([float(v) for v in val])
    show('Pearson coeff')
    show(correlation(np.asarray(map1.contact_group)[:,0], eig1), True)
    show('Spearman rank coeff')
    show(correlation(np.asarray(map1.contact_group)[:,0], eig1, rank=True), True)

def plot_replicates():
    show('''
    Check the trends of membership values in H*S
    ''', True)
    genome = '../data/hg19_chr_len.txt'
    files = [['SRR071233'             ], ## Hi-C Kalhor 1
             [             'SRR071234'], ## Hi-C Kalhor 2
             ['SRR071231', 'SRR071232'], ## TCC Kalhor
             ['SRR027962', 'SRR027963'], ## K562
            ]
    for r in [8e6, 4e6, 2e6, 1e6]:
        maps = []
        for f in files:
            map1 = ContactMap()
            map1.genome_info(genome)
            map1.create_densemap(['../work/%s_map-res1000k.hdf5'%d for d in f], reso=1e6)
            map1.get_interactions()
            map1.create_binnedmap(r)
            map1.mask_diag()
            maps.append(map1)
        for i in xrange(50,100,50):
            show(r)
            show(i)
            if True:
                map1 = maps[0]
                map1.decompose('NND', dim_num=i)
                map1.decompose('NMF-PoissonManifoldEqual', dim_num=i)
            for map2 in maps[1:]:
                map2.decompose('NND', dim_num=i)
                map2.decompose('NMF-PoissonManifoldEqual', dim_num=i)
                show(map1.compare(map2, raw=True, metric='SPC'))
                show(map1.compare(map2, raw=False, metric='SPC'))
            show()

def main(para):
    pdf = PdfPages(para['ExeFile']+'plot.pdf')
#    plot_maps_1(para, 'SRR027956', pdf)
#    plot_maps_2(pdf)
    plot_replicates()
    pdf.close()


if __name__ == '__main__': main_fun(main)
