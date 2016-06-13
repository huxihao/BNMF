from tools import *
from contact_map import *
import numpy as np
import scipy.spatial.distance as dist
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def pdb2map(pdbfile):
    pdb = open(pdbfile, 'r')
    pp = []
    for line in pdb:
        if line.startswith('ATOM'):
            pp.append([float(a) for a in line[30:55].split()])
    pdb.close()
    pos = np.array(pp)
    return dist.squareform(dist.pdist(pos, 'euclidean'))

def get_contact_map(path, cutoff=60.0, maxnum=-1):
    import os
    cc = 0
    cm = 0
    allfiles = os.listdir(path)
    for f in allfiles:
        cc += 1
        print cc, f
        m = pdb2map(path+'/'+f)
        if cutoff is None:
            if cc == 1:
                cm = m
            else:
                idx = (cm > m)
                cm[idx] = 0
                m[idx==False] = 0
                cm += m
        elif cutoff <= 0:
            cm += m
        else:
            cm += np.array(m<=cutoff, dtype='int')
        if cc == maxnum:
            break
    print cm.shape
    if cutoff is None:
        return cm
    elif cutoff <= 0:
        return cm/float(cc)
    else:
        return cm

def hic_vs_syn(para, links, reso=32000):
    show('''
    Compare HiC map and syntheic map
    ''', True)
    pdf = PdfPages(para['ExeFile']+'plot.pdf')

    map1 = ContactMap('YeastHiC')
    map1.genome_info(para['DataPath']+'/Tjong2012GR/yeast_chr_len-Tjong.txt')
    datafiles = [para['DataPath']+'/Duan2010N/interactions_HindIII_fdr0.01_inter.txt',
                 para['DataPath']+'/Duan2010N/interactions_HindIII_fdr0.01_intra.txt'] 
    for datafile in datafiles:
        map1.add_interactions(datafile)
    map1.create_binnedmap(3200)
    map2 = map1.duplicate('YeastSyn')
    refm = map1.duplicate('YeastRef')
    map2.contact_map = np.matrix(links, dtype='float')
    map2.get_interactions()

    map1.remove_small(1)
    map1.create_binnedmap(reso, throw=20e3)
    map1.mask_diag()
    map1.mask_short()
    map2.remove_small(1)
    map2.create_binnedmap(reso, throw=20e3)
    map2.mask_diag()
    map2.mask_short()

    map1.plot_map(title='Hi-C contact map')
    pdf.savefig(); plt.clf();
    map2.plot_map(title='Synthetic contact map')
    pdf.savefig(); plt.clf();

    ## Show Correlation
    n = map1.contact_map.shape[0]
    corr = np.zeros(n)
    for i in xrange(n):
        x = np.asarray(np.nan_to_num(map1.contact_map))[i,:]
        y = np.asarray(np.nan_to_num(map2.contact_map))[i,:]
        corr[i] = correlation(x,y)
    plt.hist(corr[np.logical_not(np.isnan(corr))], bins=30, color='white')
    plt.xlabel('Correlations')
    plt.ylabel('Frequency')
    plt.title("Histogram of row-based Pearson's coefficients")
    pdf.savefig(); plt.clf();

#    map1.clear()
#    map2.clear()

    ## Decompose them
    map1.decompose_auto(plot=None)
    map1.sort_groups()
    map1.plot_submap()
    pdf.savefig(); plt.clf();

    map2.decompose_auto(plot=None)
    map2.sort_groups()
    map2.plot_submap()
    pdf.savefig(); plt.clf();

    show(map1.contact_group.shape, True)
    show(map2.contact_group.shape, True)
    map1.plot_map(map1.contact_group, log=False, title='Clusters from Hi-C contact map')
    pdf.savefig(); plt.clf();
    map2.plot_map(map2.contact_group, log=False, title='Clusters from Synthetic map')
    pdf.savefig(); plt.clf();

    R1 = map1.contact_group * map1.group_map * map1.contact_group.T
    R2 = map2.contact_group * map2.group_map * map2.contact_group.T
    x = np.asarray(R1).reshape(-1).tolist()
    y = np.asarray(R2).reshape(-1).tolist()
    corr = correlation(x,y)
    p = np.poly1d(np.polyfit(x,y,1))
#    plt.plot(x[::20], y[::20], 'b.')
#    plt.plot(x[::20], p(x[::20]), 'r-')
#    plt.xlabel('Recovered read counts in %s'%map1.name)
#    plt.ylabel('Recovered read counts in %s'%map2.name)
#    plt.title('Overall correlation is %.3f'%corr)
#    pdf.savefig(); plt.clf();

    plt.imshow(map_cor(map1.contact_group, map2.contact_group),
               interpolation='none', aspect='auto', cmap='OrRd')
    plt.colorbar()
    plt.ylabel('Clusters in %s'%map1.name)
    plt.xlabel('Clusters in %s'%map2.name)
    c1 = range(0,map1.contact_group.shape[1],5)
    c2 = range(0,map2.contact_group.shape[1],5)
    plt.yticks(c1, ['C'+str(i+1) for i in c1])
    plt.xticks(c2, ['C'+str(i+1) for i in c2])
    plt.title('Pearson corr. coeff. of cluster memberships')
    pdf.savefig(); plt.clf();
    map1.save()
    map2.save()
    refm.save()
    test_enrichment(para, map1, map2, refm, pdf)
    pdf.close()

def test_enrichment(para, map1, map2, refm, pdf):
    Colo = 'AvgCCD' ## test method
    outfile = open('loci_idx.txt', 'w')

    ## check centromer
    plt.title('Check centromers')
    idx, names = map1.get_locations(para['DataPath']+'/Tjong2012GR/yeast_centromers-Tjong.txt', st=1, ch=0, po=1, nm=2)
    idx0, names = refm.get_locations(para['DataPath']+'/Tjong2012GR/yeast_centromers-Tjong.txt', st=1, ch=0, po=1, nm=2)

    srt1, val1, pval1 = map1.test_enrichment(idx, method=Colo, title='YeastHiC Centromers', plot=None)
    srt2, val2, pval2 = map2.test_enrichment(idx, method=Colo, title='YeastSyn Centromers', plot=None)

    plt.plot(np.asarray(map1.contact_group[:,srt1[0]]).reshape(-1),
             np.asarray(map2.contact_group[:,srt2[0]]).reshape(-1), 'b.')
    plt.xlabel('C%s from %s'%(srt1[0], map1.name))
    plt.ylabel('C%s from %s'%(srt2[0], map2.name))
    pdf.savefig(); plt.clf();

    outfile.write('Centromer-All\t')
    for j in xrange(len(idx)):
        outfile.write(str(idx0[j])+'\t')
    outfile.write('\n')
    for i in srt2: ## clusters in syn map
        if pval2[i] > 0.01:
            continue
        outfile.write('Centromer-C%s\t'%i)
        for j in xrange(len(idx)):
            if map2.contact_group[idx[j],i] > 1:
                outfile.write(str(idx0[j])+'\t')
        outfile.write('\n')

    ## check telomerers
    plt.title('Check telomeres')
    idx, names = map1.get_locations(para['DataPath']+'/Tjong2012GR/yeast_telomeres-Tjong.txt', st=1, ch=0, po=1, nm=2)
    idx0, names = refm.get_locations(para['DataPath']+'/Tjong2012GR/yeast_telomeres-Tjong.txt', st=1, ch=0, po=1, nm=2)

    srt1, val1, pval1 = map1.test_enrichment(idx, method=Colo, title='YeastHiC Telomeres', plot=None)
    srt2, val2, pval2 = map2.test_enrichment(idx, method=Colo, title='YeastSyn Telomeres', plot=None)

    plt.plot(np.asarray(map1.contact_group[:,srt1[0]]).reshape(-1),
             np.asarray(map2.contact_group[:,srt2[0]]).reshape(-1), 'b.')
    plt.xlabel('C%s from %s'%(srt1[0], map1.name))
    plt.ylabel('C%s from %s'%(srt2[0], map2.name))
    pdf.savefig(); plt.clf();

    outfile.write('Telomere-All\t')
    for j in xrange(len(idx)):
        outfile.write(str(idx0[j])+'\t')
    outfile.write('\n')
    for i in srt2: ## clusters in syn map
        if pval2[i] > 0.01:
            continue
        outfile.write('Telomere-C%s\t'%i)
        for j in xrange(len(idx)):
            if map2.contact_group[idx[j],i] > 1:
                outfile.write(str(idx0[j])+'\t')
        outfile.write('\n')

    ## check tRNA
    plt.title('Check tRNAs')
    idx, names = map1.get_locations(para['DataPath']+'/Duan2010N/trna_positions_and_clusters.txt', st=1, ch=2, po=3, nm=5)
    idx0, names = refm.get_locations(para['DataPath']+'/Duan2010N/trna_positions_and_clusters.txt', st=1, ch=2, po=3, nm=5)

    srt1, val1, pval1 = map1.test_enrichment(idx, method=Colo, title='YeastHiC tRNAs', plot=None)
    srt2, val2, pval2 = map2.test_enrichment(idx, method=Colo, title='YeastSyn tRNAs', plot=None)

    plt.plot(np.asarray(map1.contact_group[:,srt1[0]]).reshape(-1),
             np.asarray(map2.contact_group[:,srt2[0]]).reshape(-1), 'b.')
    plt.xlabel('C%s from %s'%(srt1[0], map1.name))
    plt.ylabel('C%s from %s'%(srt2[0], map2.name))
    pdf.savefig(); plt.clf();

    outfile.write('tRNA-All\t')
    for j in xrange(len(idx)):
        outfile.write(str(idx0[j])+'\t')
    outfile.write('\n')

    outfile.write('tRNA-DimNucleolus\t')
    for j in xrange(len(idx)):
        if names[j] == 'cluster_1_dim_nucleolus':
            outfile.write(str(idx0[j])+'\t')
    outfile.write('\n')
    outfile.write('tRNA-Centromeres\t')
    for j in xrange(len(idx)):
        if names[j] == 'cluster_2_bright_centromeres':
            outfile.write(str(idx0[j])+'\t')
    outfile.write('\n')
    outfile.write('tRNA-Other\t')
    for j in xrange(len(idx)):
        if names[j] == 'cluster_3_other':
            outfile.write(str(idx0[j])+'\t')
    outfile.write('\n')

    idx_pre, n_pre = refm.get_locations(para['DataPath']+'/genome/yeast_trna-C77.txt', st=1, ch=0, po=1, nm=3)
    outfile.write('tRNA-BNMF-Centromeres\t')
    for j in xrange(len(idx_pre)):
        outfile.write(str(idx_pre[j])+'\t')
    outfile.write('\n')

    for i in srt1: ## clusters in hi-c map
        if pval2[i] > 0.05:
            continue
        outfile.write('tRNA-check-C%s\t'%(i+1))
        for j in xrange(len(idx)):
            if map1.contact_group[idx[j],i] >= 1:
                outfile.write(str(idx0[j])+'\t')
        outfile.write('\n')

    ## check early replicating time
    plt.title('Check early replicate')
    idx, names = map1.get_locations(para['DataPath']+'/Duan2010N/origins_unchecked_early.txt', st=0, ch=0, po=1, nm=0)
    idx0, names = refm.get_locations(para['DataPath']+'/Duan2010N/origins_unchecked_early.txt', st=0, ch=0, po=1, nm=0)

    srt1, val1, pval1 = map1.test_enrichment(idx, method=Colo, title='YeastHiC Early Rep.', plot=None)
    srt2, val2, pval2 = map2.test_enrichment(idx, method=Colo, title='YeastSyn Early Rep.', plot=None)

    plt.plot(np.asarray(map1.contact_group[:,srt1[0]]).reshape(-1),
             np.asarray(map2.contact_group[:,srt2[0]]).reshape(-1), 'b.')
    plt.xlabel('C%s from %s'%(srt1[0], map1.name))
    plt.ylabel('C%s from %s'%(srt2[0], map2.name))
    pdf.savefig(); plt.clf();

    outfile.write('EarlyRep-All\t')
    for j in xrange(len(idx)):
        outfile.write(str(idx0[j])+'\t')
    outfile.write('\n')
    for i in srt2: ## clusters in syn map
        if pval2[i] > 0.01:
            continue
        outfile.write('EarlyRep-C%s\t'%i)
        for j in xrange(len(idx)):
            if map2.contact_group[idx[j],i] > 1:
                outfile.write(str(idx0[j])+'\t')
        outfile.write('\n')

    ## check replicating time
    plt.title('Check late replicate')
    idx, names = map1.get_locations(para['DataPath']+'/Duan2010N/origins_checked_late.txt', st=0, ch=0, po=1, nm=0)
    idx0, names = refm.get_locations(para['DataPath']+'/Duan2010N/origins_checked_late.txt', st=0, ch=0, po=1, nm=0)

    srt1, val1, pval1 = map1.test_enrichment(idx, method=Colo, title='YeastHiC Late Rep.', plot=None)
    srt2, val2, pval2 = map2.test_enrichment(idx, method=Colo, title='YeastSyn Late Rep.', plot=None)

    plt.plot(np.asarray(map1.contact_group[:,srt1[0]]).reshape(-1),
             np.asarray(map2.contact_group[:,srt2[0]]).reshape(-1), 'b.')
    plt.xlabel('C%s from %s'%(srt1[0], map1.name))
    plt.ylabel('C%s from %s'%(srt2[0], map2.name))
    pdf.savefig(); plt.clf();

    outfile.write('LateRep-All\t')
    for j in xrange(len(idx)):
        outfile.write(str(idx0[j])+'\t')
    outfile.write('\n')
    for i in srt2: ## clusters in syn map
        if pval2[i] > 0.01:
            continue
        outfile.write('LateRep-C%s\t'%i)
        for j in xrange(len(idx)):
            if map2.contact_group[idx[j],i] > 1:
                outfile.write(str(idx0[j])+'\t')
        outfile.write('\n')

    ## check all clusters in syn map
    file0 = map2.output_groups() ## syn
    idx, names = map2.get_locations(file0, st=0, ch=0, po=1, nm=0)
    idx0, names = refm.get_locations(file0, st=0, ch=0, po=1, nm=0)
    for i in xrange(map2.contact_group.shape[1]):
        sub_map = map2.contact_group[:,i] * map2.group_map[i,i] * map2.contact_group[:,i].T
        #map2.plot_map(sub_map, title='Cluster-C%s'%i, log=False)
        #pdf.savefig(); plt.clf();
        outfile.write('Cluster-C%s\t'%i)
        for j in xrange(len(idx)):
            if map2.contact_group[idx[j],i] > 1:
                outfile.write(str(idx0[j])+'\t')
        outfile.write('\n')

def check_dist(para, m):
    show('''
    Check the distributions of average distances of clusters
    ''', True)
    bins = range(0, 2001, 50)
    show(['', 'Mean', 'STD'])
    show(bins, True)
    show('Background')
    show(mean_std(m.reshape(-1).tolist()))
    show(histogram(m.reshape(-1).tolist(), bins, prob=False), True)
    with open('loci_idx.txt', 'r') as tempfile: 
        for line in tempfile:
            ele = line.strip().split('\t')
            idx = np.array([int(i) for i in ele[1:]])
            if len(idx) < 2:
                continue
            show(ele[0], False)
            val = [m[i,j] for i in idx for j in idx if i<j]
            show(mean_std(val))
            show(histogram(val, bins, prob=False), True)

def main(para):
    if not os.path.exists('syn_link.npy'):
        np.save('syn_link.npy', np.array(get_contact_map(para['DataPath']+'/Tjong2012GR/simulation', cutoff=60)))
        np.save('syn_dist.npy', np.array(get_contact_map(para['DataPath']+'/Tjong2012GR/simulation', cutoff=-1)))
        np.save('syn_mind.npy', np.array(get_contact_map(para['DataPath']+'/Tjong2012GR/simulation', cutoff=None)))
    hic_vs_syn(para, np.load('syn_link.npy'))
    check_dist(para, np.load('syn_dist.npy'))

if __name__ == '__main__': main_fun(main)
