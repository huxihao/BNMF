from tools import *
from contact_map import *

def example1(para):
    datapath = para['DataPath']
    dataset = {'name': 'cmp-Yeast',
         'genomefile': datapath + '/genome/yeast_chr_len.txt',
         'linkfiles': [datapath + '/Duan2010N/interactions_HindIII_fdr0.01_inter.txt', 
                       datapath + '/Duan2010N/interactions_HindIII_fdr0.01_intra.txt'],
         'resolution': 10e3,
         'trnafile': datapath + '/Duan2010N/trna_positions_and_clusters.txt',
         'refloci': {
             'Centromere':['chrI:151465', 'chrII:238207', 'chrIII:114385', 'chrIV:449711', 
                           'chrV:151987', 'chrVI:148510', 'chrVII:497038', 'chrVIII:105703', 
                           'chrIX:355629', 'chrX:436425', 'chrXI:440246', 'chrXII:150947', 
                           'chrXIII:268031', 'chrXIV:628758', 'chrXV:326702', 'chrXVI:555957']}
        }
    pdf = PdfPages(para['ExeFile']+'plot.pdf')
    map1 = ContactMap(dataset['name'])
#    map1.clear()
    if not map1.load():
        map1.genome_info(dataset['genomefile'])
        map1.create_densemap(dataset['linkfiles'], dataset['resolution'])
        map1.mask_diag()
        map1.mask_short()
        map1.mask_low()
        map1.save()
    map1.decompose_auto()
    map1.sort_groups()
    map1.plot_map()
    pdf.savefig(); plt.clf();
    map1.plot_submap()
    pdf.savefig(); plt.clf();

    H = np.asarray(map1.contact_group)
    n,r = H.shape
    target = {}
    if 'refloci' in dataset:
        marks = dataset['refloci']
        for m in marks:
            pos = map1.choose_map_loc(marks[m])
            for p in pos:
                target[p] = m
    index0, name0 = map1.get_locations(dataset['trnafile'], st=1, ch=2, po=3, nm=5)
    show([len(index0), len(set(index0))], True)
    index1 = [i for i,j in zip(index0, name0) if j == 'cluster_1_dim_nucleolus']
    index2 = [i for i,j in zip(index0, name0) if j == 'cluster_2_bright_centromeres']
    index3 = [i for i,j in zip(index0, name0) if j == 'cluster_3_other']
    show([len(index0), len(set(index0))], True)
    show([len(set(index0)), len(set(index1)), len(set(index0) & set(index1))], True)
    show([len(set(index0)), len(set(index2)), len(set(index0) & set(index2))], True)
    show([len(set(index1)), len(set(index2)), len(set(index1) & set(index2))], True)
    sel = [i for i in index0 if H[i,:].sum()>0]
    sel1 = [i for i in index1 if i in sel]
    sel2 = [i for i in index2 if i in sel]
    sel3 = [i for i in index3 if i in sel]
    idx1, val1, pval1 = map1.test_enrichment(sel, method='AUC', correct=True, plot=pdf, pages=4)
    idx2, val2, pval2 = map1.test_enrichment(sel, method='AvgCCD', correct=True, plot=pdf, pages=4)
    for k in idx1:
        if pval1[k] > 0.05:
            break
        show(k+1)
        show(val1[k])
        show(pval1[k])
        show(val2[k])
        show(pval2[k])
        if pval2[k] < 0.05: # output file
            index0, val0 = map1.get_locations(dataset['trnafile'], st=1, ch=2, po=3, nm=2)
            index0, val1 = map1.get_locations(dataset['trnafile'], st=1, ch=2, po=3, nm=3)
            index0, val2 = map1.get_locations(dataset['trnafile'], st=1, ch=2, po=3, nm=4)
            index0, val3 = map1.get_locations(dataset['trnafile'], st=1, ch=2, po=3, nm=0)
            outfile = open(datapath+'/genome/yeast_trna-C%s.txt'%(k+1),'w')
            outfile.write('Chromosome\tStart\tEnd\tName\tBin\tValue\n')
            cluster = set()
            for i in xrange(len(index0)):
                if H[index0[i],k] >= 1:
                    outfile.write('%s\t%s\t%s\t%s\t%s\t%.3f\n'%(val0[i], val1[i], val2[i], val3[i], index0[i], H[index0[i],k]))
                    cluster.add(index0[i])
            show('\nOverlaps\n')
            show(len(cluster), True)
            show(len(cluster & set(index1)))
            show(len(set(index1)), True)
            show(len(cluster & set(index2)))
            show(len(set(index2)), True)
            show(len(cluster & set(index3)))
            show(len(set(index3)), True)
            show()
            outfile.close()
        bin_idx = np.argsort(H[:,k])[::-1]
        show([map1.get_bin_str(j)+'(%s %s)'%(H[j,k], target[j]) for j in bin_idx if H[j,k]>=1 and j in target], True)
        fig = plt.figure()
        axis = fig.add_subplot('311')
        plt.title('Enrichment of tRNA clusters on C%s'%(k+1))
        axis.plot(H[:,k], 'b')
        axis.plot(sel1, H[sel1,k], 'r.', label='Cluster 1')
        axis.plot([0,n], [cal_mean(H[sel1,k])]*2, 'k--')
        axis.set_xlim([0,n])
        axis.legend()
#        axis.set_yscale('log')
        axis = fig.add_subplot('312')
        axis.plot(H[:,k], 'b')
        axis.plot(sel2, H[sel2,k], 'r.', label='Cluster 2')
        axis.plot([0,n], [cal_mean(H[sel2,k])]*2, 'k--')
        axis.set_xlim([0,n])
        axis.legend()
#        axis.set_yscale('log')
        axis = fig.add_subplot('313')
        axis.plot(H[:,k], 'b')
        axis.plot(sel3, H[sel3,k], 'r.', label='Cluster 3')
        axis.plot([0,n], [cal_mean(H[sel3,k])]*2, 'k--')
        axis.set_xlim([0,n])
        axis.legend()
#        axis.set_yscale('log')
        fig.savefig(pdf, format='pdf')
    show('\n\n')
    ################################################
    ## Show distance distribution in synthetic structures
    map2 = map1.duplicate()
    map2.create_binnedmap(3200, lazy=True)
    index0, name0 = map2.get_locations(dataset['trnafile'], st=1, ch=2, po=3, nm=5)
    show([len(index0), len(set(index0))], True)
    index1 = [i for i,j in zip(index0, name0) if j == 'cluster_1_dim_nucleolus']
    index2 = [i for i,j in zip(index0, name0) if j == 'cluster_2_bright_centromeres']
    index3 = [i for i,j in zip(index0, name0) if j == 'cluster_3_other']
    idx_new, name_new = map2.get_locations(datapath+'/genome/yeast_trna-C77.txt', st=1)

    idx_back, name_back = map2.get_locations(map1.output_groups(), st=0, ch=0, po=1, nm=79)
    background = [i for i,j in zip(idx_back, name_back) if float(j) >= 1]
    cen_idx = map2.choose_map_loc(dataset['refloci']['Centromere'])
    idx_fake = []
    for i in idx_new:
        one_dim_dist = []
        for p in cen_idx:
            one_dim_dist.append((abs(i-p), p))
        one_dim_dist.sort()
        A,B = one_dim_dist[0]
        if A == 0: ## no change
            j = i
        else:
            j = 2*B-i
            if j >= len(map2.frag_chr) or map2.frag_chr[i] != map2.frag_chr[j]:
                j = i ## set back
            else:
                print j, '<-', i
        idx_fake.append(j)

    dist = np.load('syn_dist.npy')
    fig = plt.figure()
    groups = [('All tRNAs', index0),
              #('Duan et al - Nucleolus', index1),
              ('tRNAs in Duan-cluster 2', index2),
              ('tRNAs in BNMF-C77', idx_new),
              ('All Bins in BNMF-C77', background)]
    idx_dist = []
    for i in xrange(len(groups)):
        axis = fig.add_subplot(len(groups)*100+11+i)
        name, idx = groups[i]
        D1 = dist[idx,:][:,idx]
        d1 = D1[np.triu_indices(D1.shape[0], k=1)]
        idx_dist.append(d1)
        bins = range(100, 1300, 50)
        axis.hist(d1, bins)
        axis.set_xlim([min(bins), max(bins)])
        axis.set_title('%s #=%d (Avg.=%.1fnm)'%(name, len(idx), d1.mean()))
        axis.tick_params(labelsize=8)
#        show(name)
#        show(histogram(bins, d1.tolist(), prob=False), True)

    fig.tight_layout()
    fig.subplots_adjust(left=0.2, right=0.8)
    fig.savefig(pdf, format='pdf')
    from scipy.stats import ttest_ind, ttest_rel
    show('Test1')
    show(ttest_ind(idx_dist[1], idx_dist[2]), True)
    show('Test2')
    show(ttest_ind(idx_dist[2], idx_dist[3]), True)
    pdf.close()

def cal_mean(a):
    s = 0.0
    n = 0.0
    for i in a:
        if i > 0:
            s += i
            n += 1
    m = s/n
    return m

def run1(para):
    datapath = para['DataPath']
    datasets1 = [ ## for tRNAs
        {'name': 'cmp-Yeast',
         'genomefile': datapath + '/genome/yeast_chr_len.txt',
         'linkfiles': [datapath + '/Duan2010N/interactions_HindIII_fdr0.01_inter.txt', 
                       datapath + '/Duan2010N/interactions_HindIII_fdr0.01_intra.txt'],
         'resolution': 10e3,
         'testfiles': ['yeast_centromers', 'yeast_telomeres', 'yeast_trna', 'yeast_early', 'yeast_late']
         }, 
        {'name': 'cmp-Fly',
         'genomefile': datapath + '/genome/dm3_chr_len.txt',
         'linkfiles': [datapath + '/Sexton2012C/wtb_s0.mat.links1.gz'],
         'resolution': 1e5,
         'testfiles': ['dm3_trna', 'dm3_snoRNA'],
         },
        {'name': 'cmp-Mouse',
         'genomefile': datapath + '/genome/mm9_chr_len.txt',
         'linkfiles': [datapath + '/Zhang2012C/GSM870040_ATM-67-R1.validPair.txt.gz.links.gz'], 
         'resolution': 1e6,
         'testfiles': ['mm9_trna'],
         },
        {'name': 'cmp-Human',
         'genomefile': datapath + '/genome/hg19_chr_len.txt',
         'linkfiles': [datapath + '/Imakaev2012NM/SRR027956_map-res1000k.hdf5'],
         'resolution': 1e6,
         'testfiles': ['hg19_trna'],
         }]
    datasets2 = [
        {'name': 'cmp-Parasite0',
         'genomefile': datapath + '/genome/plasmo_chr_len.txt',
         'linkfiles': [datapath + '/Ay2014GR/GSM1215592_rawContactCounts-10kb-RINGS.txt.gz'],
         'resolution': 20e3,
         'testfiles': ['plasmo_trna']
         },
        {'name': 'cmp-Parasite18',
         'genomefile': datapath + '/genome/plasmo_chr_len.txt',
         'linkfiles': [datapath + '/Ay2014GR/GSM1215593_rawContactCounts-10kb-TROPHOZOITES-XL.txt.gz'],
         'resolution': 20e3,
         'testfiles': ['plasmo_trna']
         },
        {'name': 'cmp-Parasite36',
         'genomefile': datapath + '/genome/plasmo_chr_len.txt',
         'linkfiles': [datapath + '/Ay2014GR/GSM1215594_rawContactCounts-10kb-SCHIZONTS.txt.gz'],
         'resolution': 20e3,
         'testfiles': ['plasmo_trna']
         }]
    pdf = PdfPages(para['ExeFile']+'plot.pdf')
    for dataset in datasets1:
        show(dataset['name'])
        map1 = ContactMap(dataset['name'])
#        map1.clear()
        if not map1.load():
            map1.genome_info(dataset['genomefile'])
            map1.create_densemap(dataset['linkfiles'], dataset['resolution'])
            map1.mask_diag()
            map1.mask_short()
            map1.mask_low()
            map1.save()
        map1.decompose_auto(update=False)
        map1.sort_groups()
        map1.output_groups()
        map1.plot_map()
        pdf.savefig(); plt.clf();
        map1.plot_submap()
        pdf.savefig(); plt.clf();

        H = np.asarray(map1.contact_group)
        n,r = H.shape
        show(n)
        show(r, True)
        for test in dataset['testfiles']:
            show(test)
            index, name = map1.get_locations(datapath+'/genome/%s.txt'%test, st=1, ch=0, po=1)
            sel = [i for i in index if H[i,:].sum()>0]
            idx1, val1, pval1 = map1.test_enrichment(sel, method='AUC', correct=True, plot=pdf, pages=4)
            idx2, val2, pval2 = map1.test_enrichment(sel, method='AvgCCD', correct=True, plot=pdf, pages=3)
            idx3, val3, pval3 = map1.test_enrichment(sel, method='GSEA', correct=True, plot=pdf, pages=3)
            show('P-values')
            show(idx1[0])
            show(val1[idx1[0]])
            show(pval1[idx1[0]])
            show(idx2[0])
            show(val2[idx2[0]])
            show(pval2[idx2[0]])
            show(idx3[0])
            show(val3[idx3[0]])
            show(pval3[idx3[0]])
            show()
    pdf.close()

def main(para):
#    example1(para)
    run1(para)

if __name__ == "__main__": main_fun(main)
