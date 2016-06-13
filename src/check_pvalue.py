from tools import *
from contact_map import *
from numpy.random import choice
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def run1(para):
    show('''
    Check the distribution of P-values for randomly sampled bins
    ''', True)
    pdf = PdfPages(para['ExeFile']+'plot.pdf')
    map1 = ContactMap()
    assert map1.load('demo')
    map1.decompose_auto()
    show('Now we check')
    show(map1.name, True)
    H = np.array(map1.contact_group)
    n,r = H.shape
#    mlist1 = ['AvgRnd', 'AvgUni', 'AvgCCD', 'AUC', 'GSEA']
#    mlist2 = ['Tau', 'PCC', 'SPC']
    mlist1 = ['AUC', 'AvgCCD', 'GSEA']
    mlist2 = ['SPC']
    mlist = mlist1 + mlist2
    N = 500 ## how many times to repeat
    T = n/r  ## how many bins in the query set
    print n, r, N, T

    sample = np.random.random(n)
    ridx = (sample < T/float(n))
    ridx1 = np.random.choice(n, T, replace=True)

    pvals = {}
    for m in mlist:
        pvals[m] = []
        if m in mlist1:
            pvals['#Bin-'+m] = []
    for i in xrange(N):
        print i+1, '/', N
        map1.randomize_map()
        map1.decompose('NND', dim_num=r)
        map1.decompose('NMF-PoissonManifoldEqual', dim_num=r)
        for m in mlist1:
            srt, val, pval = map1.test_enrichment(ridx, m, correct=False, normp=False)
            srt1, val1, pval1 = map1.test_enrichment(ridx1, m, correct=False, normp=False)
            pvals[m] += pval.tolist()
            pvals['#Bin-'+m] += pval1.tolist()
        for m in mlist2:
            srt, val, pval = map1.test_enrichment(sample, m, correct=False, normp=False)
            pvals[m] += pval.tolist()
    bins = [i/10.0 for i in xrange(10)]
    show('')
    show(bins, True)
    for m in sorted(pvals.keys()):
        show(m)
        show(histogram(pvals[m], bins), True)
    fig = plt.figure()
    for i in xrange(len(mlist)):
        axis = fig.add_subplot(221+i)
        axis.hist(pvals[mlist[i]], np.linspace(0,1,10), normed=1, color='white')
        axis.set_xlim([0, 1])
        axis.set_title('P-value distribution of '+mlist[i].replace('AvgCCD','CCD'))
        axis.tick_params(labelsize=10)
    fig.tight_layout()
    fig.savefig(pdf, format='pdf')
    pdf.close()

def run2():
    show('''
    Check cluster enrichments of bins enriched in each cluster
    ''', True)
    map1 = ContactMap()
    assert map1.load('YeastHiC')
    show('Now we check')
    show(map1.name, True)
    gfile = map1.output_groups()
    n,r = map1.contact_group.shape
    mlist = ['AvgRnd', 'AvgUni', 'AvgCCD', 'AUC']
    mlist2 = ['Tau', 'PCC', 'SPC']
    show('')
    show(range(r), True)
    for m in mlist+mlist2:
        ci = []
        show(m)
        for i in xrange(r):
            idx, mrk = map1.get_locations(gfile, st=0, ch=0, po=1, nm=i+3)
            avg, std = mean_std([float(j) for j in mrk])
            ridx = [i for i,j in zip(idx,mrk) if float(j) > 1]
            if m in mlist2:
                srt, val, pval = map1.test_enrichment([float(j) for j in mrk], m)
            else:
                srt, val, pval = map1.test_enrichment(ridx, m)
            ci.append(srt[0])
        show(ci)
        show()

def run3():
    show('''
    Check the cluster relationship among two maps
    ''', True)
    map1 = ContactMap()
    assert map1.load('YeastHiC')
    map2 = ContactMap()
    assert map2.load('YeastSyn')
    show('Now we compare')
    show([map1.name, 'and', map2.name], True)

    ab1 = set(); ab2 = set();

    n,r1 = map1.contact_group.shape
    for i in xrange(r1):
        srt, val, pval = map2.test_enrichment(map1.contact_group[:,i], 'PCC')
        ab1.add((i,srt[0]))

    n,r2 = map2.contact_group.shape
    for i in xrange(r2):
        srt, val, pval = map1.test_enrichment(map2.contact_group[:,i], 'PCC')
        ab2.add((srt[0],i))
    ab = ab1 & ab2
    a,b = zip(*ab)

    show('Mapped Clusters\n')
    show(map1.name)
    show([i for i,j in ab], True)
    show(map2.name)
    show([j for i,j in ab], True)
    show('Unmapped Clusters\n')
    show(map1.name)
    show([i for i in xrange(r1) if i not in a], True)
    show(map2.name)
    show([i for i in xrange(r2) if i not in b], True)

def main(para):
    run1(para)
#    run2()
#    run3()

if __name__ == '__main__': main_fun(main)

