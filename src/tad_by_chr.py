from tools import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from contact_map import ContactMap

plot_left = -5
plot_right = 6
resolution = 200
domain = 'total.combined.domain.txt'

def one_chr(path, cell, genome, ch, pdf=None):
    map1 = ContactMap('tad-%s-in-%s'%(ch,cell))
#    map1.clear()
    if not map1.load():
        map1.genome_info(path+'/%s_chr_len.txt'%genome)
        map1.focus_chromosome(ch)
        map1.create_binnedmap(40e3, lazy=True)
        map1.contact_map = np.loadtxt(path+'/'+cell+'/uij.'+ch)
        print cell, ch, map1.frag_sta.shape[0], map1.contact_map.shape[0]
        assert map1.frag_sta.shape[0] == map1.contact_map.shape[0]
        map1.get_sparse_interactions()
        map1.create_binnedmap(resolution*1000)
        map1.mask_diag()
        map1.mask_short()
        map1.mask_low()
        map1.save()
    show(cell)
    show(ch)
    if pdf is not None:
        map1.plot_map()
        pdf.savefig(); plt.clf()
    map1.decompose_auto(update=False)
    map1.sort_groups()
    show(map1.contact_group.shape)
    if pdf is not None:
        map1.plot_submap()
        pdf.savefig(); plt.clf()
    TAD_st, _ = map1.get_locations(path+'/'+cell+'/'+domain, st=0, ch=0, po=1, add=0)
    TAD_ed, _ = map1.get_locations(path+'/'+cell+'/'+domain, st=0, ch=0, po=2, add=-1)
    TAD = zip(TAD_st, TAD_ed)

    W = np.asarray(map1.contact_group * map1.group_map)
    wm = W.sum(1)
    W /= np.mean(wm[wm>0])
    group = np.argmax(W,1)+1
    group[wm==0] = -1 ## masked regions

    gini = 1-np.power(W,2).sum(1)
    gini[wm==0] = -1 ## masked regions
    log2W = np.log2(W)
    log2W[W==0] = 0
    entropy = (-W*log2W).sum(1)
    entropy[wm==0] = 0

    score = gini
    score[score<0] = 0 ## for ploting

    for i in [1,np.argmax(entropy)/50]:
        sel = np.arange(i*50, min(W.shape[0],(i+1)*50))
        pos = ['%.fM'%(j*resolution*1e-3) for j in sel]
        if pdf is not None:
            fig = plt.figure()
            axis = fig.add_subplot(211)
#            axis.plot(sel, score[sel], '--k')
        for i in xrange(W.shape[1]):
#            if W[sel,i].max() > 0.1:
                if pdf is not None:
                    axis.plot(sel, W[sel,i], label='C%s'%i)
        for i,j in TAD:
            if i in sel and j in sel:
                if pdf is not None:
                    axis.plot([i,j], [1,1], 'k-', linewidth=2)
        if pdf is not None:
            plt.ylim([0,1.2])
            plt.xticks(sel[::int(len(sel)/5)], pos[::int(len(sel)/5)])
            axis = fig.add_subplot(212)
            from matplotlib.colors import LogNorm
            axis.imshow(map1.contact_map[sel,:][:,sel], interpolation='none', norm=LogNorm(), aspect='equal', cmap='OrRd')
            axis.legend()
            fig.savefig(pdf, format='pdf')
            plt.clf()
    
    tad = np.zeros_like(gini)
    tadlen = []
    for i,j in TAD:
        for k in xrange(i+1, j-1):
            tad[k] = (i+j+1)/2 ## regions in the domain
        tadlen.append(j-i)

    tadtype = []
    for i in np.unique(tad):
        if i > 0:
            tadtype.append(len(np.unique(group[tad==i])))
    grptype = []
    for i in np.unique(group):
        if i > 0:
            grptype.append(len(np.unique(tad[group==i])))

    show(np.sum(np.logical_and(tad==0,gini>=0))) ## TADs
    for cut in [0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]:
        show(np.sum(gini>=cut)) ## clusters
        show(np.sum(np.logical_and(tad==0,gini>=cut))) ## TADs % clusters

    cci = 0; ccj = 0 ## scores around the TAD boundaries
    ni = 0; nj = 0
    for i,j in TAD:
        if i+plot_left >=0 and i+plot_right < len(score):
            cci += score[(i+plot_left):(i+plot_right)]
            ni += 1
        if j+plot_left >=0 and j+plot_right < len(score):
            ccj += score[(j+plot_left):(j+plot_right)]
            nj += 1
    show()
    return cci,ccj,ni,nj,tadlen,tadtype,grptype,gini,entropy

def one_cell(path, pdf, cell, genome):
    map1 = ContactMap()
    map1.genome_info(path+'/%s_chr_len.txt'%genome)

    cci = 0; ccj = 0
    ni = 0; nj = 0
    tadlen = []
    tadtype = []
    grptype=[]
    gini = []
    entropy = []
    for i in sorted(map1.idx2chr.keys()):
        info = one_chr(path=path, cell=cell, genome=genome, ch=map1.idx2chr[i], pdf=pdf)
        CCI,CCJ,NI,NJ,TADLEN,TADTYPE,GRPTYPE,GINI,ENTRO = info
        cci += CCI
        ccj += CCJ
        ni += NI
        nj += NJ
        tadlen += TADLEN
        tadtype += TADTYPE
        grptype += GRPTYPE
        gini += GINI[GINI>0].tolist()
        entropy += ENTRO[ENTRO>0].tolist()

    tadlen = np.array(tadlen)
    plt.hist(tadlen*resolution, np.arange(tadlen.max())*resolution)
    plt.title('Distribution of TAD sizes in %s'%cell)
    pdf.savefig(); plt.clf()

    bins = range(max(tadtype)+1)
    count = histogram(tadtype, bins, False)
    show(bins, True)
    show(count, True)
    tadtype = np.array(tadtype)
    plt.hist(tadtype, np.arange(tadtype.max()+1))
    plt.title('Distribution of covered clusters in %s'%cell)
    pdf.savefig(); plt.clf()

    bins = range(max(grptype)+1)
    count = histogram(grptype, bins, False)
    show(bins, True)
    show(count, True)
    grptype = np.array(grptype)
    plt.hist(grptype, np.arange(grptype.max()+1))
    plt.title('Distribution of covered TADs in %s'%cell)
    pdf.savefig(); plt.clf()

#    plt.plot(np.arange(plot_left, plot_right), cci/ni, '.-r', label='TAD start')
#    plt.plot(np.arange(plot_left, plot_right), ccj/nj, '.-b', label='TAD end')
    plt.plot(np.arange(plot_left, plot_right), (cci+ccj)/(ni+nj), '.-k', label='TAD boundary')
    plt.xlabel('Genomic distances (kb)')
    plt.ylabel('Gini impurity score')
    plt.xticks(np.arange(plot_left, plot_right), np.arange(plot_left, plot_right)*resolution)
    plt.xlim([plot_left, plot_right])
    plt.title('Average scores around TAD in %s'%cell)
    plt.legend()
    pdf.savefig(); plt.clf()

    show(mean_std(gini))
    plt.hist(gini, np.arange(0,1.001,0.05))
    plt.title('Distribution of Gini impurity scores')
    plt.xlabel('Gini impurity scores')
    plt.ylabel('Frequency')
    pdf.savefig(); plt.clf()
    show(mean_std(entropy))

    plt.hist(entropy, np.arange(0,6,0.2))
    plt.title('Distribution of entropy at TAD boundaries')
    plt.xlabel('Entropy')
    plt.ylabel('Frequency')
    pdf.savefig(); plt.clf()

def chr_detail(path, cell, genome, ch, loci, st=0, ed=None, pdf=None):
    map1 = ContactMap('tad-detail-%s-in-%s'%(ch,cell))
#    map1.clear()
    if not map1.load():
        map1.genome_info(path+'/%s_chr_len.txt'%genome)
        map1.focus_chromosome(ch)
        map1.create_binnedmap(40e3, lazy=True)
#        map1.contact_map = np.asmatrix(np.loadtxt(path+'/'+cell+'/uij.'+ch))
        map1.contact_map = np.asmatrix(np.loadtxt(path+'/'+cell+'/nij/nij.'+ch))
        print cell, ch, map1.frag_sta.shape[0], map1.contact_map.shape[0]
        assert map1.frag_sta.shape[0] == map1.contact_map.shape[0]
        map1.get_sparse_interactions()
        map1.focus_chromosome(ch, st=st, ed=ed)
        map1.create_contactmap(throw=0)
        map1.save()
    show(cell)
    show(ch)
    show(map1.contact_map.shape)
    map1.mask_diag()
    map1.mask_short()
    map1.decompose_auto(par_lam=1, beta=3, update=False)
    map1.sort_groups()
#    map1.add_bias_back()
    show(map1.contact_group.shape)
    show()
    if pdf is not None:
        map1.plot_map()
        pdf.savefig(); plt.clf()
        map1.plot_map(map1.contact_group*map1.group_map*map1.contact_group.T, vmin=0.01, title='H*S*H.T')
        pdf.savefig(); plt.clf()
        map1.plot_submap()
        pdf.savefig(); plt.clf()
    TAD_st, _ = map1.get_locations(path+'/'+cell+'/'+domain, st=0, ch=0, po=1, add=0)
    TAD_ed, _ = map1.get_locations(path+'/'+cell+'/'+domain, st=0, ch=0, po=2, add=-1)
    TAD = zip(TAD_st, TAD_ed)

    W = np.asarray(map1.contact_group * map1.group_map)
    Wsum = W.sum(1)
    W /= Wsum[Wsum>0].mean()
    gini = 1-np.power(W,2).sum(1)
    gini[Wsum==0] = 0

    if loci is not None:
        loc = map1.choose_map_loc(loci)
    else:
        loc = []
    grps = W[loc,:].sum(0) > 0
    map1.output_groups()
    show(loc, True)
    sel = np.arange(0, 40)
#    if pdf is not None:
#        plt.plot(sel, gini[sel], 'k.--')
    for i in xrange(W.shape[1]):
        if grps[i]:
            if pdf is not None:
                plt.plot(sel, W[sel,i], label='C%s'%i)
    tad = []
    for i,j in TAD:
        tad.append(j-i)
        if i in sel and j in sel:
            if pdf is not None:
                plt.plot([i,j], [1.1,1.1], 'k-', linewidth=2)
    if pdf is not None:
        plt.plot(loc, [1]*len(loc), 'r.')
        xt = sel[::(len(sel)/5)]
        plt.xticks(xt, ['%sM'%(X*0.04+st*1e-6) for X in xt])
        plt.ylim([0,1.2])
        plt.xlim([sel.min(), sel.max()])
        pdf.savefig(); plt.clf()
    if pdf is not None:
        map1.plot_map(map1.contact_map[sel,:][:,sel])
        pdf.savefig(); plt.clf()
        map1.plot_map(map1.contact_group[sel,:]*\
                      map1.group_map*\
                      map1.contact_group[sel,:].T,
                      vmin=0.01, title='H*S*H.T')
        pdf.savefig(); plt.clf()
    return map1

def example1(path, pdf):
    super1 = '''
chr7	3193004
chr7	4772296
chr7	13599334
chr7	30982397
chr7	31248315
chr7	38812914
chr7	52806853
chr7	56592909
chr7	71092246
chr7	86355826
chr7	87159908
chr7	87274999
chr7	87333420
chr7	91027196
chr7	119831735
chr7	140304156
chr7	147131117
chr7	152036872
    '''
    super2 = '''
chr1	13049615
chr1	13103509
chr1	13117753
chr1	13445648
chr1	13650476
chr1	13928434
chr1	13990867
chr1	14179691
chr1	14215552
chr1	14282354
chr1	14302611
chr1	14410966
'''
    list1 = [line.replace('\t',':').strip() for line in super2.split('\n') if len(line)>5]
    chr_detail(path=path, cell='mESC', genome='mm9', ch='chr1', loci=list1, st=10e6, ed=20e6, pdf=pdf)
    chr_detail(path=path, cell='mCortex', genome='mm9', ch='chr1', loci=list1, st=10e6, ed=20e6, pdf=pdf)
#    w1,w2 = map1.best_cor(map2)
#    pdf.savefig(); plt.clf();
#    plt.plot(w1[:520,105], 'b.-')
#    plt.plot(w2[:520,105], 'r.-')
#    pdf.savefig(); plt.clf();

def example2(path, pdf):
    list1 = []
    chr_detail(path=path, cell='mESC', genome='mm9', ch='chr1', loci=None, pdf=pdf)
    chr_detail(path=path, cell='mCortex', genome='mm9', ch='chr1', loci=None, pdf=pdf)

def main(para):
    pdf = PdfPages(para['ExeFile']+'plot.pdf')
    path = para['DataPath']+'/Dixon2012N'
#    example1(path, pdf)
#    example2(path, pdf)
    one_cell(path, pdf, 'IMR90', 'hg18',)
#    one_cell(path, pdf, 'hESC', 'hg18',)
#    one_cell(path, pdf, 'mESC', 'mm9',)
#    one_cell(path, pdf, 'mCortex', 'mm9',)
    pdf.close()

if __name__ == "__main__": main_fun(main)
