from tools import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from contact_map import ContactMap

def one_region(path, cell, genome, ch, bi, loci, st=0, ed=None, pdf=None):
    if bi.endswith('kb'):
        reso = int(bi.replace('kb',''))*1000
    elif bi.endswith('mb'):
        reso = int(bi.replace('mb',''))*1000000
    else:
        raise ValueError('Unknow unit %s'%bi)
    map1 = ContactMap('loop-%s-in-%s'%(ch, cell))
    map1.clear()
    if not map1.load():
        map1.genome_info(path+'/%s_chr_len.txt'%genome)
        map1.create_binnedmap(reso, lazy=True)
        map1.focus_chromosome(ch, st=st, ed=ed)
        if True: ## read files
            norm = []
            infile = open(path+'/'+cell+'/'+bi+'_resolution_intrachromosomal/'+ch+'/MAPQGE30/'+ch+'_'+bi+'.KRnorm')
            for line in infile:
                norm.append(float(line))
            infile.close()
            expect = []
            infile = open(path+'/'+cell+'/'+bi+'_resolution_intrachromosomal/'+ch+'/MAPQGE30/'+ch+'_'+bi+'.KRexpected')
            for line in infile:
                expect.append(float(line))
            infile.close()
            expect.append(1.0)
            print len(norm), len(expect)
            infile = open(path+'/'+cell+'/'+bi+'_resolution_intrachromosomal/'+ch+'/MAPQGE30/'+ch+'_'+bi+'.RAWobserved', 'r', 2<<9)
            p1 = []; p2 = []; val = []
            for line in infile:
                P1,P2,Val = line.split()
                pos1 = int(P1)
                pos2 = int(P2)
                if pos1 < st or (ed is not None and pos1 >= ed):
                    continue
                if pos2 < st or (ed is not None and pos2 >= ed): 
                    continue
                p1.append(pos1)
                p2.append(pos2)
                I = pos1 / reso
                J = pos2 / reso
                IJ = abs(pos1-pos2) / reso
#                val.append(float(Val))
                val.append(float(Val)/(norm[I]*norm[J]))
#                val.append(float(Val)/(norm[I]*norm[J]*expect[IJ]))
            map1.inter_loc1 = np.array(p1, dtype='int')
            map1.inter_loc2 = np.array(p2, dtype='int')
            map1.inter_freq = np.array(val, dtype='float')
            chidx = map1.chr2idx[ch]
            map1.inter_chr1 = chidx*np.ones(len(p1), dtype='int')
            map1.inter_chr2 = chidx*np.ones(len(p2), dtype='int')
            infile.close()
        map1.create_binnedmap(reso)
        map1.save()
    show(cell)
    show(ch)
    if pdf is not None:
        map1.plot_map()
        pdf.savefig(); plt.clf()
    map1.decompose_auto()
    map1.sort_groups()
    show(map1.contact_group.shape)
    show()
    bins = map1.choose_map_loc(loci)
    
    W = np.asarray(map1.contact_group * map1.group_map)
    n,r = W.shape
    wm = W.sum(1)
    W /= np.mean(wm[wm>0])
    gini = 1-np.power(W,2).sum(1)
    gini[wm==0] = 0

    outfile = open('loop-%s-in-%s_groups.wig'%(ch, cell), 'w')
#    outfile.write('track type=wiggle_0 name="Overall" description="BNMF" visibility=full autoScale=off viewLimits=800:1000 color=0,0,0 maxHeightPixels=100:50:20 graphType=bar priority=20\nfixedStep chrom='+ch+' start=%d'%st+' step=%d'%reso+' span=%d\n'%reso)
#    for i in xrange(n):
#        outfile.write('%d\n'%int(1000*gini[i]))
    jj = []; ww = 0
#    for j in xrange(r):
#        if W[bins,j].max() < 0.1:
#            continue
    for j in W[bins,:].argmax(1):
        ww += W[:,j]
        outfile.write('track type=wiggle_0 name="C%s'%(j+1)+'" description="BNMF" visibility=full autoScale=off viewLimits=0:200 color=0,0,0 maxHeightPixels=100:50:20 graphType=bar priority=20\nfixedStep chrom='+ch+' start=%d'%st+' step=%d'%reso+' span=%d\n'%reso)
        for i in xrange(n):
            outfile.write('%d\n'%int(1000*W[i,j]))
        jj.append(j)
#    outfile.write('track type=wiggle_0 name="Overall" description="BNMF" visibility=full autoScale=off viewLimits=0:200 color=0,0,0 maxHeightPixels=100:50:20 graphType=bar priority=20\nfixedStep chrom='+ch+' start=%d'%st+' step=%d'%reso+' span=%d\n'%reso)
#    for i in xrange(n):
#        outfile.write('%d\n'%int(1000*ww[i]))
#    outfile.close()

    sel = range(n)
    lab = ['%dk'%((i*reso+st)/1000) for i in sel]
    five = np.arange(0, len(sel), len(sel)/5)
    if pdf is not None:
        map1.plot_map(map1.contact_group*map1.group_map*map1.contact_group.T, log=False)
        pdf.savefig(); plt.clf()
        map1.plot_map(map1.contact_group[:,jj]*map1.group_map[jj,:][:,jj]*map1.contact_group[:,jj].T, log=False)
        pdf.savefig(); plt.clf();
        map1.plot_submap()
        pdf.savefig(); plt.clf()
        plt.plot(sel, gini[sel], '--k')
        for j in jj:
            plt.plot(sel, W[sel,j], '-', label='C%s'%(j+1))
#        plt.plot(sel, ww[sel], '-', label='Combined')
        plt.plot(bins, [1.1]*len(bins), 'ro')
        plt.legend()
        plt.xticks([sel[j] for j in five], [lab[j] for j in five])
        pdf.savefig(); plt.clf()
    return

def one_cell(path, pdf, cell, genome):
    map1 = ContactMap()
    map1.genome_info(path+'/%s_chr_len.txt'%genome)
    for i in sorted(map1.idx2chr.keys()):
        info = one_chr(path=path, cell=cell, genome=genome, ch=map1.idx2chr[i], pdf=pdf)

def main(para):
    pdf = PdfPages(para['ExeFile']+'plot.pdf')
    path = para['DataPath']+'/Rao2014C'
#    one_region(path, cell='IMR90', genome='hg19', ch='chr1', bi='1mb', pdf=pdf)
#    one_region(path, cell='GM12878_combined', genome='hg19', ch='chr1', bi='1kb', st=111e6, ed=112e6, pdf=pdf)
#    one_region(path, cell='GM12878_combined', genome='hg19', ch='chr1', bi='1kb', st=160e6, ed=161e6, pdf=pdf)
    one_region(path, cell='K562', genome='hg19', ch='chr12', bi='5kb', loci=['chr12:54367698', 'chr12:54378108'], st=53e6, ed=56e6, pdf=pdf)
#    one_region(path, cell='GM12878_combined', genome='hg19', ch='chr12', bi='1kb', loci=['chr12:54360413', 'chr12:54367698'], st=54e6, ed=55e6, pdf=pdf)
#    one_cell(path, pdf, 'GM12878_combined', 'hg19')
    pdf.close()

if __name__ == "__main__": main_fun(main)
