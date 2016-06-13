from tools import *
from contact_map import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def process(name, genome, files, binsize, beta=1, pdf=None):
    map1 = ContactMap(name)
#    map1.clear()
    if not map1.load():
        map1.genome_info(genome)
        map1.create_densemap(files, reso=binsize)
        map1.mask_diag()
        map1.mask_short()
    map1.decompose_auto(beta=beta, plot=None)
    map1.sort_groups()
    map1.save()
    if pdf != None:
        map1.plot_submap()
        pdf.savefig(); plt.clf();
    return map1

def get_yeast(para):
    path = para['DataPath'] + '/Duan2010N'
    genome = para['DataPath'] + '/yeast_chr_len.txt'
    binsize = 1e4
    data = {}
    for filename in os.listdir(path):
        if filename.startswith('interactions_'):
            if filename.find('_fdr') >= 0:
                continue
            name = filename.replace('interactions_','')
            if filename.endswith('_inter.txt'):
                name = name.replace('_inter.txt','')
            elif filename.endswith('_intra.txt'):
                name = name.replace('_intra.txt','')
            else: ## skip intra_all
                continue
            files = data.get(name, [])
            files.append(path+'/'+filename)
            data[name] = files
    pdf = PdfPages(para['ExeFile']+'_yeast.pdf')
    maps = []
    for name in sorted(data.keys()):
        show(name)
        map1 = process(name.replace('_beforeFDR','').replace('_noFDR','').replace('uncrosslinked','control'), genome, data[name], binsize, pdf=pdf)
        show(map1.contact_group.shape)
        idx, genes = map1.get_locations(para['DataPath']+'/Duan2010N/trna_positions_and_clusters.txt', st=1, ch=2, po=3, nm=0)
        srt, val, pval = map1.test_enrichment(idx, method='AvgCCD', normp=False, plot=None)
        show(len(idx))
        show(srt[0])
        show(pval[srt[0]])
        maps.append(map1)
        show()
    cmp_maps(maps, pdf)
    pdf.close()
    return maps

def get_plasmo(para):
    path = para['DataPath'] + '/Ay2014GR'
    genome = para['DataPath'] + '/plasmo_chr_len.txt'
    binsize = 1e4
    data = {}
    for filename in os.listdir(path):
        if filename.startswith('GSM'):
            ## need to contain the following keywords
            if filename.find('raw') < 0:
                continue
            if filename.find('1RE') < 0:
                continue
            name = filename.split('_')[1].replace('.txt.gz','')
            name = name.replace('rawContactCounts','Raw')
            name = name.replace('normalizedContactCounts','ICE')
            name = 'plasmo-'+name
            files = data.get(name, [])
            files.append(path+'/'+filename)
            data[name] = files
    pdf = PdfPages(para['ExeFile']+'_plasmo.pdf')
    maps = []
    for name in sorted(data.keys()):
        show(name)
        map1 = process(name.replace('Raw-1RE-',''), genome, data[name], binsize, pdf=pdf)
        show(map1.contact_group.shape)
        idx, names = map1.get_locations(path + '/rDNA-geneIDs.txt', st=1, ch=1, po=2, nm=4)
        print idx
        srt, val, pval = map1.test_enrichment(idx, method='AUC', plot=None)
        show(len(idx))
        show(srt[0])
        show(pval[srt[0]])
        maps.append(map1)
        show()
    cmp_maps(maps, pdf)
    pdf.close()
    return maps

def get_mouse(para):
    path = para['DataPath'] + '/Nagano2013N'
    genome = para['DataPath'] + '/mm9_chr_len.txt'
    binsize = 1e6
    data = {}
    import gzip
    for filename in os.listdir(path):
        if filename.startswith('GSM') and not filename.endswith('-new.gz'):
            name = 'mouse-'+filename.split('_')[1].replace('.txt.gz','')
            if name.find('cell') < 0:
                continue
            newfile = filename+'-new.gz'
            if not os.path.exists(path+'/'+newfile):
                infile = gzip.open(path+'/'+filename, 'rb')
                tempfile = gzip.open(path+'/'+newfile, 'wb')
                tempfile.write(infile.readline().strip()+'\tfreq\n')
                for line in infile:
                    ch1, po1, ch2, po2 = line.strip().split('\t')
                    if ch1 == 'M' or ch2 == 'M':
                        continue
                    Ch1 = int(ch1.replace('X','20'))-1
                    Ch2 = int(ch2.replace('X','20'))-1
                    tempfile.write('%s\t%s\t%s\t%s\t1\n'%(Ch1, po1, Ch2, po2))
                infile.close()
                tempfile.close()
            files = data.get(name, [])
            files.append(path+'/'+newfile)
            data[name] = files
    pdf = PdfPages(para['ExeFile']+'_mouse.pdf')
    maps = []
    for name in sorted(data.keys()):
        show(name)
        map1 = process(name, genome, data[name], binsize, pdf=pdf)
        show(map1.contact_group.shape)
        idx, names = map1.get_locations(path+'/tRNAs.txt.gz', st=0, ch=1, po=2, nm=4)
        srt, val, pval = map1.test_enrichment(idx, method='AUC', plot=None)
        show(len(idx))
        show(srt[0])
        show(pval[srt[0]])
        maps.append(map1)
        show()
    cmp_maps(maps, pdf)
    pdf.close()
    return maps

def cmp_maps(maps, plot):
    m = len(maps)
    cor_map = np.zeros((2*m, 2*m))
    from scipy.spatial.distance import squareform
    for i in xrange(2*m):
        if i < m:
            show(maps[i].name + ' Before')
        else:
            show(maps[i-m].name + ' After')
        for j in xrange(2*m):
            if i < m: 
                mapA = maps[i].contact_map
            else: 
                mapA = maps[i-m].contact_group * maps[i-m].group_map * maps[i-m].contact_group.T
            if j < m: 
                mapB = maps[j].contact_map
            else: 
                mapB = maps[j-m].contact_group * maps[j-m].group_map * maps[j-m].contact_group.T
            cor_map[i,j] = metric_js_div(mapA, mapB)
            show(cor_map[i,j])
        show()

    plt.subplots_adjust(left=0.2, right=1)
    vm = 0.4
    plt.imshow(cor_map, interpolation='none', vmin=0, vmax=vm)
    cbar = plt.colorbar(ticks=[0, vm/4, vm/2, vm*3/4, vm, vm+0.1])
    cbar.set_ticklabels(['0', '%.2f'%(vm/4), '%.2f'%(vm/2), '%.2f'%(vm*3/4), '>%.2f'%vm, ''])
    plt.yticks(range(2*m), [i.name+':$A$' for i in maps]+[i.name+':$R$' for i in maps])
    plot.savefig(); plt.clf();
    plt.subplots_adjust(left=0.125, right=0.9)

def main(para):
#    maps1 = get_yeast(para)
#    maps2 = get_plasmo(para)
    maps3 = get_mouse(para)

if __name__ == "__main__": main_fun(main)
