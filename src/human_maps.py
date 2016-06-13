from tools import *
from contact_map import *

def check_loci(map1, datapath='../data', dataset='dbgap', method='AvgCCD', dims=None, plot=None):
    ## check SNPs
    chrs = []; locs = []; names = []; groups = [];
    if dataset == 'dbgap':
        with open(datapath+'/dbgap_web_20130413.txt','r') as infile:
            for line in infile:
                database, meshterm, snpid, pvalue, \
                    chromosome, position = line.strip().split('\t')[:6]
                if int(position) < 10: ## invalid position
                    continue 
                if float(pvalue) > 1e-2:
                    continue
                chrs.append(map1.chr2idx['chr'+chromosome])
                locs.append(int(position))
                groups.append(meshterm)
                names.append('Unknown')
    elif dataset == 'gwasdb':
        with open(datapath+'/GWASdb_snp_v2.txt','r') as infile:
            for line in infile:
                snpid, pvalue, pubmed, meshterm, source, \
                    chromosome, position = line.strip().split('\t')
                if float(pvalue) > 0.01:
                    continue
                try:
                    chrs.append(map1.chr2idx['chr'+chromosome])
                    locs.append(int(position))
                    groups.append(meshterm)
                    names.append('Unknown')
                except:
                    print 'Skip data', chromosome, position
                    pass
    elif dataset == 'HouseKeeping':
        with open(datapath+'/Eisenberg2013TG/HK_genes_strong.txt','r') as infile:
            infile.readline()
            for line in infile:
                ele = line.split('\t')
                if ele[3] not in map1.chr2idx:
                    continue
                chrs.append(map1.chr2idx[ele[3]])
                locs.append(int(ele[4]))
                groups.append('HouseKeeping')
                names.append(ele[0])
    elif dataset.startswith('sup-'):
        target = dataset[4:]
        infile = open(datapath+'/Hnisz2013C/'+target+'.csv', 'r')
        infile.readline() ## head
        for line in infile:
            ele = line.strip().split(',')
            is_sup = ele[6]
            if int(is_sup) == 0: ## common enhancer
                continue
            if ele[1] in map1.chr2idx:
                chrs.append(map1.chr2idx[ele[1]])
                locs.append(int(ele[2]))
                names.append(ele[4])
                groups.append('super-'+target)
    elif dataset.startswith('gene-'):
        import gzip
        target = dataset[5:]
        if target == 'rRNA' or target == 'HLA-':
            infile = gzip.open(datapath+'/ENCODE/gencode.v7.annotation.gtf.gz', 'rb')
        elif target == 'tRNA':
            infile = gzip.open(datapath+'/ENCODE/gencode.v7.tRNAs.gtf.gz', 'rb')
        else:
            raise ValueError('Unknown dataset %s'%dataset)
        all_groups = set()
        known = set()
        for line in infile:
            if line[0] == '#':
                continue
            ele = line.split('\t')
            if ele[0] in map1.chr2idx:
                key_values = ele[8]
                ent = key_values.split()
                gene_type = ent[5].replace('"','').replace(';','')
                gene_name = ent[9].replace('"','').replace(';','')
                all_groups.add(gene_type)
                if (gene_name.find(target) >= 0 or gene_type.find(target) >= 0) and gene_type.find('pseudogene') < 0:
                    ch = map1.chr2idx[ele[0]]
                    po = int(ele[3])
                    if (ch,po) not in known: ## overlapped annotations
                        chrs.append(ch)
                        locs.append(po)
                        names.append(gene_name)
                        groups.append(gene_type)
                        known.add((ch,po))
        print 'All groups:', all_groups
    elif dataset.startswith('exp-'):
        target = dataset[4:]
        infile = open(datapath+'/ENCODE/expression_table.txt', 'r')
        head = infile.readline().split('\t')
        cidx = -1
        for i in xrange(len(head)):
            if head[i].strip() == target:
                cidx = i
                break
        assert cidx >= 0
        for line in infile:
            ele = line.strip().split('\t')
            if ele[1] in map1.chr2idx:
                chrs.append(map1.chr2idx[ele[1]])
                locs.append(int(ele[2]))
                names.append(ele[0])
                groups.append(ele[cidx])
        infile.close()
    elif dataset.startswith('highexp-'):
        target = dataset[8:]
        infile = open(datapath+'/ENCODE/expressed_genes.txt', 'r')
        infile.readline() ## header
        for line in infile:
            ele = line.strip().split('\t')
            if ele[4] != target:
                continue
            if ele[1] in map1.chr2idx:
                chrs.append(map1.chr2idx[ele[1]])
                locs.append(int(ele[2]))
                names.append(ele[0])
                groups.append(ele[4])
        infile.close()
    elif dataset.startswith('gap-'):
        target = dataset[4:]
        infile = open(datapath + '/hg19_gap.txt', 'r')
        for line in infile:
            ele = line.split('\t')
            ch = ele[1]
            po = ele[2]
            group = ele[7]
            if group != target:
                continue
            if ch in map1.chr2idx:
                chrs.append(map1.chr2idx[ch])
                locs.append(int(po))
                names.append(ch+':'+po)
                groups.append(target)
        infile.close()
    elif dataset == 'SNP-rtQTL':
        infile = open(datapath+'/Koren2014C/TableS4.rtQTL.txt', 'r')
        infile.readline()
        for line in infile:
            ele = line.split()
            chrs.append(map1.chr2idx['chr'+ele[2].replace('23','X')])
            locs.append(int(ele[3]))
            names.append(ele[0])
            groups.append('rtQTL')
        infile.close()
    elif dataset == 'SNP-RepVar':
        infile = open(datapath+'/Koren2014C/TableS3.RepVar.txt', 'r')
        infile.readline()
        for line in infile:
            ele = line.split()
            chrs.append(map1.chr2idx['chr'+ele[1].replace('23','X')])
            locs.append(int(ele[2]))
            names.append(ele[0])
            groups.append('RepVar')
        infile.close()
    elif dataset.startswith('COSMIC-'):
        target = dataset[7:]
        import gzip
        infile = gzip.open(datapath+'/COSMIC/CosmicStructExport.tsv.gz', 'rb')
        for line in infile:
            if line.strip() == 'INDIVIDUAL BREAKPOINTS':
                break
        infile.readline()
        infile.readline()
        cc = 0
        for line in infile:
            ele = line.split('\t')
            cc += 1
            if line.find(target) <0:
                continue
            ch1 = 'chr'+ele[8].replace('23','X')
            ch2 = 'chr'+ele[12].replace('23','X')
            if ch1 not in map1.chr2idx or ch2 not in map1.chr2idx:
                continue
            if ele[5].find('unknown') >= 0 and ele[8] != '24':
                chrs.append(map1.chr2idx[ch1])
                locs.append(int(ele[9]))
                names.append('SV'+str(cc))
                groups.append(ele[1])
            if ele[5].find('unknown') >= 0 and ele[12] != '24':
                chrs.append(map1.chr2idx[ch2])
                locs.append(int(ele[13]))
                names.append('SV'+str(cc))
                groups.append(ele[1])
        infile.close()
    else:
        raise ValueError('Unknown dataset %s'%dataset)
    if len(chrs) == 0:
        return -1, -1, -1
    ## save to a file
    outfile = open('loci_for_%s.txt'%dataset, 'w')
    outfile.write('Bin\tChr\tPos\tName\tGroup\n')
    idx = map1.choose_map_loc(np.array(chrs), np.array(locs))
    output = zip(idx, chrs, locs, names, groups)
    output.sort()
    for i,c,p,n,g in output:
        outfile.write('%s\t%s\t%s\t%s\t%s\n'%(i,map1.idx2chr[c],p,n,g))
    outfile.close()
    #sub = np.array([i for i in xrange(len(groups)) if groups[i]==g])
    #idx = map1.choose_map_loc(np.array(chrs)[sub], np.array(locs)[sub])
    idx = map1.choose_map_loc(np.array(chrs), np.array(locs))
    if dims is None:
        srt, val, pval = map1.test_enrichment(idx=idx, method=method, title='%s on %s'%(g, map1.name), plot=plot)
    else:
        srt, val, pval = map1.test_enrichment_dims(dims=dims, idx=idx, method=method, title='%s on %s'%(g, map1.name), plot=plot)
#   show([map1.name, g, len(idx)])
#   show(pval[srt[0]])
#   show(val[srt[0]])
#   show(srt[:3].tolist(), True)
    return srt[0], val[srt[0]], pval[srt[0]]

def check_vector(map1, datapath='../data', dataset='gene-tss', method='PCC', dims=None, plot=None):
    v = np.empty(map1.contact_group.shape[0])
    if dataset == 'TSS':
        import gzip
        known = set()
        infile = gzip.open(datapath+'/gencode.v19.annotation.gtf.gz', 'rb')
        chrs = []; locs = [];
        for line in infile:
            if line[0] == '#':
                continue
            ele = line.split('\t')
            if ele[0] in map1.chr2idx:
                key_values = ele[8]
                ent = key_values.split()
                gene_type = ent[5].replace('"','').replace(';','')
                gene_name = ent[9].replace('"','').replace(';','')
                ch = map1.chr2idx[ele[0]]
                po = int(ele[3])
                if (ch,po) not in known: ## overlapped annotations
                    known.add((ch,po))
                    chrs.append(ch)
                    locs.append(po)
        idx = map1.choose_map_loc(np.array(chrs), np.array(locs))
        v = np.zeros(v.shape)
        for i in idx:
            if i >= 0:
                v[i] += 1
    elif dataset.startswith('RE-'): 
        cell = dataset[3:]
        infile = open(datapath+'/Ryba2010GR/RT_'+cell+'_All.txt', 'r')
        for i in xrange(16):
            infile.readline() ## skip
        chrs = []; locs = []; vals = [];
        for line in infile:
            Id, Name, Chr, Sta, End, Val = line.split()
            if Chr in map1.chr2idx:
                ch = map1.chr2idx[Chr]
                po = (int(Sta)+int(End))/2
                chrs.append(ch)
                locs.append(po)
                vals.append(float(Val)*(int(End)-int(Sta)))
        infile.close()
        idx = map1.choose_map_loc(np.array(chrs), np.array(locs))
        v = np.zeros(v.shape)
        c = np.ones(v.shape)
        for i,j in zip(idx, vals):
            if i >= 0:
                v[i] += j
    elif dataset.startswith('DN-'):
        cell = dataset[3:]
        import gzip
        infile = gzip.open(datapath+'/ENCODE/wgEncodeAwgDnase'+cell+'UniPk.narrowPeak.gz', 'rb')
        chrs = []; locs = []; vals = [];
        for line in infile:
            Ch, st, ed, na, sc, sd, sv, pv, sq, pk = line.split('\t')
            if Ch in map1.chr2idx:
                ch = map1.chr2idx[Ch]
                po = (int(st)+int(ed))/2
                chrs.append(ch)
                locs.append(po)
                vals.append(float(sv)*(int(ed)-int(st)))
        infile.close()
        idx = map1.choose_map_loc(np.array(chrs), np.array(locs))
        v = np.zeros(v.shape)
        for i,j in zip(idx, vals):
            if i >= 0:
                v[i] += j
    elif dataset == 'GC-content':
        v = map1.get_gc_content(datapath+'/hg19_count.txt')
        idx = v
    else:
        raise ValueError('Unknown dataset %s'%dataset)
    ## save to a file
    outfile = open('loci_for_%s.txt'%dataset, 'w')
    outfile.write('Bin\tValue\n')
    for i in xrange(len(v)):
        outfile.write('%s\t%s\n'%(i,v[i]))
    outfile.close()
    if dims is None:
        srt, val, pval = map1.test_enrichment(idx=v, method=method, title='%s on %s'%(dataset, map1.name), plot=plot)
    else:
        srt, val, pval = map1.test_enrichment_dims(dims=dims, idx=v, method=method, title='%s on %s'%(dataset, map1.name), plot=plot)
#    show(map1.name)
#    show(dataset)
#    show(len(idx))
#    show(pval[srt[0]])
#    show(val[srt[0]])
#    show(srt[:3].tolist(), True)
    return srt[0], val[srt[0]], pval[srt[0]]

def compare_dims(names, plot=None, k=4):
    v = []; u = [];
    for i in names:
        map1 = ContactMap(i)
        map1.load()
#        map1.sort_groups('diagnal')
        for j in xrange(k):
            v.append(np.asarray(map1.contact_group[:,j]*map1.group_map[j,j]))
            u.append(i+':D%s'%(j+1))
    c = np.corrcoef(np.hstack(v).T)
    if plot is not None:
        plt.imshow(c, interpolation='none')
        plt.yticks(range(len(u)), u)
        plt.colorbar()
        plot.savefig(); plt.clf();
        plt.imshow(np.abs(c), interpolation='none', cmap='hot')
        plt.yticks(range(len(u)), u)
        plt.colorbar()
        plot.savefig(); plt.clf();

def compare_maps(names, plot=None):
    m = len(names)
    map1 = ContactMap(names[0])
    map1.load()
    n,r = map1.contact_group.shape
    rel1 = np.zeros((m,m))
    rel2 = np.zeros((m,m))
    rel3 = np.zeros((m,m))
    for i in xrange(m):
        map1 = ContactMap(names[i])
        map1.load()
        for j in xrange(m):
            map2 = ContactMap(names[j])
            map2.load()
            rel1[i,j] = map1.compare(map2, raw=True, metric='JS_D')
            rel2[i,j] = map1.compare(map2, raw=False, metric='JS_D')
            rel3[i,j] = correlation(map1.contact_group[:,0], map2.contact_group[:,0])
    if plot is not None:
        plt.subplots_adjust(left=0.28, right=1)
        plt.imshow(rel1, interpolation='none', vmin=0, vmax=0.1)
        plt.yticks(range(len(names)), names)
        plt.colorbar()
        plot.savefig(); plt.clf();
        plt.imshow(rel2, interpolation='none', vmin=0, vmax=0.1)
        plt.yticks(range(len(names)), names)
        plt.colorbar()
        plot.savefig(); plt.clf();
        plt.imshow(rel3, interpolation='none', vmin=0, vmax=1)
        plt.yticks(range(len(names)), names)
        plt.colorbar()
        plot.savefig(); plt.clf();
        plt.subplots_adjust(left=0.125, right=0.9) ## reset

def match_dims(names):
    maps = {}
    for i in xrange(len(names)):
        for j in xrange(i,len(names)):
            map1 = ContactMap(names[i])
            map2 = ContactMap(names[j])
            map1.load()
            map2.load()
            map1.group_map /= map1.group_map.sum() ## normalize
            map2.group_map /= map2.group_map.sum() ## normalize
            dims = map1.min_rmsd(map2, dims=True)
            for ii,jj in dims:
                d = maps.get((i,ii),set())
                d.add((j,jj))
                maps[(i,ii)] = d
                d = maps.get((j,jj),set())
                d.add((i,ii))
                maps[(j,jj)] = d
    match = []
    for key in maps:
        d = maps[key]
        contain_all = True
        if len(d) < len(names):
            contain_all = False
        for kk in d:
            if len(maps[kk]) < len(names):
                contain_all = False
        if contain_all:
            d = sorted(list(d))
            if d not in match:
                match.append(d)
    match.sort()
    show('''
    Find conserved cluster members
    ''', True)
    show(len(match))
    show('is the number of matched clusters', True)
    return match

def diff_mul(names, match, plot=None):
    overall = []
    for d in match:
        vs = []
        for i,ii in d:
            map1 = ContactMap(names[i])
            map1.load()
            map1.group_map /= map1.group_map.sum() ## normalize
            v = np.array(map1.contact_group[:,ii] * map1.group_map[ii,ii]).reshape(-1)
            vs.append(v)
        high = 0
        low  = 0
        for v in vs:
            high += np.array(v >= 0.3, dtype='int')
            low  += np.array(v <= 0.1, dtype='int')
        sele = np.logical_and(high==1, low==(len(names)-1))
        diff = np.flatnonzero(sele).tolist()
        overall.append(diff)
        if len(diff) < 1:
            continue
        if plot is not None:
            for i,ii in d:
                plt.plot(vs[i], label='%s:%s'%(names[i],ii))
            plt.plot(diff, [0.5]*len(diff), 'k.')
            sub_diff = diff
            plt.xlim([min(sub_diff)-100, max(sub_diff)+100])
            plt.legend()
            plot.savefig(); plt.clf();
        show([ii for i,ii in d])
        show(':')
        show(diff)
        show()

def same_mul(names, match, plot=None):
    overall = []
    cutoff = 0.1
    for d in match:
        vs = []
        for i,ii in d:
            map1 = ContactMap(names[i])
            map1.load()
            map1.group_map /= map1.group_map.sum() ## normalize
            v = np.array(map1.contact_group[:,ii] * map1.group_map[ii,ii]).reshape(-1)
            vs.append(v)
        same = None
        for i in xrange(len(d)-1):
            for j in xrange(i+1,len(d)):
                idx = np.flatnonzero(np.abs(vs[i]-vs[j]) <= cutoff)
                mk1 = np.flatnonzero(vs[i] > cutoff)
                mk2 = np.flatnonzero(vs[j] > cutoff)
                sh = set(idx.tolist()) & set(mk1.tolist()) & set(mk2.tolist())
                #sh = set(mk1.tolist()) & set(mk2.tolist())
                if same is None:
                    same = sh
                else:
                    same = same & sh
        same = sorted(list(same))
        overall.append(same)
        if len(same) < 1:
            continue
        if plot is not None:
            for i,ii in d:
                plt.plot(vs[i], label='%s:%s'%(names[i],ii))
            plt.plot(same, [cutoff]*len(same), 'k.')
            sub_same = same
            if len(same) > 5:
                sub_same = same[:5]
            plt.xlim([min(sub_same)-10, max(sub_same)+10])
            plt.legend()
            plot.savefig(); plt.clf();
        show([ii for i,ii in d])
        show(':')
        show(same)
        show()

def show_mul(names, step=20, plot=None):
    for name in names:
        map1 = ContactMap(name)
        map1.load()
        n = map1.contact_group.shape[0]
        w = np.array(map1.group_map * map1.contact_group.T)
        w /= float(w.sum()/n)
        v = np.power(w,2).sum(0)
        #v = w.max(0)
        idx = np.argsort(v)[::-1][:int(n*0.05)]
        show(name)
        show(idx.tolist())
        show()
        if plot == None:
            continue
        plt.hist(v, 50)
        plt.title(name)
        plot.savefig(); plt.clf()
        map1.plot_map(w[:,idx], log=False)
        plot.savefig(); plt.clf()
    total = 0
    for i in xrange(0,n,step):
        j = min(i+step, n)
        mark = np.ones(j-i)
        count = np.zeros(j-i)
        for name in names:
            map1 = ContactMap(name)
            map1.load()
            if plot == None:
                continue
            w = np.array(map1.group_map * map1.contact_group.T)
            w /= float(w.sum()/n)
            wm = w[:,i:j].mean(1)
            srtw = np.sort(w,0)
            plt.plot(np.arange(i,j), srtw[-1,i:j], '--.', label='%s'%name.replace('hg19-',''))
            #plt.plot(np.arange(i,j), 1-np.power(w,2).sum(0)[i:j], '--.', label='%s'%name.replace('hg19-',''))
            mark = np.logical_and(mark, srtw[-1,i:j] + srtw[-2,i:j] > 0.5)
            count += (srtw[-1,i:j] > 0.5)
        mark = np.logical_and(mark, count == len(names)-1)
        plt.plot(np.arange(i,j)[mark], [0.5]*mark.sum(), 'ro')
        total += mark.sum()
        plt.title('%s:%s - %s:%s'%(
            map1.idx2chr[map1.frag_chr[i]],
            map1.frag_sta[i],
            map1.idx2chr[map1.frag_chr[j-1]],
            map1.frag_end[j-1]))
        plt.ylim([0,1.01])
        plt.legend()
        plot.savefig(); plt.clf();
    print 'In total, we found', total

def get_known(names, match, locfile, plot=None):
    loci = []
    if True:
        infile = open(locfile, 'r')
        infile.readline()
        for line in infile:
            loci.append(int(line.split('\t')[0]))
        infile.close()
    for d in match:
        vs = []
        mk = 0
        for i,ii in d:
            map1 = ContactMap(names[i])
            map1.load()
            map1.group_map /= map1.group_map.sum() ## normalize
            v = np.array(map1.contact_group[:,ii] * map1.group_map[ii,ii]).reshape(-1)
            vs.append(v)
            mk += np.array(v>0.1, dtype='int')
        mark = np.flatnonzero(mk).tolist()
        sub_loci = list(set(loci) & set(mark))
        if len(sub_loci) < 2:
            continue
        if plot is not None:
            for i,ii in d:
                plt.plot(vs[i], label='%s:%s'%(names[i],ii))
            plt.plot(loci, [0.5]*len(loci), 'k.')
            plt.xlim([min(sub_loci)-10, max(sub_loci)+10])
            plt.ylim([0,1])
            plt.legend()
            plt.title(locfile.replace('loci_for_','').replace('.txt',''))
            plot.savefig(); plt.clf();
    for i in xrange(len(names)):
        map1 = ContactMap(names[i])
        map1.load()
        map1.group_map /= map1.group_map.sum() ## normalize
        v = np.array(map1.contact_group * map1.group_map)
        has = []
        for d in match:
            for I,J in d:
                if i == I:
                    has.append(J)
        if plot is None:
            continue
        sumv = np.zeros(len(loci))
        for j in xrange(v.shape[1]):
            if j in has:
                continue
            sumv += v[loci,j]
            mark = np.flatnonzero(v[:,j]>0.1).tolist()
            sub_loci = list(set(loci) & set(mark))
            if len(sub_loci) < 1:
                continue
            plt.plot(v[:,j], label='%s:%s'%(names[i],j))
        if True:
            plt.plot(loci, sumv, 'k.')
            plt.xlim([min(loci)-10, max(loci)+10])
            plt.ylim([0,1])
            plt.legend()
            plt.title(locfile.replace('loci_for_','').replace('.txt','')+' in '+names[i])
            plot.savefig(); plt.clf();            

def check_features(names, dims, plot=None):
    feature1 = [('Genome Features',   ['GC-content', 'TSS']),
                ('Replicating time',  ['RE-GM12878_Lymphoblastoid', 
                                       'RE-H1_ESC',
                                       'RE-IMR90_FIbroblast',
                                       'RE-K562_Leukemia', 
                                       ]),
                ('DNAseI sensitivity',['DN-UwdukeGm12878', 
                                       'DN-UwdukeH1hesc', 
                                       'DN-DukeImr90',
                                       'DN-UwdukeK562', 
                                       ])]
    feature2 = [('Genome Features', ['gap-telomere', 'gene-rRNA', 'gene-tRNA', 'HouseKeeping']),
                ('Cancer Structure Variances', ['COSMIC-lung', 'COSMIC-lymphoid', 'COSMIC-breast']),
                ('Gene Expression', ['highexp-Gm12878', 'highexp-H1hesc', 'highexp-Imr90', 'highexp-K562']),
                ('Super Enhancer',  ['sup-GM12878', 'sup-H1', 'sup-IMR90', 'sup-K562'])]
    for group, cases in feature1 + feature2:
#        show(['Dataset', 'Groups', '#Bins', 'P-value', 'Value', 'Top Clusters'], True)
        outs = []
        for name in names:
            map1 = ContactMap(name)
            map1.load()
            #dims = [map1.contact_group.shape[1]]
            ss = []; vv = []; pp = []
            for case in cases:
                if (group, cases) in feature1:
                    s,v,p = check_vector(map1, dataset=case, method='SPC', dims=dims, plot=None)
                if (group, cases) in feature2:
                    s,v,p = check_loci(map1, dataset=case, method='AvgCCD', dims=dims, plot=None)
                ss.append(s)
                vv.append(v)
                pp.append(p)
            outs.append([name]+ss+['']+vv+['']+pp)
        show(group, True)
        show('Clusters')
        show(cases)
        show('')
        show(cases)
        show('')
        show(cases)
        show()
        for out in outs:
            show(out, True)
        show()

def add_ice_result(ref, k=10):
    data = [('hg19-hic-ice', '../data/Imakaev2012NM/SRR027956_map-res1000k-ic-eig.txt'),
            ('hg19-tcc-ice', '../data/Imakaev2012NM/SRR071231_map-res1000k-ic-eig.txt')]
    map1 = ContactMap(ref)
    assert map1.load()
    names = []
    for name, eigfile in data:
        map2 = map1.duplicate(name)
        map2.contact_group = np.zeros((map2.contact_map.shape[0],k))
        for i in xrange(k):
            idx, val = map2.get_locations(eigfile, st=0, ch=0, po=1, nm=i+3, add=0)
            map2.contact_group[idx,i] = np.array([float(v) for v in val])
        map2.contact_group = np.matrix(map2.contact_group)
        map2.group_map = np.matrix(np.eye(k))
        map2.save()
        names.append(map2.name)
    return names

def run1(datalist, pdf):
    ## Construct Maps
    names = []
    genome = '../data/hg19_chr_len.txt'
    for name, files in datalist:
        files = [f+'_map-res1000k.hdf5' for f in files]
        map1 = ContactMap(name)
#        map1.clear()
        show(name)
        if not map1.load():
            map1.genome_info(genome)
            map1.create_densemap(files, reso=1e6)
            map1.mask_diag()
            map1.mask_short()
            map1.mask_low()
            if map1.name.find('K562') >= 0:
                map1.trim_high(0.05) ## this one is very noisy
            else:
                map1.trim_high(0.005)
        map1.decompose_auto(beta=1, plot=None)
        map1.save()
        if os.path.exists('loci_for_TSS.txt'):
            gene_density = [float(v) for v in map1.fread('loci_for_TSS.txt', column=1, start=1, end=-1)]
            idx, val, pval = map1.test_enrichment(gene_density, method='SPC-mrk')
        map1.sort_groups()
        map1.output_groups()
        show(map1.contact_group.shape)
        show(gini_impurity(np.diag(map1.group_map)))
        show()
        map1.plot_submap()
        pdf.savefig(); plt.clf();
        names.append(name)
    ##########################################################
#    names += add_ice_result(names[0])

    check_features(names, None, plot=pdf)
#    compare_dims(names, plot=pdf)

#    m = match_dims(names)
#    same_mul(names, m, plot=pdf)
#    same_mul(names, m, plot=pdf)
#    show_mul(names, plot=pdf)
#    get_known(names, m, 'loci_for_gene-HLA-.txt', plot=pdf)
#    get_known(names, m, 'loci_for_exp-GM12878.txt', plot=pdf)
#    get_known(names, m, 'loci_for_sup-GM12878.txt', plot=pdf)
#    get_known(names, m, 'loci_for_exp-Imr90.txt', plot=pdf)
#    get_known(names, m, 'loci_for_sup-Imr90.txt', plot=pdf)

def run2(datalist, pdf):
    names = []
    genome = '../data/hg19_chr_len.txt'
    for name, files in datalist:
        for fname in files:
            new_name = name+'-%s'%(fname)
            map1 = ContactMap(new_name)
#            map1.clear()
            print new_name
            if True:
                map1.genome_info(genome)
                map1.create_densemap(['human/'+fname+'_map-res1000k.hdf5'], 1e6)
                map1.mask_diag()
            map1.decompose_auto(plot=None)
            map1.plot_submap()
            pdf.savefig(); plt.clf();
            map1.save()
            names.append(new_name)
#    check_features(names, plot=pdf)
#    show_mul(names, plot=pdf) 
#    same_mul(names, plot=pdf)

def main(para):
    show('''
    Compare HiC heatmap from Human genomes
    ''', True)
    datalist0 = [("hg19-hic-nmf", ["SRR027956"]), 
                 ("hg19-tcc-nmf", ["SRR071231"])]
    datalist1 = [
        ("hg19-GM12878", ["SRR071231", "SRR071232"]),
        ("hg19-hESC", ["SRR400261", "SRR400262", "SRR400263",
                       "SRR442155", "SRR442156", "SRR442157"]),
        ("hg19-IMR90", ["SRR400264", "SRR400265", "SRR400266", "SRR400267",
                      "SRR400268", "SRR442158", "SRR442159", "SRR442160"]),
        ("hg19-K562", ["SRR027962", "SRR027963"]),
                ]
    datalist2 = [
        ("Kalhor2012HiC", ["SRR071233", "SRR071234"]),
        ("Kalhor2012TCC", ["SRR071231", "SRR071232"]),
        ("Lieber2009GM", ["SRR027956", "SRR027957", "SRR027958", "SRR027959"]),
        ("Lieber2009K5", ["SRR027962", "SRR027963"]),
        ("Dixon2012hESC-rep1", ["SRR400260", "SRR400261", "SRR400262", "SRR400263"]),
        ("Dixon2012hESC-rep2", ["SRR442155", "SRR442156", "SRR442157"]),
        ("Dixon2012IMR90-rep1", ["SRR400264", "SRR400265", "SRR400266", "SRR400267", "SRR400268"]),
        ("Dixon2012IMR90-rep2", ["SRR442158", "SRR442159", "SRR442160"])]
    if 'SetStart' not in para:
        para['SetStart'] = '0'
    if 'SetEnd' not in para:
        para['SetEnd'] = 'None'
    pdf = PdfPages(para['ExeFile']+'plot%s%s.pdf'%(para['SetStart'], para['SetEnd']))
    run1(datalist1, pdf)
#    run2(datalist2, pdf)
    pdf.close()

if __name__ == '__main__': main_fun(main)
