from tools import *
from contact_map import ContactMap

def gene_to_term(filename):
    infile = open(filename, 'r')
    cc = 0
    for line in infile:
        orf, gene, sgdid, golink, term, goid, ftype = line.split('\t')
        if sgdid == '' or goid == '':
            continue
        yield sgdid, goid
        cc += 1
    infile.close()
    print 'Read in', cc, 'pairs from', filename

def read_map(mapfile):
    mp = {}
    for gene, term in gene_to_term(mapfile):
        gs = mp.get(term, [])
        gs.append(gene)
        mp[term] = gs
    print 'There are', len(mp), 'unique GO IDs'
    return mp

def read_complex(comfile):
    co = {}
    with open(comfile, 'r') as tempfile:
        for line in tempfile:
            cc, ge = line.split('\t')
            term = 'GO:'+cc.split('/')[-1].zfill(7)
            gene = [g for g in ge.split('/') if len(g)==10 and g[0]=='S']
            co[term] = gene
    print 'There are', len(co), 'protein complex'
    return co

def read_go(gofile):
    go = {}
    with open(gofile, 'r') as tempfile:
        for line in tempfile:
            ele = line.split('\t')
            go[int(ele[0])] = (ele[1], ele[2], ele[3].strip())
    print 'Read in', len(go), 'GO terms from', gofile
    return go

def read_gene(genefile):
    gene = {}
    with open(genefile, 'r') as tempfile:
        for line in tempfile:
            ele = line.split('\t')
            gid = ele[0]
            Chr = ele[8]
            Sta = ele[9]
            End = ele[10]
            if Chr == '' or Sta == '' or End == '':
                continue
            gene[gid] = (Chr, int(Sta), int(End))
    print 'Read in', len(gene), 'gene locations from', genefile
    return gene

def save_gogene(go, gene, loci):
    fname = '%s.txt'%go.replace(':','_')
    outfile = open(fname, 'w')
    outfile.write('GeneID\tChr\tStart\tEnd\n')
    for ge in gene[go]:
        Chr, Sta, End = loci[ge]
        if Chr.find('micron') >= 0:
            continue
        outfile.write('%s\t%s\t%s\t%s\n'%(ge, Chr, Sta, End))
    outfile.close()
    return fname

def gene_bin(para, data='YeastSyn'):
    show("""
    Check gene distributions among the bins in heatmap
    """, True)
    map1 = ContactMap()
    assert map1.load(data, False)
    GE = read_gene(para['DataPath']+'/SGD/SGD_features.tab')
    bin_idx = []
    for gid in GE:
        ch, st, ed = GE[gid]
        try:
            idx = map1.choose_map_loc([int(ch)], [int(st)-1])
            for i in idx:
                if i >= 0:
                    bin_idx.append(i)
        except:
            print 'Skip', ch, st, gid
    n = map1.contact_map.shape[0]
    cout = histogram(bin_idx, range(n), False)
    for i in xrange(n):
        show(map1.idx2chr[map1.frag_chr[i]])
        show(map1.frag_sta[i])
        show(map1.frag_end[i])
        show(cout[i])
        show()

def go_test(para,  data='YeastHiC'):
    show("""
    Check gene groups relating to the same GO term
    """, True)
    MP = read_map(para['DataPath']+'/SGD/go_slim_mapping.tab')
    CO = read_complex(para['DataPath']+'/SGD/go_protein_complex_slim.tab')
    GO = read_go(para['DataPath']+'/SGD/go_terms.tab')
    GE = read_gene(para['DataPath']+'/SGD/SGD_features.tab')
    cc = 0
    go2gene = MP.copy()
    go2gene.update(CO)

    map1 = ContactMap()
    assert map1.load(data)
    map1.output_groups()
    for go in go2gene:
        go = go.strip()
        fname = save_gogene(go, go2gene, GE)
        idx, names = map1.get_locations(fname, st=1, ch=1, po=2, nm=0)
        os.remove(fname)
        if len(idx) < 1:
            continue
        srt, val, pval = map1.test_enrichment(idx, method='AvgCCD')
        cc += 1
        sign = []
        for i in srt:
            if pval[i] < 0.01:
                sign.append(i)
        if len(sign) > 0:
            show(go)
            show(len(idx))
            show(GO[int(go.split(':')[1])])
            show(pval[sign[0]])
            show(pval[sign[-1]])
            show(sign)
            show()
    show('We tested %s GO terms for %s.\n'%(cc,data))

def main(para):
    go_test(para, 'YeastHiC')
    go_test(para, 'YeastSyn')
#    gene_bin(para)

if __name__ == '__main__': main_fun(main)
