from tools import *

def read_gene(genefile='../data/SGD/SGD_features.tab', i=[1,2]):
    with open(genefile, 'r') as tempfile:
        for line in tempfile:
            ele = line.split('\t')
            gid = ele[0]
            typ = ele[1]
            nam = ele[3]
            Chr = ele[8]
            Sta = ele[9]
            End = ele[10]
            yield [ele[a] for a in i]

def get_centromere(para):
    genes = read_gene(para['DataPath']+'/SGD/SGD_features.tab', i=[1,3,8,9,10])
    outfile = 'geneset_centromere.txt'
    out = open(outfile, 'w')
    out.write('Chro\tStart\tEnd\tName\n')
    types = set()
    for Type, Name, Chro, Start, End in genes:
        types.add(Type)
        if Chro == '' or Start == '':
            continue
        if Type != 'centromere':
            continue
        out.write('%s\t%s\t%s\t%s\n'%(Chro, Start, End, Name))
    print types
    return outfile

def get_telomere(para):
    genes = read_gene(para['DataPath']+'/SGD/SGD_features.tab', i=[1,3,8,9,10])
    outfile = 'geneset_telomere.txt'
    out = open(outfile, 'w')
    out.write('Chro\tStart\tEnd\tName\n')
    for Type, Name, Chro, Start, End in genes:
        if Chro == '' or Start == '':
            continue
        if Type != 'telomere':
            continue
        out.write('%s\t%s\t%s\t%s\n'%(Chro, Start, End, Name))
    return outfile

def get_tRNA(para):
    genes = read_gene(para['DataPath']+'/SGD/SGD_features.tab', i=[1,3,8,9,10])
    outfile = 'geneset_tRNA.txt'
    out = open(outfile, 'w')
    out.write('Chro\tStart\tEnd\tName\n')
    for Type, Name, Chro, Start, End in genes:
        if Chro == '' or Start == '':
            continue
        if Type != 'tRNA':
            continue
        out.write('%s\t%s\t%s\t%s\n'%(Chro, Start, End, Name))
    return outfile

def get_rRNA(para):
    genes = read_gene(para['DataPath']+'/SGD/SGD_features.tab', i=[1,3,8,9,10])
    outfile = 'geneset_rRNA.txt'
    out = open(outfile, 'w')
    out.write('Chro\tStart\tEnd\tName\n')
    for Type, Name, Chro, Start, End in genes:
        if Chro == '' or Start == '':
            continue
        if Type != 'rRNA':
            continue
        out.write('%s\t%s\t%s\t%s\n'%(Chro, Start, End, Name))
    return outfile

def get_expressed(para, time='0 min'):
    infile = open(para['DataPath']+'/Pramila2006GD/GSE4987_setA_family.pcl', 'r')
    ele = infile.readline().split('\t')
    idx = -1
    for i in xrange(len(ele)):
        if ele[i].find(time) >= 0:
            print ele[i]
            idx = i
            break
    if idx == -1: return ''
    select = set()
    infile.readline()
    for line in infile:
        ele = line.split('\t')
        gene = ele[0]
        score = float(ele[idx])
        if score > 0.2:
            select.add(gene)
    genes = read_gene(para['DataPath']+'/SGD/SGD_features.tab', i=[1,3,8,9,10])
    outfile = 'geneset_%s.txt'%time.replace(' ', '')
    out = open(outfile, 'w')
    out.write('Chro\tStart\tEnd\tName\n')
    for Type, Name, Chro, Start, End in genes:
        if Chro == '' or Start == '':
            continue
        if Name not in select:
            continue
        out.write('%s\t%s\t%s\t%s\n'%(Chro, Start, End, Name))
    return outfile

def get_paralogs(para):
    infile = open(para['DataPath']+'/Chursov2011B/Paralogs_RNA_Alignments.txt', 'r')
    select = set()
    for line in infile:
        if line.startswith('>'):
            select.add(line.replace('>','').strip())
    genes = read_gene(para['DataPath']+'/SGD/SGD_features.tab', i=[1,3,8,9,10])
    outfile = 'geneset_paralogs.txt'
    out = open(outfile, 'w')
    out.write('Chro\tStart\tEnd\tName\n')
    for Type, Name, Chro, Start, End in genes:
        if Chro == '' or Start == '':
            continue
        if Name not in select:
            continue
        out.write('%s\t%s\t%s\t%s\n'%(Chro, Start, End, Name))
    return outfile

def bin_count(infile, pdf):
    if not os.path.exists(infile): return
    from contact_map import ContactMap
    map1 = ContactMap()
    map1.load('YeastHiC')
    n = map1.contact_map.shape[0]
    idx, name = map1.get_locations(infile)
    srt, val, pval = map1.test_enrichment(idx, 'AUC', title=infile, plot=pdf, pages=9)
    show(infile)
    show(len(idx))
    show(len(set(idx)))
    show(pval[srt[0]])
    sign = [i for i in srt if pval[i] < 0.01]
    show(sign)
    show()
    if infile.find('telomere') > -1 or infile.find('tRNA') > -1 or infile.find('paralogs') > -1:
        outfile = open(infile+'_val.csv', 'w')
        outfile.write('Name,Bin Idx,Membership\n')
        for Name, Idx in zip(name, idx):
            outfile.write('%s,%s'%(Name, Idx))
            for i in sign:
                outfile.write(',%s'%map1.contact_group[Idx,i])
            outfile.write('\n')
        outfile.close()
    os.remove(infile)

def main(para):
    files = []
    files.append(get_centromere(para))
    files.append(get_telomere(para))
    files.append(get_tRNA(para))
    files.append(get_rRNA(para))
    files.append(get_paralogs(para))
    for t in xrange(0, 121, 5):
        files.append(get_expressed(para, '%s min  2001-08'%t))
    for t in xrange(0, 121, 5):
        if t == 0:
            files.append(get_expressed(para, '%s min  2001-10'%t))
            continue
        files.append(get_expressed(para, '%s min  2001-11'%t))
    show('Dataset\t#Bins\tUnique Bins\tMin. P-value\tEnriched Clusters\n')
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages(para['ExeFile']+'plot.pdf')
    for f in files:
        bin_count(f, pdf)
    pdf.close()

if __name__ == '__main__': main_fun(main)

