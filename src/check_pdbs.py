from tools import *
from contact_map import ContactMap

def run1(para):
    show('''
    Decompose headmap and show clusters in PDB format
    ''', True)
    path = para['DataPath']+'/Duan2010N'
    map1 = ContactMap('PDBMAP')
    map1.genome_info(path+'/restriction_fragments_mappable_HindIII.txt', i2=3, i3=0)
    map1.add_interactions(path+'/interactions_HindIII_fdr0.01_inter.txt')
    map1.add_interactions(path+'/interactions_HindIII_fdr0.01_intra.txt')
    map1.create_binnedmap(10000)
    map1.decompose_auto(dims=range(5,51,5))
    map_idx, pdb_idx = map1.get_locations(path+'/3d_model_of_yeast_with_genomic_positions.txt', st=1, ch=0, po=1, nm=-1)
    H = map1.contact_group
    n,r = H.shape
    for i in xrange(3):
        members = set()
        for j in xrange(n):
            if H[j,i] > 1:
                members.add(j)
        mark_idx = [ip for im,ip in zip(map_idx, pdb_idx) if im in members]
        output_pdb('Yeast3D-C%s.pdb'%i, path+'/3d_model_of_yeast_genome.pdb', mark_idx)

def output_pdb(out_pdb, ref_pdb, mark_idx):
    outfile = open(out_pdb, 'w')
    infile = open(ref_pdb, 'r')
    cc = 1
    for line in infile:
        if cc in mark_idx:
            line = line[:12]+'C'+line[13:]
        outfile.write(line)
        cc += 1
    infile.close()
    outfile.close()

def run2(para, maxi=50):
    show('''
    Combine simulated structure population
    ''', True)
    import numpy as np
    import matplotlib.pyplot as plt
    import scipy.spatial.distance as dist
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages(para['ExeFile']+'plot.pdf')
    path = para['DataPath']+'/Tjong2012GR/simulation/'
    pdbfile = os.listdir(path)
    cc = 0
    avg_pos = 0
    avg_map = 0
    ## REF: http://cupnet.net/pdb_format/
    for f in pdbfile:
        cc += 1
        print cc, f
        pdb = open(path+f, 'r')
        ref = []
        pp = []
        for line in pdb:
            ref.append(line)
            pp.append([float(a) for a in line[30:54].split()])
        pdb.close()
        pos1 = np.array(pp)
        map1 = 2000-dist.squareform(dist.pdist(pos1, 'euclidean'))
        avg_pos += pos1
        avg_map += map1
        if cc < 5:
            plt.imshow(map1[::5,::5])
            plt.colorbar()
            pdf.savefig(); plt.clf()
        if cc == maxi: break
    avg_pos /= cc
    avg_map /= cc
    plt.imshow(avg_map[::5,::5])
    plt.colorbar()
    pdf.savefig(); plt.clf()
    outfile = open('average.pdb', 'w')
    pos = avg_pos.tolist()
    for line,xyz in zip(ref,pos):
        line = line[:30] + '%8.3f%8.3f%8.3f'%tuple(xyz) + line[54:]
        outfile.write(line)
    outfile.close()
    pdf.close()

def main(para):
    run1(para)
    run2(para)

if __name__ == '__main__': main_fun(main)

