import numpy as np
from tools import *
from contact_map import ContactMap
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def plot1(para):
    pdf = PdfPages(para['ExeFile']+'plot1.pdf')
    ## initalization
    map1 = ContactMap('plot1')
    map1.clear()
    ## read chromosome sizes
    if not map1.load():
        map1.genome_info('../data/yeast_chr_len.txt')
        datafiles = ['../data/Duan2010N/interactions_HindIII_fdr0.01_inter.txt',
                     '../data/Duan2010N/interactions_HindIII_fdr0.01_intra.txt'] 
        for datafile in datafiles:
            map1.add_interactions(datafile)
        map1.create_binnedmap(binsize=20e3)
        map1.mask_diag()
        map1.mask_short()
        map1.mask_low()

    map1.decompose_auto(plot=pdf)
    map1.sort_groups()
    map1.save()

    map1.plot_map(vmin=1, vmax=1000, title='$X$')
    pdf.savefig(); plt.clf();
    map1.plot_map(np.diag(map1.bias_vector), log=False, title='$B$')
    pdf.savefig(); plt.clf();
    map1.plot_map(map1.contact_group, log=False, title='$H$')
    pdf.savefig(); plt.clf();
    map1.plot_map(map1.group_map, log=False, title='$S$')
    pdf.savefig(); plt.clf();
    map1.plot_map(map1.group_map * map1.contact_group.T, log=False, title='$W=SH^T$')
    pdf.savefig(); plt.clf();
    map1.plot_map(map1.contact_group * map1.group_map * map1.contact_group.T, vmin=1, vmax=1000, title='$R=HSH^T$')
    pdf.savefig(); plt.clf();
    grps = map1.label_groups(plot=pdf)
    r = map1.contact_group.shape[1]
    for i in [0,r-2,r-1]:
        map1.plot_map(map1.contact_group[:,i] * map1.contact_group[:,i].T, vmin=1, title=str(i+1))
        pdf.savefig(); plt.clf();
    map1.plot_map(np.outer(map1.bias_vector, map1.bias_vector), log=False)
    pdf.savefig(); plt.clf();
    map1.add_bias_back()
    map1.plot_map(map1.contact_group * map1.group_map * map1.contact_group.T, vmin=1, vmax=1000, title='$R=HSH^T$')
    pdf.savefig(); plt.clf();
    pdf.close()

def plot2(para):
    pdf = PdfPages(para['ExeFile']+'plot2.pdf')
    ## initalization
    map1 = ContactMap('plot2')
    if True:
        map1.genome_info('../data/yeast_chr_len.txt')
        datafiles = ['../data/Duan2010N/interactions_HindIII_fdr0.01_inter.txt',
                     '../data/Duan2010N/interactions_HindIII_fdr0.01_intra.txt'] 
        for datafile in datafiles:
            map1.add_interactions(datafile)
        map1.create_binnedmap(binsize=10e3)
        map1.mask_diag()
        map1.mask_short()
        map1.mask_low()
    map1.plot_map(map1.contact_map, log=True, vmin=1, vmax=100)
    pdf.savefig(); plt.clf();
    sel = np.arange(200,400)
    map1.plot_map(map1.contact_map[sel,:][:,sel], log=True, vmin=1, vmax=100)
    pdf.savefig(); plt.clf();
    for l in [0, 0.1, 1, 10]:
        map1.reset_solution()
        map1.decompose('NMF-PoissonManifoldEqual', dim_num=55, par_lam=l)
        R = map1.contact_group * map1.group_map * map1.contact_group.T
        map1.plot_map(R[sel,:][:,sel], vmin=1, vmax=100, title=str(l))
        pdf.savefig(); plt.clf();
    pdf.close()

def main(para):
    plot1(para)
    plot2(para)

if __name__ == "__main__": main_fun(main)
