from tools import *
from contact_map import ContactMap
import gzip
import numpy as np

class Armatus(ContactMap):
    def __init__(self, exepath, name='DataName', enzyme='Enzyme', sparse=False):
        ContactMap.__init__(self, name, enzyme, sparse)
        path = os.path.abspath(exepath)
        if not os.path.exists(path):
            print 'Please install Armatus to', path
            print 'The project is at https://github.com/kingsfordgroup/armatus'
            exit(0)
        self.exepath = path

    def decompose(self):
        ## chr by chr
        r = self.frag_end[0] - self.frag_sta[0] ## resolution
        tad = []
        for c in self.chr2idx:
            ## save the contact map in the format
            mrk = (self.frag_chr == self.chr2idx[c])
            dat = np.nan_to_num(self.contact_map[mrk,:][:,mrk])
            out = gzip.open('armatus_chr.gz', 'wb')
            for i in xrange(dat.shape[0]):
                for j in xrange(dat.shape[1]):
                    out.write(str(dat[i,j])+'\t')
                out.write('\n')
            out.close()
            ## run the program
            os.system(self.exepath+' -g .5 -m -n 20 -c '+c+' -r '+str(r)+' -i armatus_chr.gz -o armatus_out')
            ## read domains
            p1, p1v = self.get_locations('armatus_out.consensus.txt', st=0, ch=0, po=1, nm=1, add=1, skip=False)
            p2, p2v = self.get_locations('armatus_out.consensus.txt', st=0, ch=0, po=2, nm=2, add=-r, skip=False)
            for i,j,I,J in zip(p1, p2, p1v, p2v):
                if i > j or i <0 or j<0:
                    print '!! Some problem in reading the domain', I,J, i,j
                    exit(0)
            tad += zip(p1,p2)
        self.contact_group = np.zeros((self.contact_map.shape[0], len(tad)+1))
        for k in xrange(len(tad)):
            p1,p2 = tad[k]
            for i in xrange(p1,p2+1):
                self.contact_group[i,k] = 1
        self.contact_group[self.contact_group.sum(1)==0,-1] = 2 ## remaining ones
        return tad 

def main(para):
    map1 = Armatus(para['WorkPath']+'/../tools/armatus2.1/armatus')
    map1.genome_info(para['DataPath']+'/yeast_chr_len.txt')
    map1.add_interactions(para['DataPath']+'/Duan2010N/interactions_HindIII_fdr0.01_intra.txt')
    map1.create_binnedmap(binsize=4000)
    map1.decompose()
    print (map1.contact_group.sum(1) >= 1).sum(), map1.contact_group.shape

if __name__ == '__main__': main_fun(main)

