from tools import *
from contact_map import ContactMap
import numpy as np

class DomainCall(ContactMap):
    def __init__(self, exepath, name='DataName', enzyme='Enzyme', sparse=False):
        ContactMap.__init__(self, name, enzyme, sparse)
        path = os.path.abspath(exepath)
        if not os.path.exists(path):
            print 'Please install domaincall to', path
            print 'The project is at https://github.com/kingsfordgroup/armatus'
            exit(0)
        self.exepath = path

    def decompose(self):
        oldpath = os.path.abspath('.')
        os.chdir(self.exepath)
        ## save chrom size
        with open('genome_size.txt', 'w') as tmpfile:
            for c in self.idx2chr:
                m = (self.frag_end[self.frag_chr==c])[-1]
                tmpfile.write('chr%s\t%s\t0\t0\t0\n'%(c, m))
        ## chr by chr
        r = self.frag_end[0] - self.frag_sta[0] ## resolution
        tad = []
        os.system('> DI_file.txt') ## whole genome DI
        for c in self.idx2chr:
            ## save the contact map in the format
            mrk = (self.frag_chr == c)
            dat = np.nan_to_num(self.contact_map[mrk,:][:,mrk])
            out = open('matrix_file.txt', 'w')
            for i in xrange(dat.shape[0]):
                out.write('chr%s\t%s\t%s'%(c, i*r, i*r+r))
                for j in xrange(dat.shape[1]):
                    out.write('\t'+str(dat[i,j]))
                out.write('\n')
            out.close()
            ## run the program
            os.system('perl perl_scripts/DI_from_matrix.pl matrix_file.txt %s %s genome_size.txt >> DI_file.txt'%(r, r*50))
        ## run HMM
        os.system('''matlab -nojvm -nodisplay -nodesktop -nosplash -r "run('HMM_calls');exit;" ''')
        os.system('perl perl_scripts/file_ends_cleaner.pl hmm_output.txt DI_file.txt | perl perl_scripts/converter_7col.pl > hmm_7colfile')
        for c in self.idx2chr:
            os.system('cat hmm_7colfile | grep "^chr%s\t" > hmm_7colfile_chr'%c)
            os.system('perl perl_scripts/hmm_probablity_correcter.pl hmm_7colfile_chr 2 0.99 '+str(r)+' | perl perl_scripts/hmm-state_caller.pl genome_size.txt chr'+str(c)+' | perl perl_scripts/hmm-state_domains.pl > tad_domains.txt')
            ## read domains
            p1, p1v = self.get_locations('tad_domains.txt', st=0, ch=0, po=1, nm=1, add=1, skip=False)
            p2, p2v = self.get_locations('tad_domains.txt', st=0, ch=0, po=2, nm=2, add=-r, skip=False)
            for i,j,I,J in zip(p1, p2, p1v, p2v):
                if i > j or i <0 or j<0:
                    print '!! Some problem in reading the domain', I,J, i,j
                    exit(0)
            tad += zip(p1,p2)
        os.chdir(oldpath)
        self.contact_group = np.zeros((self.contact_map.shape[0], len(tad)+1))
        for k in xrange(len(tad)):
            p1,p2 = tad[k]
            for i in xrange(p1,p2+1):
                self.contact_group[i,k] = 1
        self.contact_group[self.contact_group.sum(1)==0,-1] = 2 ## remaining ones
        return tad 

def main(para):
    map1 = DomainCall(para['WorkPath']+'/../tools/domaincall/')
    map1.genome_info(para['DataPath']+'/yeast_chr_len.txt')
    map1.add_interactions(para['DataPath']+'/Duan2010N/interactions_HindIII_fdr0.01_intra.txt')
    map1.create_binnedmap(binsize=4000)
    map1.decompose()
    print (map1.contact_group.sum(1) >= 1).sum(), map1.contact_group.shape

if __name__ == '__main__': main_fun(main)

