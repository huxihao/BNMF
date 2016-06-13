from tools import *

def main(para):
    i = 0
    methods = ['PCA', 'ICE', 'K-means', 'NMF', 'BNMF']
    show('#Clusters')
    show([m+'\tintra\tinter' for m in methods], True)
    for r in range(10,201,10):
        show(r)
        for method in methods:
            i += 1
            if not os.path.exists('temp_par.log%s'%i):
                os.system('mkdir temp_par'+str(i))
                os.system('cp syn_link.npy temp_par'+str(i))
                os.system('cp syn_dist.npy temp_par'+str(i))
                os.system('cp pdb.txt temp_par'+str(i))
                os.system('python ../src/compare_cluster.py DataPath=%s Method=%s Cluster=%s WorkPath=temp_par%s LogFile=../temp_par.log%s &'%(para['DataPath'], method, r, i, i))
            else:
                if os.path.exists('temp_par'+str(i)):
                    os.system('rm -rf temp_par'+str(i))
                with open('temp_par.log%s'%i, 'r') as inp:
                    show(inp.readline().strip())
        show()

if __name__ == '__main__': main_fun(main)
