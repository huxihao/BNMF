## This script trys to replicate the results in the paper:
##      Imakaev, M.; Fudenberg, G.; McCord, R. P.; Naumova, N.; Goloborodko,
##      A.; Lajoie, B. R.; Dekker, J. & Mirny, L. A. 
##      Iterative correction of Hi-C data reveals hallmarks of chromosome
##      organization. Nat Methods, ., 2012 
##
## Need to install all required packages as indicated by:
##      http://mirnylab.bitbucket.org/hiclib/index.html
##
## Author is Xihao <huxihao@gmail.com>
## Most of codes are adopted from the tutorial at:
##      http://mirnylab.bitbucket.org/hiclib/tutorial.html

from tools import *

def step1(hiclib_path, ## the path of hiclib folder on machine
          dataset='Kalhor2012NB', 
          sraid = 'SRR071231', 
          readlen = 40): ## each read with length 40
    ''' 1. Map reads to the genome
        http://mirnylab.bitbucket.org/hiclib/tutorial/01_iterative_mapping.html
    '''

    ## Adopted from hiclib tutorial
    import os
    import logging
    from hiclib import mapping
    from mirnylib import h5dict, genome

    logging.basicConfig(level=logging.DEBUG)

    # A. Map the reads iteratively.
    mapping.iterative_mapping(
        bowtie_path=hiclib_path+'/bin/bowtie2/bowtie2',
        bowtie_index_path=hiclib_path+'/bin/bowtie2/index/hg19',
        fastq_path='../data/SRA/'+dataset+'/'+sraid+'/'+sraid+'.sra',
        out_sam_path='../data/SRA/'+sraid+'_1.bam',
        min_seq_len=25,
        len_step=5,
        seq_start=0,
        seq_end=readlen,
        nthreads=12, # on intel corei7 CPUs 4 threads are as fast as
                     # 8, but leave some room for you other applications
        #max_reads_per_chunk = 10000000,  #optional, on low-memory machines
        temp_dir='../data/SRA/',  # optional, keep temporary files here
        bowtie_flags='--very-sensitive',
        bash_reader=hiclib_path+'/bin/sra/bin/fastq-dump -Z')

    mapping.iterative_mapping(
        bowtie_path=hiclib_path+'/bin/bowtie2/bowtie2',
        bowtie_index_path=hiclib_path+'/bin/bowtie2/index/hg19',
        fastq_path='../data/SRA/'+dataset+'/'+sraid+'/'+sraid+'.sra',
        out_sam_path='../data/SRA/'+sraid+'_2.bam',
        min_seq_len=25,
        len_step=5,
        seq_start=readlen,
        seq_end=2*readlen,
        nthreads=12,  
        #max_reads_per_chunk = 10000000, 
        temp_dir='../data/SRA/',  
        bowtie_flags='--very-sensitive',
        bash_reader=hiclib_path+'/bin/sra/bin/fastq-dump -Z')

    # B. Parse the mapped sequences into a Python data structure,
    #    assign the ultra-sonic fragments to restriction fragments.
    mapped_reads = h5dict.h5dict(sraid + '_mapped_reads.hdf5') ## to local folder
    genome_db    = genome.Genome(hiclib_path+'/fasta/hg19', readChrms=['#', 'X'])

    mapping.parse_sam(
        sam_basename1='../data/SRA/'+sraid+'_1.bam',
        sam_basename2='../data/SRA/'+sraid+'_2.bam',
        out_dict=mapped_reads,
        genome_db=genome_db, 
        enzyme_name='HindIII')

def step2(hiclib_path, sraid, res=1000000):
    ''' 2. Filter the dataset at the restriction fragment level.
        http://mirnylab.bitbucket.org/hiclib/tutorial/02_fragment_filtering.html
    '''
    from mirnylib import genome
    from hiclib import fragmentHiC 

    # Create a HiCdataset object.
    genome_db = genome.Genome(hiclib_path+'/fasta/hg19', readChrms=['#', 'X'])
    fragments = fragmentHiC.HiCdataset(
        filename=sraid+'_fragment_dataset.hdf5',
        genome=genome_db,
        maximumMoleculeLength=500,
        enzymeName='HindIII',
        mode='w')

    # Load the parsed reads into the HiCdataset. The dangling-end filter is applied
    # at this stage, with maximumMoleculeLength specified at the initiation of the 
    # object.
    fragments.parseInputData(dictLike=sraid+'_mapped_reads.hdf5')

    fragments.filterRsiteStart(offset=5)
    fragments.filterDuplicates()
    fragments.filterLarge()
    if sraid in ["SRR071231", "SRR071232"]: ## set to 0.1% for TCC
        fragments.filterExtreme(cutH=0.001, cutL=0) 
    else: ## default for Hi-C is 0.5%
        fragments.filterExtreme(cutH=0.005, cutL=0) 

#    fragments.saveFragments()
    fragments.saveHeatmap(sraid+'_map-res%sk.hdf5'%(res/1000), resolution=res)

def step3(hiclib_path, sraid, res=1000000):
    ''' 3. Filter and iteratively correct heatmaps.
        http://mirnylab.bitbucket.org/hiclib/tutorial/03_heatmap_processing.html
    '''
    import matplotlib.pyplot as plt
    import numpy as np

    from mirnylib import genome
    from mirnylib import h5dict
    from mirnylib import plotting
    from hiclib import binnedData

    genome_db = genome.Genome(hiclib_path+'/fasta/hg19', readChrms=['#', 'X'])

    # Read resolution from the dataset.
    raw_heatmap = h5dict.h5dict(sraid+'_map-res%sk.hdf5'%(res/1000), mode='r') 
    resolution = int(raw_heatmap['resolution'])

    # Create a binnedData object, load the data.
    BD = binnedData.binnedData(resolution, genome_db)
    BD.simpleLoad(sraid+'_map-res%sk.hdf5'%(res/1000), 'DataName')

    # Plot the heatmap directly.
    plotting.plot_matrix(np.log(BD.dataDict['DataName']))
    plt.savefig(sraid+'_map-res%sk.pdf'%(res/1000))
    plt.clf()

    # Remove the contacts between loci located within the same bin.
    BD.removeDiagonal()

    # Remove bins with less than half of a bin sequenced.
    BD.removeBySequencedCount(0.5)

    # Remove 1% of regions with low coverage.
    BD.removePoorRegions(cutoff=1)

    # Truncate top 0.05% of interchromosomal counts (possibly, PCR blowouts).
    BD.truncTrans(high=0.0005)

    # Perform iterative correction.
    BD.iterativeCorrectWithoutSS()

    # Save the iteratively corrected heatmap.
    BD.export('DataName', sraid+'_map-res%sk-ic.hdf5'%(res/1000))

    # Plot the heatmap directly.
    plotting.plot_matrix(np.log(BD.dataDict['DataName']))
    plt.savefig(sraid+'_map-res%sk-ic.pdf'%(res/1000))
    plt.clf()

    # Save Bias
    outfile = open(sraid+"_map-res%sk-ic-bias.txt"%(res/1000), "w")
    for i in xrange(len(BD.chromosomeIndex)):
        chro = BD.genome.idx2label[BD.chromosomeIndex[i]]
        posi = BD.positionIndex[i]
        outfile.write("chr%s\t%s\t%s"%(chro, posi, posi+res))
        outfile.write("\t%s"%BD.biasDict['DataName'][i])
        outfile.write("\n")
    outfile.close()

def step4(hiclib_path, sraid, res=1000000):
    ''' 4. Eigen vector decomposition
    /examples/iterativeCorrectionEigenvectorExpansion/eigenvectorAnalysis.py
    '''
    import matplotlib.pyplot as plt
    import numpy as np
    from mirnylib import genome
    from mirnylib import h5dict
    from mirnylib import plotting
    from hiclib import binnedData

    genome_db = genome.Genome(hiclib_path+'/fasta/hg19', readChrms=['#', 'X'])

    # Read resolution from the dataset.
    raw_heatmap = h5dict.h5dict(sraid+'_map-res%sk.hdf5'%(res/1000), mode='r')  
    resolution = int(raw_heatmap['resolution'])

    # Create a binnedData object, load the data.
    BD = binnedData.binnedData(resolution, genome_db)
    BD.simpleLoad(sraid+'_map-res%sk.hdf5'%(res/1000), 'DataName')
    
    # Do eigen decomposition
    BD.removeDiagonal()
    BD.removeBySequencedCount(0.5)
    BD.removeCis()
    BD.truncTrans(high=0.0005)
    BD.removePoorRegions(cutoff=1)
    BD.fakeCis()
    BD.removeZeros()
    BD.doEig(numPCs=30, force=True) ## First 30 EIGs
    BD.restoreZeros(value=0)

    eig = BD.eigEigenvalueDict['DataName']
    eig_v = BD.EigDict['DataName']

    # Plot the heatmap directly.
    plotting.plot_matrix(np.log(np.dot(np.dot(eig_v.T, np.diag(eig)), eig_v)))
    plt.savefig(sraid+'_map-res%sk-eig.pdf'%(res/1000))
    plt.clf()

    outfile = open(sraid+"_map-res%sk-ic-eig.txt"%(res/1000), "w")
    for i in xrange(len(BD.chromosomeIndex)):
        chro = BD.genome.idx2label[BD.chromosomeIndex[i]]
        posi = BD.positionIndex[i]
        outfile.write("chr%s\t%s\t%s"%(chro, posi, posi+res))
        for eigenvector in eig_v:
            outfile.write("\t%s"%eigenvector[i])
        outfile.write("\n")
    outfile.close()
        
def tcc(hiclib_path, binsize=1000000):
    for sraid in ["SRR071233", "SRR071234", "SRR071231", "SRR071232"]:
        show(sraid+'\n')
        if not os.path.exists(sraid + '_mapped_reads.hdf5'): ## to save time
            step1(hiclib_path=hiclib_path, dataset='Kalhor2012NB', sraid=sraid, readlen=40)
        step2(hiclib_path=hiclib_path, sraid=sraid, res=binsize)
        step3(hiclib_path=hiclib_path, sraid=sraid, res=binsize)
        step4(hiclib_path=hiclib_path, sraid=sraid, res=binsize)

def hic(hiclib_path, binsize=1000000):
    for sraid in ["SRR027956", "SRR027957", "SRR027958", "SRR027959", "SRR027962", "SRR027963"]:
        show(sraid+'\n')
        if not os.path.exists(sraid + '_mapped_reads.hdf5'): ## to save time
            step1(hiclib_path=hiclib_path, dataset='Lieberman-Aiden2009S', sraid=sraid, readlen=76)
        step2(hiclib_path=hiclib_path, sraid=sraid, res=binsize)
        step3(hiclib_path=hiclib_path, sraid=sraid, res=binsize)
        step4(hiclib_path=hiclib_path, sraid=sraid, res=binsize)

def ren(hiclib_path, binsize=1000000):
    for sraid in ['SRR400260','SRR400261','SRR400262','SRR400263',
        'SRR400264','SRR400265','SRR400266','SRR400267','SRR400268',
        'SRR442155','SRR442156','SRR442157',
        'SRR442158','SRR442159','SRR442160']:
        show(sraid+'\n')
        if not os.path.exists(sraid + '_mapped_reads.hdf5'): ## to save time
            if sraid == 'SRR400260': ## only one with 100 read length
                step1(hiclib_path=hiclib_path, dataset='Dixon2012N', sraid=sraid, readlen=100)
            else: ## all remainings have 36 read length
                step1(hiclib_path=hiclib_path, dataset='Dixon2012N', sraid=sraid, readlen=36)
        step2(hiclib_path=hiclib_path, sraid=sraid, res=binsize)
        step3(hiclib_path=hiclib_path, sraid=sraid, res=binsize)
        step4(hiclib_path=hiclib_path, sraid=sraid, res=binsize)

def main(para):
    reso = int(1e6)
    hiclib_path = "../../hiclib/hiclib"
    hic(hiclib_path, reso)
    tcc(hiclib_path, reso)
    ren(hiclib_path, reso)

if __name__ == "__main__": main_fun(main)
