import argparse, sys, os, gzip
import cPickle as pickle
import gc

teflonBase = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.insert(1, teflonBase)

from teflon_scripts import mean_stats as ms
from teflon_scripts import genotyper_poolType as pt
from teflon_scripts import pseudo2refConvert as p2rC

def load_pickle(pFILE):
    tmpFILE=pFILE.replace(os.path.basename(pFILE),os.path.basename(pFILE)+".tmp")
    cmd="gunzip -c %s > %s" %(pFILE,tmpFILE)
    print "cdm:",cmd
    os.system(cmd)
    print "loading pickle:",tmpFILE
    print "NOTE: this step can be time and memory intesive for large reference genomes"
    inFILE = open(tmpFILE, "rb")
    gc.disable()
    pDICT=pickle.load(inFILE)
    gc.enable()
    inFILE.close()
    os.system("rm %s" %(tmpFILE))
    print "pickle loaded!"
    return pDICT

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-wd',dest='wd',help='full path to working directory', default=-1)
    parser.add_argument('-d',dest='DIR',help='full path to prep_TF directory')
    parser.add_argument('-s',dest='samples',help='tab delimited text file with full paths to indexed bamFILEs and sorted te positions')
    parser.add_argument('-f',dest='filt',help='sites genotyped as -9 if adjusted read counts greater than this threshold (default=mean_coverage + 4*STDEV)', type=int, default=-1)
    parser.add_argument('-dt',dest='dataType',help='haploid, diploid, or pooled')
    args = parser.parse_args()

    # identify current working directory
    if args.wd == -1:
        cwd=os.getcwd()
    else:
        cwd=os.path.realpath(args.wd)

    # import options
    hardFilt=args.filt
    prep_TF=args.DIR
    prefix=os.path.dirname(prep_TF).split("/")[-1].split(".prep_TF")[0]
    dataType=args.dataType
    if dataType not in "haploid, diploid, or pooled":
        return "Error datatype must be either haploid, diploid, or pooled"
        sys.exit()

    # read samples and stats
    samples=[]
    with open(args.samples, 'r') as fIN:
        for line in fIN:
            bamFILE=line.split()[0].replace(".bam",".subsmpl.bam")
            statsFile = bamFILE.replace(".bam", ".stats.txt")
            with open(statsFile, 'r') as fIN:
                for l in fIN:
                    if 'average length' in l:
                        readLen=int(float(l.split()[-1]))
                    if 'insert size average' in l:
                        insz=int(float(l.split()[-1]))
                    if 'insert size standard deviation' in l:
                        sd=int(float(l.split()[-1]))
            covFILE = bamFILE.replace(".bam", ".cov.txt")
            with open(covFILE, "r") as fIN:
                for l in fIN:
                    if l.startswith("Av"):
                        cov = float(l.split()[-1])
                    if l.startswith("St"):
                        cov_sd = float(l.split()[-1])
            samples.append([bamFILE, line.split()[1], [readLen, insz, sd, cov, cov_sd]])

    # average the stats for each sample
    stats=ms.mean_stats_portal(samples)

    # create the genotype directory
    genoDir = os.path.join(cwd,"finalPos")

    # define coverage threshold
    if hardFilt == -1:
        c_thresh=[]
        for sample in samples:
            c_thresh.append(int(sample[2][3]+ (4*sample[2][4])))
    else:
        c_thresh=[]
        for sample in samples:
            c_thresh.append(hardFilt)
    print "Coverage threshold filters (for each sample): %s" %(c_thresh)
    print "NOTE: all sites with adjusted read counts > coverage threshold will be marked -9"

    #load pseudo2ref.pickle
    pickleFILE=os.path.join(prep_TF,prefix+".pseudo2ref.pickle.gz")
    #posMap=pickle.load(gzip.open(pickleFILE, "rb"))
    posMap=load_pickle(pickleFILE)

    # genotype samples
    if dataType == "pooled":
        pt.pt_portal(genoDir,samples, posMap, stats, p2rC, c_thresh)
    else:
        print """coming soon...use "pooled" to obtain presence and absence read counts"""
    print "TEFLON GENOTYPE FINISHED!"

if __name__ == "__main__":
	main()
