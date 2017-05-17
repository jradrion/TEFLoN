import argparse, sys, os, gzip
import cPickle as pickle


teflonBase = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.insert(1, teflonBase)

from teflon_scripts import mean_stats as ms
from teflon_scripts import genotyper_poolType as pt
from teflon_scripts import pseudo2refConvert as p2rC

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-wd',dest='wd',help='full path to working directory')
    parser.add_argument('-s',dest='samples',help='tab delimited text file with full paths to indexed bamFILEs and sorted te positions')
    parser.add_argument('-p',dest='pickle',help='pseudo2ref.pickle')
    parser.add_argument('-t',dest='dataType',help='haploid, diploid, or pooled')
    args = parser.parse_args()

    wd=os.path.realpath(args.wd)
    dataType=args.dataType
    if dataType not in "haploid, diploid, or pooled":
        return "Error datatype must be either haploid, diploid, or pooled"
        sys.exit()

    #load pseudo2ref.pickle
    print "loading",args.pickle,"..."
    posMap=pickle.load(gzip.open(args.pickle, "rb"))
    print args.pickle,"loaded"

    samples=[]
    #each sample will be [path to sorted_position.txt, path to bamFILE, uniqueID, stats]
    with open(args.samples, 'r') as fIN:
        for line in fIN:
            statsFile = line.split()[0].replace(".bam", ".stats.txt")
            with open(statsFile, 'r') as fIN:
                for l in fIN:
                    if 'average length' in l:
                        readLen=int(float(l.split()[-1]))
                    if 'insert size average' in l:
                        insz=int(float(l.split()[-1]))
                    if 'insert size standard deviation' in l:
                        sd=int(float(l.split()[-1]))
            covFILE = line.split()[0].replace(".bam", ".cov.txt")
            with open(covFILE, "r") as fIN:
                for l in fIN:
                    if l.startswith("Av"):
                        cov = float(l.split()[-1])
                    if l.startswith("St"):
                        cov_sd = float(l.split()[-1])
            samples.append([line.split()[0], line.split()[1], [readLen, insz, sd, cov, cov_sd]])


    #average the stats for each sample
    stats=ms.mean_stats_portal(samples)

    #create the genotype directory
    genoDir = os.path.join(wd,"finalPos")

    print "genotyping"
    if dataType == "pooled":
        pt.pt_portal(genoDir,samples, posMap, stats, p2rC)
    else:
        print """coming soon... Use "pooled" for temporary read counts"""
    print "finished genotyping"

if __name__ == "__main__":
	main()
