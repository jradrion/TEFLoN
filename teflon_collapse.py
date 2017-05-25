import argparse, sys, os
import subprocess as sp
import shlex

teflonBase = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.insert(1, teflonBase)

from teflon_scripts import sort_positions as sortp
from teflon_scripts import collapse_union as cu
from teflon_scripts import mean_stats as ms
from teflon_scripts import subsample_alignments as sa


def mkdir_if_not_exist(*dirs):
    for dir in dirs:
        if not os.path.exists(dir):
            os.makedirs(dir)
            print "creating directory: %s" %(dir)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-wd',dest='wd',help='full path to working directory',default=-1)
    parser.add_argument('-s',dest='samples',help='tab delimited text file with full paths to indexed bamFILEs and sorted te positions')
    parser.add_argument('-es',dest='exeSAM',help='full path to samtools executable')
    parser.add_argument('-n',dest='thresh',help='TEs must be supported by >= n reads in at least one sample', type=int, default=1)
    parser.add_argument('-cov',dest='cov',help='subsample to coverage override', type=float, default=-1)
    parser.add_argument('-q',dest='qual',help='map quality threshold', type=int, default=20)
    parser.add_argument('-t',dest='nProc',help='number of processors', type=int, default=1)
    args = parser.parse_args()

    # identify current working directory
    if args.wd == -1:
        cwd=os.getcwd()
    else:
        cwd=os.path.realpath(args.wd)

    # import options
    exeSAM=args.exeSAM
    thresh=args.thresh
    qual=args.qual
    nProc=args.nProc
    covOverride=args.cov

    # read the samples and stats
    samples=[]
    # each sample will be formatted [path to bamFILE, uniqueID, [stats]]
    with open(args.samples, 'r') as fIN:
        for line in fIN:
            statsOutFile = line.split()[0].replace(".bam", ".stats.txt")
            with open(statsOutFile, 'r') as fIN:
                for l in fIN:
                    if "raw total sequences:" in l:
                        total_n=int(l.split()[-1])
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
            samples.append([line.split()[0], line.split()[1], [readLen, insz, sd, total_n,cov,cov_sd]])

    # generate subsampled alignments for use in teflon_count
    sa.subsample_alignments_portal(samples, exeSAM, nProc, qual, covOverride)

    # average the stats for each sample
    stats=ms.mean_stats_portal(samples)
    readLen,insz,sd=stats[0],stats[1],stats[2]

    # create the genotype directory
    genoDir = os.path.join(cwd,"initialPos")
    mkdir_if_not_exist(genoDir)

    # concatonate position estimates for each sample
    catFile = os.path.join(genoDir, "union.txt")
    try:
        files = ""
        for sample in samples:
            files += os.path.join(genoDir, sample[1]+".all_positions_sorted.txt" + " ")
        cmd = "cat %s" %(files)
        #print "cmd:", cmd  #p = sp.Popen(shlex.split(cmd), stdout=open(catFile, 'w'), stderr=sp.PIPE)
        p = sp.Popen(shlex.split(cmd), stdout=open(catFile, 'w'), stderr=sp.PIPE)
        perr = p.communicate()[1] # communicate returns a tuple (stdout, stderr)
        #print perr
        if p.returncode != 0:
            print "error concatenating positions"
            sys.exit(1)
    except OSError:
        print "Cannot concatenate positions"
        sys.exit(1)

    # sort positions
    print "Sorting positions"
    catFILE_sorted = sortp.sort_portal(catFile)

    # collapse union of all samples
    print "Collapse union of all samples"
    cu.collapse_union_portal(catFILE_sorted, readLen, insz, sd, thresh)
    print "TEFLON COLLAPSE FINISHED!"


if __name__ == "__main__":
	main()
