import argparse, sys, os
import multiprocessing as mp
import subprocess as sp
import shlex

teflonBase = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.insert(1, teflonBase)

from teflon_scripts import mean_stats as ms
from teflon_scripts import genotyper_writeBed as wb
from teflon_scripts import genotyper_countReads as cr
from teflon_scripts import genotyper_poolType as pt

def check_samtools(exePATH):
    try:
        cmd = "%s" %(exePATH)
        p = sp.Popen(shlex.split(cmd), stdout=sp.PIPE, stderr=sp.PIPE)
    except OSError:
        print "Cannot find samtools"
        sys.exit(1)

def mkdir_if_not_exist(*dirs):
    for dir in dirs:
        if not os.path.exists(dir):
            os.makedirs(dir)
            print "creating directory: %s" %(dir)

def assign_task(siteID, task_q, nProcs):
    c,i,nth_job=0,0,1
    while (i+1)*nProcs <= len(siteID):
        i+=1
    nP1=nProcs-(len(siteID)%nProcs)
    for j in range(nP1):
        #print "put", len(siteID[c:c+i]), nth_job
        task_q.put((siteID[c:c+i], nth_job))
        nth_job += 1
        c=c+i
    for j in range(nProcs-nP1):
        #print "put", len(siteID[c:c+i+1]), nth_job
        task_q.put((siteID[c:c+i+1], nth_job))
        nth_job += 1
        c=c+i+1

def create_procs(nProcs, task_q, result_q, params):
    for _ in range(nProcs):
        p = mp.Process(target=worker, args=(task_q, result_q, params))
        p.daemon = True
        p.start()

def worker(task_q, result_q, params):
    while True:
        try:
            siteID, nth_job = task_q.get()
            #unpack parameters
            genoDir, tmpDir, exePATH, union, samples, qual, annotation, dataType = params
            tmp=[]
            for ID in siteID:
                #print union[ID]
                # 1. write sam files for all samples
                wb.wb_portal(ID, union, samples, tmpDir, exePATH, qual)
                # 2. count support and non support reads for each event
                #print "counting support reads: siteID", ID
                counts = cr.countReads_portal(ID, union, samples, tmpDir, annotation) # counts = [[#p(sample1),#a(sample1),#NA(sample1)],[#p(sampleN),#a(sampleN),#NA(sampleN)]]
                print ID, counts, union[ID]
                tmp.append([ID,counts])
            #print "result put",tmp
            result_q.put(tmp)
        finally:
            task_q.task_done()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-wd',dest='wd',help='full path to working directory')
    parser.add_argument('-ex',dest='exe',help='full path to samtools executable')
    parser.add_argument('-g',dest='genome',help='genomeSize file')
    parser.add_argument('-s',dest='samples',help='tab delimited text file with full paths to indexed bamFILEs and sorted te positions')
    parser.add_argument('-a',dest='annotation',help='full path to pseudo-transformed bedFILE')
    parser.add_argument('-t',dest='te_hierarchy',help='full path to te hierarchy file')
    parser.add_argument('-pm',dest='posMap',help='map of positions from pseudo ref to ref')
    parser.add_argument('-dt',dest='dataType',help='haploid, diploid, or pooled')
    parser.add_argument("-x", dest="nProc", type=int, default=4, help="Specify number of processes")
    args = parser.parse_args()

    wd=os.path.realpath(args.wd)
    exePATH=args.exe
    posMap=args.posMap
    dataType=args.dataType
    if dataType not in "haploid, diploid, or pooled":
        return "Error datatype must be either haploid, diploid, or pooled"
        sys.exit()
    nProc=args.nProc

    #check samtools executabe for function
    check_samtools(exePATH)

    #groups=[id,family,superfamily,suborder,order,class]
    hierarchy,label={},[]
    with open(args.te_hierarchy, 'r') as fIN:
        for line in fIN:
            if line.startswith('id'):
                label=line.split()[1:]
            if not line.startswith('id'):
                hierarchy[line.split()[0]] = line.split()[1:]

    annotation=[]
    with open(args.annotation, 'r') as fIN:
        for line in fIN:
            arr=line.split()
            annotation.append([arr[0],int(arr[1]),int(arr[2]),arr[3]])

    chromosomes,lengths=[],[]
    with open(args.genome, 'r') as fIN:
        for line in fIN:
            arr=line.split()
            chromosomes.append(arr[0])
            lengths.append(int(arr[2]))

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
        #print outCounts
        pt.pt_portal(genoDir,samples, posMap, stats)
    else:
        print "coming soon... Use "pooled" for temporary read counts"
    print "finished genotyping"

if __name__ == "__main__":
	main()
