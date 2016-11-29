import argparse, sys, os
import multiprocessing as mp
import subprocess as sp
import shlex
import shutil

teflonBase = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.insert(1, teflonBase)

from teflon_scripts import sort_positions as sortp
from teflon_scripts import collapse_union as cu
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
            genoDir, tmpDir, exePATH, union, samples, qual, annotation = params
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
    parser.add_argument('-b',dest='bam',help='bam file for focal sample')
    parser.add_argument('-s',dest='samples',help='tab delimited text file of full path to bam files')
    parser.add_argument('-a',dest='annotation',help='full path to pseudo-transformed bedFILE')
    parser.add_argument('-t',dest='te_hierarchy',help='full path to te hierarchy file')
    parser.add_argument('-q',dest='qual',help='map quality threshold',type=int)
    parser.add_argument("-x", dest="nProc", type=int, default=4, help="Specify number of processes")
    args = parser.parse_args()

    wdPath=args.wd
    exePATH=args.exe
    qual=args.qual
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
            if line.split()[0] == args.bam:
                statsFILE = line.split()[0].replace(".bam", ".stats.txt")
                pre=line.split()[1]
                with open(statsFILE, 'r') as fIN:
                    for l in fIN:
                        if 'average length' in l:
                            readLen=int(float(l.split()[-1]))
                        if 'insert size average' in l:
                            insz=int(float(l.split()[-1]))
                        if 'insert size standard deviation' in l:
                            sd=int(float(l.split()[-1]))
                samples.append([line.split()[0], line.split()[1], [readLen, insz, sd]])



    #create the genotype directory
    genoDir = wdPath + "finalPos"
    tmpDir = os.path.join(wdPath,pre+".tmp")
    mkdir_if_not_exist(genoDir, tmpDir)

    union=[]
    with open(os.path.join(wdPath,"initialPos","union_sorted.collapsed.txt"), "r") as fIN:
        for line in fIN:
            union.append(line.split())


    #partition data for multiprocessing
    siteID = range(len(union))
    task_q = mp.JoinableQueue()
    result_q = mp.Queue()
    params=[genoDir, tmpDir, exePATH, union, samples, qual, annotation]
    #do the work boyeeee!
    print "counting reads..."
    create_procs(nProc, task_q, result_q, params)
    assign_task(siteID, task_q, nProc)
    try:
        task_q.join()
    except KeyboardInterrupt:
        print "KeyboardInterrupt"
        sys.exit(0)
    else:
        outCounts=[]
        while nProc:
            x = result_q.get()
            #print "result get",x
            outCounts.extend(y for y in x)
            nProc-=1
        print "finished counting reads"
    cts=sorted(outCounts)
    #print len(cts), len(union)
    #sys.exit()
    for i in range(len(samples)):
        outFILE=os.path.join(genoDir,samples[i][1]+".counts.txt")
        with open(outFILE, "w") as fOUT:
            for j in range(len(union)):
                fOUT.write("\t".join([str(x) for x in union[j]])+"\t%s\n" %("\t".join([str(y) for y in cts[j][1][0]])))

    #remove directory with beds/sams
    shutil.rmtree(tmpDir)







if __name__ == "__main__":
	main()
