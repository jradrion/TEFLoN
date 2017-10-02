import argparse, sys, os
import multiprocessing as mp
import subprocess as sp
import shlex
import shutil

teflonBase = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.insert(1, teflonBase)

from teflon_scripts import genotyper_writeBed as wb
from teflon_scripts import genotyper_countReads as cr
from teflon_scripts import reduceSearchSpace as rss

def drawProgressBar(percent, barLen = 50):
    sys.stdout.write("\r")
    progress = ""
    for i in range(barLen):
        if i < int(barLen * percent):
            progress += "="
        else:
            progress += " "
    sys.stdout.write("[ %s ] %.2f%%" % (progress, percent * 100))
    sys.stdout.flush()

def check_dependency(exePATH):
    try:
        cmd = "%s" %(exePATH)
        sp.Popen(shlex.split(cmd), stdout=sp.PIPE, stderr=sp.PIPE)
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
            genoDir, tmpDir, exePATH, union, samples, qual, hierarchy, label, l2, annotation = params
            tmp=[]
            ct=1
            for ID in siteID:
                flag=0
                #print union[ID]
                # 1. write sam files for all samples
                samFILE = wb.wb_portal(ID, union, samples, tmpDir, exePATH, qual, "firstPass")
                # 2. count support and non support reads for each event
                counts = cr.countReads_portal(ID, union, samples, tmpDir, samFILE, hierarchy, label, l2, annotation) # counts = [[#p(sample1),#a(sample1),#NA(sample1)],[#p(sampleN),#a(sampleN),#NA(sampleN)]]
                # 3. were new clipped positions identified?
                if counts[1]!="-":
                    union[ID][1]=counts[1]
                    union[ID][7]="+"
                    flag=1
                if counts[2]!="-":
                    union[ID][2]=counts[2]
                    union[ID][8]="+"
                    flag=1
                if flag==1: #position estimates were refined, recount
                    # 1. write sam files for all samples
                    samFILE = wb.wb_portal(ID, union, samples, tmpDir, exePATH, qual, "refine")
                    # 2. count support and non support reads for each event
                    counts = cr.countReads_portal(ID, union, samples, tmpDir, samFILE, hierarchy, label, l2, annotation) # counts = [[#p(sample1),#a(sample1),#NA(sample1)],[#p(sampleN),#a(sampleN),#NA(sampleN)]]
                    union[ID].extend(ct for ct in counts[0])
                    tmp.append([ID,union[ID]])
                    if nth_job==1:
                        drawProgressBar(ct/float(len(siteID)))
                else:
                    union[ID].extend(ct for ct in counts[0])
                    tmp.append([ID,union[ID]])
                    if nth_job==1:
                        drawProgressBar(ct/float(len(siteID)))
                ct+=1
            result_q.put(tmp)
        finally:
            task_q.task_done()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-wd',dest='wd',help='full path to working directory',default=-1)
    parser.add_argument('-d',dest='DIR',help='full path to prep_TF directory')
    parser.add_argument('-s',dest='samples',help='samples file')
    parser.add_argument('-i',dest='ID',help='unique id of this sample')
    parser.add_argument('-eb',dest='exeBWA',help='full path to bwa executable')
    parser.add_argument('-es',dest='exeSAM',help='full path to samtools executable')
    parser.add_argument('-l2',dest='cLevel',help='level of hierarchy to cluster')
    parser.add_argument('-q',dest='qual',help='map quality threshold',type=int)
    parser.add_argument("-t", dest="nProc", type=int, default=4, help="Specify number of processes")
    args = parser.parse_args()

    # identify current working directory
    if args.wd == -1:
        cwd=os.getcwd()
    else:
        cwd=os.path.realpath(args.wd)

    # import options
    prep_TF=os.path.realpath(args.DIR)
    #prefix=os.path.dirname(prep_TF).split("/")[-1].split(".prep_TF")[0]
    prefix=os.path.basename(args.DIR).replace(".prep_TF","")
    exeSAM=args.exeSAM
    exeBWA=args.exeBWA
    qual=args.qual
    nProc=args.nProc

    #check dependencies for function
    check_dependency(exeSAM)
    check_dependency(exeBWA)

    # import hierarchy
    l2=args.cLevel
    hierFILE=os.path.join(prep_TF,prefix+".hier")
    hierarchy,label={},[]
    ct=0
    with open(hierFILE, 'r') as fIN:
        for line in fIN:
            if ct==0:
                label=line.split()[1:]
            else:
                hierarchy[line.split()[0]] = line.split()[1:]
            ct+=1

    # read annotation
    annotation={}
    with open(os.path.join(prep_TF,prefix+".te.pseudo.bed"), 'r') as fIN:
        for line in fIN:
            arr=line.split()
            annotation[arr[3]]=[arr[0],int(arr[1]),int(arr[2])]

    # read chromosome lengths
    chromosomes,lengths=[],[]
    genomeSizeFILE=os.path.join(prep_TF,prefix+".genomeSize.txt")
    with open(genomeSizeFILE, 'r') as fIN:
        for line in fIN:
            arr=line.split()
            chromosomes.append(arr[0])
            lengths.append(int(arr[2]))

    # read samples and stats
    samples=[]
    with open(os.path.realpath(args.samples), 'r') as fIN:
        for line in fIN:
            bamFILE = line.split()[0].replace(".bam",".subsmpl.bam")
            #print "#",line.split()[1],args.ID
            if line.split()[1] == args.ID:
                #print "hah"
                statsFILE = bamFILE.replace(".bam", ".stats.txt")
                #print statsFILE
                #sys.exit()
                pre=line.split()[1]
                with open(statsFILE, 'r') as fIN:
                    for l in fIN:
                        if 'average length' in l:
                            readLen=int(float(l.split()[-1]))
                        if 'insert size average' in l:
                            insz=int(float(l.split()[-1]))
                        if 'insert size standard deviation' in l:
                            sd=int(float(l.split()[-1]))
                samples.append([bamFILE, pre, [readLen, insz, sd]])

    # create the genotype directory
    genoDir = os.path.join(cwd,"finalPos")
    tmpDir = os.path.join(cwd,pre+".tmp")
    mkdir_if_not_exist(genoDir, tmpDir)

    # read positions to search
    union=[]
    with open(os.path.join(cwd,"initialPos","union_sorted.collapsed.txt"), "r") as fIN:
        for line in fIN:
            union.append(line.split())

    # reduce search space
    samples[0][0]=rss.reduceSearchSpace_portal(union,samples,tmpDir,exeSAM,qual,nProc)
    # index new reduced alignment
    cmd="%s index %s" %(exeBWA,samples[0][0])
    print "cmd:", cmd
    os.system(cmd)

    # partition data for multiprocessing
    siteID = range(len(union))
    task_q = mp.JoinableQueue()
    result_q = mp.Queue()
    params=[genoDir, tmpDir, exeSAM, union, samples, qual, hierarchy, label, l2, annotation]

    # do the work boyeeee!
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
    cts=sorted(outCounts)

    # write the results
    for i in range(len(samples)):
        outFILE=os.path.join(genoDir,samples[i][1]+".counts.txt")
        with open(outFILE, "w") as fOUT:
            for j in range(len(cts)):
                fOUT.write("\t".join([str(y) for y in cts[j][1]])+"\n")

    # remove directory with beds/sams
    shutil.rmtree(tmpDir)

    print "\nTEFLON COUNT FINISHED!"

if __name__ == "__main__":
	main()
