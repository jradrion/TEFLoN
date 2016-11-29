# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 15:29:00 2016
@author: jeffreyadrion
"""

import argparse, sys, os
import multiprocessing as mp
import subprocess as sp
import shlex
import glob
import shutil

teflonBase = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.insert(1, teflonBase)

from teflon_scripts import write_bed as wb
from teflon_scripts import cluster_positions as cp
from teflon_scripts import position_estimate as pe
from teflon_scripts import sort_positions as sortp

def mkdir_if_not_exist(*dirs):
    for dir in dirs:
        if not os.path.exists(dir):
            os.makedirs(dir)
            print "creating directory: %s" %(dir)

def check_samtools(exePATH):
    try:
        cmd = "%s" %(exePATH)
        p = sp.Popen(shlex.split(cmd), stdout=sp.PIPE, stderr=sp.PIPE)
    except OSError:
        print "Cannot find samtools"
        sys.exit(1)

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


def create_procs(nProcs, task_q, params):
    for _ in range(nProcs):
        p = mp.Process(target=worker, args=(task_q, params))
        p.daemon = True
        p.start()


def worker(task_q, params):
    while True:
        try:
            groups, nth_job = task_q.get()
            #unpack parameters
            annotation, bam, chromosomes, exePATH, hierarchy, insz, label, lengths, level, cLevel, qual, readLen, sd, cov, bedDir, samDir, posDir, suppDir = params
            for group in groups:
                levelIndex=label.index(level)
                clustIndex=label.index(cLevel)
                clustGroup=""
                for ID in hierarchy:
                    if hierarchy[ID][levelIndex] == group:
                        clustGroup=hierarchy[ID][clustIndex]
                        break
                # 1. write bed
                teIDs = wb.write_bed_portal(hierarchy, label, group, level, bedDir)
                # 2. samtools view > complete.sam
                try:
                    samFile_complete = os.path.join(samDir, "%s_complete.sam" %(group))
                    cmd = "%s view %s -L %s/%s_complete.bed" %(exePATH, bam, bedDir, group)
                    print "cmd:", cmd
                    p = sp.Popen(shlex.split(cmd), stdout=open(samFile_complete, 'w'), stderr=sp.PIPE)
                    perr = p.communicate()[1] # communicate returns a tuple (stdout, stderr)
                    print perr
                    if p.returncode != 0:
                        print "error running samtools"
                        sys.exit(1)
                except OSError:
                    print "Cannot run samtools"
                    sys.exit(1)
                # 3. cluster positions
                suppOutFile = os.path.join(suppDir, "%s_supplemental_alignments.txt" %(group))
                print "clustering:", group
                cp.cluster_positions_portal(samFile_complete, group, chromosomes, lengths, readLen, insz, sd, bedDir, qual, suppOutFile)
                #4. samtools view > clustered.sam
                try:
                    samFile_clust = os.path.join(samDir, "%s_clustered.sam" %(group))
                    cmd = "%s view %s -q %s -L %s/%s_clustered.bed" %(exePATH, bam, qual, bedDir, group)
                    print "cmd:", cmd
                    p = sp.Popen(shlex.split(cmd), stdout=open(samFile_clust, 'w'), stderr=sp.PIPE)
                    # communicate returns a tuple (stdout, stderr)
                    perr = p.communicate()[1]
                    print perr
                    if p.returncode != 0:
                        print "error running samtools"
                        sys.exit(1)
                except OSError:
                    print "Cannot run samtools"
                    sys.exit(1)
                #5. position estimate
                print 'estimating positions:',group
                pe.position_estimate_portal(samFile_clust, suppOutFile, group, annotation, teIDs, readLen, insz, sd, cov, posDir, clustGroup)
        finally:
            task_q.task_done()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-wd',dest='wd',help='full path to working directory')
    parser.add_argument('-ex',dest='exe',help='full path to samtools executable')
    parser.add_argument('-g',dest='genome',help='genomeSize file')
    parser.add_argument('-s',dest='samples',help='samples file')
    parser.add_argument('-b',dest='bam',help='full path to indexed bam file')
    parser.add_argument('-a',dest='annotation',help='full path to pseudo-transformed bedFILE')
    parser.add_argument('-t',dest='te_hierarchy',help='full path to te hierarchy file')
    parser.add_argument('-l',dest='level',help='level of hierarchy to search')
    parser.add_argument('-cl',dest='cLevel',help='level of hierarchy to cluster')
    parser.add_argument('-e',dest='exclude',help='newline separated list of te families to exclude from analysis', default=0)
    parser.add_argument('-q',dest='qual',help='map quality threshold',type=int)
    parser.add_argument('-sd',dest='stdev',help='insert size standard deviation override',type=int, default=-1)
    parser.add_argument("-x", dest="nProc", type=int, default=4, help="Specify number of processes")
    args = parser.parse_args()

    wdPath=args.wd
    exePATH=args.exe
    bam=args.bam
    level=args.level
    cLevel=args.cLevel
    qual=args.qual

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

    pre=""
    with open(args.samples, "r") as fIN:
        for line in fIN:
            if line.split()[0] == bam:
                pre = line.split()[1]
    if pre=="":
        print ".bam path in sample different from path in options"
        sys.exit()

    #Identify the name of all tes for the specified level of the hierarchy
    groups=[]
    groupIndex=label.index(level)
    for ID in hierarchy:
        groups.append(hierarchy[ID][groupIndex])
    groups=sorted(set(groups))


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

    #run stats and read results
    statsOutFile = bam.replace(".bam", ".stats.txt")
    #print "Calculating alignment statistics"
    cmd = "%s stats %s" %(exePATH, bam)
    print cmd
    p = sp.Popen(shlex.split(cmd), stdout=open(statsOutFile, 'w'), stderr=sp.PIPE)
    perr = p.communicate()[1]
    if p.returncode != 0:
        print "samtools stats issued error: %s" %(perr)
        sys.exit(1)

    #calculate coverage
    covFILE = bam.replace(".bam", ".cov.txt")
    #print covFILE
    cmd="""%s depth -Q %s %s | awk '{sum+=$3; sumsq+=$3*$3} END {print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' > %s""" %(exePATH, str(qual), bam, covFILE)
    print cmd
    os.system(cmd)

    #read stats file
    with open(statsOutFile, 'r') as fIN:
        for line in fIN:
            if 'average length' in line:
                readLen=int(float(line.split()[-1]))
            if 'insert size average' in line:
                insz=int(float(line.split()[-1]))
            if 'insert size standard deviation' in line:
                sd=int(float(line.split()[-1]))

    if sd > 75:
        if args.stdev == -1:
            print "!!! Warning: insert size standard deviation reported as",sd,"!!!"
            print "Please ensure this is correct or override in options"
            sys.exit()
        else:
            sd=args.stdev

    cov=25
    with open(covFILE, "r") as fIN:
        for line in fIN:
            if line.startswith("Av"):
                cov=int(float(line.split()[-1]))

    bedDir = args.wd + pre+".bed_files"
    samDir = args.wd + pre+".sam_files"
    posDir = args.wd + pre+".te_positions"
    suppDir = args.wd + pre+".supplemental_alignments"
    outDir =  args.wd + "initialPos"

    if args.exclude == 0:
        excludeList=[]
    else:
        excludeList=[]
        with open(args.exclude, "r") as fIN:
            for line in fIN:
                excludeList.append(line.split()[0])

    mkdir_if_not_exist(bedDir, posDir, samDir, suppDir, outDir)
    groups = [group for group in groups if group not in excludeList]
    #groups= ["1360"] #debug single family
    print "Groups:", groups
    task_q = mp.JoinableQueue()
    params=[annotation, bam, chromosomes, exePATH, hierarchy, insz, label, lengths, level, cLevel, qual, readLen, sd, cov, bedDir, samDir, posDir, suppDir]
    create_procs(args.nProc, task_q, params)
    assign_task(groups, task_q, args.nProc)
    try:
        task_q.join()
    except KeyboardInterrupt:
        print "KeyboardInterrupt"
        sys.exit(0)
    else:
        print "finished position estimates"
    #concatonate and cluster position estimates
    catFile = os.path.join(outDir,pre + ".all_positions.txt")
    try:
        files = ""
        for file in glob.glob(os.path.join(posDir, "*.txt")):
            files += file + " "
        cmd = "cat %s" %(files)
        print "cmd:", cmd  #p = sp.Popen(shlex.split(cmd), stdout=open(catFile, 'w'), stderr=sp.PIPE)
        p = sp.Popen(shlex.split(cmd), stdout=open(catFile, 'w'), stderr=sp.PIPE)
        perr = p.communicate()[1] # communicate returns a tuple (stdout, stderr)
        print perr
        if p.returncode != 0:
            print "error concatenating positions"
            sys.exit(1)
    except OSError:
        print "Cannot concatenate positions"
        sys.exit(1)

    #sort all positions
    print "Sorting positions"
    sortp.sort_portal(catFile)

    shutil.rmtree(bedDir)
    shutil.rmtree(samDir)
    shutil.rmtree(posDir)
    shutil.rmtree(suppDir)


if __name__ == "__main__":
	main()
