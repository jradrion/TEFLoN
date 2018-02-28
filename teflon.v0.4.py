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

teflonBase = os.path.dirname(os.path.abspath(sys.argv[0]))
sys.path.insert(1, teflonBase)

from teflon_scripts import write_bed as wb
from teflon_scripts import cluster_positions as cp
from teflon_scripts import position_estimate as pe
from teflon_scripts import sort_positions as sortp

def progress_bar(percent, barLen = 50):
    sys.stdout.write("\r")
    progress = ""
    for i in range(barLen):
        if i < int(barLen * percent):
            progress += "="
        else:
            progress += " "
    sys.stdout.write("[ %s ] %.2f%%" % (progress, percent * 100))
    sys.stdout.flush()

def mkdir_if_not_exist(*dirs):
    for dir in dirs:
        if not os.path.exists(dir):
            os.makedirs(dir)
            print "creating directory: %s" %(dir)

def check_dependency(exeSAM):
    try:
        cmd = "%s" %(exeSAM)
        sp.Popen(shlex.split(cmd), stdout=sp.PIPE, stderr=sp.PIPE)
    except OSError:
        print "Cannot find %s" %(exeSAM)
        sys.exit(1)

def assign_task(siteID, task_q, nProcs):
    c,i,nth_job=0,0,1
    while (i+1)*nProcs <= len(siteID):
        i+=1
    nP1=nProcs-(len(siteID)%nProcs)
    for j in range(nP1):
        task_q.put((siteID[c:c+i], nth_job))
        nth_job += 1
        c=c+i
    for j in range(nProcs-nP1):
        task_q.put((siteID[c:c+i+1], nth_job))
        nth_job += 1
        c=c+i+1

def create_proc1(nProcs, task_q, params):
    for _ in range(nProcs):
        p = mp.Process(target=worker1, args=(task_q, params))
        p.daemon = True
        p.start()

def create_proc2(nProcs, task_q, params):
    for _ in range(nProcs):
        p = mp.Process(target=worker2, args=(task_q, params))
        p.daemon = True
        p.start()

def create_proc3(nProcs, task_q, params):
    for _ in range(nProcs):
        p = mp.Process(target=worker3, args=(task_q, params))
        p.daemon = True
        p.start()

def worker1(task_q, params):
    while True:
        try:
            groups, nth_job = task_q.get()
            #unpack parameters
            annotation, bam, chromosomes, exeSAM, hierarchy, insz, label, lengths, level, cLevel, qual, readLen, sd, cov, bedDir, samDir, posDir, suppDir = params
            for group in groups:
                #print "group:",group
                wb.write_bed_portal(hierarchy, label, group, level, bedDir)
        finally:
            task_q.task_done()

def worker2(task_q, params):
    while True:
        try:
            groups, nth_job = task_q.get()
            #unpack parameters
            annotation, bam, chromosomes, exeSAM, hierarchy, insz, label, lengths, level, cLevel, qual, readLen, sd, cov, bedDir, samDir, posDir, suppDir = params
            ct=1
            for group in groups:
                # 2. samtools view > complete.sam
                try:
                    samFile_complete = os.path.join(samDir, "%s_complete.sam" %(group))
                    cmd = "%s view %s -L %s/%s_complete.bed" %(exeSAM, bam, bedDir, group)
                    #print "cmd:", cmd
                    p = sp.Popen(shlex.split(cmd), stdout=open(samFile_complete, 'w'), stderr=sp.PIPE)
                    perr = p.communicate()[1] # communicate returns a tuple (stdout, stderr)
                    #print perr
                    if p.returncode != 0:
                        print "error running samtools"
                        sys.exit(1)
                except OSError:
                    print "Cannot run samtools"
                    sys.exit(1)
                # 3. cluster positions
                suppOutFile = os.path.join(suppDir, "%s_supplemental_alignments.txt" %(group))
                cp.cluster_positions_portal(samFile_complete, group, chromosomes, lengths, readLen, insz, sd, bedDir, qual, suppOutFile)
                if nth_job==1:
                    progress_bar(ct/float(len(groups)))
                ct+=1
        finally:
            task_q.task_done()

def worker3(task_q, params):
    while True:
        try:
            groups, nth_job = task_q.get()
            #unpack parameters
            annotation, bam, chromosomes, exeSAM, hierarchy, insz, label, lengths, level, cLevel, qual, readLen, sd, cov, bedDir, samDir, posDir, suppDir = params
            ct=1
            for group in groups:
                teIDs=[]
                bed=os.path.join(bedDir,group+"_complete.bed")
                with open(bed, "r") as fIN:
                    for line in fIN:
                        teIDs.append(line.split()[0])
                levelIndex=label.index(level)
                clustIndex=label.index(cLevel)
                clustGroup=""
                for ID in hierarchy:
                    if hierarchy[ID][levelIndex] == group:
                        clustGroup=hierarchy[ID][clustIndex]
                        break
                suppOutFile = os.path.join(suppDir, "%s_supplemental_alignments.txt" %(group))
                #4. samtools view > clustered.sam
                try:
                    samFile_clust = os.path.join(samDir, "%s_clustered.sam" %(group))
                    cmd = "%s view %s -L %s/%s_clustered.bed" %(exeSAM, bam, bedDir, group)
                    #print "cmd:", cmd
                    p = sp.Popen(shlex.split(cmd), stdout=open(samFile_clust, 'w'), stderr=sp.PIPE)
                    # communicate returns a tuple (stdout, stderr)
                    perr = p.communicate()[1]
                    #print perr
                    if p.returncode != 0:
                        print "error running samtools"
                        sys.exit(1)
                except OSError:
                    print "Cannot run samtools"
                    sys.exit(1)
                #5. position estimate
                pe.position_estimate_portal(samFile_clust, suppOutFile, group, annotation, teIDs, readLen, insz, sd, cov, posDir, clustGroup, hierarchy)
                if nth_job==1:
                    progress_bar(ct/float(len(groups)))
                ct+=1
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
    parser.add_argument('-l1',dest='level',help='level of hierarchy to guide initial search')
    parser.add_argument('-l2',dest='cLevel',help='level of hierarchy to cluster')
    parser.add_argument('-q',dest='qual',help='map quality threshold',type=int)
    parser.add_argument('-exclude',dest='exclude',help='newline separated list of te families to exclude from analysis', default=-1)
    parser.add_argument('-sd',dest='stdev',help='insert size standard deviation override',type=int, default=-1)
    parser.add_argument('-cov',dest='cov',help='manual coverage override',type=int, default=-1)
    parser.add_argument("-t", dest="nProc", type=int, default=1, help="Specify number of processes")
    args = parser.parse_args()

    # identify current working directory
    if args.wd == -1:
        cwd=os.getcwd()
    else:
        cwd=os.path.abspath(args.wd)

    # import options
    prep_TF=os.path.abspath(args.DIR)
    prefix=os.path.abspath(args.DIR).split("/")[-1].replace(".prep_TF","")
    exeSAM=args.exeSAM
    exeBWA=args.exeBWA
    level=args.level
    cLevel=args.cLevel
    qual=args.qual
    nProc=args.nProc

    # check dependencies for function
    check_dependency(exeSAM)
    check_dependency(exeBWA)

    # import hierarchy
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
    bam,pre="",""
    with open(os.path.abspath(args.samples), "r") as fIN:
        for line in fIN:
            if line.split()[1] == args.ID:
                pre = line.split()[1]
                bam = line.split()[0]
    if pre=="" or bam=="":
        print "Warning: prefix in samples file different from path in options"
        sys.exit()

    # identify the group-name of all TEs for the specified level of the hierarchy
    groups=[]
    groupIndex=label.index(level)
    for ID in hierarchy:
        groups.append(hierarchy[ID][groupIndex])
    groups=sorted(set(groups))

    # import the TE annotation
    annotation=[]
    with open(os.path.join(prep_TF,prefix+".te.pseudo.bed"), 'r') as fIN:
        for line in fIN:
            arr=line.split()
            annotation.append([arr[0],int(arr[1]),int(arr[2]),arr[3]])

    # import the chromosome lengths
    chromosomes,lengths=[],[]
    genomeSizeFILE=os.path.join(prep_TF,prefix+".genomeSize.txt")
    with open(genomeSizeFILE, 'r') as fIN:
        for line in fIN:
            arr=line.split()
            chromosomes.append(arr[0])
            lengths.append(int(arr[2]))

    # run samtools stats
    statsOutFile = bam.replace(".bam", ".stats.txt")
    print "Calculating alignment statistics"
    cmd = "%s stats -t %s %s" %(exeSAM, genomeSizeFILE, bam)
    print "cmd:",cmd
    p = sp.Popen(shlex.split(cmd), stdout=open(statsOutFile, 'w'), stderr=sp.PIPE)
    perr = p.communicate()[1]
    if p.returncode != 0:
        print "samtools stats issued error: %s" %(perr)
        sys.exit(1)

    # calculate coverage
    covFILE = bam.replace(".bam", ".cov.txt")
    cmd="""%s depth -Q %s %s | awk '{sum+=$3; sumsq+=$3*$3} END {print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' > %s""" %(exeSAM, str(qual), bam, covFILE)
    print "cmd:",cmd
    os.system(cmd)

    # read samtools stats file
    with open(statsOutFile, 'r') as fIN:
        for line in fIN:
            if 'average length' in line:
                readLen=int(float(line.split()[-1]))
            if 'insert size average' in line:
                insz=int(float(line.split()[-1]))
            if 'insert size standard deviation' in line:
                sd=int(float(line.split()[-1]))

    if args.stdev == -1:
        print "Insert size standard deviation estimated as %s. Use the override option if you suspect this is incorrect!" %(sd)
        if sd > 100:
            print "!!! Warning: insert size standard deviation reported as",sd,"!!!"
            print "Please ensure this is correct and use the override option!"
            sys.exit()
    else:
        sd=args.stdev

    # read coverage file
    cov=args.cov
    with open(covFILE, "r") as fIN:
        for line in fIN:
            if line.startswith("Av"):
                cov=int(float(line.split()[-1]))
    if cov == -1:
        print "Warning: coverage could not be estimated, enter coverage manually"
        sys.exit()

    # read list of TE groups to exclude from analysis
    if args.exclude == -1:
        excludeList=[]
    else:
        excludeList=[]
        with open(args.exclude, "r") as fIN:
            for line in fIN:
                excludeList.append(line.split()[0])

    # define and create subdirectories
    bedDir = os.path.join(cwd,pre+".bed_files")
    samDir = os.path.join(cwd,pre+".sam_files")
    posDir = os.path.join(cwd,pre+".te_positions")
    suppDir = os.path.join(cwd,pre+".supplemental_alignments")
    outDir =  os.path.join(cwd,"countPos")

    mkdir_if_not_exist(bedDir, posDir, samDir, suppDir, outDir)
    groups = [group for group in groups if group not in excludeList]
    #groups= ["doc3"] #debug single family
    print "Groups to search:", groups

    # run multiprocess 1
    print "writing TE bed files..."
    task_q = mp.JoinableQueue()
    params=[annotation, bam, chromosomes, exeSAM, hierarchy, insz, label, lengths, level, cLevel, qual, readLen, sd, cov, bedDir, samDir, posDir, suppDir]
    create_proc1(nProc, task_q, params)
    assign_task(groups, task_q, nProc)
    try:
        task_q.join()
    except KeyboardInterrupt:
        print "KeyboardInterrupt"
        sys.exit(0)
    else:
        print "writing TE bed files completed!"

    # reduce search-space 1
    print "reducing search space..."
    try:
        bedFILE = os.path.join(bedDir,"mega_complete.bed")
        bamFILE = os.path.join(samDir,"mega_complete.bam")
        cmd = "%s view -@ %s -L %s %s -b" %(exeSAM, str(nProc), bedFILE, bam)
        print "cmd:", cmd
        p = sp.Popen(shlex.split(cmd), stdout=open(bamFILE, 'w'), stderr=sp.PIPE)
        perr = p.communicate()[1] # communicate returns a tuple (stdout, stderr)
        #print perr
        if p.returncode != 0:
            print "Error running samtools: p.returncode =",p.returncode
            sys.exit(1)
    except OSError:
        print "Cannot run samtools"
    print "search space succesfully reduced..."
    print "new reduced bam file:",bamFILE


    # index new reduced alignment 1
    cmd="%s index %s" %(exeBWA,bamFILE)
    print "cmd:", cmd
    os.system(cmd)

    # run multiprocess 2
    print "clustering TE positions..."
    task_q = mp.JoinableQueue()
    params=[annotation, bamFILE, chromosomes, exeSAM, hierarchy, insz, label, lengths, level, cLevel, qual, readLen, sd, cov, bedDir, samDir, posDir, suppDir]
    create_proc2(nProc, task_q, params)
    assign_task(groups, task_q, nProc)
    try:
        task_q.join()
    except KeyboardInterrupt:
        print "KeyboardInterrupt"
        sys.exit(0)
    else:
        print "\nclustering TE positions completed!"

    # reduce search-space 2
    print "final reduction of search space..."
    try:
        bedFILE = os.path.join(bedDir,"mega_clustered.bed")
        bamFILE = os.path.join(samDir,"mega_clustered.bam")
        cmd = "%s view -@ %s -q %s -L %s %s -b" %(exeSAM, str(nProc), str(qual), bedFILE, bam)
        print "cmd:", cmd
        p = sp.Popen(shlex.split(cmd), stdout=open(bamFILE, 'w'), stderr=sp.PIPE)
        perr = p.communicate()[1] # communicate returns a tuple (stdout, stderr)
        #print perr
        if p.returncode != 0:
            print "Error running samtools: p.returncode =",p.returncode
            sys.exit(1)
    except OSError:
        print "Cannot run samtools"
    print "search space succesfully reduced..."
    print "new reduced bam file:",bamFILE

    # index new reduced alignment 2
    cmd="%s index %s" %(exeBWA,bamFILE)
    print "cmd:", cmd
    os.system(cmd)


    # run multiprocess 3
    print "estimating TE breakpoints..."
    bamFILE = os.path.join(samDir,"mega_clustered.bam")
    task_q = mp.JoinableQueue()
    params=[annotation, bamFILE, chromosomes, exeSAM, hierarchy, insz, label, lengths, level, cLevel, qual, readLen, sd, cov, bedDir, samDir, posDir, suppDir]
    create_proc3(nProc, task_q, params)
    assign_task(groups, task_q, nProc)
    try:
        task_q.join()
    except KeyboardInterrupt:
        print "KeyboardInterrupt"
        sys.exit(0)
    else:
        print "\nestimating TE breakpoints completed!"

    # concatonate position estimates
    catFile = os.path.join(outDir,pre + ".all_positions.txt")
    try:
        files = ""
        for file in glob.glob(os.path.join(posDir, "*.txt")):
            files += file + " "
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

    # sort position estimates
    print "Sorting positions..."
    sortp.sort_portal(catFile)

    # remove temporary directories
    shutil.rmtree(bedDir)
    shutil.rmtree(samDir)
    shutil.rmtree(posDir)
    shutil.rmtree(suppDir)

    print "TEFLON DISCOVERY FINISHED!"


if __name__ == "__main__":
	main()
