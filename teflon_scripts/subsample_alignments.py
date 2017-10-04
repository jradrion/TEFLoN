import sys, os
import multiprocessing as mp
import subprocess as sp
import shlex

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

def create_proc(nProcs, task_q, params):
    for _ in range(nProcs):
        p = mp.Process(target=worker, args=(task_q, params))
        p.daemon = True
        p.start()

def worker(task_q, params):
    while True:
        try:
            samples, nth_job = task_q.get()
            #unpack parameters
            exePATH, qual, min_seqs, DIR = params
            prep_TF=DIR
            #prefix=os.path.dirname(prep_TF).split("/")[-1].split(".prep_TF")[0]
            prefix=os.path.basename(DIR).replace(".prep_TF","")
            for sample in samples:
                bamFILE = sample[0].replace(".bam",".subsmpl.bam")
                if not os.path.exists(bamFILE):
                    print "Creating",bamFILE
                    try:
                        fract=round(1/float(int(sample[2][4])/float(min_seqs)),3)
                        if fract < 1.0:
                            cmd = "%s view -b -s %s %s" %(exePATH, fract, sample[0])
                            print "cmd:", cmd
                            p = sp.Popen(shlex.split(cmd), stdout=open(bamFILE, 'w'), stderr=sp.PIPE)
                            perr = p.communicate()[1] # communicate returns a tuple (stdout, stderr)
                            #print perr
                            if p.returncode != 0:
                                print "Error running samtools: p.returncode =",p.returncode
                                continue
                        else:
                            cmd = "cp %s %s" %(sample[0], bamFILE)
                            print "cmd:", cmd
                            os.system(cmd)
                    except OSError:
                        print "Cannot run samtools"
                        sys.exit(1)
                else:
                    print bamFILE,"exists"

                # 2. index new bam file
                cmd = "%s index %s" %(exePATH, bamFILE)
                print "cmd:", cmd
                os.system(cmd)

                # 3. run samtools stats on subsampled alignments
                statsOutFile = bamFILE.replace(".bam", ".stats.txt")
                genomeSizeFILE=os.path.join(prep_TF,prefix+".genomeSize.txt")
                if not os.path.exists(statsOutFile):
                    print "Calculating alignment statistics for", bamFILE
                    cmd = "%s stats -t %s %s" %(exePATH, genomeSizeFILE, bamFILE)
                    print "cmd:",cmd
                    p = sp.Popen(shlex.split(cmd), stdout=open(statsOutFile, 'w'), stderr=sp.PIPE)
                    perr = p.communicate()[1]
                    if p.returncode != 0:
                        print "samtools stats issued error: %s" %(perr)
                        sys.exit(1)

                    #calculate coverage
                    covFILE = bamFILE.replace(".bam", ".cov.txt")
                    #print covFILE
                    cmd="""%s depth -Q %s %s | awk '{sum+=$3; sumsq+=$3*$3} END {print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' > %s""" %(exePATH, str(qual), bamFILE, covFILE)
                    print "cmd:",cmd
                    os.system(cmd)
                else:
                    print statsOutFile, "exists"
        finally:
            task_q.task_done()

def subsample_alignments_portal(samples, exePATH, nProc, qual, cov, DIR):
    if cov == -1:
        min_seqs=int(min(x[2][4] for x in samples))
        print "Subsampling to the minimum read depth from a single sample:", min_seqs
    else:
        min_seqs=min(x[2][4] for x in samples)
        if cov > int(min_seqs):
            print "Coverage override ,%sX, is greater than minimum coverage in the group, %sX" %(cov,min_seqs)
            print "Reduce coverage override and rerun!"
            sys.exit()
        else:
            min_seqs=int(cov)
            print "Subsampling to the minimum read depth from a single sample:", min_seqs

    # multithread
    task_q = mp.JoinableQueue()
    params=[exePATH, qual, min_seqs, DIR]
    create_proc(nProc, task_q, params)
    assign_task(samples, task_q, nProc)
    try:
        task_q.join()
    except KeyboardInterrupt:
        print "KeyboardInterrupt"
        sys.exit(0)
    else:
        print "finished standardaizing sample depth"
