#Takes a bam file and a list of sites (ie. the union file) and creates a new, smaller, bam file to search
import sys, os
import subprocess as sp
import shlex

def reduceSearchSpace_portal(union, samples, tmpDir, exePATH, qual, nProc):
    print "reducing search space..."
    # 1. write bed
    beds=[]
    for ID in range(len(union)):
        bed=[]
        if union[ID][7] == "+":
            bed.append([union[ID][0], int(union[ID][1])-1, int(union[ID][1])+1])
            beds.append([union[ID][0], int(union[ID][1])-1, int(union[ID][1])+1])
        if union[ID][8] == "+":
            bed.append([union[ID][0], int(union[ID][2])-1, int(union[ID][2])+1])
            beds.append([union[ID][0], int(union[ID][2])-1, int(union[ID][2])+1])
        if union[ID][7] == "-" and union[ID][8] == "-" and union[ID][1] != "-" and union[ID][2] != "-":
            bed.append([union[ID][0], min(int(union[ID][1]),int(union[ID][2]))-1, max(int(union[ID][1]),int(union[ID][2]))+1])
            beds.append([union[ID][0], min(int(union[ID][1]),int(union[ID][2]))-1, max(int(union[ID][1]),int(union[ID][2]))+1])
        outBED = os.path.join(tmpDir, "union_%s.bed" %(str(ID)))
        with open(outBED, 'w') as fOUT:
            for line in bed:
                try:
                    fOUT.write("\t".join([str(x) for x in line]) + "\n")
                except IOError:
                    continue

    megaBed=os.path.join(tmpDir,"megaBed.bed")
    with open(megaBed, 'w') as fOUT:
        for line in beds:
            try:
                fOUT.write("\t".join([str(x) for x in line]) + "\n")
            except IOError:
                continue

    try:
        bamFILE = megaBed.replace(".bed", ".bam")
        cmd = "%s view -@ %s -q %s -L %s %s -b" %(exePATH, str(nProc), str(qual), megaBed, samples[0][0])
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
    return bamFILE

