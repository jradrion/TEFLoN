import sys, os
import subprocess as sp
import shlex

def wb_portal(ID, union, samples, tmpDir, exePATH, qual):
    # 1. write bed
    bed=[]
    if union[ID][7] == "+":
        bed.append([union[ID][0], int(union[ID][1])-1, int(union[ID][1])+1])
    if union[ID][8] == "+":
        bed.append([union[ID][0], int(union[ID][2])-1, int(union[ID][2])+1])
    if union[ID][7] == "-" and union[ID][8] == "-" and union[ID][1] != "-" and union[ID][2] != "-":
        bed.append([union[ID][0], min(int(union[ID][1]),int(union[ID][2]))-1, max(int(union[ID][1]),int(union[ID][2]))+1])
    outBED = os.path.join(tmpDir, "union_%s.bed" %(str(ID)))
    with open(outBED, 'w') as fOUT:
        for line in bed:
            fOUT.write("\t".join([str(x) for x in line]) + "\n")

    # 2. samtools view > sample.sam
    for sample in samples:
        try:
            samFILE = outBED.replace(".bed", ".%s.sam" %(sample[1]))
            cmd = "%s view -q %s -L %s %s" %(exePATH, str(qual), outBED, sample[0])
            #print "cmd:", cmd
            p = sp.Popen(shlex.split(cmd), stdout=open(samFILE, 'w'), stderr=sp.PIPE)
            perr = p.communicate()[1] # communicate returns a tuple (stdout, stderr)
            #print perr
            if p.returncode != 0:
                print "error running samtools"
                sys.exit(1)
        except OSError:
            print "Cannot run samtools"
            sys.exit(1)
    os.remove(outBED)
