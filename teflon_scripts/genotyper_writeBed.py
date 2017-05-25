import sys, os
import subprocess as sp
import shlex

def wb_portal(ID, union, samples, tmpDir, exePATH, qual):
    # 1. write bed
    outBED = os.path.join(tmpDir, "union_%s.bed" %(str(ID)))

    # 2. samtools view > sample.sam
    for sample in samples:
        try:
            samFILE = outBED.replace(".bed", ".%s.sam" %(sample[1]))
            cmd = "%s view -L %s %s" %(exePATH, outBED, sample[0])
            #print "cmd:", cmd
            p = sp.Popen(shlex.split(cmd), stdout=sp.PIPE)
            OUT = p.communicate()[0] # communicate returns a tuple (stdout, stderr)
            #print perr
            #os.remove(outBED)
            with open(samFILE, "w") as fOUT:
                for z in OUT:
                    fOUT.write(z)
            return OUT
            if p.returncode != 0:
                print "Error running samtools: p.returncode =",p.returncode
                #os.remove(outBED)
                return ""
        except OSError:
            print "Cannot run samtools"
            #os.remove(outBED)
            return ""
