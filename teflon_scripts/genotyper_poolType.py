import os, sys
def avgReads(line,readLen):
    F,R="",""
    if line[1].isdigit():
        F=int(line[1])
    if line[2].isdigit():
        R=int(line[2])
    totalReads = int(line[9]) + int(line[10]) + int(line[11])
    if F != "" and R != "" and max([F,R])-min([F,R]) != 0:
        if max([F,R])-min([F,R]) <= readLen:
            #print "reads",totalReads,1, totalReads/1, line
            return totalReads/1
        else:
            #print "reads", totalReads,float(max([F,R])-min([F,R])/readLen), int(totalReads/float(max([F,R])-min([F,R])/readLen)), line
            return int(totalReads/float(max([F,R])-min([F,R])/readLen))
    else:
        #print "reads",totalReads,1, totalReads/1, line
        return totalReads/1

def fq(line):
    if int(line[9]) + int(line[10]) == 0:
        return 0.00
    else:
        return round(int(line[9])/float(int(line[9]) + int(line[10])), 3)


def pt_portal(genoDir, samples, posMap, stats, p2rC):
    readLen=stats[0]
    cts=[]
    for i in range(len(samples)):
        tmp=[]
        inFILE=os.path.join(genoDir,samples[i][1]+".counts.txt")
        with open(inFILE, "r") as fIN:
            for line in fIN:
                tmp.append(line.split())
        cts.append(tmp)
    skip=[]
    c_thresh=[]
    for sample in samples:
        c_thresh.append(sample[2][3] + (4*sample[2][4]))
    c_thresh=sum(c_thresh)/float(len(c_thresh))
    print "coverage threshold (mean cov across all samples + 4*stdDev) :", c_thresh
    print "all TEs with read counts > coverage threshold will be reported with allele frequency = -1"
    outCall=[]
    for i in range(len(cts[0])):
        tmp=[]
        for j in range(len(cts)):
            if avgReads(cts[j][i],readLen) > c_thresh:
                tmp.append(-1)
            else:
                tmp.append(fq(cts[j][i]))
        if -1 in tmp and all(x==tmp[0] for x in tmp):
            #print tmp
            #skip.append(i)
            pass
        elif 0.00 in tmp and all(x==tmp[0] for x in tmp):
            #print tmp
            #skip.append(i)
            pass
        outCall.append(tmp)
    for i in range(len(cts)):
        outFILE1=os.path.join(genoDir,samples[i][1]+".pseudoSpace.genotypes.txt")
        with open(outFILE1, "w") as fOUT:
            for j in range(len(cts[i])):
                if j not in skip:
                    fOUT.write("\t".join([str(x) for x in cts[i][j]])+"\t%s\n" %(str(outCall[j][i])))
        outFILE2=os.path.join(genoDir,samples[i][1]+".refSpace.genotypes.txt")
        p2rC.pseudo2refConvert_portal(outFILE1,posMap,outFILE2)








