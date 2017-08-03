import os, sys
def avgReads(line,readLen):
    F,R="",""
    if line[1].isdigit():
        F=int(line[1])
    if line[2].isdigit():
        R=int(line[2])
    totalReads = int(line[9]) + int(line[10]) + int(line[11])
    if F != "" and R != "" and max([F,R])-min([F,R]) != 0:
        if line[7] == "-" and line[8] == "-":
            return int(totalReads/float(int((max([F,R])-min([F,R]))/float(readLen))+1))+1
        else:
            return int(totalReads/2.0)+1
    else:
        return totalReads

def fq(line):
    if int(line[9]) + int(line[10]) == 0:
        return -9
    else:
        return round(int(line[9])/float(int(line[9]) + int(line[10])), 3)


def pt_portal(genoDir, samples, posMap, stats, p2rC, c_thresh):
    readLen=stats[0]
    cts=[]
    for i in range(len(samples)):
        tmp=[]
        inFILE=os.path.join(genoDir,samples[i][1]+".counts.txt")
        with open(inFILE, "r") as fIN:
            for line in fIN:
                tmp.append(line.split())
        cts.append(tmp)

    outCall=[]
    for i in range(len(cts[0])):
        tmp=[]
        for j in range(len(cts)):
            if avgReads(cts[j][i],readLen) > c_thresh[j]:
                tmp.append(-9)
            else:
                tmp.append(fq(cts[j][i]))
        outCall.append(tmp)
    for i in range(len(cts)):
        outFILE1=os.path.join(genoDir,samples[i][1]+".pseudoSpace.genotypes.txt")
        with open(outFILE1, "w") as fOUT:
            for j in range(len(cts[i])):
                fOUT.write("\t".join([str(x) for x in cts[i][j]])+"\t%s\n" %(str(outCall[j][i])))
        outFILE2=os.path.join(genoDir,samples[i][1]+".refSpace.genotypes.txt")
        p2rC.pseudo2refConvert_portal(outFILE1,posMap,outFILE2)








