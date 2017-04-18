import os
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
    if len(samples) == 1:
        c_thresh=samples[0][2][3]+(10*samples[0][2][4])
        print "c_thresh:",c_thresh
        outFILE=os.path.join(genoDir,samples[0][1]+".genotypes.txt")
        with open(outFILE, "w") as fOUT:
            for j in range(len(cts[0])):
                if avgReads(cts[0][j],readLen) <= c_thresh and fq(cts[0][j]) != 0.00:
                    fOUT.write("\t".join([str(x) for x in cts[0][j]])+"\t%s\n" %(str(fq(cts[0][j]))))
    else:
        skip=[]
        c_thresh=[]
        for sample in samples:
            c_thresh.append(sample[2][3] + (10*sample[2][4]))
        print "c_thresh", c_thresh
        outCall=[]
        for i in range(len(cts[0])):
            tmp=[]
            for j in range(len(cts)):
                if avgReads(cts[j][i],readLen) > c_thresh[j]:
                    tmp.append(-1)
                else:
                    tmp.append(fq(cts[j][i]))
            if -1 in tmp and all(x==tmp[0] for x in tmp):
                #print tmp
                skip.append(i)
            elif 0.00 in tmp and all(x==tmp[0] for x in tmp):
                #print tmp
                skip.append(i)
            outCall.append(tmp)
        for i in range(len(cts)):
            outFILE1=os.path.join(genoDir,samples[i][1]+".pseudoSpace.genotypes.txt")
            with open(outFILE1, "w") as fOUT:
                for j in range(len(cts[i])):
                    if j not in skip:
                        fOUT.write("\t".join([str(x) for x in cts[i][j]])+"\t%s\n" %(str(outCall[j][i])))
            outFILE2=os.path.join(genoDir,samples[i][1]+".refSpace.genotypes.txt")
            p2rC.pseudo2refConvert_portal(outFILE1,posMap,outFILE2)








