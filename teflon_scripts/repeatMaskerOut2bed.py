def cluster(fam,prev,readLen):
    if fam[0] == prev[0]:
        if int(fam[1]) <= int(prev[2])+(readLen) and fam[3] == prev[3]:
            if int(fam[4]) >= int(prev[4]):
                carry_name=fam[-1]
            else:
                carry_name=prev[-1]
            return [fam[0],min([int(x) for x in [fam[1],fam[2],prev[1],prev[2]]]),max([int(x) for x in [fam[1],fam[2],prev[1],prev[2]]]), fam[3], max(int(fam[4]),int(prev[4])), carry_name]
        else:
            return 0
    else:
        return 0

def combine(fam,readLen):
    new_fam=[]
    for f in fam:
        tmp=[]
        for i in range(len(fam[f])):
            if not tmp:
                tmp.append(fam[f][i])
            else:
                clustered_elements=cluster(fam[f][i],tmp[-1],readLen)
                if clustered_elements == 0:
                    tmp.append(fam[f][i])
                else:
                    tmp[-1]=clustered_elements
        if tmp:
            for i in range(len(tmp)):
                if tmp[i][3]=="+":
                    strand = "+"
                else:
                    strand = "-"
                new_fam.append([tmp[i][0],int(tmp[i][1]),int(tmp[i][2]),tmp[i][-1]+"#"+f,strand])
    return new_fam


##Converts RepeatMasker out file to a bed file to be used with TEFLoN
def rep2bed_portal(inFILE,outFILE,readLen,minLen):
    print "Writing annotation file:",outFILE
    repMasked=[]
    with open(inFILE, "r") as fIN:
        for line in fIN:
            ar=line.split()
            if ar:
                if ar[0] not in "SWscore":
                    repMasked.append(ar)

    fam={}
    for i in range(len(repMasked)):
        name_pre=repMasked[i][9]
        name=repMasked[i][10]
        try:
            fam[name].append([repMasked[i][j] for j in [4,5,6,8,0]])
            fam[name][-1].append(name_pre)
        except KeyError:
            fam[name] = [[repMasked[i][j] for j in [4,5,6,8,0]]]
            fam[name][-1].append(name_pre)

    combined_fam=sorted(combine(fam,readLen), key= lambda x: [x[0],x[1]])

    with open(outFILE, "w") as fOUT:
        ct=1
        for line in combined_fam:
            if int(line[2])-int(line[1]) >= minLen:
                if ct<10:
                    ct="0"+str(ct)
                fOUT.write("%s\t%s\t%s\t%s\t.\t%s\n" %(line[0],line[1],line[2],"teid"+str(ct)+"#"+line[3],line[4]))
                ct=int(ct)
                ct+=1

