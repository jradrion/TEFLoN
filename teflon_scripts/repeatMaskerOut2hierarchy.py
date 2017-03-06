##Converts RepeatMasker out file to a bed file to be used with TEFLoN
def rep2hier_portal(inFILE,outFILE):
    print "Writing hierarchy file:",outFILE
    ID,fam,superfam=[],[],[]
    with open(inFILE, "r") as fIN:
        for line in fIN:
            ar = line.split()
            ID.append(ar[3])
            fam.append(ar[3].split("//")[1].split("|")[0])
            superfam.append(ar[3].split("|")[1])

    with open(outFILE, "w") as fOUT:
        fOUT.write("%s\t%s\t%s\n" %("id","family","superfamily"))
        for i in range(len(ID)):
            fOUT.write("%s\t%s\t%s\n" %(ID[i],ID[i].split("//")[1].split("|")[0],ID[i].split("|")[1]))


