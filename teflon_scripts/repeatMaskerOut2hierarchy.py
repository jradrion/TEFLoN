##Converts RepeatMasker out file to a bed file to be used with TEFLoN
def rep2hier_portal(repBase,inFILE,outFILE):
    print "Writing hierarchy file:",outFILE
    ID=[]
    with open(inFILE, "r") as fIN:
        for line in fIN:
            ar = line.split()
            ID.append(ar[3])

    with open(repBase, "r") as fIN:
        for line in fIN:
            if line.startswith(">"):
                ID.append(line[1:].replace("\n",""))

    with open(outFILE, "w") as fOUT:
        fOUT.write("%s\t%s\t%s\n" %("id","family","superfamily"))
        for i in range(len(ID)):
            if "//" in ID[i]:
                fOUT.write("%s\t%s\t%s\n" %(ID[i],ID[i].split("//")[1].split("|")[0],ID[i].split("|")[1]))
            else:
                fOUT.write("%s\t%s\t%s\n" %(ID[i],ID[i].split("|")[0],ID[i].split("|")[1]))


