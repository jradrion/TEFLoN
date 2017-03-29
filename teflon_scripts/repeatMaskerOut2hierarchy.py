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
        fOUT.write("%s\t%s\t%s\n" %("id","hier_level_1","hier_level_2"))
        for i in range(len(ID)):
            if len(ID[i].split(".")) > 2:
                fOUT.write("%s\t%s\t%s\n" %(ID[i], ID[i].split(".")[1], ID[i].split(".")[2]))
            else:
                fOUT.write("%s\t%s\t%s\n" %(ID[i], ID[i].split(".")[0], ID[i].split(".")[1]))
