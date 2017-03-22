import sys

##Converts RepeatMasker out file to a bed file to be used with TEFLoN
def rep2bed_portal(inFILE,outFILE):
    print "Writing annotation file:",outFILE
    repMasked=[]
    with open(inFILE, "r") as fIN:
        for line in fIN:
            ar=line.split()
            if ar:
                if ar[0] not in "SWscore":
                    repMasked.append(ar)

    for i in range(10):
        print repMasked[i]
    sys.exit()

    with open(outFILE, "w") as fOUT:
        for line in repMasked:
            if line[8] == "C":
                strand = "-"
            else:
                strand = "+"
            fOUT.write("%s\t%s\t%s\t%s\t.\t%s\n" %(line[4],line[5],line[6],"TEid"+line[14]+"."+line[9],strand))

