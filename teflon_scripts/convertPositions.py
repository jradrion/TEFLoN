'''convert a bed file with coordinates from a true reference genome to the pseudo genome'''

def convertPositions_portal(bedFILE,mapFILE,outFILE):
    print "reading map"
    pseudoMap={}
    with open(mapFILE, 'r') as fIN:
        for line in fIN:
            if line.startswith(">"):
                chrom=line.replace("\n","").replace(">","")
            else:
                pseudoMap[chrom]=[int(x) for x in line.split(",")]
    print "converting bed"
    with open(bedFILE, 'r') as fIN, open(outFILE, 'w') as fOUT:
        for line in fIN:
            ls=line.split()
            chrom,start,stop=ls[0],int(ls[1]),int(ls[2])
            for ch in pseudoMap:
                if ch == chrom:
                    ls[1]=pseudoMap[chrom][start-1]
                    ls[2]=pseudoMap[chrom][stop-1]
                    fOUT.write("\t".join([str(x) for x in ls])+"\n")

