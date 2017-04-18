'''convert a bed file with coordinates from a true reference genome to the pseudo genome'''
def ref2pseudoConvert_portal(bedFILE,mapFILE,outFILE):
    pseudoMap=mapFILE
    print "Converting "+bedFILE+" to pseudospace ..."
    with open(bedFILE, 'r') as fIN, open(outFILE, 'w') as fOUT:
        for line in fIN:
            ls=line.split()
            chrom,start,stop=ls[0],int(ls[1]),int(ls[2])
            for ch in pseudoMap:
                if ch == chrom:
                    ls[1]=pseudoMap[chrom][start-1]
                    ls[2]=pseudoMap[chrom][stop-1]
                    fOUT.write("\t".join([str(x) for x in ls])+"\n")

