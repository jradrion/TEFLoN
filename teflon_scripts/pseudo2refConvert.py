import sys
def gapFlanks(ch,pos,psmap,FR):
    for i in range(int(pos)-4,int(pos)+5):
        if psmap[ch][i+1] - psmap[ch][i] > 1:
            if FR == "F":
                return psmap[ch][i]
            else:
                return psmap[ch][i+1]
    print "Error: no gap found", psmap[ch][int(pos)-4], psmap[ch][int(pos)+5]
    sys.exit()

def isGap(ch,pos,psmap):
    """if there is a gap near the position to be transformed
    return the left most position of the gap for a forward read and the rightmost for a reverse read
    else, return the normal transformation of the position"""
    if int(pos)+5 > len(psmap[ch]) or int(pos)-4 < 0:
        return 0
    elif abs(psmap[ch][int(pos)+5] - psmap[ch][int(pos)-4]) > 10:
        return 1
    else:
        return 0

def pseudo2refConvert_portal(bedFILE,pseudo2refMap,outFILE):
    pseudoMap=pseudo2refMap
    print "converting bed"
    with open(bedFILE, 'r') as fIN, open(outFILE, 'w') as fOUT:
        for line in fIN:
            ls=line.split()
            chrom=ls[0]
            for ch in pseudoMap:
                if ch == chrom:
                    if ls[1].isdigit() and not ls[2].isdigit():
                        if isGap(chrom,ls[1],pseudoMap) == 0:
                            ls[1]=pseudoMap[chrom][int(ls[1])-1]
                        else:
                            ls[1]=gapFlanks(chrom,ls[1],pseudoMap,"F")
                    elif not ls[1].isdigit() and ls[2].isdigit():
                        if isGap(chrom,ls[2],pseudoMap) == 0:
                            ls[2]=pseudoMap[chrom][int(ls[2])]
                        else:
                            ls[2]=gapFlanks(chrom,ls[2],pseudoMap,"R")
                    elif ls[1].isdigit() and ls[2].isdigit():
                        if isGap(chrom,ls[1],pseudoMap) == 0:
                            ls[1]=pseudoMap[chrom][int(ls[1])-1]
                        else:
                            ls[1]=gapFlanks(chrom,ls[1],pseudoMap,"F")
                        if isGap(chrom,ls[2],pseudoMap) == 0:
                            ls[2]=pseudoMap[chrom][int(ls[2])]
                        else:
                            ls[2]=gapFlanks(chrom,ls[2],pseudoMap,"R")
                    fOUT.write("\t".join([str(x) for x in ls])+"\n")

