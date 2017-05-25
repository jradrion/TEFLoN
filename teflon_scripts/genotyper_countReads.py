import os, sys

def mean(lst):
    return sum(lst)/float(len(lst))

def mode(lst):
    return max(set(lst), key=lst.count)

def bitFlag(flag):
    '''parse the bitwise flag, returning a list of explanations with 0 if the bit is absent and 1 if present'''
    '''0 read paired
       1 read mapped in proper pair
       2 read unmapped
       3 mate unmapped
       4 read reverse strand
       5 mate reverse strand
       6 first in pair
       7 second in pair
       8 not primary alignment
       9 read fails platform/vendor quality checks
       10 read is PCR or optical duplicate
       11 supplementary alignment'''
    explanations=[0,0,0,0,0,0,0,0,0,0,0,0]
    for i in xrange(len(explanations)):
        if flag & 2**i:
            explanations[i] =1
    return explanations

def cigarParse(cig):
    operation = ['M','I','D','N','S','H','P','=','X']
    opValue = [0,0,0,0,0,0,0,0,0]
    for i in xrange(1):
        arr = cig.split(operation[i])
        for j in xrange(len(arr)):
            if arr[j].isdigit():
                opValue[i] += int(arr[j])
                arr[j]=''
            else:
                for k in xrange(len(operation)):
                    switch =1
                    arr2 = arr[j].split(operation[k])
                    for l in xrange(len(arr2)):
                        if arr2[l].isdigit() and switch ==1:
                            opValue[k] += int(arr2[l])
                            switch =0
                        elif arr2[l].isdigit() and switch ==0:
                            opValue[i] += int(arr2[l])
    return opValue

def rightMostPos(leftMostPos,opValue):
    '''add M=match/mismatch and D=deletions to left-most position'''
    return int(leftMostPos) + opValue[0]+opValue[2] - 1

def side_clipped(cigar):
    '''which side of the read is clipped?'''
    temp=[]
    for i in cigar:
        if i in 'MS':
            temp.append(i)
    if temp[0] == 'S' and temp[-1] != 'S':
        return 'L'
    if temp[0] != 'S' and temp[-1] == 'S':
        return 'R'
    if temp[0] == 'S' and temp[-1] == 'S':
        return 'LR'
    else:
        return ''

def isPrimary(read):
    if bitFlag(int(read[1]))[8]==0 and bitFlag(int(read[1]))[11]==0:
        return 1
    else:
        return 0


def ct_type0(reads, union, TEID, annotation):
    '''Type 0 events are characterized by non-clipped reads flanking both sides of an insertion'''
    """try to find the breakpoints from clipped reads not aligning to tes"""
    F_anno,R_anno="",""
    F_flag,R_flag=0,0
    try:
        F_anno=annotation[union[6]][1]
        R_anno=annotation[union[6]][2]
    except KeyError:
        F_anno,R_anno="",""
    F,R,f_clip,r_clip="","",[],[]
    for read in reads:
        if "R" in side_clipped(read[5]):
            f_clip.append(rightMostPos(int(read[3]),cigarParse(read[5])))
        if "L" in side_clipped(read[5]):
            r_clip.append(int(read[3]))
    if f_clip:
        F=mode(f_clip)
    if r_clip:
        R=mode(r_clip)
    if F == "":
        if F_anno:
            F=F_anno
        else:
            F_flag=1
            F=int(union[1])
    if R == "":
        if R_anno:
            R=R_anno
        else:
            F_flag=1
            R=int(union[2])
    if F_flag == 0 and R_flag == 0:
        cts=ct_type1FR(reads, F, R, union, TEID)
        cts[-2]=F
        cts[-1]=R
        return cts
    elif F_flag == 0 and R_flag == 1:
        cts=ct_type1F(reads, F, R, union, TEID)
        cts[-2]=F
        return cts
    elif F_flag == 1 and R_flag == 0:
        cts=ct_type1R(reads, F, R, union, TEID)
        cts[-1]=R
        return cts
    else:
        return [[0,0,0],"-","-"]

def ct_type1FR(reads, F, R, union, TEID):
    '''Type 1 events are characterized by the F support breakpoint following the R support breakpoint (ie TSD)'''
    #print "type1FR"
    cts=[0,0,0]
    overShoot=3
    for read in reads:
        if isPrimary(read) == 1:
            try:
                TEID[read[6]]
                if bitFlag(int(read[1]))[4] == 0 and rightMostPos(int(read[3]),cigarParse(read[5])) < F+overShoot:
                    #print read[0], "\t", "p"
                    cts[0]+=1
                elif bitFlag(int(read[1]))[4] == 1 and int(read[3]) > R-overShoot:
                    #print read[0], "\t", "p"
                    cts[0]+=1
                else:
                    TEID["EXCEPT!"]
            except KeyError:
                if int(read[3]) <= R-overShoot and rightMostPos(int(read[3]),cigarParse(read[5])) >= F+overShoot:
                    #print read[0], "\t", "a"
                    cts[1]+=1
                elif bitFlag(int(read[1]))[4] == 0 and rightMostPos(int(read[3]),cigarParse(read[5])) < F+overShoot and "R" in side_clipped(read[5]):
                    #print read[0], "\t", "p"
                    cts[0]+=1
                elif bitFlag(int(read[1]))[4] == 0 and int(read[3])  > R-overShoot and  "L" in side_clipped(read[5]):
                    #print read[0], "\t", "p"
                    cts[0]+=1
                elif bitFlag(int(read[1]))[4] == 1 and int(read[3])  > R-overShoot and  "L" in side_clipped(read[5]):
                    #print read[0], "\t", "p"
                    cts[0]+=1
                elif bitFlag(int(read[1]))[4] == 1 and rightMostPos(int(read[3]),cigarParse(read[5])) < F+overShoot and "R" in side_clipped(read[5]):
                    #print read[0], "\t", "p"
                    cts[0]+=1
                else:
                    #print read[0], "\t", "na"
                    cts[2]+=1
        else:
            #print read[0], "\t", "na"
            cts[2]+=1
    return [cts,"-","-"]

def ct_type1F(reads, F, R, union, TEID):
    '''Type 1 events are characterized by the F support breakpoint following the R support breakpoint (ie TSD)'''
    #print "type1F"
    if union[6] == "-":
        R=F-12
    else:
        R=F-1
    cts=[0,0,0]
    overShoot=3
    for read in reads:
        if isPrimary(read) == 1:
            try:
                TEID[read[6]]
                if bitFlag(int(read[1]))[4] == 0 and rightMostPos(int(read[3]),cigarParse(read[5])) < F+overShoot:
                    #print read[0], "\t", "p"
                    cts[0]+=1
                elif bitFlag(int(read[1]))[4] == 1 and int(read[3]) > R-overShoot:
                    #print read[0], "\t", "p"
                    cts[0]+=1
                else:
                    TEID["EXCEPT!"]
            except KeyError:
                if int(read[3]) <= R-overShoot and rightMostPos(int(read[3]),cigarParse(read[5])) >= F+overShoot:
                    #print read[0], "\t", "a"
                    cts[1]+=1
                elif bitFlag(int(read[1]))[4] == 0 and rightMostPos(int(read[3]),cigarParse(read[5])) < F+overShoot and "R" in side_clipped(read[5]):
                    #print read[0], "\t", "p"
                    cts[0]+=1
                elif bitFlag(int(read[1]))[4] == 0 and int(read[3])  > R-overShoot and  "L" in side_clipped(read[5]):
                    #print read[0], "\t", "p"
                    cts[0]+=1
                elif bitFlag(int(read[1]))[4] == 1 and int(read[3])  > R-overShoot and  "L" in side_clipped(read[5]):
                    #print read[0], "\t", "p"
                    cts[0]+=1
                elif bitFlag(int(read[1]))[4] == 1 and rightMostPos(int(read[3]),cigarParse(read[5])) < F+overShoot and "R" in side_clipped(read[5]):
                    #print read[0], "\t", "p"
                    cts[0]+=1
                else:
                    #print read[0], "\t", "na"
                    cts[2]+=1
        else:
            #print read[0], "\t", "na"
            cts[2]+=1
    return [cts,"-","-"]

def ct_type1R(reads, F, R, union, TEID):
    #print "type1R"
    """Same with R only"""
    if union[6] == "-":
        F=R+12
    else:
        F=R+1
    '''Type 1 events are characterized by the F support breakpoint following the R support breakpoint (ie TSD)'''
    cts=[0,0,0]
    overShoot=3
    for read in reads:
        if isPrimary(read) == 1:
            try:
                TEID[read[6]]
                if bitFlag(int(read[1]))[4] == 0 and rightMostPos(int(read[3]),cigarParse(read[5])) < F+overShoot:
                    #print read[0], "\t", "p"
                    cts[0]+=1
                elif bitFlag(int(read[1]))[4] == 1 and int(read[3]) > R-overShoot:
                    #print read[0], "\t", "p"
                    cts[0]+=1
                else:
                    TEID["EXCEPT!"]
            except KeyError:
                if int(read[3]) <= R-overShoot and rightMostPos(int(read[3]),cigarParse(read[5])) >= F+overShoot:
                    #print read[0], "\t", "a"
                    cts[1]+=1
                elif bitFlag(int(read[1]))[4] == 0 and rightMostPos(int(read[3]),cigarParse(read[5])) < F+overShoot and "R" in side_clipped(read[5]):
                    #print read[0], "\t", "p"
                    cts[0]+=1
                elif bitFlag(int(read[1]))[4] == 0 and int(read[3])  > R-overShoot and  "L" in side_clipped(read[5]):
                    #print read[0], "\t", "p"
                    cts[0]+=1
                elif bitFlag(int(read[1]))[4] == 1 and int(read[3])  > R-overShoot and  "L" in side_clipped(read[5]):
                    #print read[0], "\t", "p"
                    cts[0]+=1
                elif bitFlag(int(read[1]))[4] == 1 and rightMostPos(int(read[3]),cigarParse(read[5])) < F+overShoot and "R" in side_clipped(read[5]):
                    #print read[0], "\t", "p"
                    cts[0]+=1
                else:
                    #print read[0], "\t", "na"
                    cts[2]+=1
        else:
            #print read[0], "\t", "na"
            cts[2]+=1
    return [cts,"-","-"]

def ct_type2(reads, union, TEID):
    '''Type 2 events are characterized by any other event where either one or both sides have evidence for a clipped position'''
    if union[7] == "-" and union[8] == "+":
        #return ct_type1FR(reads, int(union[2]), int(union[2]))
        return ct_type1R(reads, int(union[2]), int(union[2]), union, TEID)
    elif union[7] == "+" and union[8] == "-":
        #return ct_type1FR(reads, int(union[1]), int(union[1]))
        return ct_type1F(reads, int(union[1]), int(union[1]), union, TEID)
    else:
        print "ct_type2 error"

def countReads(reads, union, unionType, TEID, annotation):
    if unionType == 0:
        return ct_type0(reads, union, TEID, annotation)
    elif unionType == 1:
        return  ct_type1FR(reads, int(union[1]), int(union[2]), union, TEID)
    elif unionType == 2:
        return ct_type2(reads, union, TEID)

def unionType(union):
    if union[7] == "-" and union[8] == "-":
        return 0
    elif union[7] == "+" and union[8] == "+":
        return 1
    else:
        return 2

def countReads_portal(ID, union, samples, tmpDir, samFILE, hierarchy, label, l2, annotation):
    # make a dict of all the TEids that match at the cluster level
    TEID ={}
    groupIndex=label.index(l2)
    for h in hierarchy:
        if hierarchy[h][groupIndex]==union[ID][4]:
            TEID[h]=""
    # count the reads
    if not samFILE:
        return [[0,0,0],"-","-"]
    else:
        reads=[x.split("\t") for x in samFILE.split("\n")]
        del reads[-1]
        return countReads(reads, union[ID], unionType(union[ID]), TEID, annotation)

