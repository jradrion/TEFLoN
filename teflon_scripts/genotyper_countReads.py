import os, sys
import multiprocessing as mp

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

def isReference(ch, pos, side, focal_annotation, readLen, insz, sd):
    '''does insertion support a reference te'''
    if side == 'F':
        nearby_tes=[]
        for i in xrange(len(focal_annotation)):
            if focal_annotation[i][0] == ch and pos-(readLen/2) <= focal_annotation[i][1] < pos+(insz+(2*sd)):
                nearby_tes.append([focal_annotation[i][3], focal_annotation[i][1]])
        sorted_nearby_tes=sorted(nearby_tes, key = lambda x: x[1])
    elif side == 'R':
        nearby_tes=[]
        focal_annotation.reverse()
        for i in xrange(len(focal_annotation)):
            if focal_annotation[i][0] == ch and pos+(readLen/2) >= focal_annotation[i][2] > pos-(insz+(2*sd)):
                nearby_tes.append([focal_annotation[i][3],focal_annotation[i][2]])
        sorted_nearby_tes=sorted(nearby_tes, key = lambda x: x[1], reverse=True)
        focal_annotation.reverse()
    if sorted_nearby_tes:
        return sorted_nearby_tes[0][0]
    else:
        return '-'

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

def idOrientation(flag):
    #orientation not accurate
    #this will effectively exclude orientaton until fixed
    return "."
    if bitFlag(int(flag))[4]==0 and bitFlag(int(flag))[5]==0:
        return '-'
    elif bitFlag(int(flag))[4]==1 and bitFlag(int(flag))[5]==1:
        return '-'
    else:
        return '+'

def isPrimary(read):
    if bitFlag(int(read[1]))[8]==0 and bitFlag(int(read[1]))[11]==0:
        return 1
    else:
        return 0


def ct_type0(reads, union, annotation):
    '''Type 0 events are characterized by non-clipped reads flanking both sides of an insertion'''
    #try to find the breakpoints from clipped reads not aligning to tes
    F_anno,R_anno="",""
    for te in annotation:
        if te[3] == union[6]:
            F_anno=int(te[1])
            R_anno=int(te[2])
    F,R,f_clip,r_clip="","",[],[]
    for read in reads:
        if "R" in side_clipped(read[5]):
            f_clip.append(rightMostPos(int(read[3]),cigarParse(read[5])))
        if "L" in side_clipped(read[5]):
            r_clip.append(int(read[3]))
    if f_clip:
        F=mode(f_clip)
        redefinedF=F
    if r_clip:
        R=mode(r_clip)
        redefinedR=R
    if F == "":
        if F_anno:
            F=F_anno
        else:
            F=int(union[1])
        #elif union[1].isdigit():
        #    F=int(union[1])
        #else:
        #    F=int(union[2])
    if R == "":
        if R_anno:
            R=R_anno
        else:
            R=int(union[2])
        #elif union[2].isdigit():
        #    R=int(union[2])
        #else:
        #    R=int(union[1])
    return ct_type1FR(reads, F, R)

def ct_type1FR(reads, F, R):
    '''Type 1 events are characterized by the F support breakpoint following the R support breakpoint (ie TSD)'''
    cts=[0,0,0]
    overShoot=3
    for read in reads:
        if isPrimary(read) == 1:
            if bitFlag(int(read[1]))[4] == 0 and rightMostPos(int(read[3]),cigarParse(read[5])) < F+overShoot and "FB" in read[6]: # change to incorporate given te chromosomes names from database
                #print read[0], "\t", "p"
                cts[0]+=1
            elif bitFlag(int(read[1]))[4] == 1 and int(read[3]) > R-overShoot and "FB" in read[6]: # change to incorporate given te chromosomes names from database
                #print read[0], "\t", "p"
                cts[0]+=1
            elif int(read[3]) <= R-overShoot and rightMostPos(int(read[3]),cigarParse(read[5])) >= F+overShoot:
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
    #return [cts, F, R]
    return cts

def ct_type1F(reads, F, R):
    """Since only the F pos is a known clip you should really only count F reads, this might be a shortcut accomplish this"""
    R=F-20 #20 used as a conservative estimate for the largest TSD
    '''Type 1 events are characterized by the F support breakpoint following the R support breakpoint (ie TSD)'''
    cts=[0,0,0]
    overShoot=3
    for read in reads:
        if isPrimary(read) == 1:
            if bitFlag(int(read[1]))[4] == 0 and rightMostPos(int(read[3]),cigarParse(read[5])) < F+overShoot and "FB" in read[6]: # change to incorporate given te chromosomes names from database
                #print read[0], "\t", "p"
                cts[0]+=1
            elif bitFlag(int(read[1]))[4] == 1 and int(read[3]) > R-overShoot and "FB" in read[6]: # change to incorporate given te chromosomes names from database
                #print read[0], "\t", "p"
                cts[0]+=1
            elif int(read[3]) <= R-overShoot and rightMostPos(int(read[3]),cigarParse(read[5])) >= F+overShoot:
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
    #return [cts, F, R]
    return cts


def ct_type1R(reads, F, R):
    """Same with R only"""
    F=R+20
    '''Type 1 events are characterized by the F support breakpoint following the R support breakpoint (ie TSD)'''
    cts=[0,0,0]
    overShoot=3
    for read in reads:
        if isPrimary(read) == 1:
            if bitFlag(int(read[1]))[4] == 0 and rightMostPos(int(read[3]),cigarParse(read[5])) < F+overShoot and "FB" in read[6]: # change to incorporate given te chromosomes names from database
                #print read[0], "\t", "p"
                cts[0]+=1
            elif bitFlag(int(read[1]))[4] == 1 and int(read[3]) > R-overShoot and "FB" in read[6]: # change to incorporate given te chromosomes names from database
                #print read[0], "\t", "p"
                cts[0]+=1
            elif int(read[3]) <= R-overShoot and rightMostPos(int(read[3]),cigarParse(read[5])) >= F+overShoot:
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
    #return [cts, F, R]
    return cts

def ct_type2(reads, union):
    """test using type1FR, type1F, and type1R
    if not producing better estimates go back to all statements returning type1FR"""
    '''Type 2 events are characterized by any other event where either one or both sides have evidence for a clipped position'''
    #if union[7] == "+" and union[8] == "+":
    #    return ct_type1FR(reads, int(union[1]), int(union[2]))
    if union[7] == "-" and union[8] == "+":
        #return ct_type1FR(reads, int(union[2]), int(union[2]))
        return ct_type1R(reads, int(union[2]), int(union[2]))
    elif union[7] == "+" and union[8] == "-":
        #return ct_type1FR(reads, int(union[1]), int(union[1]))
        return ct_type1F(reads, int(union[1]), int(union[1]))
    else:
        print "ct_type2 error"


def countReads(reads, union, unionType, annotation):
    if unionType == 0:
        return ct_type0(reads, union, annotation)
    elif unionType == 1:
        return  ct_type1FR(reads, int(union[1]), int(union[2]))
    elif unionType == 2:
        return ct_type2(reads, union)


def unionType(union):
    if union[7] == "-" and union[8] == "-":
        return 0
    elif union[7] == "+" and union[8] == "+":
        return 1
    else:
        return 2

def countReads_portal(ID, union, samples, tmpDir, annotation):
    #print "element:", union[ID], unionType(union[ID])
    counts=[]
    for i in range(len(samples)):
        reads=[]
        samFILE = os.path.join(tmpDir, "union_%s.%s.sam" %(str(ID),str(samples[i][1])))
        try:
            with open(samFILE, 'r') as fIN:
                try:
                    for line in fIN:
                        reads.append(line.split())
                except IOError:
                    continue
        except IOError:
            continue
        if not reads:
            counts.append([0,0,0]) #[# presence support reads, # absence support reads, # non-informative reads]
        else:
            counts.append(countReads(reads, union[ID], unionType(union[ID]), annotation))
    return counts


