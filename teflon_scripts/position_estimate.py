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

def isReference(ch, pos, side, focal_annotation, readLen, insz, sd, clustGroup, hierarchy):
    '''does insertion support a reference te'''
    if side == 'F':
        nearby_tes=[]
        for i in xrange(len(focal_annotation)):
            if focal_annotation[i][0] == ch and pos-(readLen/4) <= focal_annotation[i][1] < pos+(insz+(2*sd)):
                nearby_tes.append([focal_annotation[i][3], focal_annotation[i][1]])
        sorted_nearby_tes=sorted(nearby_tes, key = lambda x: x[1])
    elif side == 'R':
        nearby_tes=[]
        focal_annotation.reverse()
        for i in xrange(len(focal_annotation)):
            if focal_annotation[i][0] == ch and pos+(readLen/4) >= focal_annotation[i][2] > pos-(insz+(2*sd)):
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

def idOrientation(flag):
    if bitFlag(int(flag))[4]==0 and bitFlag(int(flag))[5]==0:
        return '-'
    elif bitFlag(int(flag))[4]==1 and bitFlag(int(flag))[5]==1:
        return '-'
    else:
        return '+'

def te_support(clustered_sam,focal_annotation,readLen,insz,sd, clustGroup, hierarchy):
    '''return support for tes'''
    arr=clustered_sam.split()
    ch=arr[2]
    if bitFlag(int(arr[1]))[4]==0:
        strand='+'
    else:
        strand='-'
    if strand== '+' and 'S' in arr[5] and 'R' in side_clipped(arr[5]):
        clip='+'
    elif strand== '-' and 'S' in arr[5] and 'L' in side_clipped(arr[5]):
        clip='+'
    else:
        clip='-'
    if strand=='+':
        support_side='F'
    if strand=='-':
        support_side='R'
    if support_side=='F':
        pos=rightMostPos(int(arr[3]),cigarParse(arr[5]))
    if support_side=='R':
        pos=int(arr[3])
    ref_tes=isReference(ch, pos, support_side, focal_annotation, readLen, insz, sd, clustGroup, hierarchy)
    orientation=idOrientation(arr[1])
    return [ch,pos,support_side,clip,ref_tes,orientation]

def cluster(data, maxgap):
    '''Arrange data into groups where successive elements differ by no more than maxgap'''
    data.sort()
    groups = [[data[0]]]
    for x in data[1:]:
        if abs(x - groups[-1][-1]) <= maxgap:
            groups[-1].append(x)
        else:
            groups.append([x])
    return groups

def sharedTEs(te_cluster,potential):
    cluster_tes=[]
    for i in xrange(len(te_cluster)):
        cluster_tes.append(te_cluster[i][4])
    if potential[4] in cluster_tes:
        return 1
    else:
        return


def combineClusters(events, group, clustGroup):
    combinedEvents=[]
    for i in xrange(len(events)):
        combinedEvent=['-','-','-',group, clustGroup, "." , '-','-','-','-','-']
        if len(events[i]) <= 2:
            ref=[]
            for j in xrange(len(events[i])):
                if events[i][j][6]=="+":
                    combinedEvent[5]="+"
                if events[i][j][6]=="-":
                    combinedEvent[5]="-"
                if events[i][j][2] == 'F':
                    ref.append(events[i][j][5])
                    combinedEvent[0]=events[i][j][0]
                    combinedEvent[1]=events[i][j][1]
                    combinedEvent[7]=events[i][j][3]
                    combinedEvent[9]=events[i][j][4]
                if events[i][j][2] == 'R':
                    ref.append(events[i][j][5])
                    combinedEvent[0]=events[i][j][0]
                    combinedEvent[2]=events[i][j][1]
                    combinedEvent[8]=events[i][j][3]
                    combinedEvent[10]=events[i][j][4]
            if ref:
                if len(set(ref))>1:
                    combinedEvent[6]=ref[0]+','+ref[1]
                else:
                    combinedEvent[6]=ref[0]
            else:
                break
        else:
            F,R=[],[]
            ref=[]
            for j in xrange(len(events[i])):
                if events[i][j][6]=="+":
                    combinedEvent[5]="+"
                if events[i][j][6]=="-":
                    combinedEvent[5]="-"
                if events[i][j][2] == 'F':
                    F.append(events[i][j])
                if events[i][j][2] == 'R':
                    R.append(events[i][j])
            if F:
                clipped,notClipped=[],[]
                for j in xrange(len(F)):
                    if F[j][3] == '+':
                        clipped.append(F[j])
                    else:
                        notClipped.append(F[j])
                if clipped:
                    sorted_clipped=sorted(clipped, key = lambda x: x[1])
                    ref.append(sorted_clipped[-1][5])
                    combinedEvent[0]=sorted_clipped[-1][0]
                    combinedEvent[1]=sorted_clipped[-1][1]
                    combinedEvent[7]=sorted_clipped[-1][3]
                    combinedEvent[9]=sorted_clipped[-1][4]
                else:
                    sorted_notClipped=sorted(notClipped, key = lambda x: x[1])
                    ref.append(sorted_notClipped[-1][5])
                    combinedEvent[0]=sorted_notClipped[-1][0]
                    combinedEvent[1]=sorted_notClipped[-1][1]
                    combinedEvent[7]=sorted_notClipped[-1][3]
                    combinedEvent[9]=sorted_notClipped[-1][4]
            if R:
                clipped,notClipped=[],[]
                for j in xrange(len(R)):
                    if R[j][3] == '+':
                        clipped.append(R[j])
                    else:
                        notClipped.append(R[j])
                if clipped:
                    sorted_clipped=sorted(clipped, key = lambda x: x[1])
                    ref.append(sorted_clipped[0][5])
                    combinedEvent[0]=sorted_clipped[0][0]
                    combinedEvent[2]=sorted_clipped[0][1]
                    combinedEvent[8]=sorted_clipped[0][3]
                    combinedEvent[10]=sorted_clipped[0][4]
                else:
                    sorted_notClipped=sorted(notClipped, key = lambda x: x[1])
                    ref.append(sorted_notClipped[0][5])
                    combinedEvent[0]=sorted_notClipped[0][0]
                    combinedEvent[2]=sorted_notClipped[0][1]
                    combinedEvent[8]=sorted_notClipped[0][3]
                    combinedEvent[10]=sorted_notClipped[0][4]
            if len(set(ref))>1:
                combinedEvent[6]=ref[0]+','+ref[1]
            else:
                combinedEvent[6]=ref[0]
        combinedEvents.append(combinedEvent)
    return combinedEvents

def matchOrient(old, new):
    if old == new:
        return 1
    else:
        if old == "." or new == ".":
            return 1
        else:
            return 0

def clusterClusters(old, readLen):
    old=sorted(old, key = lambda x: (x[0],x[1]))
    new_cluster=[]
    temp=[]
    for i in xrange(len(old)):
        if temp:
            if old[i][5] != '-' and old[i][5] == temp[-1][5] and matchOrient(old[i][6],temp[-1][6]) == 1 :
                temp.append(old[i])
            elif old[i][0] == temp[-1][0] and old[i][2] == temp[-1][2] and old[i][1] <= temp[-1][1]+readLen and matchOrient(old[i][6],temp[-1][6]) == 1:
                if old[i][5]=="-" or temp[-1][5]=="-":
                    temp.append(old[i])
            elif old[i][0] == temp[-1][0] and old[i][5] == '-' and temp[-1][5] == '-'and old[i][2] == 'R' and temp[-1][2] =='F' and old[i][1] <= temp[-1][1]+(2*readLen) and matchOrient(old[i][6],temp[-1][6]) == 1:
                temp.append(old[i])
            elif old[i][0] == temp[-1][0] and old[i][5] == '-' and temp[-1][5] == '-'and old[i][2] == 'F' and temp[-1][2] == 'R' and old[i][1] <= temp[-1][1]+20 and matchOrient(old[i][6],temp[-1][6]) == 1:
                temp.append(old[i])
            else:
                new_cluster.append(temp)
                temp=[old[i]]
        else:
            temp.append(old[i])
    new_cluster.append(temp)
    return new_cluster

def orient(cluster):
    plus,minus=0,0
    for x in cluster:
        if len(x) == 6:
            if x[5] == "+":
                plus+=1
            elif x[5] == "-":
                minus+=1
    if plus>minus:
        return "+"
    elif minus>plus:
        return "-"
    else:
        return "."

def clusterReads(calls, readLen, insz, sd, cov):
    clustered_calls,clustered_index,new_cluster,start=[],[],[],0
    while len(clustered_index) < len(calls):
        count=[]
        if start + (5*cov) <= len(calls):
            end=start + (5*cov)
        else:
            end=len(calls)
        for i in xrange(start,end):
            if new_cluster and not i in clustered_index:
                if calls[i][0] == new_cluster[-1][0] and calls[i][2] == new_cluster[-1][2] and calls[i][1] < new_cluster[-1][1]+(readLen+sd) and sharedTEs(new_cluster,calls[i]): #new_cluster[-1][1]+(insz+(2*sd)-readLen) and sharedTEs(new_cluster,calls[i]):
                    new_cluster.append(calls[i])
                    clustered_index.append(i)
                else:
                    count.append(i)
            if not new_cluster and not i in clustered_index:
                new_cluster=[calls[i]]
                clustered_index.append(i)
        clustered_calls.append(new_cluster)
        new_cluster=[]
        if count:
            start=min(count)
        else:
            start=end

    clusters=[]
    for i in xrange(len(clustered_calls)):
        clustered,clip_pos,est_pos=[],[],[]
        for j in xrange(len(clustered_calls[i])):
            if clustered_calls[i][j][3] == '+':
                clip_pos.append(clustered_calls[i][j][1])
            else:
                est_pos.append(clustered_calls[i][j][1])
        if clustered_calls[i][0][2] == 'F':
            if clip_pos and est_pos:
                ct=0
                for x in est_pos:
                    if x > mode(clip_pos):
                        ct+=1
                if mode(clip_pos) >= max(est_pos) or max(est_pos) > mode(clip_pos) >= max(est_pos)-(readLen/2):
                    p=mode(clip_pos)
                    clustered=[clustered_calls[i][0][0],p,clustered_calls[i][0][2],'+',len(clustered_calls[i]),clustered_calls[i][0][4],orient(clustered_calls[i])]
                else:
                    p=max(est_pos)
                    clustered=[clustered_calls[i][0][0],p,clustered_calls[i][0][2],'-',len(clustered_calls[i]),clustered_calls[i][0][4],orient(clustered_calls[i])]
            elif clip_pos and not est_pos:
                p=mode(clip_pos)
                clustered=[clustered_calls[i][0][0],p,clustered_calls[i][0][2],'+',len(clustered_calls[i]),clustered_calls[i][0][4],orient(clustered_calls[i])]
            else:
                p=max(est_pos)
                clustered=[clustered_calls[i][0][0],p,clustered_calls[i][0][2],'-',len(clustered_calls[i]),clustered_calls[i][0][4],orient(clustered_calls[i])]
        else:
            if clip_pos and est_pos:
                if mode(clip_pos) <= min(est_pos) or min(est_pos) < mode(clip_pos) <= min(est_pos)+(readLen/2):
                    p=mode(clip_pos)
                    clustered=[clustered_calls[i][0][0],p,clustered_calls[i][0][2],'+',len(clustered_calls[i]),clustered_calls[i][0][4],orient(clustered_calls[i])]
                else:
                    p=min(est_pos)
                    clustered=[clustered_calls[i][0][0],p,clustered_calls[i][0][2],'-',len(clustered_calls[i]),clustered_calls[i][0][4],orient(clustered_calls[i])]
            elif clip_pos and not est_pos:
                p=mode(clip_pos)
                clustered=[clustered_calls[i][0][0],p,clustered_calls[i][0][2],'+',len(clustered_calls[i]),clustered_calls[i][0][4],orient(clustered_calls[i])]
            else:
                p=min(est_pos)
                clustered=[clustered_calls[i][0][0],p,clustered_calls[i][0][2],'-',len(clustered_calls[i]),clustered_calls[i][0][4],orient(clustered_calls[i])]
        clusters.append(clustered)
    return clusters


def position_estimate_portal(samFile_clust, suppFile, group, annotation, teIDs, readLen, insz, sd, cov, posDir, clustGroup, hierarchy):
    # read annotation
    focal_annotation=[]
    for i in xrange(len(annotation)):
        if annotation[i][3] in teIDs:
            focal_annotation.append(annotation[i])
    # read primary support
    primary_support=[]
    with open(samFile_clust, 'r') as fIN:
        for line in fIN:
            arr=line.split()
            bFlags=bitFlag(int(arr[1]))
            if arr[6] in teIDs and bFlags[8]==0 and bFlags[9]==0 and bFlags[10]==0 and bFlags[11]==0: #count only primary alignments whose mate is a te as support
                primary_support.append(te_support(line,focal_annotation,readLen,insz,sd, clustGroup, hierarchy))
    # read supplementary support
    supp_support=[]
    with open(suppFile, 'r') as fIN:
        for line in fIN:
            supp_support.append(line.split())
    for x in supp_support:
        primary_support.append([x[0],int(x[1]),x[2],"+",isReference(x[0], int(x[1]), x[2], focal_annotation, readLen, insz, sd, clustGroup, hierarchy),"."])
    # sort by position
    sorted_reads=sorted(primary_support, key = lambda x: (x[0],x[1]))
    # combine and cluster
    clusteredCalls=combineClusters(clusterClusters(clusterReads(sorted_reads, readLen, insz, sd, cov), readLen), group, clustGroup)
    # write position estimates
    outFile = os.path.join(posDir, "%s_positions.txt" %(group))
    with open(outFile, 'w') as fOUT:
        for line in clusteredCalls:
            fOUT.write('\t'.join([str(x) for x in line])+'\n')




