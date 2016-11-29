import sys,random

def refTEset(refs):
    return ','.join(set(refs.split(',')))

def isGenotypable(lst):
    if lst[2] == "-" and lst[7] == "-":
        return 1
    if lst[1] == "-" and lst[8] == "-":
        return 1
    else:
        return lst

def totalReads(lst):
    reads=0
    if lst[9].isdigit():
        reads+=int(lst[9])
    if lst[10].isdigit():
        reads+=int(lst[10])
    return reads

def readThresh(lst, thresh):
    ct=0
    for x in lst:
        if x[9].isdigit():
            ct+=int(x[9])
        if x[10].isdigit():
            ct+=int(x[10])
    if ct >= thresh:
        return 0
    else:
        return 1

def randomFamilySelect(families):
    return families[random.randint(0,len(families))-1]

def collapseElements(lst, thresh):
    #print "========================================================================="
    #for x in lst:
    #    print x
    #return
    chrom,Fpos,Rpos,fam,cLev,orient,ref,Fclip,Rclip=lst[0][0],"-","-","-","-","-","-","-","-"
    if readThresh(lst, thresh) == 0:
        pass
    else:
        return ["-","-","-","-","-","-","-","-","-"]
    if len(lst) == 1:
        return [chrom,lst[0][1],lst[0][2],lst[0][3],lst[0][4],lst[0][5],lst[0][6],lst[0][7],lst[0][8]]
    else:
        Fpos_clp,Fpos_noClp = [],[]
        Rpos_clp,Rpos_noClp = [],[]
        famReads={}
        for i in range(len(lst)):
            if lst[i][7] != "-":
                Fpos_clp.append(int(lst[i][1]))
            elif lst[i][7] == "-" and lst[i][1].isdigit():
                Fpos_noClp.append(int(lst[i][1]))
            if lst[i][8] != "-":
                Rpos_clp.append(int(lst[i][2]))
            elif lst[i][8] == "-" and lst[i][2].isdigit():
                Rpos_noClp.append(int(lst[i][2]))
            if lst[i][6] != "-":
                fam,cLev,orient,ref=lst[i][3],lst[i][4],lst[i][5],lst[i][6]
            if lst[i][3] in famReads:
                famReads[lst[i][3]][0]+=totalReads(lst[i])
            else:
                famReads[lst[i][3]]=[totalReads(lst[i]),lst[i][4],lst[i][5]]
        if ref == "-":
            famSelect=[]
            readCounts=[]
            for x in famReads:
                readCounts.append(famReads[x][0])
            for x in famReads:
                if famReads[x][0] == max(readCounts):
                    famSelect.append(x)
            fam=randomFamilySelect(famSelect)
            cLev,orient=famReads[fam][1],famReads[fam][2]
        if Fpos_clp:
            Fpos=max(Fpos_clp)
            Fclip="+"
        else:
            if Fpos_noClp:
                Fpos=max(Fpos_noClp)
        if Rpos_clp:
            Rpos=min(Rpos_clp)
            Rclip="+"
        else:
            if Rpos_noClp:
                Rpos=min(Rpos_noClp)
        return [chrom,str(Fpos),str(Rpos),fam,cLev,orient,ref,Fclip,Rclip]


def addToCluster(elem, cluster, readLen, insz, sd):
    for i in xrange(len(cluster)):
        if cluster[i][6] != "-" and elem[6] != "-" and cluster[i][6] != elem[6]:
            return 0
        elif cluster[i][6] != "-" and elem[6] != "-" and cluster[i][6] == elem[6]:
            return 1
    for i in xrange(len(cluster)):
        if cluster[i][0] == elem[0] and cluster[i][1] != "-" and elem[1] != "-" and int(elem[1]) <= int(cluster[i][1])+readLen:
            return 1
        if cluster[i][0] == elem[0] and cluster[i][2] != "-" and elem[2] != "-" and int(elem[2]) <= int(cluster[i][2])+readLen:
            return 1
        if cluster[i][0] == elem[0] and cluster[i][2] != "-" and elem[1] != "-" and int(elem[1]) <= int(cluster[i][2])+20: #20 should be longer than nearly all TSDs
            return 1
    else:
        return 0


def passReadThresh(tmp,thresh):
    ct=0
    for i in range(len(tmp)):
        if tmp[i][9].isdigit() and tmp[i][10].isdigit():
            ct+= (int(tmp[i][9])+int(tmp[i][10]))
        elif tmp[i][9].isdigit() and not tmp[i][10].isdigit():
            ct+=int(tmp[i][9])
        elif tmp[i][10].isdigit() and not tmp[i][9].isdigit():
            ct+=int(tmp[i][10])
    if ct >= thresh:
        return 1
    else:
        return 0

def collapse(oldlst, readLen, insz, sd, thresh):
    #filter out noisy calls. This method is crude and should be refined
    lst=[]
    for x in oldlst:
        if passReadThresh([x],thresh) == 1:
            lst.append(x)
    #cluster calls
    genotypableElements=[]
    nonGenotypableElements=[]
    collapsedElements=[]
    preCollapse=[]
    while lst:
        print "Clustering iterations remaining:", len(lst)
        newLst=[]
        tmp=[]
        for i in xrange(len(lst)):
            #print lst[i]
            if not tmp:
                tmp.append(lst[i])
            else:
                if i < 500:
                    #newLst.append(lst[i])
                    if addToCluster(lst[i], tmp, readLen, insz, sd) == 1:
                        tmp.append(lst[i])
                    else:
                        newLst.append(lst[i])
                else:
                    newLst.append(lst[i])
        if passReadThresh(tmp,thresh) == 1:
            preCollapse.append(tmp)
            lst=newLst
        else:
            lst=newLst

    for x in preCollapse:
        collapsedElements.append(collapseElements(x, thresh))
    for x in collapsedElements:
        if isGenotypable(x) != 1:
            genotypableElements.append(x)
        else:
            nonGenotypableElements.append(x)
    return [genotypableElements,nonGenotypableElements]

def collapse_union_portal(inFILE, readLen, insz, sd, thresh):
    elements, collapsed = [],[]
    with open(inFILE, "r") as fIN:
        for line in fIN:
            elements.append(line.split())
    collapsed = collapse(elements, readLen, insz, sd, thresh)[0]
    with open(inFILE.replace(".txt",".collapsed.txt"), "w") as fOUT:
        for line in collapsed:
            fOUT.write("\t".join([str(x) for x in line])+"\n")

    return collapsed


