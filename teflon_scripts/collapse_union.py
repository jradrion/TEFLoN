import sys,random

def progress_bar(percent, barLen = 50):
    sys.stdout.write("\r")
    progress = ""
    for i in range(barLen):
        if i < int(barLen * percent):
            progress += "="
        else:
            progress += " "
    sys.stdout.write("[ %s ] %.2f%%" % (progress, percent * 100))
    sys.stdout.flush()

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
    chrom,Fpos,Rpos,fam,cLev,orient,ref,Fclip,Rclip,Fct,Rct=lst[0][0],"-","-","-","-","-","-","-","-","-","-"
    if readThresh(lst, thresh) == 0:
        pass
    else:
        return ["-","-","-","-","-","-","-","-","-","-","-"]
    if len(lst) == 1:
        return [chrom,lst[0][1],lst[0][2],lst[0][3],lst[0][4],lst[0][5],lst[0][6],lst[0][7],lst[0][8],lst[0][9],lst[0][10]]
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
        Fct_tmp,Rct_tmp=0,0
        for i in range(len(lst)):
            if lst[i][3] == fam:
                if lst[i][9].isdigit():
                    Fct_tmp+=int(lst[i][9])
                if lst[i][10].isdigit():
                    Rct_tmp+=int(lst[i][10])
        if Fct_tmp>0:
            Fct=Fct_tmp
        if Rct_tmp>0:
            Rct=Rct_tmp
        return [chrom,str(Fpos),str(Rpos),fam,cLev,orient,ref,Fclip,Rclip,Fct,Rct]


#def addToGroup(elem, cluster, readLen, insz):
#    print "elem:",elem,":elem"
#    print "readLen;",readLen,":readLen"
#    sys.exit()
#    maxTSDLen=20
#    for i in xrange(len(cluster)):
#        if cluster[i][6] != "-" and elem[6] != "-" and cluster[i][6] != elem[6]:
#            return 0
#        if cluster[i][6] != "-" and elem[6] != "-" and cluster[i][6] == elem[6]:
#            return 1
#        if cluster[i][0] == elem[0] and cluster[i][1] != "-" and elem[1] != "-" and int(elem[1]) <= int(cluster[i][1])+readLen and cluster[i][4]==elem[4]:
#            return 1
#        if cluster[i][0] == elem[0] and cluster[i][2] != "-" and elem[2] != "-" and int(elem[2]) <= int(cluster[i][2])+readLen and cluster[i][4]==elem[4]:
#            return 1
#        if cluster[i][0] == elem[0] and cluster[i][2] != "-" and elem[1] != "-" and int(elem[1]) <= int(cluster[i][2])+maxTSDLen and cluster[i][4]==elem[4]:
#            return 1
#    else:
#        return 0

def addToGroup(elem, cluster, readLen, insz):
    #print "cluster:",cluster,":cluster"
    #print "elem:",elem,":elem"
    #print "readLen;",readLen,":readLen"
    #sys.exit()
    maxTSDLen=20
    for i in xrange(len(cluster)):
        if cluster[i][6] != "-" and elem[6] != "-" and cluster[i][6] != elem[6]:
            #print "if0"
            return 0
        if cluster[i][6] != "-" and elem[6] != "-" and cluster[i][6] == elem[6]:
            #print "if1"
            return 1
        if cluster[i][0] == elem[0] and cluster[i][1] != "-" and elem[1] != "-" and int(elem[1]) <= int(cluster[i][1])+readLen and cluster[i][4]==elem[4]:
            #print "if2"
            return 1
        if cluster[i][0] == elem[0] and cluster[i][2] != "-" and elem[1] == "-" and  elem[2] != "-" and int(elem[2]) <= int(cluster[i][2])+readLen and cluster[i][4]==elem[4]:
            #print "if3"
            return 1
        if cluster[i][0] == elem[0] and cluster[i][2] != "-" and elem[1] != "-" and int(elem[1]) <= int(cluster[i][2])+maxTSDLen and cluster[i][4]==elem[4] and elem[6] == "-" and cluster[i][6] == "-":
            #print "if4"
            return 1
    else:
        return 0

def grouper(lst,readLen,insz,uniqueGroups):
    group=[]
    ct=1
    for item in lst:
        if not group:
            group.append([item])
        else:
            groupsBack=min(uniqueGroups,len(group))
            flag=0
            for i in range(groupsBack):
                index=(i+1)*-1
                if addToGroup(item,group[index],readLen,insz)==1:
                    group[index].append(item)
                    flag=1
                    break
            if flag==0:
                group.append([item])
        progress_bar(ct/float(len(lst)))
        ct+=1
    #print "group:",group,":group"
    return group


def collapse(lst, readLen, insz, thresh, uniqueGroups):
    #cluster calls
    genotypableElements=[]
    nonGenotypableElements=[]
    collapsedElements=[]
    preCollapse=dict(enumerate(grouper(lst,readLen,insz,uniqueGroups),0))
    for i in range(len(preCollapse)):
        collapsedElements.append(collapseElements(preCollapse[i], thresh))

    #print collapsedElements
    #sys.exit()

    for x in collapsedElements:
        if isGenotypable(x) != 1:
            genotypableElements.append(x)
        else:
            nonGenotypableElements.append(x)
    return [genotypableElements,nonGenotypableElements]

def collapse_union_portal(inFILE, readLen, insz, thresh):
    groups=[]
    elements, collapsed = [],[]
    with open(inFILE, "r") as fIN:
        for line in fIN:
            ar=line.split()
            elements.append(ar)
            groups.append(ar[3])
    uniqueGroups=len(set(groups))
    collapsed = collapse(elements, readLen, insz, thresh, uniqueGroups)[0]
    collapsedFILE=inFILE.replace(".txt",".collapsed.txt")
    with open(collapsedFILE, "w") as fOUT:
        for line in collapsed:
            fOUT.write("\t".join([str(x) for x in line])+"\n")

    return collapsedFILE


