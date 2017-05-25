import os

def write_bed_portal(hierarchy, label, group, level, bedDir):
    teIDs=[]
    groupIndex=label.index(level)
    for ID in hierarchy:
        if hierarchy[ID][groupIndex] == group:
            teIDs.append(ID)
    bedOutFile = os.path.join(bedDir, "%s_complete.bed" %(group))
    megaOutFile = os.path.join(bedDir, "mega_complete.bed")
    with open(bedOutFile, 'w') as fOUT:
        for x in sorted(teIDs):
            fOUT.write(x+'\t1\t1000000'+'\n')
    with open(megaOutFile, 'a') as fOUT:
        for x in sorted(teIDs):
            fOUT.write(x+'\t1\t1000000'+'\n')
