# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 11:01:34 2016

@author: jeffreyadrion

convert the full bed of all TEs to a bed of TEs for a group of interest
"""
import os

def write_bed_portal(hierarchy, label, group, level, bedDir):
    teIDs=[]
    groupIndex=label.index(level)
    for ID in hierarchy:
        if hierarchy[ID][groupIndex] == group:
            teIDs.append(ID)
    bedOutFile = os.path.join(bedDir, "%s_complete.bed" %(group))
    with open(bedOutFile, 'w') as fOUT:
        for x in sorted(teIDs):
            fOUT.write(x+'\t1\t1000000'+'\n')
    return sorted(teIDs)
