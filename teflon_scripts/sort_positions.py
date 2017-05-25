def sort(lines):
    sorted_lines,idTags=[],[]
    for i in xrange(len(lines)):
        if lines[i][1] == '-':
            pos=int(lines[i][2])+1
        else:
            pos=int(lines[i][1])
        idTags.append([lines[i][0],pos,i])
    sorted_idTags=sorted(idTags,key = lambda x: (x[0],x[1]))
    for i in xrange(len(sorted_idTags)):
        sorted_lines.append(lines[sorted_idTags[i][2]])
    return sorted_lines


def sort_portal(positions):
    elements=[]
    with open(positions, 'r') as fIN:
        for line in fIN:
            elements.append(line.split())
    sorted_elements=sort(elements)
    fileNAME=positions.replace('.txt','_sorted.txt')
    with open(fileNAME, 'w') as fOUT:
        for line in sorted_elements:
            fOUT.write('\t'.join([str(x) for x in line])+'\n')
    return fileNAME
