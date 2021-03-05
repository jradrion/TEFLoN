import os

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
        if int(flag) & 2**i:
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

def cluster(data, maxgap):
    '''Arrange data into groups where successive elements differ by no more than maxgap'''
    sortedData=sorted(data)
    groups = [[sortedData[0]]]
    for x in sortedData[1:]:
        if abs(x - groups[-1][-1]) <= maxgap:
            groups[-1].append(x)
        else:
            groups.append([x])
    return groups

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

def parseSuppAlign(line, chromosomes, qual):
	if not 'SA:Z:' in line:
		return ['-','-','-']
	if line.split()[6] in chromosomes:
		return ['-','-','-']
	cig=line.split()[5]
	arr=line.split('SA:Z:')[-1].split(';')[0].split(',')
	support_side = ""
	if int(arr[4]) < qual:
		return ['-','-','-']
	if side_clipped(arr[3]) == 'R':
		support_side='F'
	if side_clipped(arr[3]) == 'L':
		support_side='R'
	if side_clipped(arr[3]) == 'LR':
		if side_clipped(cig) == 'R':
			support_side='F'
		elif side_clipped(cig) == 'L':
			support_side='R'
		else:
			return ['-','-','-']
	pos = 0
	if support_side=='F':
		pos=rightMostPos(arr[1],cigarParse(arr[3]))
	if support_side=='R':
		pos=int(arr[1])
	return[arr[0],pos,support_side,line.split()[0]]

def cluster_positions_portal(sam, group, chromosomes, lengths, readLen, insz, sd, bedDir, qual, suppOutFile):
    # record position of all paired ends aligning to TEs
    leftpos=[]
    with open(sam, 'r') as fIN, open(suppOutFile, 'w') as fOUT:
        for line in fIN:
            arr=line.split()
            bFlags=bitFlag(int(arr[1]))
            if bFlags[8] == 1 or bFlags[11] == 1:
                if bFlags[9] == 0 and bFlags[10] == 0 and parseSuppAlign(line, chromosomes, qual)[0] in chromosomes:
                    fOUT.write('\t'.join([str(x) for x in parseSuppAlign(line, chromosomes, qual)])+'\n')
            if arr[6] in chromosomes:
                leftpos.append([arr[6],int(arr[7])])
        lpos=sorted(leftpos, key = lambda x: (x[0],x[1]))
        # cluster positions
        clusters=[]
        for i in range(len(chromosomes)):
            temp=[]
            for j in xrange(len(lpos)):
                if lpos[j][0] == chromosomes[i]:
                    temp.append(lpos[j][1])
            if temp:
                clusters.append(cluster(temp, (insz+(2*sd))-readLen))
            else:
                clusters.append([])

    # write positions
    bedOutFile = os.path.join(bedDir, "%s_clustered.bed" %(group))
    with open(bedOutFile, 'w') as fOUT:
        for i in xrange(len(clusters)):
            for j in xrange(len(clusters[i])):
                fOUT.write(chromosomes[i]+'\t'+str(clusters[i][j][0])+'\t'+str(clusters[i][j][-1]+readLen)+'\n')
