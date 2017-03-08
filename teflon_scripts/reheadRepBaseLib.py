def reheadRepBaseLib_portal(libIN,libOUT):
    print "Reheading repBase library:",libIN
    dupCheck=[]
    with open(libIN, "r") as fIN, open(libOUT, "w") as fOUT:
        for line in fIN:
            if line[0] == '>':
                 arr = line.split("\t")
                 if len(arr[1].split()) > 1:
                     tmp="_".join([str(x) for x in arr[1].split()])
                 else:
                     tmp=arr[1]
                 name = arr[0][1:] + "|" + tmp
                 dupCheck.append(name)
                 fOUT.write(">"+name+"\n")
            else:
                fOUT.write(line)


    if len(dupCheck) != len(set(dupCheck)):
        print "duplicate found!"

