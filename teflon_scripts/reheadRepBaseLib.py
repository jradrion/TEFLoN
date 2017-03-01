def reheadRepBaseLib_portal(lib):
    dupCheck=[]
    with open(lib, "r") as fIN, open(lib+".rehead", "w") as fOUT:
        for line in fIN:
            if line[0] == '>':
                 arr = line.split()
                 name = arr[0][1:] + "|" + arr[1]
                 dupCheck.append(name)
                 fOUT.write(">"+name+"\n")
            else:
                fOUT.write(line)


    if len(dupCheck) != len(set(dupCheck)):
        print "duplicate found!"

