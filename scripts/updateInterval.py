import sys
import glob

# folder = sys.argv[1]
originalH = sys.argv[1]
originalS = sys.argv[2]
intersected = sys.argv[3]
outputFile = sys.argv[4]

# fileList = glob.glob(folder + "/*.dat")

humanHandler = open(originalH)
speciesHandler = open(originalS)
commonHandler = open(intersected)

outputHandler = open(outputFile, "w")

# commonDict = dict()
# old_file = ""

# TODO !!!ATTENTION!!! "-" : intervalle dans un sens sur l'humain & dans l'autre sur le chimpanzé : corriger problèmes
# if - : startH = stopS // stopH = startS

diffInt = dict() # will be ID : diffStart, diffStop
commonInt = dict() # will be chrom : list(start, stop)

print("File " + originalS + " : getting common ints")
for line in commonHandler:
    line = line.strip().split("\t") # line is chr : start : stop : uselessID : . : strand
    line[1:3] = [int(value) for value in line[1:3]]
    if line[0] not in commonInt:
        commonInt[line[0]] = list()
    commonInt[line[0]].append((line[1], line[2]))

for key in commonInt:
    commonInt[key].sort(key = lambda x:x[0])

print("File " + originalS + " : getting reductions")

for line in humanHandler: # line is chr : start : stop : ID : . : strand
    line = line.strip().split("\t")
    line[1:3] = [int(value) for value in line[1:3]]

    A = 0
    B = len(commonInt[line[0]])
    found = False
    stop = False

    while(found == False and stop == False):
        index = int((A + B)/2)
        inter = commonInt[line[0]][index] # is START : STOP
        if inter[0] >= line[1]: # START smaller or equal to current
            if inter[1] <= line[2]:
                if line[3] in diffInt:
                    print(line[3])
                    raise Exception("IDs are not unique !!!")
                diffInt[line[3]] = (inter[0] - line[1], line[2] - inter[1])
                found = True
            else: # START is smaller AND STOP is smaller
                A = index
        else: # START bigger than current
            B = index
        if not B-A > 1:
            stop = True

    # for inter in commonInt[line[0]]:
    #     if inter[0] <= line[1]:
    #         if inter[1] >= line[2]:
    #             if line[3] in diffInt:
    #                 raise("IDs are not unique !!!")
    #             diffInt[line[3]] = (line[1] - inter[0], inter[1] - line[2])
    #             break

print("File " + originalS + " : reading diffs and writing results")
for line in speciesHandler: # line is chr : start : stop : ID : . : strand
    line = line.strip().split("\t")
    line[1:3] = [int(value) for value in line[1:3]]
    if line[3] in diffInt:
        if line[5] == "+":
            print(diffInt[line[3]])
            start = int(line[1]) + diffInt[line[3]][0]
            end = int(line[2]) - diffInt[line[3]][1]
        else:
            start = int(line[1]) + diffInt[line[3]][1]
            end = int(line[2]) - diffInt[line[3]][0]
        outputHandler.write(line[0] + "\t" + str(start) + "\t" + str(end) + "\t" + line[3] + "\t" + "\t".join(line[-2:]) + "\n")