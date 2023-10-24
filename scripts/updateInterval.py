import sys

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

originalInt = dict() # will be dict of lists
print("File " + originalS + " : getting original ints")

for line in humanHandler:
    line = line.strip().split("\t") # line is chr : start : stop : ID : . : strand
    line[1:3] = [int(value) for value in line[1:3]]
    if line[0] not in originalInt:
        originalInt[line[0]] = list()
    originalInt[line[0]].append((line[1], line[2], line[3]))

diffInt = dict() # will be ID : diffStart, diffStop
count = 0

print("File " + originalS + " : getting reductions")

for line in commonHandler: # line is chr : start : stop : . : . : strand
    line = line.strip().split("\t")
    count += 1
    line[1:3] = [int(value) for value in line[1:3]]

    A = 0
    B = len(originalInt[line[0]])
    found = False
    stop = False

    while(found == False and stop == False):
        index = int((A + B)/2)
        inter = originalInt[line[0]][index] # is original START : STOP
        # original start is <= than new start
        # original stop is >= than new stop
        if inter[0] <= line[1]: # START smaller or equal to current
            if inter[1] >= line[2]:
                new_id = str(inter[2]) + "|" + str(count)
                if new_id in diffInt:
                    print(line[3])
                    raise Exception("IDs are not unique !!!")
                diffInt[new_id] = (line[1] - inter[0],  inter[1] - line[2], count)
                found = True
            else: # START is smaller AND STOP is smaller
                A = index + 1
        else: # START bigger than current
            B = index - 1
        if (B - A) < 0:
            stop = True

originalInt = dict()
lines = list()
print("File " + originalS + " : getting species int")
for line in speciesHandler: # line is chr : start : stop : ID : . : strand
    line = line.strip().split("\t")
    line[1:3] = [int(value) for value in line[1:3]]
    if line[3] in originalInt:
        print(line[3])
        raise Exception("Duplicated IDs !!!")
    originalInt[line[3]] = line

print("File " + originalS + " : getting and writing diffs")
for key in diffInt:
    originalKey = key.split("|")[0]
    line = originalInt[originalKey]
    if line[5] == "+":
        start = int(line[1]) + diffInt[key][0]
        end = int(line[2]) - diffInt[key][1]
    else:
        start = int(line[1]) + diffInt[key][1]
        end = int(line[2]) - diffInt[key][0]
    lines.append((line[0] + "\t" + str(start) + "\t" + str(end) + "\t" + key + "\t" + "\t".join(line[-2:]) + "\n", diffInt[key][2]))

lines.sort(key = lambda x: x[1]) # sort following line order in common bed
lines = [value[0] for value in lines] # delete line index now that the line is sorted
outputHandler.writelines(lines)
