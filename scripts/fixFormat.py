import sys

wideFile = open(sys.argv[1])
longOutput = open(sys.argv[2], "w")

header = wideFile.readline()
header = header.strip().split('\t')
for line in wideFile:
    line = line.strip().split('\t')
    for base, count in zip(header[1:], line[1:]):
        if count != "0":
            longOutput.writelines([line[0], '\t', base, '\t', count, '\n']) # pos, base, count
