import sys

header = ""
count = 0
oH = open(sys.argv[2], "w")

for line in open(sys.argv[1]):
    if not line.startswith(">"):
        count += len(line)
    else:
        if header != "":
            oH.write("\t".join((header, str(count) + "\n")))
            count = 0
        header = line[1:].strip()
