import sys
import re

input = open(sys.argv[1])
CPGoutput = open(sys.argv[2], "w")
nCPGoutput = open(sys.argv[3], "w")
header_CPG = set()
header_nCPG = set()

for line in input:
    if line.startswith(">"):
        header = line.strip()
        reg = re.findall(">(\\w+):(\\d+)-\\d+", line)[0]
        chrom = reg[0]
        start = int(reg[1])
        before_char = ""
    else:
        countChar = 0 # count nb of chars in line
        startCPG = 0
        startnCPG = 0
        for char in line:
            char = char.upper()
            countChar += 1
            if before_char in ("J", "K", "L"):
                if char in ("F", "H", "I"):
                    if header in header_CPG:
                        gapSize = countChar - startCPG
                        for i in range(0, gapSize):
                            CPGoutput.write("-")
                            
                        CPGoutput.write(before_char + char)
                    else:
                        header_CPG.add(header)
                        startCPG = start + countChar - 1 # bc we have a "before" char
                        CPGoutput.write(">" + chrom + ":" + str(startCPG) + "\n")
                        CPGoutput.write(before_char + char)
            else:
                if header in header_nCPG:
                    gapSize = countChar - startnCPG
                    for i in range(0, gapSize):
                        nCPGoutput.write("-")
                            
                    nCPGoutput.write(before_char + char)
                else:
                    header_nCPG.add(header)
                    startnCPG = start + countChar - 1
                    CPGoutput.write(">" + chrom + ":" + str(startnCPG) + "\n")
                    CPGoutput.write(before_char + char)
            before_char = char
