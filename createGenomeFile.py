#!/usr/bin/env python3
import sys
import collections
header = ""
count = 0
oH = open(sys.argv[2], "w")
output_content = dict()

for line in open(sys.argv[1]):
    if not line.startswith(">"):
        # count += len(line)
        count += len(line.strip())
    else:
        if header != "":
            output_content[header] = count
            count = 0
        header = line[1:].strip().split(" ")[0]

output_content[header] = count
output_content = collections.OrderedDict(sorted(output_content.items()))
oH.write("\n".join(["\t".join(x) for x in [(key, str(output_content[key])) for key in output_content]]))