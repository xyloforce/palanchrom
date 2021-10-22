import glob
import sys

folder = sys.argv[1]
output = sys.argv[2]

print("conversion started, output will be : " + output)

fileList = glob.glob(folder + "/*.dat")
outHandler = open(output, "w")

for file in fileList:
    chromH = file[:-4].split("/")[-1] # keep only filename
    datHandler = open(file)
    for line in datHandler.readlines():
        line = line.split("\t")
        startH = line[0]
        endH = line[1]
        outHandler.writelines(chromH + "\t" + startH +  "\t" + endH + "\n")
