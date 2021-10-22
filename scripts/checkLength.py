import glob
import sys

folder = sys.argv[1]
# fileList = glob.glob(folder + "/common*.fna")
fileList = sys.argv[1:]

fileHandlerList = list()

maxLine = 0
for line in open(fileList[0]):
    maxLine += 1

for file in fileList:
    fileHandlerList.append(open(file))

for i in range(1, maxLine):
    lineContent = [file.readline().strip() for file in fileHandlerList]
    if lineContent[0][0] != ">":
        if all(len(element) == len(lineContent[0]) for element in lineContent):
            # nothing to do
            pass
        else:
            print("Line length not ok")
            print(lineContent)
            raise Exception("Length not equal between files")
