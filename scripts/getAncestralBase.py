import sys
import glob
import re
import tqdm

# input is ref1 ref2 outgroups
fileList = sys.argv[1:-1]

fileHandlerList = list()
fileHandlerList = [open(file) for file in fileList]

outputHandler = open(sys.argv[-1], "w")

maxLine = 0
for line in open(fileList[0]):
    maxLine += 1
    
for i in tqdm.tqdm(range(maxLine), desc = "Working..."):
    
    lineContent = [file.readline().strip() for file in fileHandlerList] # lineContent is list of lines
    
    if lineContent[0].startswith(">"):
        reg = re.findall(">(\w+):(\d+)-\d+", lineContent[0]) # get chrom, start, stop of seq
        reg = reg[0]
    else:
        maxChar = len(lineContent[0])
        for j in range(0, maxChar):
            char = [char.upper() for char in lineContent] # for char in line for each line:
            if char[0] == char[1]:
                if all(modChar == char[0] for modChar in char[2:]):
                    outputHandler.writelines(reg[0] + "\t" + str(int(reg[1]) + j) + "\t" + char[0] + char[0] + char[1]) 
                else:
                    outputHandler.writelines(reg[0] + "\t" + str(int(reg[1]) + j) + "\t" + "°" + char[0] + char[1])
                
            elif all(modChar == char[0] for modChar in char[2:]): # si ts les derniers char sont égaux à ref 1
                outputHandler.writelines(reg[0] + "\t" + str(int(reg[1]) + j) + "\t" + char[0] + char[0] + char[1])
                    
            elif all(modChar == char[1] for modChar in char[2:]):
                outputHandler.writelines(reg[0] + "\t" + str(int(reg[1]) + j) + "\t" + char[1] + char[0] + char[1])
            else:
                outputHandler.writelines(reg[0] + "\t" + str(int(reg[1]) + j) + "\t" + "*" + char[0] + char[1])
