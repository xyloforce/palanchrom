import sys
import glob

folder = sys.argv[1]
intersected = sys.argv[2]
outputFile = sys.argv[3]

fileList = glob.glob(folder + "/*.dat")

outputHandler = open(outputFile, "w")

commonHandler = open(intersected)
commonDict = dict()
old_file = ""

# TODO !!!ATTENTION!!! "-" : intervalle dans un sens sur l'humain & dans l'autre sur le chimpanzé : corriger problèmes
# if - : startH = stopS // stopH = startS

for intervalR in commonHandler.readlines(): # intervalR is line in intersected file
    intervalR = intervalR.strip().split("\t") # intervalR is chr:start:stop
    if old_file != folder + "/" + intervalR[0] + ".dat":
        fileHandler = open(folder + "/" + intervalR[0] + ".dat")
        print("\nOpening file " + intervalR[0] + ".dat")
        intervalList = list()
        for line2 in fileHandler.readlines(): # is a list
            intervalList.append(line2.strip().split("\t")) # intervalList is list of startH:stopH:chrS:startS:stopS:strand
        print("Loading ended, size of list : " + str(len(intervalList)))
        old_file = folder + "/" + intervalR[0] + ".dat"
        intervalList.sort(key=lambda x: int(x[0]))
        windowBegin = 0
        windowEnd = len(intervalList)
    else:
        windowEnd = len(intervalList)

    while windowEnd-windowBegin > 0:
        intervalO = intervalList[int((windowEnd + windowBegin)/2)]
        if int(intervalO[0]) > int(intervalR[1]): # reduced start is smaller than original start
            windowEnd = int((windowEnd + windowBegin)/2)
        else: # reduced start bigger
            if int(intervalO[1]) >= int(intervalR[2]): # is reduced end smaller than original ?
                # you found it !!!
                startSlide = int(intervalR[1]) - int(intervalO[0]) # feature starts later...
                endSlide = int(intervalO[1]) - int(intervalR[2]) # or end earlier
                if intervalO[5] == "+":
                    start = int(intervalO[3]) + startSlide
                    end = int(intervalO[4]) - endSlide
                else:
                    start = int(intervalO[3]) + endSlide
                    end = int(intervalO[4]) - startSlide

                break # no way yo u can find 2 overlaping elements
            else:
                windowBegin = int((windowEnd + windowBegin)/2)

    outputHandler.writelines(intervalO[2] + "\t" + str(start) + "\t" + str(end) + "\t" + intervalO[5] + "\n")
