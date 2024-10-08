import gzip
import sys

def load_regions(input_file) :
    """
    This function is made to load the liftOver chain file and translate them into a set of corresponding intervals in species 1 and 2
    """
    f = gzip.open(input_file, 'rt')
    # Data will be stored in a dictionnary ordered by species 1 chromosomes 
    regions = {}
    identifier = 0 
    for line in f :
        # print(line)
        # For some reason, there is an empty line at the end of each chain, skip it
        # SKIP ALSO COMMENT LINES
        if line.strip() != "" and not line.startswith("#") :
            # Except for chain lines, chain files are in tab separated format 
            line = line.strip().split("\t")
            # New chain 
            if "chain" in line[0] :
                # Get data from chain line
                chain = line[0].split()
                score = int(chain[1])
                chrom_1 = chain[2]
                size_1 = int(chain[3])
                strand_1 = chain[4]
                start_1 = int(chain[5])
                end_1 = int(chain[6])
                chrom_2 = chain[7]
                size_2 = int(chain[8])
                strand_2 = chain[9]
                start_2 = int(chain[10])
                end_2 = int(chain[11])
                chain_id = int(chain[12])
                # Check that chrom 1 and 2 are well sequenced autosome (we didnt calculate NIEBs on the others). If so, put the flag to True to keep working with these data. If not, put the flag to False to discard these data.
                # DELETED FLAG BC NOT PORTABLE TO OTHER SPECIES

                # For now, - strand in species 1 is not managed. We'll see later if adding it is necessary (hope not)
                if strand_1 == "-" :
                    print("ERROR : Minus strand for species 1 in following chain : ")
                    print(chain)
                    print("Exiting")
                    sys.exit()
                # Add the chrom from species 1 to dictionnary if it is not already in it, and init a list to put the region objects
                if chrom_1 not in regions.keys() :
                    regions[chrom_1] = []
                # Establish start species 1 and species 2 from line data
                # These variables will next be updated with each region
                # If both strand are +, we are in a classical case, both chromosome must be browsed from start to end 
                if strand_1 == "+" and strand_2 == "+" : 
                    s1 = start_1
                    s2 = start_2
                    # s1 and s2 variables will then be updated by adding the size of each region to the variable (and of gaps if necessary). In the end, we should obtain a value corresponding to the end variable for each species 
                    
                # If strand for species 2 is -, the region of species 1 is aligned on the reverse complemented chromosome of species 2
                # We must then browse the chromosome like in following schema :
                # + : 1 -------> 10  Classical case for species 1
                # - : 10 <------- 1  Reverse complemented case for species 2
                # Start (corresponding to 1 in previous schema) can be calculated like this : size_chromosome - start_position.
                elif strand_1 == "+" and strand_2 == "-" :
                    s1 = start_1
                    s2 = size_2 - start_2
                    # In this case, s1 variable will be updated as in classical case while s2 variable will be updated by removing the size of each region, not adding it. In the end, we should obtain a value corresponding to size_chromosome - end_position
                else :
                    print(chain)
            # New interval 	
            else :
                c1 = chrom_1
                identifier += 1 
                # Last entry seems to have no gap column
                if len(line) == 3 : 
                    # Get data from line
                    size = int(line[0])
                    gap_1 = int(line[1])
                    gap_2 = int(line[2])
                else :
                    size = int(line[0])
                    gap_1 = 0
                    gap_2 = 0
                # Calculate end position of the interval in genome 1 and 2 
                # genome 1 
                e1 = s1 + size
                # genome 2
                c2 = chrom_2
                # If regions is + in genome 2, calculate end position as in genome 1
                if strand_2 == "+" : 
                    e2 = s2 + size
                else :
                    e2 = s2
                    s2 = s2-size
                # Init a new inter object and put it in the species 1 chromosome list 
                inter = (c1, s1, e1, chain_id, c2, s2, e2, strand_2, identifier)
                regions[chrom_1].append(inter)
                # Update s1 and s2 variable by adding (or removing if - strand) size of the interval + gap to next interval for the corresponding species 
                s1 = e1 + gap_1
                if strand_2 == "+" : 
                    s2 = e2 + gap_2
                else :
                    s2 = s2 - gap_2

    f.close()
    # Sort each list of object by start position, to put interval in ascending order on chromosome 
    for c in sorted(regions.keys()) :
        # #print(c)
        regions[c].sort(key=lambda x:x[2])

    # At the end, we have a species 1 chromosome-as-keys dictionnary, with, for each chromosome, a list of interval objects sorted by coordinates in ascending order. Each interval object represents a corresponding interval in species 1 and 2
    return regions

def writeIntervals(outputNameH, outputNameS, regions):
    outputH = open(outputNameH, "w")
    outputS = open(outputNameS, "w")
    for key in regions:
        for intervalList in regions[key]:
            intervalList = [str(value) for value in intervalList]
            outputH.write("\t".join(intervalList[:3]) + "\t" + intervalList[3] + "-" + intervalList[-1] + "\t.\t+\n")
            outputS.write("\t".join(intervalList[4:7]) + "\t" + intervalList[3] + "-" + intervalList[-1] + "\t.\t" + intervalList[7] + "\n")

print("Started conversion to bed...")
inputFile = sys.argv[1]
outputH = sys.argv[2]
outputS = sys.argv[3]

writeIntervals(outputH, outputS, load_regions(inputFile))
