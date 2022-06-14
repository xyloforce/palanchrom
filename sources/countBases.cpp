#include <iostream>
#include "fasta_tools.h"
#include "bed_tools.h"
#include "bio_tools.h"

int main(int argc, char* argv[]) {
    bool lowMem = false;
    bool restart = false;
    int block_quantity(100000);

    if(argc < 5) {
        std::cout << "Doesnt have enough args, need fasta AOE bed and output names. Optionnal : flag TRUE if you need low-mem, again TRUE if you want to restart from dump, int to set number of blocks to load" << std::endl;
        exit(1);
    } else if(argc == 5) {
      std::cout << "normal start" << std::endl;
    } else if(argc == 6) {
        std::cout << "starting in lowmem" << std::endl;
        if(std::string(argv[5]) == "TRUE") {
            lowMem = true;
        }
    }  else if(argc == 7) {
        std::cout << "restarting from dump" << std::endl;
        if(std::string(argv[6]) == "TRUE") {
            restart = true;
        }
    } else if(argc == 8) {
        std::cout << "starting with differing block value" << std::endl;
        if(std::string(argv[5]) == "TRUE") {
            lowMem = true;
        }
        if(std::string(argv[6]) == "TRUE") {
            restart = true;
        }
        try {
            block_quantity = std::stoi(argv[7]);
        } catch(std::invalid_argument) {
            block_quantity = 100000;
        }
    } else {
        std::cout << "too many args, exiting." << std::endl;
        exit(1);
    }

    if(!restart) {
        std::cout << "Loading AOEs..." << std::endl;
        AOEbed intsOfInterest(argv[2]);

        if(lowMem) {
            bed mask(argv[3], openType::read_line);

            std::cout << "Intersecting... " << std::endl;
            intsOfInterest.cutToMask(mask);
        } else {
            std::cout << "Loading bed... " << std::endl;
            sorted_bed mask(argv[3]);

            std::cout << "Intersecting... " << std::endl;
            intsOfInterest.cutToMask(mask);
        }

    //     intsOfInterest.writeToFile(".savestate.tmp");
        intsOfInterest.dumpAOE(intsOfInterest.size());
    }

    std::map <int, std::map<char, std::map <char, int>>> counts; // pos on NIEB : base : type of int : count
    fasta source(argv[1], read, standard);
    AOEbed inputFile("dump.AOE", read_line);

    std::cout << "Loading input block by block" << std::endl;

    while(!inputFile.isEOF()) {
        inputFile.loadBlock(block_quantity);
        std::vector <fasta_entry> toCount = source.getSeqFromInts(inputFile);
        for(int i(0); i < toCount.size(); i ++) {
//             std::cout << toCount[i].getHeader() << std::endl;
            std::string sequence = toCount[i].getUppercaseSequence();
            for(int j(0); j < sequence.size(); j++) {
                counts[inputFile.getEntryByIndex(i).getRelativePos(toCount[i].getPos(j))][sequence[j]][inputFile.getEntryByIndex(i).getType()] ++;
            }
            if(i / toCount.size() % 10 == 0.0) {
                std::cout << i << "\r";
            }
        }
    }

    std::cout << "Writing results ... " << std::endl;
    std::ofstream resultFile(argv[4]);
    for(const auto &intToMap: counts) {
        for(const auto &charToMap: intToMap.second) {
            for(const auto &charToInt: charToMap.second) {
                resultFile << intToMap.first << '\t' << charToMap.first << '\t'  << charToInt.first << '\t'  << charToInt.second << '\n';
            }
        }
    }
    return 0;
}
