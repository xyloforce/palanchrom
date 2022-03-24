#include <iostream>
#include "fasta_tools.h"
#include "bed_tools.h"
#include "bio_tools.h"

int main(int argc, char* argv[]) {
    bool lowMem = false;
    if(argc < 5) {
        throw std::domain_error("Doesnt have enough args, need fasta AOE bed and output names. Optionnal : flag TRUE if you need low-mem");
    } else if(argc == 6) {
        if(std::string(argv[5]) == "TRUE") {
            lowMem = true;
        }
    }

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

    std::map <int, std::map<char, std::map <char, int>>> counts; // pos on NIEB : base : type of int : count

    fasta source(argv[1], read, standard);
    AOEbed inputFile("dump.AOE", read_line);
    std::cout << "Loading input block by block" << std::endl;

    while(!inputFile.isEOF()) {
        inputFile.loadBlock(100000);
        std::vector <fasta_entry> toCount = source.getSeqFromInts(inputFile);
        for(int i(0); i < toCount.size(); i ++) {
            std::string sequence = toCount[i].getUppercaseSequence();
            for(int j(0); j < sequence.size(); j++) {
                counts[inputFile.getEntryByIndex(i).getRelativePos(toCount[i].getPos(j))][sequence[j]][inputFile.getEntryByIndex(i).getType()] ++;
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
