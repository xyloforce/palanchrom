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
        std::cout << "Loading bed... " << std::endl;
        bed mask(argv[3], openType::read_line);

        std::cout << "Intersecting... " << std::endl;
        intsOfInterest.cutToMask(mask);
    } else {
        std::cout << "Loading bed... " << std::endl;
        sorted_bed mask(argv[3]);

        std::cout << "Intersecting... " << std::endl;
        intsOfInterest.cutToMask(mask);
    }

    intsOfInterest.dumpAOE(intsOfInterest.size()/2);
    std::map <int, std::map<char, std::map <char, int>>> counts;
    std::cout << "Loading fasta... " << std::endl;
    fasta source(argv[1], read_line, standard);

    for(int i(0); i < 2; i ++) {
        std::cout << "Getting seqs... " << std::endl;
        source.subsetFromInts(intsOfInterest); // change this to intra-fasta code
        // even better : load entry / cut / load next entry
         std::cout << source.size() << std::endl;

        std::cout << "Counting seqs..." << std::endl;
        for(int j(0); j < source.size(); j++) {
            std::string sequence = source.getFastaByIndex(j).getUppercaseSequence();
            for(int i(0); i < sequence.size(); i++) {
                char base = sequence[i]; // FIXME TO COUNT IS USELESS !!!!
                counts[intsOfInterest.getEntryByIndex(j).getRelativePos(source.getFastaByIndex(j).getPos(i))][base][intsOfInterest.getEntryByIndex(j).getType()] ++;
            }
        }
        if(i == 0) {
            std::cout << "Loading dump.." << std::endl;
            intsOfInterest = AOEbed("dump.AOE");
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
