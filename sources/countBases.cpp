#include <iostream>
#include "fasta_tools.h"
#include "bed_tools.h"
#include "bio_tools.h"

int main(int argc, char* argv[]) {
    if(argc < 5) {
        throw std::domain_error("Doesnt have enough args, need fasta AOE bed and output names");
    }

    std::cout << "Loading fasta... " << std::endl;
    fasta source(argv[1], read, standard);
    std::cout << "Loading AOEs..." << std::endl;
    AOEbed intsOfInterest(argv[2]);
    std::vector <AOE_entry> toCount;

    if(argc == 6) {
        std::cout << "Loading bed... " << std::endl;
        bed mask(argv[3], openType::read_line);

        std::cout << "Intersecting... " << std::endl;
        toCount = intsOfInterest.getIntersects(mask);
    } else {
        std::cout << "Loading bed... " << std::endl;
        sorted_bed mask(argv[3]);

        std::cout << "Intersecting... " << std::endl;
        toCount = intsOfInterest.getIntersects(mask);
    }

    sort(toCount.begin(), toCount.end());

    std::cout << "Getting seqs... " << std::endl;
    std::vector <fasta_entry> seqs = source.getSeqFromInts(intsOfInterest.convertToBed(toCount));
    std::cout << seqs.size() << std::endl;

    std::map <int, std::map<char, std::map <char, int>>> counts;
    std::cout << "Counting seqs..." << std::endl;
    for(int j(0); j < seqs.size(); j++) {
        std::string sequence = seqs[j].getUppercaseSequence();
        for(int i(0); i < sequence.size(); i++) {
            char base = sequence[i];
            counts[toCount[j].getRelativePos(seqs[j].getPos(i))][base][toCount[j].getType()] ++;
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
