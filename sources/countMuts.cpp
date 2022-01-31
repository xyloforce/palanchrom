#include <iostream>
#include <map>
#include "bed_tools.h"
#include "vcf_tools.h"
#include "bio_tools.h"

int main(int argc, char* argv[]) {
    if(argc < 4) {
        throw std::logic_error("Not enough args were given : needs AOE, bed, vcf, output file");
    }
    std::map <int, std::map <std::string, std::map <char, int>>> counts;

    std::cout << "Loading AOEs..." << std::endl;
    AOEbed intsOfInterest(argv[1]);
    std::cout << "Loading bed..." << std::endl;
    sorted_bed mask(argv[2]);
    
    std::cout << "Intersecting ints..." << std::endl;
    AOEbed virtualF(intsOfInterest.getIntersects(mask));
    
    std::cout << "Loading mutations..." << std::endl;
    vcf muts(argv[3], read);

    std::cout << "Getting muts in ints" << std::endl;
    std::map <bed_entry, std::vector <AOE_entry>> overlaps = virtualF.getOverlap(muts);

    std::cout << "Getting muts by pos" << std::endl;
    for(const auto &pair: overlaps) {
        // pair.first is converted vcf & pair.second is a vector of AOE entry
        if(pair.second.size() > 1) {
            throw std::logic_error("More than one overlap");
        }
        vcf_entry entry(pair.first);
        std::string str_mut("");
        str_mut += toupper(entry.getAlternate()[0][0]);
        str_mut += toupper(entry.getRef()[0]);
        counts[pair.second[0].getRelativePos(entry.getPos()-1)][str_mut][pair.second[0].getType()] ++;
    }

    std::ofstream outputFile(argv[4]);
    for(const auto &pair: counts) {
        for(const auto &pair2: pair.second) {
            for(const auto &pair3: pair2.second)
            outputFile << pair.first << '\t' << pair2.first << '\t' << pair3.first << '\t' << pair3.second << '\n';
        }
    }
    return 0;
}
