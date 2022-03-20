#include <iostream>
#include <map>
#include "bed_tools.h"
#include "vcf_tools.h"
#include "bio_tools.h"

int main(int argc, char* argv[]) {
    bool lowMem = false;

    if(argc < 5) {
        throw std::logic_error("Not enough args were given : needs AOE, bed, vcf, output file. Optionnal : flag TRUE if you need low-mem");
    }  else if(argc == 6) {
        if(std::string(argv[5]) == "TRUE") {
            lowMem = true;
        }
    }

    std::map <int, std::map <std::string, std::map <char, int>>> counts;

    std::cout << "Loading AOEs..." << std::endl;
    AOEbed intsOfInterest(argv[1]);
    std::cout << "Loading bed..." << std::endl;

    std::vector <AOE_entry> intersects;
    if(lowMem) {
        bed mask(argv[2], openType::read_line);
        intsOfInterest.cutToMask(mask);
    } else {
        sorted_bed mask(argv[2]);
        intsOfInterest.cutToMask(mask);
    }
    std::cout << "Intersecting finished, dumping..." << std::endl;
    intsOfInterest.dumpAOE(intsOfInterest.size()/2);

    for(int i(0); i < 2; i++) {
        std::map <bed_entry, std::vector <AOE_entry>> overlaps;
        std::cout << "Loading mutations..." << std::endl;
        vcf muts(argv[3], read);
        std::cout << "Getting muts in ints" << std::endl;
        overlaps = intsOfInterest.getOverlap(muts);
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
        if(i == 0) {
            std::cout << "Loading dump..." << std::endl;
            intsOfInterest = AOEbed("dump.AOE");
        }
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
