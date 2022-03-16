#include <iostream>
#include <map>
#include "bed_tools.h"
#include "vcf_tools.h"
#include "bio_tools.h"

std::vector <AOE_entry> intersect (std::string AOEfilename, std::string bedFilename, int argc) {
    std::cout << "Loading AOEs..." << std::endl;
    AOEbed intsOfInterest(AOEfilename);
    std::cout << "Loading bed..." << std::endl;

    std::vector <AOE_entry> intersects;
    if(argc == 5) {
        sorted_bed mask(bedFilename);
        intersects = intsOfInterest.getIntersects(mask);
    } else {
        bed mask(bedFilename, openType::read_line);
        intersects = intsOfInterest.getIntersects(mask);
    }
    return intersects;
}

std::map <bed_entry, std::vector <AOE_entry>> overlap (std::string vcfFilename, std::string AOEfilename, std::string bedFilename, int argc) {
    AOEbed virtualF(intersect (AOEfilename, bedFilename, argc));
    std::map <bed_entry, std::vector <AOE_entry>> overlaps;

    if(argc == 5) {
        std::cout << "Loading mutations..." << std::endl;
        vcf muts(vcfFilename, read);
        std::cout << "Getting muts in ints" << std::endl;
        overlaps = virtualF.getOverlap(muts);
    } else {
        vcf muts(vcfFilename, read_line);
        std::cout << "Getting muts in ints" << std::endl;
        overlaps = virtualF.getOverlapLowMem(muts);
    }
    return overlaps;
}

int main(int argc, char* argv[]) {
    if(argc < 5) {
        throw std::logic_error("Not enough args were given : needs AOE, bed, vcf, output file");
    }
    std::map <int, std::map <std::string, std::map <char, int>>> counts;
    std::map <bed_entry, std::vector <AOE_entry>> overlaps = overlap(std::string(argv[3]), std::string(argv[1]), std::string(argv[2]), argc);

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
