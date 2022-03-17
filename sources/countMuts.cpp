#include <iostream>
#include <map>
#include "bed_tools.h"
#include "vcf_tools.h"
#include "bio_tools.h"

void dump(std::vector <AOE_entry> &data, int limit) {
    std::ofstream dumpH("dump.AOE");
    int count (0);
    for(const auto &entry: data) {
        dumpH << entry.to_string() << "\n";
        count ++;
        if(count == limit) {
            break;
        }
    }
    data.erase(data.begin(), data.begin()+limit);
}

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
    std::cout << "Intersecting finished, dumping..." << std::endl;
    dump(intersects, intersects.size()/2);
    return intersects;
}

std::map <bed_entry, std::vector <AOE_entry>> overlap (AOEbed &virtualF, std::string vcfFilename, int argc) {

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
    AOEbed virtualF(intersect (std::string(argv[1]), std::string(argv[2]), argc));

    for(int i(0); i < 2; i++) {
        std::map <bed_entry, std::vector <AOE_entry>> overlaps = overlap(virtualF, std::string(argv[3]), argc);
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
            virtualF = AOEbed("dump.AOE");
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
