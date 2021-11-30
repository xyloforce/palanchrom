#include <iostream>
#include <map>
#include "bed_tools.h"
#include "vcf_tools.h"

int baseToIndex(char base) {
    switch(toupper(base)) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
        case 'N':
            return 4;
        default:
            std::cout << "Given base : " << base << std::endl;
            throw std::domain_error("Unexpected base");
    }
}

int main(int argc, char* argv[]) {
// in each CG is first and not CG second
    std::array<std::map <int, std::map<char, std::array<int, 5>>>, 2> ACGTbyType;
    std::array<std::vector <vcf_entry>, 2> mutsByType;

    // files
    if(argc < 6) {
        throw std::domain_error("Not enough args were given : need vcf, CG bed, AOE file, outputs name");
    }
    std::cout << "Loading vcf ... " << std::endl;
    vcf muts(argv[1], true);

    std::cout << "Loading bed ... " << std::endl;
    sorted_bed CGints (argv[2]);

    std::cout << "Loading AOEs ... " << std::endl;
    AOEbed AOEs (argv[3]);

    // 1 find pos that match with CG entries
    std::cout << "Starting analysis" << std::endl;
    for(const auto &entry: muts.getVCFEntries()) {
        if(CGints.isInside(bed_entry(entry.getChrom(), entry.getPos(), entry.getPos() + 1))) {
            mutsByType[0].push_back(entry);
        } else if (!(entry == vcf_entry())) {
            mutsByType[1].push_back(entry);
        }
    }

    std::cout << "Found " << mutsByType[0].size() << " mutations in CPG zones" << std::endl;
    std::cout << "Found " << mutsByType[1].size() << " mutations in nCPG zones" << std::endl;

    // 2 for the two lists of pos find matching AOE & translate pos relative to AOE
    for(unsigned int i(0); i < mutsByType.size(); i++) {
        int count(0);
        for(const auto &entry : mutsByType[i]) {
            std::vector <bed_entry> convertedInts;
            bed_entry convert(entry.getChrom(), entry.getPos()-1, entry.getPos());
            convertedInts.push_back(convert);
            std::map <bed_entry, std::vector <AOE_entry>> matching_AOEs = AOEs.getOverlap(entry.getChrom(), convertedInts);
            if(matching_AOEs[convert].size() == 1) {
                AOE_entry result = matching_AOEs[convert][0];
                int posAOE(result.getRelativePos(entry.getPos()));
                ACGTbyType[i][posAOE][entry.getRef()[0]][baseToIndex(entry.getAlternate())] ++;
            } else if (matching_AOEs[convert].size() > 1) {
                std::cout<< matching_AOEs[convert][0].getStringEntry() << std::endl;
                std::cout<< matching_AOEs[convert][1].getStringEntry() << std::endl;
                throw std::domain_error("Overlapping ints in AOE file");
            }
            count ++;
            if(count % 1000 == 0) {
                std::cout << count << "                \r";
            }
        }
    }

    // 3 output in the right file
    std::array<std::ofstream, 2> files;
    files[0] = std::ofstream(argv[4]);
    files[1] = std::ofstream(argv[5]);

    for(unsigned int i(0); i < ACGTbyType.size(); i++) {
        for(const auto &pair: ACGTbyType[i]) { // pair is pos : map(char, int[5])
            for(const auto &pair2 : pair.second) { // pair2 is char: int[5]
                files[i] << pair.first << '\t' << pair2.first << '\t' << 'A' << '\t' << pair2.second[0] << '\n';
                files[i] << pair.first << '\t' << pair2.first << '\t' << 'C' << '\t' << pair2.second[1] << '\n';
                files[i] << pair.first << '\t' << pair2.first << '\t' << 'G' << '\t' << pair2.second[2] << '\n';
                files[i] << pair.first << '\t' << pair2.first << '\t' << 'T' << '\t' << pair2.second[3] << '\n';
                files[i] << pair.first << '\t' << pair2.first << '\t' << 'N' << '\t' << pair2.second[4] << '\n';
            }
        }
    }
}
