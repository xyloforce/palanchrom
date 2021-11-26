#include <iostream>
#include <map>
#include "bed_tools.h"
#include "vcf_tools.h"

int baseToIndex(char base) {
    switch(base) {
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
    // vectors of muts
    std::map <int, std::map<char, std::array<int, 5>>> ACGTMutByPosCG;
    std::map <int, std::map<char, std::array<int, 5>>> ACGTMutByPosNotCG;
    // vectors of vcf entries
    std::vector <vcf_entry> mutsForCG;
    std::vector <vcf_entry> mutsForNotCG;
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
            mutsForCG.push_back(entry);
        } else if (!(entry == vcf_entry())) {
            mutsForNotCG.push_back(entry);
        }
    }

    std::cout << "Found " << mutsForNotCG.size() << " mutations in nCPG zones" << std::endl;
    std::cout << "Found " << mutsForCG.size() << " mutations in CPG zones" << std::endl;

    // 2 for the two lists of pos find matching AOE & translate pos relative to AOE
    for(const auto &entry : mutsForCG) {
        std::vector <bed_entry> convertedInts;
        bed_entry convert(entry.getChrom(), entry.getPos(), entry.getPos() + 1);
        convertedInts.push_back(convert);
        std::map <bed_entry, std::vector <AOE_entry>> matching_AOEs = AOEs.getOverlap(entry.getChrom(), convertedInts);
        if(matching_AOEs[convert].size() == 1) {
            AOE_entry result = matching_AOEs[convert][0];
            int posAOE(result.getRelativePos(entry.getPos()));
            ACGTMutByPosCG[posAOE][entry.getRef()[0]][baseToIndex(entry.getAlternate())] ++;
        } else if (matching_AOEs[convert].size() > 1) {
            std::cout<< matching_AOEs[convert][0].getStringEntry() << std::endl;
            std::cout<< matching_AOEs[convert][1].getStringEntry() << std::endl;
            throw std::domain_error("Overlapping ints in AOE file");
        }
    }

    for(const auto &entry : mutsForNotCG) {
        std::vector <bed_entry> convertedInts;
        bed_entry convert(entry.getChrom(), entry.getPos(), entry.getPos() + 1);
        convertedInts.push_back(convert);
        std::map <bed_entry, std::vector <AOE_entry>> matching_AOEs = AOEs.getOverlap(entry.getChrom(), convertedInts);
        if(matching_AOEs[convert].size() == 1) {
            AOE_entry result = matching_AOEs[convert][0];
            int posAOE(result.getRelativePos(entry.getPos()));
            ACGTMutByPosNotCG[posAOE][entry.getRef()[0]][baseToIndex(entry.getAlternate())] ++;
        } else if (matching_AOEs[convert].size() > 1) {
            throw std::domain_error("Overlapping ints in AOE file");
        }
    }
    // 3 output in the right file
    std::ofstream CG_output(argv[4]);
    for(const auto &pair: ACGTMutByPosCG) { // pair is pos : map(char, int[5])
        for(const auto &pair2 : pair.second) { // pair2 is char: int[5]
            CG_output << pair.first << '\t' << pair2.first << '\t' << 'A' << '\t' << pair2.second[0] << '\n';
            CG_output << pair.first << '\t' << pair2.first << '\t' << 'C' << '\t' << pair2.second[1] << '\n';
            CG_output << pair.first << '\t' << pair2.first << '\t' << 'G' << '\t' << pair2.second[2] << '\n';
            CG_output << pair.first << '\t' << pair2.first << '\t' << 'T' << '\t' << pair2.second[3] << '\n';
            CG_output << pair.first << '\t' << pair2.first << '\t' << 'N' << '\t' << pair2.second[4] << '\n';
        }
    }
    std::ofstream notCG_output(argv[5]);
    for(const auto &pair: ACGTMutByPosNotCG) { // pair is pos : map(char, int[5])
        for(const auto &pair2 : pair.second) { // pair2 is char: int[5]
            notCG_output << pair.first << '\t' << pair2.first << '\t' << 'A' << '\t' << pair2.second[0] << '\n';
            notCG_output << pair.first << '\t' << pair2.first << '\t' << 'C' << '\t' << pair2.second[1] << '\n';
            notCG_output << pair.first << '\t' << pair2.first << '\t' << 'G' << '\t' << pair2.second[2] << '\n';
            notCG_output << pair.first << '\t' << pair2.first << '\t' << 'T' << '\t' << pair2.second[3] << '\n';
            notCG_output << pair.first << '\t' << pair2.first << '\t' << 'N' << '\t' << pair2.second[4] << '\n';
        }
    }
}
