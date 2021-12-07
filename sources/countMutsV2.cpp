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
    std::array<std::map <std::string, std::vector <int>>, 2> mutsByType;

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
    std::cout << std::endl;
    std::cout << "Starting analysis" << std::endl;
    int count(0);
    std::cout << "Searching muts in CGs ..." << std::endl;
    std::vector <bed_entry> convertedInts;
    for(const auto &entry: muts.getVCFEntries()) {
        convertedInts.push_back(bed_entry(entry));
    }
    std::vector <bool> areInside = CGints.areInside(convertedInts);
    std::cout << "Finished searching, sorting" << std::endl;
    int nCPG(0), CPG(0);
    for(unsigned int i(0); i < areInside.size(); i ++) {
        if(areInside[i]) {
            mutsByType[0][muts.getVCFEntry(i).getChrom()].push_back(i);
            CPG ++;
        } else {
            mutsByType[1][muts.getVCFEntry(i).getChrom()].push_back(i);
            nCPG ++;
        }
    }
    std::cout << "Found " << CPG << " mutations in CPG zones" << std::endl;
    std::cout << "Found " << nCPG << " mutations in nCPG zones" << std::endl;
    std::cout << std::endl;

    // 2 for the two lists of pos find matching AOE & translate pos relative to AOE
    for(unsigned int i(0); i < mutsByType.size(); i++) {
        count = 0;
        for(const auto &chromIndex: mutsByType[i]) { // chromIndex is std::string : std::vector <int>
            std::vector <bed_entry> convertedInts;
            for(const auto &index : chromIndex.second) { // index is on muts.getVCFentries
                bed_entry tmp(muts.getVCFEntry(index));
                tmp.setName(std::to_string(index));
                convertedInts.push_back(tmp); // convert vcf to bed and put it in vector
            }
            std::map <bed_entry, std::vector <AOE_entry>> matching_AOEs = AOEs.getOverlap(chromIndex.first, convertedInts);
            for(const auto &result: matching_AOEs) { // bed_entry : one AOE_entry normally
                if(result.second.size() > 1) {
                    std::cout<< result.second[0].getStringEntry() << std::endl;
                    std::cout<< result.second[1].getStringEntry() << std::endl;
                    throw std::domain_error("Overlapping ints in AOE file");
                }
                vcf_entry entry = muts.getVCFEntry(stoi(result.first.getName()));
                int posAOE(result.second[0].getRelativePos(entry.getPos()-1)); // one based vs zero based
                ACGTbyType[i][posAOE][entry.getAlternate()][baseToIndex(entry.getRef()[0])] ++;
            }
            count ++;
            std::cout << "Treated " << count << "chromosomes                 \r";
        }
    }
    std::cout << std::endl;
    std::cout << "Writing results to file" << std::endl;

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
