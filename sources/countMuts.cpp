#include <iostream>
#include <map>
#include "bed_tools.h"
#include "vcf_tools.h"
#include "bio_tools.h"

int main(int argc, char* argv[]) {
    bool lowMem = false;
    bool restart = false;

    if(argc < 5) {
        std::cout << "Not enough args were given : needs AOE, bed, vcf, output file. Optionnal : flag TRUE if you need low-mem, TRUE again if you want to restart from dump" << std::endl;
        exit(1);
    }  else if(argc == 5) {
        std::cout << "Normal start" << std::endl;
    }   else if(argc == 6) {
        if(std::string(argv[5]) == "TRUE") {
            std::cout << "Starting in low mem" << std::endl;
            lowMem = true;
        }
    }  else if(argc == 7) {
        if(std::string(argv[6]) == "TRUE") {
            std::cout << "Restarting from dump" << std::endl;
            restart = true;
        }
    } else {
        std::cout << "Too much args" << std::endl;
        exit(1);
    }

    if(!restart) {
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
        intsOfInterest.dumpAOE(intsOfInterest.size());
    }

    AOEbed inputFile("dump.AOE", read_line);

    std::cout << "Loading mutations..." << std::endl;
    vcf muts(argv[3], read);

    std::map <int, std::map <std::string, std::map <char, int>>> counts;
    int count_lines(0);

    std::cout << "overlapping..." << std::endl;
    while(!inputFile.isEOF()) {
        inputFile.loadBlock(1000000);
        for(const auto &pair: inputFile.getOverlap(muts)) {
            // pair.first is converted vcf & pair.second is a vector of AOE entry
            if(pair.second.size() > 1) {
                std::cout << "Warning : more than one overlap" << std::endl;
            }
            vcf_entry entry(pair.first);
            for(const auto &a_entry: pair.second) {
                std::string str_mut {static_cast<char>(toupper(entry.getAlternate()[0][0])), static_cast<char>(toupper(entry.getRef()[0]))};
                counts[a_entry.getRelativePos(entry.getPos()-1)][str_mut][a_entry.getType()] ++;
                count_lines ++;
            }
        }
        std::cout << count_lines << "\r";
    }

    std::cout << "Writing results ... " << std::endl;
    std::ofstream outputFile(argv[4]);
    for(const auto &pair: counts) {
        for(const auto &pair2: pair.second) {
            for(const auto &pair3: pair2.second)
            outputFile << pair.first << '\t' << pair2.first << '\t' << pair3.first << '\t' << pair3.second << '\n';
        }
    }
    return 0;
}
