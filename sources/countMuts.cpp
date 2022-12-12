#include <iostream>
#include <map>
#include "bed_tools.h"
#include "vcf_tools.h"
#include "bio_tools.h"

int main(int argc, char* argv[]) {
    bool lowMem = false;
    bool restart = false;
    bool strand = false;

    std::map <char, std::string> args = getArgs(std::vector<std::string>(argv, argv + argc));

    std::string AOEfilename, bedFilename, vcfFilename, outputFilename;

    try {
        AOEfilename = args.at('a');
        bedFilename = args.at('b');
        vcfFilename = args.at('v');
        outputFilename = args.at('o');
    } catch(std::out_of_range) {
        std::cout << "Missing obligatory parameters. Parameters are : \n";
        std::cout << "\t-a AOE file \n";
        std::cout << "\t-b bed file \n";
        std::cout << "\t-v vcf file \n";
        std::cout << "\t-o output file \n";
        std::cout << "Optionnal stuff includes : \n";
        std::cout << "\t-l low mem mode \n";
        std::cout << "\t-d restart from dump \n";
        std::cout << "\t-s use strand information in intervals \n";
        exit(1);
    }

    if(args.find('l') != args.end()) {
        lowMem = true;
    }

    if(args.find('d') != args.end()) {
        restart = true;
    }

    if(args.find('s') != args.end()) {
        strand = true;
    }

    if(!restart) {
        std::cout << "Loading AOEs..." << std::endl;
        AOEbed intsOfInterest(AOEfilename);
        std::cout << "Loading bed..." << std::endl;

        std::vector <AOE_entry> intersects;
        if(lowMem) {
            bed mask(bedFilename, openType::read_line);
            intsOfInterest.cutToMask(mask, false, false, strand);
        } else {
            sorted_bed mask(bedFilename);
            intsOfInterest.cutToMask(mask, false, false, strand);
        }
        std::cout << "Intersecting finished, dumping..." << std::endl;
        intsOfInterest.dumpAOE(intsOfInterest.size());
    }

    AOEbed inputFile("dump.AOE", read_line);

    std::cout << "Loading mutations..." << std::endl;
    vcf muts(vcfFilename, read);

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
    std::ofstream outputFile(outputFilename);
    for(const auto &pair: counts) {
        for(const auto &pair2: pair.second) {
            for(const auto &pair3: pair2.second)
            outputFile << pair.first << '\t' << pair2.first << '\t' << pair3.first << '\t' << pair3.second << '\n';
        }
    }
    return 0;
}
