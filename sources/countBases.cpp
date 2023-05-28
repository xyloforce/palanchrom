#include <iostream>
#include "fasta_tools.h"
#include "bed_tools.h"
#include "bio_tools.h"

int main(int argc, char* argv[]) {
    bool lowMem = false;
    bool restart = false;
    bool no_int = false;
    int block_quantity(50000);
    bool strand = false;

    std::map <char, std::string> args = getArgs(std::vector<std::string>(argv, argv + argc));

    std::string AOEfilename, bedFilename, fastaFilename, outputFilename;

    try {
        AOEfilename = args.at('a');
        fastaFilename = args.at('f');
        outputFilename = args.at('o');
    } catch(std::out_of_range) {
        std::cout << "Missing obligatory parameters. Parameters are : \n";
        std::cout << "\t-a AOE file \n";
        std::cout << "\t-b bed file \n";
        std::cout << "\t-f fasta file \n";
        std::cout << "\t-o output file \n";
        std::cout << "Optionnal stuff includes : \n";
        std::cout << "\t-l low mem mode \n";
        std::cout << "\t-d restart from dump \n";
        std::cout << "\t-s use strand information in intervals \n";
        std::cout << "\t-n ignore bed \n";
        exit(1);
    }

    bedFilename = "";
    if(args.find('b') != args.end()) {
        bedFilename = args.at('b');
    } else if(args.find('n') != args.end()) {
        no_int = true;
    } else {
        std::cout << "Need to set either bed or n" << std::endl;
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
        std::cout << "Strand set to true" << std::endl;
    }

    if(!restart) {
        std::cout << "Loading AOEs..." << std::endl;
        AOEbed intsOfInterest(AOEfilename);

        if(lowMem) {
            bed mask(bedFilename, openType::read_line);

            std::cout << "Intersecting... " << std::endl;
            intsOfInterest.cutToMask(mask, false, false, strand);
        } else if(no_int) {
            std::cout << "Skipping intersection as asked by n flag..." << std::endl;
        } else {
            std::cout << "Loading bed... " << std::endl;
            sorted_bed mask(bedFilename);

            std::cout << "Intersecting... " << std::endl;
            intsOfInterest.cutToMask(mask, false, false, strand);
        }
    //     intsOfInterest.writeToFile(".savestate.tmp");
        intsOfInterest.dumpAOE();
    }

    std::map <int, std::map<char, std::map <char, int>>> counts; // pos on NIEB : base : type of int : count
    fasta source(fastaFilename, read, standard);
    AOEbed inputFile("dump.AOE", read_line);

    std::cout << "Loading input block by block" << std::endl;

    while(!inputFile.isEOF()) {
        inputFile.loadBlock(block_quantity);
        std::vector <fasta_entry> toCount = source.getSeqFromInts(inputFile);
        std::vector <AOE_entry> sorted_entries = inputFile.getSortedEntries();
        for(int i(0); i < toCount.size(); i ++) {
//             std::cout << toCount[i].getHeader() << std::endl;
            std::string sequence = toCount[i].getUppercaseSequence();
            for(int j(0); j < sequence.size(); j++) {
//                 std::cout << inputFile.getSortedEntries()[i].getStringEntry() << "  " << inputFile.getSortedEntries()[i].getZero() << "  " << j << "  "<< sequence[j] << "  " << inputFile.getSortedEntries()[i].getRelativePos(toCount[i].getPos(j)) << std::endl;
//                 std::cout << toCount[i].getHeader() << std::endl;
                counts[sorted_entries[i].getRelativePos(toCount[i].getPos(j, 0))][sequence[j]][sorted_entries[i].getType()] ++;
            }
            if(i / toCount.size() % 10 == 0.0) {
                std::cout << i << "\r";
            }
        }
    }

    std::cout << "Writing results ... " << std::endl;
    std::ofstream resultFile(outputFilename);
    for(const auto &intToMap: counts) {
        for(const auto &charToMap: intToMap.second) {
            for(const auto &charToInt: charToMap.second) {
                resultFile << intToMap.first << '\t' << charToMap.first << '\t'  << charToInt.first << '\t'  << charToInt.second << '\n';
            }
        }
    }
    return 0;
}
