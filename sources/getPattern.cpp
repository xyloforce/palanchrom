#include <iostream>
#include <fstream>
#include <regex>

#include "fasta_tools.h"
#include "bed_tools.h"
#include "bio_tools.h"

int main(int argc, char *argv[])
{
    // parse pattern
    // create regex
    // merge pattern
    std::string regex, Nregex;
    std::string fasta_file, output_1, output_2;
    bool noN = false;
    bool capturingGroups = false;
    bool needToBeReversed = false; // pattern does not need to be reverse complemented, is the same in the two directions

    std::map <char, std::string> args = getArgs(std::vector<std::string>(argv, argv + argc));

    try {
        regex = args.at('p');
        fasta_file = args.at('f');
        output_1 = args.at('1');
    } catch(std::out_of_range) {
        std::cout << "Missing obligatory parameters. Parameters are : \n";
        std::cout << "\t-p pattern \n";
        std::cout << "\t-f fasta \n";
        std::cout << "\t-1 output 1 \n";
        std::cout << "Optionnal stuff includes : \n";
        std::cout << "\t-2 output 2 (non matches) \n";
        std::cout << "\t-n pattern that must NOT be matched (set -2 before) \n";
        std::cout << "\t-c has capture groups ? (T if true) \n";
        std::cout << "\t-r count both strands. Can't be used with -2/-n or with regex as a first value" << std::endl;
        exit(1);
    }

    // then : either we want both AND we setted the -n or not

    if(args.find('2') != args.end()) {
        output_2 = args.at('2');
        if(args.find('n') != args.end()) {
            Nregex = args.at('n');
        } else {
            Nregex = constructRegex(regex, true);
        }
    } else { // no output for 2 : no 2
        noN = true;
    }

    if(args.find('c') != args.end()) {
        capturingGroups = true;
    }

    if(args.find('r') != args.end()) {
        if(reverseComp(regex) != regex) {
            needToBeReversed = true;
        } else {
            std::cout << "Pattern is identical when reversed" << std::endl;
        }
    }

    std::cout << "Reading input..." << std::endl;
    fasta inputFile(fasta_file, read_line, standard);

    std::cout << "Creating outputs..." << std::endl;
    bed outputFile(output_1, openType::write);
    bed outputConvert(output_2, openType::write);

    std::cout << "Starting analysis..." << std::endl;

    while(!inputFile.isEOF()) {
        fasta_entry entry(inputFile.readFastaLine());
        std::cout << entry.getHeader() << "        \r";
        std::vector <bed_entry> matchs = entry.matchPatterns(regex, capturingGroups);
        for(auto bed_line: matchs) {
            if(needToBeReversed) {
                bed_line.setStrand('+');
            }
            outputFile.writeBedLine(bed_line);
        }
        if(needToBeReversed) {
            regex = reverseComp(regex);
            std::vector <bed_entry> matchs = entry.matchPatterns(regex, capturingGroups);
            for(auto bed_line: matchs) {
                if(needToBeReversed) {
                    bed_line.setStrand('-');
                }
                outputFile.writeBedLine(bed_line);
            }
        }
        if(!noN) {
            std::vector <bed_entry> tmp = entry.matchPatterns(Nregex, capturingGroups);
            std::vector <bed_entry> convert = entry.reverseInts(tmp);
            for(const auto &bed_line: convert) {
                outputConvert.writeBedLine(bed_line);
            }
        }
    }
    return 0;
}
