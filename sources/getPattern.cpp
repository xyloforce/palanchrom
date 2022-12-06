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
    std::string regex;
    std::string Nregex;
    bool noN = false;
    bool capturingGroups = false;

    if(argc == 5) {
        std::string pattern(argv[1]);
        regex = constructRegex(pattern);
        Nregex = constructRegex(pattern, true);
        std::cout << "regex : " << regex << std::endl;
        std::cout << "n regex : " << Nregex << std::endl;
    } else if(argc >= 6) {
        regex = argv[1];
        std::cout << "matching : " << regex << std::endl;
        if (std::string(argv[5]) == "FALSE") {
            noN = true;
        } else {
            Nregex = argv[5];
            std::cout << "avoiding : " << Nregex << std::endl;
        }
    }
    if(argc == 7 && std::string(argv[6]) == "TRUE") {
        std::cout << "Using capturing groups" << std::endl;
        capturingGroups = true;
    }
    if(argc < 5) {
        std::cout << "Incorrect number of args : need pattern, source fasta and name of bed outputs (2 names). Optionnal : nregex (but first arg is also regex then) or FALSE + regex in first pos. Add TRUE at the end if regex contains capture groups" << std::endl;
        exit(1);
    }

    std::cout << "Reading input..." << std::endl;
    fasta inputFile(argv[2], read_line, standard);

    std::cout << "Creating outputs..." << std::endl;
    bed outputFile(argv[3], openType::write);
    bed outputConvert(argv[4], openType::write);

    std::cout << "Starting analysis..." << std::endl;

    while(!inputFile.isEOF()) {
        fasta_entry entry(inputFile.readFastaLine());
        std::cout << entry.getHeader() << "        \r";
        std::vector <bed_entry> matchs = entry.matchPatterns(regex, capturingGroups);
        for(const auto &bed_line: matchs) {
            outputFile.writeBedLine(bed_line);
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
