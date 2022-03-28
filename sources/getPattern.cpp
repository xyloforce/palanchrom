#include <iostream>
#include <fstream>
#include <regex>

#include "fasta_tools.h"
#include "bed_tools.h"
#include "bio_tools.h"

int main(int argc, char *argv[])
{
    if (argc < 5)
    {
        throw std::domain_error("Unsufficient number of args : need pattern, source fasta and name of bed outputs (2 names). Optionnal : nregex (but first arg is also regex then)");
    }
    // parse pattern
    // create regex
    // merge pattern
    std::string regex;
    std::string Nregex;

    if(argc == 5) {
        std::string pattern(argv[1]);
        regex = constructRegex(pattern);
        Nregex = constructRegex(pattern, true);
        std::cout << "regex : " << regex << std::endl;
        std::cout << "n regex : " << Nregex << std::endl;
    } else {
        regex = argv[1];
        Nregex = argv[5];
    }


    std::cout << "Reading input..." << std::endl;
    fasta inputFile(argv[2], read_line, standard);

    std::cout << "Creating outputs..." << std::endl;
    bed outputFile(argv[3], openType::write);
    bed outputConvert(argv[4], openType::write);

    std::cout << "Starting analysis..." << std::endl;

    while(!inputFile.isEOF()) {
        fasta_entry entry(inputFile.readFastaLine());
        std::cout << entry.getHeader() << std::endl;
        std::vector <bed_entry> matchs = entry.matchPatterns(regex);
        for(const auto &bed_line: matchs) {
            outputFile.writeBedLine(bed_line);
        }

        std::vector <bed_entry> tmp = entry.matchPatterns(Nregex);
        std::vector <bed_entry> convert = entry.reverseInts(tmp);
        for(const auto &bed_line: convert) {
            outputConvert.writeBedLine(bed_line);
        }
    }
    return 0;
}
