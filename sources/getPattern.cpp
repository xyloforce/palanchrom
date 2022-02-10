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
        throw std::domain_error("Unsufficient number of args : need pattern, source fasta and name of bed outputs (2 names)");
    }
    std::string pattern(argv[1]);
    std::string regex = constructRegex(addNEachPos(pattern));

    std::cout << "Reading input..." << std::endl;
    fasta inputFile(argv[2], read_line, standard);

    std::cout << "Creating outputs..." << std::endl;
    bed outputFile(argv[3], false);
    bed outputConvert(argv[4], false);

    std::cout << "Starting analysis..." << std::endl;

    while(!inputFile.isEOF()) {
        fasta_entry entry(inputFile.readFastaLine());
        std::vector <bed_entry> matchs = entry.matchPatterns(pattern);
        for(const auto &bed_line: matchs) {
            outputFile.writeBedLine(bed_line);
        }

        std::vector <bed_entry> tmp = entry.matchPatterns(regex);
        std::vector <bed_entry> convert = entry.reverseInts(tmp);
        for(const auto &bed_line: convert) {
            outputConvert.writeBedLine(bed_line);
        }
    }
    return 0;
}
