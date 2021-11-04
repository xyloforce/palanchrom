#include "fasta_tools.h"
#include "vcf_tools.h"
#include "bio_tools.h"
#include <iostream>
#include <fstream>
#include <regex>
#include <stdexcept>

int main(int argc, char* argv[])
{
    // defines array of fasta objects
    fasta files[argc -2];
    // defines array of fasta entries objects
    fasta_entry entries[argc - 2];
    vcf outputFile(argv[argc - 1], false);
    int count = 0;
    
    std::cout << "Opening files..." << std::endl;
    // open files
    for (int arg(1); arg < argc -1; arg++) {
        files[arg - 1] = fasta(argv[arg], "read_line", true);
    }
    
    std::cout << "Starting analysis..." << std::endl;
    
    while(!files[0].isEOF()) {
        // read one fasta entry for each file
        for(int i(0); i < argc -2; i ++) {
            entries[i] = files[i].read_fasta_line();
        }
        // now you need to compare each base of sequences
        // for this you need to know if outgroups is defined
        // then you need to check if first one is ok or not
        std::string consensus = entries[1].getSequence();
        std::string tstring;
        for(int i(1); i < argc - 2; i++) {
            tstring = entries[i].getSequence();
            for (int j(0); j < tstring.size(); j++) {
                if(tstring[j] != consensus[j]) {
                    consensus[j] = 'N';
                }
            }
        }
        // now we have a string with "N" if not common
        std::string ref = entries[0].getSequence();
        for(int i(0); i<ref.size(); i++) {
            if(consensus[i] != ref[i] && consensus[i] != 'N' && ref[i] != 'N') {
                // write mutation in file
                std::string refS = "";
                char commonBase = consensus[i];
                if(entries[0].getStrand() == '-') {
                    refS += reverseComp(ref[i]);
                    commonBase = reverseComp(commonBase);
                }
                std::cout << consensus[i] << "   " << ref[i] << std::endl;
                std::cout << commonBase << "   " << refS << std::endl;
                outputFile.vcf_writeline(vcf_entry(entries[0].getChrom(), entries[0].getPos(i), ".", refS, commonBase));
            }
        }
        count ++;
        if(count % 10 == 0) {
            std::cout << count << "                    \r";
        }
    }
    return 0;
}
