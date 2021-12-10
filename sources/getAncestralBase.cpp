#include "fasta_tools.h"
#include "vcf_tools.h"
#include "bio_tools.h"
#include <iostream>
#include <fstream>
#include <regex>
#include <stdexcept>

int main(int argc, char* argv[])
{
    if(argc < 4) {
        throw std::domain_error("Need at last three args : fasta1, fasta2, vcf output");
    }
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
            entries[i] = files[i].readFastaLine();
        }
        // now you need to compare each base of sequences
        // for this you need to know if outgroups is defined
        // then you need to check if first one is ok or not
        std::string consensus = entries[1].getSequence();
        std::string tstring;
        for(unsigned int i(1); i < argc - 2; i++) {
            tstring = entries[i].getSequence();
            for (unsigned long j(0); j < tstring.size(); j++) {
                if(toupper(tstring[j]) != toupper(consensus[j])) {
                    consensus[j] = 'N';
                }
            }
        }
        // now we have a string with "N" if not common
        std::string ref = entries[0].getSequence();
        for(unsigned long i(0); i<ref.size(); i++) {
            if(toupper(consensus[i]) != toupper(ref[i]) && toupper(ref[i]) != 'N') { // no need to consider N on the ref side bc there already muted
                // write mutation in file
                std::string refS = "";
                char commonBase = consensus[i];
                if(entries[0].getStrand() == '-') {
                    refS += reverseComp(ref[i]);
                    commonBase = reverseComp(commonBase);
                } else {
                    refS += ref[i];
                }

                outputFile.vcf_writeline(vcf_entry(entries[0].getChrom(), entries[0].getPos(i), ".", refS, commonBase));
            }
        }
        count ++;
        if(count % 1000 == 0) {
            std::cout << count << "                    \r";
        }
    }
    return 0;
}
