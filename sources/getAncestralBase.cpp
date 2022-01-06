#include "fasta_tools.h"
#include "vcf_tools.h"
#include "bio_tools.h"
#include <iostream>
#include <fstream>
#include <regex>
#include <stdexcept>

int main(int argc, char* argv[])
{
    if(argc < 5) {
        throw std::domain_error("Need at last three args : ref fasta 1, ref fasta 2, outgroup fasta 3, vcf output");
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

        std::string consensus = entries[2].getUppercaseSequence();
        std::string tstring;
        for(unsigned int i(2); i < argc - 2; i++) { // second file has "half-ref" status
            tstring = entries[i].getUppercaseSequence();
            for (unsigned long j(0); j < tstring.size(); j++) {
                if(tstring[j] != consensus[j]) {
                    consensus[j] = 'N';
                }
            }
        }
        std::string ref = entries[0].getUppercaseSequence();
        std::string ref2 = entries[1].getUppercaseSequence();
        for(unsigned long i(0); i<ref.size(); i++) {
            if(consensus[i] != ref[i] && ref[i] != 'N') {
                std::string refS("");
                std::string info(".");
                char commonBase('\0');

                if(consensus[i] == ref2[i]) {
                    commonBase = consensus[i];
                    if(entries[0].getStrand() == '-') {
                        refS += reverseComp(ref[i]);
                        commonBase = reverseComp(commonBase);
                    } else {
                        refS += ref[i];
                    }
                }
                else if(consensus[i] != 'N') {
                    commonBase = 'N';
                    info = "source=non-match;ref=";
                    info += consensus[i];
                    info += ";alt=";
                    info += ref2[i];
                    info += ";";
                    if(entries[0].getStrand() == '-') {
                        refS += reverseComp(ref[i]);
                    } else {
                        refS += ref[i];
                    }
                }
                else {
                    commonBase = 'N';
                    if(entries[0].getStrand() == '-') {
                        refS += reverseComp(ref[i]);
                    } else {
                        refS += ref[i];
                    }
                }
                std::vector <std::string> alternate(1);
                alternate[0] += commonBase;
                outputFile.vcf_writeline(vcf_entry(entries[0].getChrom(), entries[0].getPos(i), ".", refS, alternate, 0, ".", info));
            }
        }
    }
}
