#include "fasta_tools.h"
#include "vcf_tools.h"
#include <iostream>
#include <fstream>
#include <regex>
#include <stdexcept>

// class header {
// public:
//     int getPos(int count) {
//         if(m_strand == '+') {
//             return (m_start + count);
//         } else if(m_strand == '-') {
//             return (m_stop - count);
//         } else {
//             std::cout << "Strand is : " << m_strand << std::endl;
//             throw std::logic_error("Undefined strand");
//         }
//     }
//     std::string getChrom() {return m_chrom;}
//     header(int start, int stop, char strand, std::string chrom) {
//         m_start = start;
//         m_chrom = chrom;
//         m_strand = strand;
//         m_stop = stop;
//     }
//     header(){}
// private:
//     int m_start;
//     int m_stop;
//     char m_strand;
//     std::string m_chrom;
// };

// void vcf_writer(std::ofstream& output, std::string chrom, int pos, char ref, char alt, bool toInit) {
//     if(toInit) {
//         output<<"##fileformat=VCFv4.2\n";
//         output<<"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
//     }
//     output<<chrom<<'\t'<<pos<<'\t'<<'.'<<'\t'<<ref<<'\t'<<alt<<'\t'<<".\t.\t.\n";
// }

int main(int argc, char* argv[])
{
    // defines array of fasta objects
    fasta files[argc -2];
    // defines array of fasta entries objects
    fasta_entry entries[argc-2];
    vcf outputFile(argv[argc], false);
    int count = 0;
    
    // open files
    for (int arg(1); arg < argc -1; arg++) {
        files[arg] = fasta(argv[arg], "read_line", true);
    }
    
    while(!files[0].isEOF()) {
        // read one fasta entry for each file
        for(int i(0); i < argc -2; i ++) {
            entries[i] = files[i].read_fasta_line(true);
        }
        // now you need to compare each base of sequences
        // for this you need to know if outgroups is defined
        // then you need to check if first one is ok or not
        std::string consensus = entries[2].getSequence();
        std::string tstring;
        for(int i(2); i < argc - 2; i++) {
            tstring = entries[i].getSequence();
            for (int j(0); j < tstring.size(); j++) {
                if(tstring[j] != consensus[j]) {
                    consensus[j] = 'N';
                }
            }
        }
        // now we have a string with "N" if not common
        std::string ref = entries[1].getSequence();
        for(int i(0); i<ref.size(); i++) {
            if(consensus[i] != ref[i] && consensus[i] != 'N' && ref[i] != 'N') {
                // write mutation in file
                outputFile.vcf_writeline(vcf_entry(entries[1].getChrom(), entries[1].getPos(i), ".", ".", ref[i], consensus[i]));
            }
        }
        count ++;
        if(count % 100 == 0) {
            std::cout << count << "                    \r";
        }
    }
    return 0;
}
