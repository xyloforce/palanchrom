#include "vcf_tools.h"
#include "bed_tools.h"
#include "bio_tools.h"
#include "fasta_tools.h"
#include <fstream>
#include <iostream>
#include <map>

int main(int argc, char* argv[]) {
    if(argc < 4) {
        std::cout << "Unsufficient number of args" << std::endl;
        throw std::length_error("need vcf, bed, fasta input and fasta output");
    }
    
    std::cout << "Loading vcf... " << std::endl;
    vcf mutations(argv[1], read);
    std::cout << "Loading bed... " << std::endl; 
    sorted_bed intervals(argv[2]);
    std::cout << "Started analysis" << std::endl;
    
    fasta inputFasta(argv[3], read_line, standard);
    fasta outputFasta(argv[4], write, standard);
    fasta_entry entry;
    
    std::string temp = "";
    
    while(!inputFasta.isEOF()) {
        entry = inputFasta.readFastaLine();
        std::string sequence;
        for(int i(0); i< entry.getSize(); i++) {
            sequence += 'n';
        }
        
        std::vector <bed_entry> currentInt = intervals.getBedByID(entry.getChrom());
        int count = 0;

        for(const auto &line : currentInt) {
            count ++;
            int size(line.getStop()-line.getStart());
            // first arg is POS of char so it's zero based
            sequence.replace(line.getStart(), size, entry.subsetEntry(line.getStart(), line.getStop()).getSequence());

            if(count % 1000 == 0) {
                std::cout << count << "         \r";
            }
        }
        count = 0;
        std::vector <vcf_entry> currentVCF = mutations.getVCFByChrom(entry.getChrom());
        for(const auto &VCFentry : currentVCF) {
            temp = "";
            count ++;
            // check old base
            int posMut = VCFentry.getPos();
            posMut --; // bc vcf are 1-based
            std::string strBase = "";
            strBase += sequence[posMut];

            if(toUpper(strBase) == toUpper(VCFentry.getRef())) {
                sequence[posMut] = VCFentry.getAlternate()[0][0]; // get alternate is vector of string but in this case we can safely assume that it's only one mut by pos;
            } else {
                std::cout << "Ref is : " << VCFentry.getRef() << " and current is : " << sequence.substr(posMut -2, 4) << " at pos : " << posMut << " or as in vcf " << VCFentry.getPos() << std::endl;
                throw std::logic_error("Probable index issue");
            }
            if(count % 100 == 0) {
                std::cout << count << "         \r";
            }
            
        }
        entry.editSeq(sequence, 0, entry.getSize());
        outputFasta.write_fasta_entry(entry);
    }
    
    return 0;
}
