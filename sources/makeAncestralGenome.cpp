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
    vcf mutations(argv[1], true);
    std::cout << "Loading bed... " << std::endl; 
    bed intervals(argv[2], true);
    std::cout << "Started analysis" << std::endl;
    
    fasta inputFasta(argv[3], "read_line", false);
    fasta outputFasta(argv[4], "write", false);
    fasta_entry entry;
    
    std::string temp = "";
    
    while(!inputFasta.isEOF()) {
        entry = inputFasta.read_fasta_line();
        // for each entry in bed
        // subset entry accordingly
        // get seq
        // fill seq
        
        fasta_entry entry2 = entry;
        std::string N = "";
        for(int i(0); i< entry2.getSize(); i++) {
            N += 'N';
        }
        
        entry2.editSeq(N, 0, entry2.getSize());
        std::cout << "Checking int..." << std::endl;
        
        std::map <int, bed_entry> currentInt = intervals.getBedByID(entry2.getChrom());
        for(const auto &pair : currentInt) {
            entry2.editSeq(
                entry.subsetEntry(pair.second.getStart(), pair.second.getStop()).getSequence(),
                           pair.second.getStart(),
                           pair.second.getStop());
        }
        std::vector <vcf_entry> currentVCF = mutations.getVCFByID(entry2.getChrom());
        std::cout << "Checking mutations..." << currentVCF.size() << " mutations left" << std::endl;
        for(const auto &VCFentry : currentVCF) {
            temp = "";
            // check old base
            int posMut = VCFentry.getPos();
            posMut --; // bc vcf are 1-based
            std::cout<< entry2.getSequence() << std::endl;
            if(entry2.subsetEntry(posMut, posMut+1).getSequence() == VCFentry.get_ref()) {
                temp += VCFentry.get_alternate();
                entry2.editSeq(temp, posMut, posMut +1);
            } else {
                std::cout << "Ref is : " << VCFentry.get_ref() << " and current is : " << entry2.subsetEntry(posMut, posMut+1).getSequence() << " at pos : " << posMut << std::endl;
                throw std::logic_error("Probable index issue");
            }
            
        }
        outputFasta.write_fasta_entry(entry2);
    }
    
    return 0;
}
