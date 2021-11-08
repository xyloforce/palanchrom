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
        
        std::string sequence;
        //fasta_entry entry2 = entry;
        // std::string N = "";
        for(int i(0); i< entry.getSize(); i++) {
            sequence += 'N';
        }
        
        // entry2.editSeq(N, 0, entry2.getSize());
        
        std::map <int, bed_entry> currentInt = intervals.getBedByID(entry.getChrom());
        std::cout << "Checking int..." << currentInt.size() << " intervals left" << std::endl;
        int count = 0;

        for(const auto &pair : currentInt) {
            count ++;
            int size(pair.second.getStop()-pair.second.getStart());
            sequence.replace(pair.second.getStart(),  size, entry.subsetEntry(pair.second.getStart(), pair.second.getStop()).getSequence());

            if(count % 1000 == 0) {
                std::cout << count << "         \r";
            }
        }

        count = 0;
        std::vector <vcf_entry> currentVCF = mutations.getVCFByID(entry.getChrom());
        std::cout << "Checking mutations..." << currentVCF.size() << " mutations left" << std::endl;
        for(const auto &VCFentry : currentVCF) {
            temp = "";
            count ++;
            // check old base
            int posMut = VCFentry.getPos();
            posMut --; // bc vcf are 1-based
            std::string strBase = "";
            strBase += sequence[posMut];

            if(strBase == VCFentry.get_ref()) {
                sequence[posMut] = VCFentry.get_alternate();
            } else {
                std::cout << "Ref is : " << VCFentry.get_ref() << " and current is : " << sequence[posMut] << " at pos : " << posMut << std::endl;
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
