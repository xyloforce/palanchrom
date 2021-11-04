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
    fasta outputFasta(argv[4], "write", true);
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
        std::map <int, vcf_entry> currentVCF = mutations.getVCFByID(entry2.getChrom());
        std::cout << "Checking mutations..." << std::endl;
        for(const auto &pair : currentVCF) {
            temp = "";
            // check old base
            int posMut = pair.second.getPos();
            posMut --; // bc vcf are 1-based
            std::cout<< entry2.getSequence() << std::endl;
            if(entry2.subsetEntry(posMut, posMut+1).getSequence() == pair.second.get_ref()) {
                temp += pair.second.get_alternate();
                entry2.editSeq(temp, pair.first -1, pair.first);
            } else {
                std::cout << "Ref is : " << pair.second.get_ref() << " and current is : " << entry2.subsetEntry(posMut, posMut+1).getSequence() << " at pos : " << posMut << std::endl;
                throw std::logic_error("Probable index issue");
            }
            
        }
        outputFasta.write_fasta_entry(entry2);
    }
    /*
    std::ifstream fasta(argv[3]);
    std::ofstream ancestralFasta(argv[4]);
    std::string header;
    char tchar;
    int start;
    
    const int size = 1000; // size of int to check for overlap with bed
    
    int count = 0;
    bed_entry temp;
    
    while(!fasta.eof()) {
        fasta.get(tchar);
        if(tchar == '>') {
            while(tchar != '\n') {
                fasta.get(tchar);
                if(tchar != '\n') {
                    header += tchar;
                }
            }
            start = 0;
            ancestralFasta << '>' << header << '\n';
        } else if(tchar != '\n') {
            start ++;
            // check if were in a bed region like on ten nt
            // yes : check each time for mutations
            // no : write N
            // overlapping : write N until you reach the border
            if(count == 0) {
                temp = intervals.inInt(header, start, size);
                if(temp == bed_entry()) {
                    ancestralFasta << 'N';
                    count = size;
                } else if(start < temp.getStart()) {
                    ancestralFasta << 'N';
                    count = temp.getStart() - (start + 1);
                } else {
                    std::string tstring = "";
                    // in all other cases i'm for now at last in the current int
                    if(temp.getStrand() == '+') {
                        tstring += tchar;
                        char base = mutations.isMuted(header, start, tstring)[0];
                        ancestralFasta << base;
                    } else if(temp.getStrand() == '-') {
                        tchar = reverseComp(tchar);
                        tstring += tchar;
                        char base = mutations.isMuted(header, start, tstring)[0];
                        ancestralFasta << reverseComp(base);
                    } else {
                        std::cout << temp.getStrand() << std::endl;
                        throw std::domain_error("Strand unset");
                    }
                }
            } else {
                ancestralFasta << 'N';
                count --;
            }
        } else {
            ancestralFasta << '\n';
        }
    }*/
    
    return 0;
}
