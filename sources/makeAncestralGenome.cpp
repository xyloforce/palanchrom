#include "vcf_tools.h"
#include "bed_tools.h"
#include "bio_tools.h"
#include <fstream>
#include <iostream>

int main(int argc, char* argv[]) {
    std::cout << "Loading vcf... " << std::endl;
    vcf mutations(argv[1], true);
    std::cout << "Loading bed... " << std::endl; 
    bed intervals(argv[2], true);
    std::cout << "Started analysis" << std::endl;
    
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
    }
    
    return 0;
}
