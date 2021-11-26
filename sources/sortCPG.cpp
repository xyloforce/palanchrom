#include <iostream>
#include <fstream>
#include <regex>

#include "fasta_tools.h"

// change it : create CPG bed then search vcf pos in it

int main(int argc, char* argv[]) {
    if(argc < 4) {
        throw std::domain_error("Unsufficient number of args");
    }
    
    std::cout << "Reading input..." << std::endl;
    fasta inputFile(argv[1], "read_line", false);
    fasta outputCPG(argv[2], "write", false);
    fasta outputNCPG(argv[3], "write", false);
    
    // read an entry
    // create two strings : one with CG only one without
    std::cout << "Searching for CG..." << std::endl;
    int count(0);
    while(!inputFile.isEOF()) {
        fasta_entry entry(inputFile.read_fasta_line());
        std::string no_cpg("");
        std::string cpg("");

        std::string sequenceData(entry.getSequence());
        char lastChar('\0');
        for(int i(0); i<sequenceData.size(); i++) {
            if(lastChar != '\0') {
                if(toupper(lastChar) == 'C' && toupper(sequenceData[i]) == 'G') {
                    cpg = cpg + lastChar + sequenceData[i];
                    no_cpg += "--";
                    lastChar = '\0';
                } else {
                    cpg += "--";
                    no_cpg = no_cpg + lastChar + sequenceData[i];
                    lastChar = '\0';
                }
            } else if((i+1) > sequenceData.size()) {
                no_cpg += sequenceData[i];
            } else {
                lastChar = sequenceData[i];
            }
        } // sequence done

        // replace seqs
        fasta_entry no_cpg_entry(entry);
        no_cpg_entry.editSeq(no_cpg, 0, no_cpg_entry.getSize());
        fasta_entry cpg_entry(entry);
        cpg_entry.editSeq(cpg, 0, cpg_entry.getSize());

        // write entries
        outputCPG.write_fasta_entry(cpg_entry);
        outputNCPG.write_fasta_entry(no_cpg_entry);
        count ++;
        if(count % 10 == 0) {
            std::cout << "Treated " << count << " lines                      \r";
        }
    }
    return 0;
}
    
//     std::ifstream globalFile(argv[1]);
//     std::ofstream CPGFile(argv[2]);
//     std::ofstream nCPGFile(argv[3]);
    
//     char tchar;
//     char beforeChar = '/';
//     // int countChars = 0;
//     bool writtenCHeader = false;
//     bool writtenNHeader = false;
//     int lastCPG = 0;
//     int lastNCPG = 0;
//     bool firstLineC = true;
//     bool firstLineN = true;
    
//     int start = 0;
//     std::string chrom = "";
    
//     while(!globalFile.eof()) {
//         globalFile.get(tchar);
//         if(tchar != '\n' && !globalFile.eof()) {
//             // countChars ++;
//             if(tchar == '>') {
//                 std::string header = "";
//                 while(tchar != '\n') {
//                     header += tchar;
//                     globalFile.get(tchar);
//                 }
//                 std::smatch m;
//                 std::regex e(">(\\w+):(\\d+)-\\d+");
//                 regex_search(header, m, e);
//                 chrom = m.str(1);
//                 start = stoi(m.str(2)) - 1;
//                 writtenCHeader = false;
//                 writtenNHeader = false;
//             } else if(tchar == '\n') {
//                 // countChars = 0;
//             } else { // not a > not a '\n' neither a letter in a header
//                 start ++;
//                 if((beforeChar == 'J' || beforeChar == 'K' || beforeChar == 'L' || beforeChar == 'G') && (tchar == 'F' || tchar == 'H' || tchar == 'I' || tchar == 'C')) {
//                     if(writtenCHeader) {
//                         for(int i(0); i< start - (lastCPG+ 1); i++) {
//                             CPGFile << "-";
//                         }
//                         CPGFile << beforeChar << tchar;
//                         lastCPG = start;
//                     } else {
//                         lastCPG = start;
//                         if(firstLineC) { CPGFile << '>' << chrom << ":" << std::to_string(lastCPG-1) << '\n'; firstLineC = false ; }
//                         else { CPGFile << "\n>" << chrom << ":" << std::to_string(lastCPG-1) << '\n'; }
//                         CPGFile << beforeChar << tchar;
//                         writtenCHeader = true;
//                     }
//                     tchar = ':';
//                 } else if(beforeChar != '\n' && beforeChar != ':') {
//                     if(writtenNHeader) {
//                         for(int i(0); i< start - (lastNCPG + 1); i++) {
//                             nCPGFile << "-";
//                         }
//                         nCPGFile << beforeChar << tchar;
//                         lastNCPG = start;
//                     } else {
//                         lastNCPG = start;
//                         if(firstLineN) { nCPGFile  << '>' << chrom << ":" << std::to_string(lastNCPG - 1) << '\n'; firstLineN = false; }
//                         else { nCPGFile << "\n>" << chrom << ":" << std::to_string(lastNCPG - 1) << '\n'; }
//                         nCPGFile << beforeChar << tchar;
//                         writtenNHeader = true;
//                     }
//                 }
//             }
//             beforeChar = tchar;
//         }
//     }
    
//     return 0;
// }
