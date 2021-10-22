#include <iostream>
#include <fstream>
#include <regex>

int main(int argc, char* argv[]) {
    std::ifstream globalFile(argv[1]);
    std::ofstream CPGFile(argv[2]);
    std::ofstream nCPGFile(argv[3]);
    
    char tchar;
    char beforeChar = '/';
    // int countChars = 0;
    bool writtenCHeader = false;
    bool writtenNHeader = false;
    int lastCPG = 0;
    int lastNCPG = 0;
    bool firstLineC = true;
    bool firstLineN = true;
    
    int start = 0;
    std::string chrom = "";
    
    while(!globalFile.eof()) {
        globalFile.get(tchar);
        if(tchar != '\n' && !globalFile.eof()) {
            // countChars ++;
            if(tchar == '>') {
                std::string header = "";
                while(tchar != '\n') {
                    header += tchar;
                    globalFile.get(tchar);
                }
                std::smatch m;
                std::regex e(">(\\w+):(\\d+)-\\d+");
                regex_search(header, m, e);
                chrom = m.str(1);
                start = stoi(m.str(2)) - 1;
                writtenCHeader = false;
                writtenNHeader = false;
            } else if(tchar == '\n') {
                // countChars = 0;
            } else { // not a > not a '\n' neither a letter in a header
                start ++;
                if((beforeChar == 'J' || beforeChar == 'K' || beforeChar == 'L' || beforeChar == 'G') && (tchar == 'F' || tchar == 'H' || tchar == 'I' || tchar == 'C')) {
                    if(writtenCHeader) {
                        for(int i(0); i< start - (lastCPG+ 1); i++) {
                            CPGFile << "-";
                        }
                        CPGFile << beforeChar << tchar;
                        lastCPG = start;
                    } else {
                        lastCPG = start;
                        if(firstLineC) { CPGFile << '>' << chrom << ":" << std::to_string(lastCPG-1) << '\n'; firstLineC = false ; }
                        else { CPGFile << "\n>" << chrom << ":" << std::to_string(lastCPG-1) << '\n'; }
                        CPGFile << beforeChar << tchar;
                        writtenCHeader = true;
                    }
                    tchar = ':';
                } else if(beforeChar != '\n' && beforeChar != ':') {
                    if(writtenNHeader) {
                        for(int i(0); i< start - (lastNCPG + 1); i++) {
                            nCPGFile << "-";
                        }
                        nCPGFile << beforeChar << tchar;
                        lastNCPG = start;
                    } else {
                        lastNCPG = start;
                        if(firstLineN) { nCPGFile  << '>' << chrom << ":" << std::to_string(lastNCPG - 1) << '\n'; firstLineN = false; }
                        else { nCPGFile << "\n>" << chrom << ":" << std::to_string(lastNCPG - 1) << '\n'; }
                        nCPGFile << beforeChar << tchar;
                        writtenNHeader = true;
                    }
                }
            }
            beforeChar = tchar;
        }
    }
    
    return 0;
}
