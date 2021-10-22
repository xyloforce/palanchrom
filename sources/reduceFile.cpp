#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <stdexcept>
#include <regex>
#include <unordered_set>

inline bool exists (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

int main(int argc, char* argv[]) {
    std::ifstream fasta(argv[1]);
    
    if(!exists(argv[2])) {
        // create output
        std::ofstream compressedFasta(argv[2]);
        
        std::unordered_set<std::string> doneChroms;
        std::string fullChromStart[3];
        std::string line;
        bool notFirst = false;
        int beginSeq;
        
        while(std::getline(fasta, line)) {
            std::smatch result;
            if(std::regex_search(line, result, std::regex(">(\\w+):(\\d+)-\\d+"))) { // if header is chrom:start:end
                for(unsigned index = 0; index < result.size(); index++) {
                    fullChromStart[index] = result.str(index);
                }
                if(doneChroms.count(fullChromStart[1]) == 0) { // first time we see this chrom
                    if(notFirst) {
                        compressedFasta << "\n";
                    }
                    doneChroms.insert(fullChromStart[1]);
                    compressedFasta << ">" << fullChromStart[1] << ":" << fullChromStart[2] << "\n";
                    beginSeq = stoi(fullChromStart[2]);
                    notFirst = true;
                }
                else { // then we're on the seq line
                    compressedFasta << ":" << stoi(fullChromStart[2]) - beginSeq << ":";
                }
            } else { // its a seq line
                compressedFasta << line;
            }
        }
    } else {
        throw std::invalid_argument("output file exists ; no overwriting allowed");
    }
}
