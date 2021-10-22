#include <iostream>
#include <fstream>
#include <regex>
#include <sys/stat.h>
#include <stdexcept>

inline bool exists (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

int main(int argc, char* argv[])
{
    // first open files
    std::ifstream files[argc-2]; // create array of file handlers ; minus one for the command and minus one for the output
    
    for (int arg(1); arg < argc - 1; arg++) { // fill array with handlers
        files[arg-1] = std::ifstream(argv[arg]);
    }
    
    // then get max line for one file (they're all equal)
    std::ifstream test(argv[1]);
    int maxChar = 0;
    int currentChar = 0;
    std::cout << "Getting max lines..."<<std::endl;
    test.seekg(0, std::ios::end);
    maxChar = test.tellg();
	std::cout << "Finished" << std::endl;
    std::cout<<"Progress = 0 %"<<std::endl;
    // open output file
    if(!exists(argv[argc - 1])) {
        std::ofstream outputFile(argv[argc - 1]);
        bool not_finished = true;
        
        while(not_finished) {
            if(currentChar/maxChar * 100 % 10 == 0 && currentChar/maxChar * 100 != 0) {
                std::cout<<"Progress = "<<currentChar/maxChar * 100<<" %"<<std::endl;
            }
            currentChar ++;
            char character[argc-2];
            for(int file(0); file < argc-2; file++) { // for each file
                char tchar;
                files[file].get(tchar);
                
                // beware of ":"
                if(tchar == ':') {
                    char digit = 'X';
                    std::string tstart;
                    while(digit != ':') {
                        files[file].get(digit);
                        if(file == 2 && digit != ':') {
                            tstart = tstart += digit;
                        }
                    }
                    if(file == 2) {
                        outputFile << ":" << tstart << ":"; // write ref2 pos
                    }
                    files[file].get(tchar); // get first char AFTER ":"
                } else if(tchar == '>' || tchar == '\n') {
                    std::string header = "";
                    header += tchar;
                    
                    do { // not at the end of the line
                        files[file].get(tchar);
                        header += tchar;
                    } while(tchar != '\n');
                    
                    if(file == 1) {
                        outputFile<<header<<'\n';
                    }
                    files[file].get(tchar); // get beyond \n
                }
                character[file] = tchar;
            } // we now have an array of chars
            
            if(character[0] != '\n' && !files[0].eof()) {
                // ok
                char ref1 = toupper(character[0]);
                char ref2 = toupper(character[1]);
                char ancestralBase = 'X';
                
                bool equalRef1 = true;
                bool equalRef2 = true;
                
                for(int charX(2); charX < argc -2; charX++) {
                    char currentChar = toupper(character[charX]);
                    if(ref2 != currentChar) {
                        equalRef2 = false;
                    }
                    if(ref1 != currentChar) {
                        equalRef1 = false;
                    }
                }
                if(equalRef2) {
                    ancestralBase = ref2;
                } else if(equalRef1) {
                    switch(ref1) {
                        case 'A':
                            switch(ref2) {
                                case 'C':
                                    ancestralBase = 'B';
                                    break;
                                case 'G':
                                    ancestralBase = 'D';
                                    break;
                                case 'T':
                                    ancestralBase = 'E';
                                    break;
                                case 'N':
                                    ancestralBase = 'F';
                                    break;
                                default:
                                    std::cout<<ref1<<"    "<<ref2<<std::endl;
                                    throw std::domain_error("base value unexpected");
                                    break;
                            }
                            break;
                                case 'C':
                                    switch(ref2) {
                                        case 'A':
                                            ancestralBase = 'H';
                                            break;
                                        case 'G':
                                            ancestralBase = 'I';
                                            break;
                                        case 'T':
                                            ancestralBase = 'J';
                                            break;
                                        case 'N':
                                            ancestralBase = 'K';
                                            break;
                                        default:
                                            std::cout<<ref1<<"    "<<ref2<<std::endl;
                                            throw std::domain_error("base value unexpected");
                                            break;
                                    }
                                    break;
                                        case 'G':
                                            switch(ref2) {
                                                case 'A':
                                                    ancestralBase = 'L';
                                                    break;
                                                case 'C':
                                                    ancestralBase = 'M';
                                                    break;
                                                case 'T':
                                                    ancestralBase = 'O';
                                                    break;
                                                case 'N':
                                                    ancestralBase = 'P';
                                                    break;
                                                default:
                                                    std::cout<<ref1<<"    "<<ref2<<std::endl;
                                                    throw std::domain_error("base value unexpected");
                                                    break;
                                            }
                                            break;
                                                case 'T':
                                                    switch(ref2) {
                                                        case 'A':
                                                            ancestralBase = 'Q';
                                                            break;
                                                        case 'C':
                                                            ancestralBase = 'R';
                                                            break;
                                                        case 'G':
                                                            ancestralBase = 'S';
                                                            break;
                                                        case 'N':
                                                            ancestralBase = 'V';
                                                            break;
                                                        default:
                                                            std::cout<<ref1<<"    "<<ref2<<std::endl;
                                                            throw std::domain_error("base value unexpected");
                                                            break;
                                                    }
                                                    break;
                                                        case 'N':
                                                            ancestralBase = '*';
                                                            break;
                                                        default:
                                                            std::cout<<ref1<<"    "<<ref2<<std::endl;
                                                            throw std::domain_error("base value unexpected");
                                                            break;
                    }
                } else {
                    ancestralBase = '*';
                }
                outputFile << ancestralBase;
            }
            else if(files[0].eof()) {
                not_finished = false;
            }
        }
    } else {
        throw std::invalid_argument("output file exists ; no overwriting allowed");
    } 
    return 0;
}
