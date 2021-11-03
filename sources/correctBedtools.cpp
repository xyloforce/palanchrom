#include "fasta_tools.h"
#include <fstream>
#include <iostream>

int main(int argc, char* argv[]) {
    fasta inputFile(argv[1], "read_line", true);
    fasta outputFile(argv[2], "write", true);
    fasta_entry entry;
    
    std::cout<< "Started trimming..." << std::endl;
    
    while(!inputFile.isEOF()) {
        entry = inputFile.read_fasta_line(true);
        entry.trimSequence(1,3);
        outputFile.write_fasta_entry(entry);
    }
    return 0;
}
