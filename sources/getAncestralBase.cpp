#include <iostream>
#include <fstream>
#include <regex>
#include <stdexcept>

class header {
public:
    int getPos(int count) {
        if(m_strand == '+') {
            return (m_start + count);
        } else if(m_strand == '-') {
            return (m_stop - count);
        } else {
            std::cout << "Strand is : " << m_strand << std::endl;
            throw std::logic_error("Undefined strand");
        }
    }
    std::string getChrom() {return m_chrom;}
    header(int start, int stop, char strand, std::string chrom) {
        m_start = start;
        m_chrom = chrom;
        m_strand = strand;
        m_stop = stop;
    }
    header(){};
private:
    int m_start;
    int m_stop;
    char m_strand;
    std::string m_chrom;
};

void vcf_writer(std::ofstream& output, std::string chrom, int pos, char ref, char alt, bool toInit) {
    if(toInit) {
        output<<"##fileformat=VCFv4.2\n";
        output<<"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    }
    output<<chrom<<'\t'<<pos<<'\t'<<'.'<<'\t'<<ref<<'\t'<<alt<<'\t'<<".\t.\t.\n";
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
    int maxLine = 0;
    int currentLine = 0;
    std::cout<<"Getting max lines..."<<std::endl;
    while(!test.eof()) {
        if(test.get() == '\n') {
            maxLine ++;
        }
    }
    std::cout<<"Finished."<<std::endl;
    // open output file
    
    std::ofstream outputFile(argv[argc - 1]);
    bool not_finished = true;
    int start = 0;
    header headers[argc-2];
    bool firstLine = true;
    char character[argc-2];
    int countChars = 0;
    char strand;
    int stop;
    
    while(not_finished) { // not at the end of the file
        for(int file(0); file < argc-2; file++) { // for each file
            char tchar;
            files[file].get(tchar);
            
            if(tchar == '>') {
                // header is >chrom:start-stop
                // we want chrom & start
                int info = 0;
                std::string chrom = "";
                std::string startS = "";
                std::string stopS = "";
                bool notStrand = true;
                
                while(tchar != '\n') {
                    files[file].get(tchar);
                    if(tchar == ':' || tchar == '(' || tchar == ')') {
                        info ++;
                    } else if(tchar == '-' && notStrand) {
                        notStrand = false;
                        info ++;
                    } else if(info == 0){
                        chrom += tchar;
                    } else if(info == 1){
                        startS += tchar;
                    } else if(info == 2) {
                        stopS += tchar;
                    } else if(info == 3) {
                        strand = tchar;
                    }
                }
                start = stoi(startS);
                stop = stoi(stopS);
                headers[file] = header(start, stop, strand, chrom);
                files[file].get(tchar); // get beyond \n
                countChars = 0;
            }
            
            character[file] = tchar;
        } // we now have an array of chars
        
        if(character[0] != '\n' && !files[0].eof()) {
            // ok
            countChars ++;
            char ref = toupper(character[0]);
            // char ref2 = toupper(character[1]);
            char outgroup1 = toupper(character[1]);
            
            bool ancestralDefined = true;
            bool equalRef = true;
            // bool equalRef2 = true;
            char tchar;
            
            if(ref != 'N') {
                for(int charX(1); charX <argc-2; charX++) {
                    tchar = toupper(character[charX]);
                    if(tchar != ref) {equalRef = false;}
                    // if(tchar != ref2) {equalRef2 = false;}
                    if(outgroup1 != tchar) {ancestralDefined = false;}
                }
                
                if(ancestralDefined) {
                    if(!equalRef) { // diff ref1 
                        vcf_writer(outputFile, headers[0].getChrom(), headers[0].getPos(countChars), ref, outgroup1, firstLine);
                        firstLine = false;
                    }
                } else {
                    vcf_writer(outputFile, headers[0].getChrom(), headers[0].getPos(countChars), ref, 'N', firstLine);
                    firstLine = false;
                }
            }
        }
        else if(files[0].eof()) {
            not_finished = false;
        } else {
            if((currentLine/maxLine *100)+1 % 10 == 0) {
                std::cout<<"Progress = "<<currentLine/maxLine *100<<"                \r";
            }
            currentLine ++;
        }
    }
    
    std::cout<<std::endl;
    return 0;
}
