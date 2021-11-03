#ifndef FASTA_TOOLS_INCLUDED
#define FASTA_TOOLS_INCLUDED

#include <string>
#include <vector>
#include <fstream>

class header {
private:
    std::string m_chrom;
    int m_start;
    int m_stop;
    char m_strand;
public:
    char getStrand();
    void setStart(int start);
    void setEnd(int stop);
    int getStart();
    int getEnd();
    std::string getID();
    std::string getID_full();
    header(std::string chrom, int start, int stop, char strand);
    header();
};

class sequence {
private:
    std::string m_sequence;
public:
    sequence(std::string sequence);
    sequence();
    std::string getReverseComplement();
    std::string getSequence();
    void setSequence(std::string sequence);
};

class fasta_entry {
private:
    sequence m_sequence;
    header m_header;
public:
    fasta_entry(std::string inputSeq, std::string id, int start, int stop, char strand);
    fasta_entry();
    std::string getSequence();
    std::string getPluStrand();
    std::string getMinusStrand();
    std::string getChrom();
    void trimSequence(int size, int end);
    void write_fasta_entry(std::ofstream& outputFile);
    int getPos(int pos_sequence);
};

class fasta {
private:
    std::vector <fasta_entry> m_content;
    std::ifstream m_input;
    std::ofstream m_output;
public:
    fasta(std::string filename, std::string read, bool bedtools_type);
    fasta();
    void read_fasta(bool bedtools_type);
    fasta_entry read_fasta_line(bool bedtools_type);
    void write_fasta_entry(fasta_entry entry);
    bool isEOF() const;
};

#endif
