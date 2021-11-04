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
    sequence subsetSequence(int begin, int end);
    int getSize();
};

class fasta_entry {
private:
    sequence m_sequence;
    header m_header;
public:
    fasta_entry(std::string inputSeq, std::string id, int start, int stop, char strand);
    fasta_entry();
    fasta_entry(sequence seq, header head);
    std::string getSequence();
    std::string getPluStrand();
    std::string getMinusStrand();
    std::string getChrom();
    char getStrand();
    void trimSequence(int size, int end);
    void write_fasta_entry(std::ofstream& outputFile);
    int getPos(int pos_sequence);
    fasta_entry subsetEntry(int begin, int end);
    int getSize();
    void editSeq(std::string edit, int start, int end);
};

class fasta {
private:
    std::vector <fasta_entry> m_content;
    std::ifstream m_input;
    std::ofstream m_output;
    bool m_bedtools_type;
public:
    fasta(std::string filename, std::string read, bool bedtools_type);
    fasta();
    void read_fasta();
    fasta_entry read_fasta_line();
    void write_fasta_entry(fasta_entry entry);
    bool isEOF() const;
};

#endif
