#ifndef FASTA_TOOLS_INCLUDED
#define FASTA_TOOLS_INCLUDED

#include <string>
#include <vector>
#include <fstream>
#include <map>
#include <unordered_set>
#include "bed_tools.h"
#include "openType.h"

class header {
private:
    std::string m_chrom;
    int m_start;
    int m_stop;
    char m_strand;
public:
    char getStrand() const;
    void setStart(int start);
    void setEnd(int stop);
    int getStart() const;
    int getEnd() const;
    std::string getID() const;
    std::string getID_full() const;
    header(std::string chrom, int start, int stop, char strand);
    header();
};

class sequence {
private:
    std::string m_sequence;
public:
    sequence(std::string sequence);
    sequence();
    std::string to_string();
    std::string getReverseComplement();
    std::string getSequence() const;
    std::string getUppercaseSequence() const;
    char getChar(int index) const;
    void setSequence(std::string sequence);
    sequence subsetSequence(int begin, int end) const;
    void editSequence(std::vector <vcf_entry> entries);
    int getSize() const;
    int searchChar(char searched, int pos) const;
};

class fasta_entry {
private:
    sequence m_sequence;
    header m_header;
    fastaType m_type;
public:
    fasta_entry(std::string inputSeq, std::string id, int start, int stop, char strand, fastaType type);
    fasta_entry();
    fasta_entry(sequence seq, header head, fastaType type);
    std::string getHeader() const;
    std::string getSequence() const;
    std::string getUppercaseSequence() const;
    std::string getPluStrand();
    std::string getMinusStrand();
    std::string getChrom();
    char getStrand();
    void trimSequence(int size, int end);
    void write_fasta_entry(std::ofstream& outputFile, fastaType type);
    long getPos(long pos_sequence) const;
    fasta_entry subsetEntry(int begin, int end) const;
    int getSize() const;
    void editSeq(std::string edit, int start, int end);
    void editSeq(std::vector <vcf_entry> entries);
    int searchChar(char searched, int pos) const;
    fasta_entry getSubset(bed_entry entry);
    fasta_entry getSubset(vcf_entry entry);
    std::vector <bed_entry> matchPattern(std::string pattern) const;
    std::vector <bed_entry> matchPatterns(std::string pattern) const;
    std::vector <bed_entry> reverseInts (std::vector <bed_entry> ints) const;
    bool isValid(vcf_entry entry);
};

class fasta {
private:
    std::vector <fasta_entry> m_content;
    std::map <std::string, int> m_indexes;
    std::ifstream m_input;
    std::ofstream m_output;
    fastaType m_type;
    openType m_read;
    bool m_warned;
public:
    fasta(std::string filename, openType type, fastaType typeH);
    fasta();
    fasta_entry readFastaLine();
    void editSeq(vcf_entry entry);
    void write_fasta_entry(fasta_entry entry);
    void write_fasta_file(fasta &fileHandler);
    std::vector <fasta_entry> getEntries();
    bool isEOF() const;
    fasta_entry getFastaById(std::string id) const;
    fasta_entry getSubset(bed_entry entry);
    fasta_entry getSubset(vcf_entry entry);
    std::vector <fasta_entry> getSeqFromInts (std::vector <bed_entry> intsOfInterest);
    bool isValid(vcf_entry entry);
    std::vector <bool> areValid(std::vector <vcf_entry> entries);
};

#endif
