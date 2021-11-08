#include "fasta_tools.h"
#include "bio_tools.h"
#include <iostream>
#include <fstream>

header::header(std::string chrom, int start = 0, int stop = 0, char strand = 'U')
{
    m_start = start;
    m_stop = stop;
    m_chrom = chrom;
    m_strand = strand;
}

char header::getStrand()
{
    return m_strand;
}

void header::setStart(int start)
{
    m_start = start;
}

void header::setEnd(int stop)
{
    m_stop = stop;
}

int header::getEnd()
{
    return m_stop;
}

int header::getStart() {
    return m_start;
}

std::string header::getID()
{
    return m_chrom;
}


std::string header::getID_full()
{
    return (m_chrom + std::to_string(m_start) + std::to_string(m_stop) + m_strand);
}

header::header() {}
sequence::sequence() {}

sequence::sequence(std::string sequence)
{
    m_sequence = sequence;
}

std::string sequence::getReverseComplement()
{
    std::string result = "";
    for(int i(0); i < m_sequence.size(); i ++) {
        result = reverseComp(m_sequence[i]) + result;
    }
    return result;
}

std::string sequence::getSequence()
{
    return m_sequence;
}

int sequence::getSize()
{
    return m_sequence.size();
}


void sequence::setSequence(std::string sequence) {
    m_sequence = sequence;
}

sequence sequence::subsetSequence ( int begin, int end )
{
    std::string tString;
    tString = m_sequence.substr(begin, end - begin);
    sequence tSeq(tString);
    return tSeq;
}



fasta_entry::fasta_entry(std::string inputSeq, std::string id, int start = 0, int stop = 0, char strand = 'U', bool bedtools_type = false)
{
    m_sequence = sequence(inputSeq);
    m_header = header(id, start, stop, strand);
    m_bedtools_type = bedtools_type;
}

fasta_entry::fasta_entry()
{
    m_sequence = sequence("");
    m_header = header("");
    m_bedtools_type = false;
}

fasta_entry::fasta_entry(sequence seq, header head, bool bedtools_type)
{
    m_sequence = seq;
    m_header = head;
    m_bedtools_type = bedtools_type;
}

std::string fasta_entry::getMinusStrand()
{
    if(m_header.getStrand() == '+') {
        return m_sequence.getReverseComplement();
    } else if(m_header.getStrand() == '-') {
        return m_sequence.getSequence();
    } else {
        if(m_bedtools_type == true) {
            std::cout << "unitialized strand, using + as default" << std::endl;
        }
        return m_sequence.getReverseComplement();
    }
}

std::string fasta_entry::getPluStrand()
{
    if(m_header.getStrand() == '+') {
        return m_sequence.getSequence();
    } else if(m_header.getStrand() == '-') {
        return m_sequence.getReverseComplement();
    } else {
        if(m_bedtools_type == true) {
            std::cout << "unitialized strand, using + as default" << std::endl;
        }
        return m_sequence.getSequence();
    }
}

std::string fasta_entry::getSequence()
{
    return m_sequence.getSequence();
}

void fasta_entry::trimSequence(int size, int end)
{
    std::string seq = m_sequence.getSequence();
    if(m_header.getStrand() == '+') {
        if(end == 5) {
            m_sequence.setSequence(seq.erase(0, size));
            m_header.setStart(m_header.getStart() + size);
        } else if(end == 3) {
            m_sequence.setSequence(seq.erase((seq.size() - size), seq.size()));
            m_header.setEnd(m_header.getEnd() - size);
        }
    } else if(m_header.getStrand() == '-') {
        if(end == 5) {
            m_sequence.setSequence(seq.erase((seq.size() - size), seq.size()));
            m_header.setEnd(m_header.getEnd() - size);     
        } else if(end == 3) {
            m_sequence.setSequence(seq.erase(0, size));
            m_header.setStart(m_header.getStart() + size);
        }
    }
}

void fasta_entry::write_fasta_entry(std::ofstream& outputFile)
{
    outputFile << m_header.getID() << ":" << m_header.getStart() << "-" << m_header.getEnd() << "(" << m_header.getStrand() << ")" << '\n';
    outputFile << m_sequence.getSequence() << '\n';
}

int fasta_entry::getPos(int pos_sequence)
{
    if(m_header.getStrand() == '+') {
        return (m_header.getStart() + pos_sequence);
    } else if(m_header.getStrand() == '-') {
        return (m_header.getEnd() - pos_sequence);
    } else {
        if(m_bedtools_type) {
            std::cout << "Undefined strand, using + as default" << std::endl;
        }
        return (m_header.getStart() + pos_sequence);
    }
}

std::string fasta_entry::getChrom()
{
    return m_header.getID();
}

fasta_entry fasta_entry::subsetEntry(int begin, int end)
{
    sequence tSequence = m_sequence.subsetSequence(begin, end);
    header tHeader = m_header;
    tHeader.setStart(m_header.getStart() + begin);
    tHeader.setEnd(m_header.getEnd() - end);
    fasta_entry tEntry(tSequence, tHeader, m_bedtools_type);
    
    return tEntry;
}

int fasta_entry::getSize()
{
    return m_sequence.getSize();
}

void fasta_entry::editSeq ( std::string edit, int start, int end )
{
    std::string seq = m_sequence.getSequence();
    seq.replace(start, end - start, edit);
    m_sequence.setSequence(seq);
    int editSize = seq.size() - (end - start);
    if(m_header.getStrand() == '+') {
        m_header.setEnd(m_header.getEnd() + editSize);
    } else if (m_header.getStrand() == '-') {
        m_header.setStart(m_header.getStart() - editSize);
    } else {
        if(m_bedtools_type) {
            std::cout << "Undefined strand, using + as default" << std::endl;
        }
        m_header.setEnd(m_header.getEnd() + editSize);
    }
}

char fasta_entry::getStrand()
{
    return m_header.getStrand();
}

fasta::fasta()
{
}

fasta::fasta(std::string filename, std::string read, bool bedtools_type)
{
    m_bedtools_type = bedtools_type;
    if(read == "read") {
        m_input = std::ifstream(filename);
        read_fasta();
    } else if(read == "read_line") {
        m_input = std::ifstream(filename);
    } else if(read == "write") {
        m_output = std::ofstream(filename);
    }
}


void fasta::read_fasta()
{
    while(!m_input.eof()) {
        m_content.push_back(read_fasta_line());
    }
}

fasta_entry fasta::read_fasta_line()
{
    char tchar = '\0';
    std::string headerF = "";
    int start = 0;
    std::string startS = "";
    int stop = 0;
    std::string stopS = "";
    char strand = 'U';
    std::string sequence = "";
    
    int info = 0;
    bool notStrand = true;
    bool continueIter = true;
    bool currentHeader = true;
    
    while(continueIter && !m_input.eof()) {
        m_input.get(tchar);
        if(tchar == '>' && currentHeader) {
            currentHeader = false;
            while(tchar != '\n') {
                if(m_bedtools_type) {
                    if(tchar == ':' || tchar == '(' || tchar == ')') {
                        info ++;
                    } else if(tchar == '-' && notStrand) {
                        notStrand = false;
                        info ++;
                    } else if(info == 0 && tchar != '>'){
                        headerF += tchar;
                    } else if(info == 1){
                        startS += tchar;
                    } else if(info == 2) {
                        stopS += tchar;
                    } else if(info == 3) {
                        strand = tchar;
                    }
                } else {
                    if(tchar != '>') {
                        headerF += tchar;
                    }
                }
                m_input.get(tchar);
            }
        } else if (tchar == '>') {
            m_input.seekg(-1, std::ios::cur);
            continueIter = false;
        } else if(tchar != '\n') {
            sequence += tchar;
        }
    }
    
    if(m_bedtools_type) {
        if(startS != "" && stopS != "") {
            start = stoi(startS);
            stop = stoi(stopS);
        } else {
            start = 0;
            stop = sequence.size();
        }
    } else {
        start = 0;
        stop = sequence.size();
    }
    return fasta_entry(sequence, headerF, start, stop, strand, m_bedtools_type);
}

void fasta::write_fasta_entry(fasta_entry entry)
{
    entry.write_fasta_entry(m_output);
}

bool fasta::isEOF() const
{
    return m_input.eof();
}

