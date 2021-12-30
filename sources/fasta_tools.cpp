#include "fasta_tools.h"
#include "bio_tools.h"
#include <iostream>
#include <fstream>
#include <regex>

header::header(std::string chrom, int start = 0, int stop = 0, char strand = 'U')
{
    m_start = start;
    m_stop = stop;
    m_chrom = chrom;
    m_strand = strand;
}

char header::getStrand() const
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

int header::getEnd() const
{
    return m_stop;
}

int header::getStart() const {
    return m_start;
}

std::string header::getID() const
{
    return m_chrom;
}


std::string header::getID_full() const
{
    return (m_chrom + " " + std::to_string(m_start) + " " + std::to_string(m_stop) + " " + m_strand);
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
    for(unsigned int i(0); i < m_sequence.size(); i ++) {
        result = reverseComp(m_sequence[i]) + result;
    }
    return result;
}

std::string sequence::getSequence() const
{
    return m_sequence;
}

int sequence::getSize() const
{
    return m_sequence.size();
}


void sequence::setSequence(std::string sequence) {
    m_sequence = sequence;
}

sequence sequence::subsetSequence ( int begin, int end ) const
{
    std::string tString;
    tString = m_sequence.substr(begin, end - begin);
    sequence tSeq(tString);
    return tSeq;
}

char sequence::getChar(int index) const {
    return m_sequence[index];
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

std::string fasta_entry::getSequence() const
{
    return m_sequence.getSequence();
}

std::string sequence::getUppercaseSequence() const
{
    return toUpper(m_sequence);
}

std::string fasta_entry::getUppercaseSequence() const
{
    return m_sequence.getUppercaseSequence();
}

std::string fasta_entry::getHeader() const
{
    return m_header.getID_full();
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

void fasta_entry::write_fasta_entry(std::ofstream& outputFile, bool bedtools_type)
{
    if(bedtools_type) {
        outputFile << ">" << m_header.getID() << ":" << m_header.getStart() << "-" << m_header.getEnd() << "(" << m_header.getStrand() << ")" << '\n';
    } else {
        outputFile << ">" << m_header.getID() << '\n';
    }
    outputFile << m_sequence.getSequence() << '\n';
}

long fasta_entry::getPos(long pos_sequence) const
{
    // beware : were 0-based and need to return 1-based pos
    // for negative strand its ok bc the last pos is not included so it's like zero
    // for plus strand need to add one
    if(m_header.getStrand() == '+') {
        return (m_header.getStart() + pos_sequence +1);
    } else if(m_header.getStrand() == '-') {
        return (m_header.getEnd() - pos_sequence);
    } else {
        if(m_bedtools_type) {
            std::cout << "Undefined strand, using + as default" << std::endl;
        }
        return (m_header.getStart() + pos_sequence +1);
    }
}

std::string fasta_entry::getChrom()
{
    return m_header.getID();
}

fasta_entry fasta_entry::subsetEntry(int begin, int end) const
{
    sequence tSequence = m_sequence.subsetSequence(begin, end);
    header tHeader = m_header;
    tHeader.setStart(m_header.getStart() + begin);
    tHeader.setEnd(m_header.getEnd() - end);
    fasta_entry tEntry(tSequence, tHeader, m_bedtools_type);
    
    return tEntry;
}

int fasta_entry::getSize() const
{
    return m_sequence.getSize();
}

void fasta_entry::editSeq ( std::string edit, int start, int end )
{
    std::string seq = m_sequence.getSequence();
    int editSize = seq.size() - (end - start);
    seq.replace(start, end - start, edit);
    m_sequence.setSequence(seq);
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
    int index(0);
    m_bedtools_type = bedtools_type;
    if(read == "read") {
        m_input = std::ifstream(filename);
        while(!m_input.eof()) {
            fasta_entry entry = readFastaLine();
            m_content.push_back(entry);
            m_indexes[entry.getChrom()] = index;
            index ++;
            std::cout << index << "     \r" << std::flush;
        }
    } else if(read == "read_line") {
        m_input = std::ifstream(filename);
    } else if(read == "write") {
        m_output = std::ofstream(filename);
    }
}

fasta_entry fasta::getFastaById(std::string id) {
    return m_content[m_indexes[id]];
}

fasta_entry fasta::getSubset(bed_entry entry) {
    fasta_entry entryF = getFastaById(entry.getChrom());
    return entryF.getSubset(entry);
}

fasta_entry fasta_entry::getSubset(bed_entry entry) {
    return subsetEntry(entry.getStart(), entry.getStop());
}

fasta_entry fasta::readFastaLine()
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
    entry.write_fasta_entry(m_output, m_bedtools_type);
}

bool fasta::isEOF() const
{
    return m_input.eof();
}

std::vector <fasta_entry> fasta::getSeqFromInts (std::vector <bed_entry> intsOfInterest) {
    std::sort(intsOfInterest.begin(), intsOfInterest.end());
    std::string lastChrom = "";
    fasta_entry entryF;
    std::vector <fasta_entry> results;
    for(const auto &entryB: intsOfInterest) {
        if(lastChrom != entryB.getChrom()) {
            entryF = getFastaById(entryB.getChrom());
            lastChrom = entryB.getChrom();
        }
        results.push_back(entryF.getSubset(entryB));
    }
    return results;
}

int sequence::searchChar(char searched, int pos) const {
    for(int i(pos); i<m_sequence.size(); i++) {
        if(m_sequence[i] == searched) {
            return i - pos + 1; // bc wait will be decremented post-settings
        }
    }
    return m_sequence.size();
}

int fasta_entry::searchChar(char searched, int pos) const {
    return m_sequence.searchChar(searched, pos);
}

std::vector <bed_entry> fasta_entry::matchPattern(std::string pattern) const {
    std::vector <bed_entry> values;
    int size(pattern.size());
    int wait(0);
    int endPos(size);
    for(int i(0); i< m_sequence.getSize(); i++) {
        if(wait == 0) {
            if(subsetEntry(i, i + size).getSequence() == pattern) {
                endPos = i + size -1;
                for(int j(1); j*size<m_sequence.getSize(); j++) {
                    if(subsetEntry(i + (size*(j-1)), i + (size*j)).getUppercaseSequence() != pattern) {
                        break;
                    } else {
                        endPos = i+size*j;
                    }
                }
                values.push_back(bed_entry(m_header.getID(), i, endPos));
                i = endPos - size;
            } else {
                endPos = i + size;
            }
            wait = m_sequence.searchChar(pattern[0], endPos - size + 1); // -2 bc pattern begin
        }
        wait --;
        // obvious solution : take substring for each pos and check
        // but safer : take substring -> check -> then search first letter of pattern -> skip search until finding it again

    }
    return values;
}

std::vector <bed_entry> fasta_entry::matchPatterns(std::string pattern) const {
    std::string sequence(m_sequence.getUppercaseSequence());
    std::regex regex(pattern);
    std::vector <bed_entry> results;

    std::regex_iterator <std::string::iterator> rit (sequence.begin(), sequence.end(), regex);
    std::regex_iterator <std::string::iterator> rend;

    while(rit != rend) {
        int position = rit -> position();
        std::string match = rit -> str();
        int length = position + match.size();
        bed_entry entry(m_header.getID(), position, length, match, 0, '.');
        results.push_back(entry);
        rit ++;
    }

    return results;
}

std::vector <bed_entry> fasta_entry::reverseInts (std::vector <bed_entry> ints) const {
    std::sort(ints.begin(), ints.end());
    std::vector <bed_entry> results;
    for(int i(0); i < ints.size(); i++) {
        if(i == 0) {
            if(ints[i].getStart() - m_header.getStart() > 0) {
                results.push_back(bed_entry(m_header.getID(), m_header.getStart(), ints[i].getStart()));
            }
        } else if(i+1 == ints.size()) {
            if(m_header.getEnd() - ints[i].getStop() > 0) {
                results.push_back(bed_entry(m_header.getID(), ints[i].getStop(), m_header.getEnd()));
            }
        } else {
            if(ints[i].getStart() - ints[i-1].getStop() > 0) {
                results.push_back(bed_entry(m_header.getID(), ints[i-1].getStop(), ints[i].getStart()));
            }
        }
    }
    return results;
}
