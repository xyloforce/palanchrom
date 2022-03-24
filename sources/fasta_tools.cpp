#include "fasta_tools.h"
#include <iostream>
#include <fstream>
#include <regex>
#include "bio_tools.h"

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
    // HALF OPEN INT BC begin and len = end - begin so if 1-2 : base 1 and 1 base
    std::string tString;
    tString = m_sequence.substr(begin, end - begin);
    sequence tSeq(tString);
    return tSeq;
}

char sequence::getChar(int index) const {
    return m_sequence[index];
}

fasta_entry::fasta_entry(std::string inputSeq, std::string id, int start = 0, int stop = 0, char strand = 'U', fastaType type = fastaType::standard)
{
    m_sequence = sequence(inputSeq);
    m_header = header(id, start, stop, strand);
    m_type = type;
}

fasta_entry::fasta_entry()
{
    m_sequence = sequence("");
    m_header = header("");
    m_type = standard;
}

fasta_entry::fasta_entry(sequence seq, header head, fastaType type)
{
    m_sequence = seq;
    m_header = head;
    m_type = type;
}

std::string fasta_entry::getMinusStrand()
{
    if(m_header.getStrand() == '+') {
        return m_sequence.getReverseComplement();
    } else if(m_header.getStrand() == '-') {
        return m_sequence.getSequence();
    } else {
        if(m_type == fastaType::bedtools) {
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
        if(m_type == fastaType::bedtools) {
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

std::string sequence::to_string() {
    std::string return_string("");
    for(int i(0); i < m_sequence.size()/80; i ++) {
        return_string += m_sequence.substr(i*80,80) + '\n';
    }
    return return_string;
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

void fasta_entry::write_fasta_entry(std::ofstream& outputFile, fastaType type)
{
    if(type == fastaType::bedtools) {
        outputFile << ">" << m_header.getID() << ":" << m_header.getStart() << "-" << m_header.getEnd() << "(" << m_header.getStrand() << ")" << '\n';
    } else {
        outputFile << ">" << m_header.getID() << '\n';
    }
    outputFile << m_sequence.to_string();
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
        if(m_type == fastaType::bedtools) {
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
    fasta_entry tEntry;
    try {
        sequence tSequence = m_sequence.subsetSequence(begin, end);
        header tHeader = m_header;
        tHeader.setStart(m_header.getStart() + begin);
        tHeader.setEnd(m_header.getEnd() - end);
        tEntry = fasta_entry(tSequence, tHeader, m_type);
    } catch (std::out_of_range) {
        std::cout << "Catched exception, return empty entry" << std::endl;
        std::cout << begin << " : " << end << std::endl;
        std::cout << getHeader() << std::endl;
    }
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
        if(m_type == fastaType::bedtools) {
            std::cout << "Undefined strand, using + as default" << std::endl;
        }
        m_header.setEnd(m_header.getEnd() + editSize);
    }
}

void fasta_entry::editSeq(std::vector <vcf_entry> entries) {
    m_sequence.editSequence(entries);
}

void sequence::editSequence(std::vector <vcf_entry> entries)
{
    for(const auto &entry: entries) {
        int editSize = entry.getAlternate()[0].size();
        if(entry.getRef() == toUpper(m_sequence.substr(entry.getPos()-1, entry.getRef().size()))) {
            m_sequence.replace(entry.getPos()-1, editSize, entry.getAlternate()[0]);
        } else {
            throw std::logic_error("non-match");
        }
    }
}


char fasta_entry::getStrand()
{
    return m_header.getStrand();
}

fasta::fasta()
{
}

fasta::fasta(std::string filename, openType type, fastaType typeH) {
    int index(0);
    m_type = typeH;
    m_read = type;
    if(type == openType::read) {
        m_input = std::ifstream(filename);
        while(!m_input.eof()) {
            fasta_entry entry = readFastaLine();
            m_content.push_back(entry);
            m_indexes[entry.getChrom()] = index;
            index ++;
            std::cout << index << "     \r" << std::flush;
        }
    } else if(type == openType::read_line) {
        m_input = std::ifstream(filename);
    } else if(type == openType::write) {
        m_output = std::ofstream(filename);
    }
}

fasta_entry fasta::getFastaById(std::string id) const {
    std::map<std::string,int>::const_iterator it = m_indexes.find(id);
    if(it != m_indexes.end()) {
        return m_content[it->second];
    } else {
        std::cout << "id : " << id << std::endl;
        throw std::logic_error("id not found");
    }
}

fasta_entry fasta::getSubset(bed_entry entry) {
    fasta_entry entryF = getFastaById(entry.getChrom());
    return entryF.getSubset(entry);
}

fasta_entry fasta_entry::getSubset(bed_entry entry) {
    return subsetEntry(entry.getStart(), entry.getStop());
}

fasta_entry fasta::getSubset(vcf_entry entry) {
    fasta_entry entryF = getFastaById(entry.getChrom());
    return entryF.subsetEntry(entry.getPos() -1, entry.getPos());
}

fasta_entry fasta_entry::getSubset(vcf_entry entry)
{
    return subsetEntry(entry.getPos() - 1, entry.getPos()); // one based vs 0 based
}

bool fasta_entry::isValid(vcf_entry entry)
{
    if(toUpper(m_sequence.subsetSequence(entry.getPos() - 1, entry.getPos()).getSequence()) != toUpper(entry.getRef())) {
        return false;
    } else {
        return true;
    }
}

bool fasta::isValid(vcf_entry entry)
{
    return getFastaById(entry.getChrom()).isValid(entry);
}

std::vector<bool> fasta::areValid(std::vector<vcf_entry> entries)
{
    std::vector <bool> result;
    std::string currentChrom;
    fasta_entry current;
    for(const auto &entry: entries) {
        if(entry.getChrom() == currentChrom) {
            result.push_back(current.isValid(entry));
        } else {
            currentChrom = entry.getChrom();
            current = getFastaById(currentChrom);
            result.push_back(current.isValid(entry));
        }
    }
    return result;
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

    bool skip = false;
    
    while(continueIter && !m_input.eof()) {
        m_input.get(tchar);
        if(tchar == '>' && currentHeader) {
            currentHeader = false;
            while(tchar != '\n') {
                if(m_type == fastaType::bedtools) {
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
                } else if(m_type == fastaType::ucsc) {
                    if(tchar != ' ' && !skip && tchar != '>') {
                        headerF += tchar;
                    } else if(tchar == ' ') {
                        skip = true;
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
    
    if(m_type == fastaType::bedtools) {
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
    return fasta_entry(sequence, headerF, start, stop, strand, m_type);
}

void fasta::write_fasta_entry(fasta_entry entry)
{
    entry.write_fasta_entry(m_output, m_type);
}

bool fasta::isEOF() const
{
    return m_input.eof();
}

void fasta::editSeq(vcf_entry entry)
{
    int index = m_indexes[entry.getChrom()];
    m_content[index].editSeq(std::vector<vcf_entry>(1, entry));
}

std::vector<fasta_entry> fasta::getEntries()
{
    return m_content;
}


void fasta::write_fasta_file(fasta &fileHandler)
{
    for(const auto &entry: fileHandler.getEntries()) {
        write_fasta_entry(entry);
    }
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

std::vector<fasta_entry> fasta::getSeqFromInts(AOEbed &fileI)
{
    return getSeqFromInts(fileI.convertToBed());
}

bool sequence::operator==(const sequence& entry) const {
    return m_sequence == entry.getSequence();
}

bool header::operator==(const header& entry) const
{
    return m_start == entry.getStart() && m_stop == entry.getEnd() && m_chrom == entry.getID() && m_strand == entry.getStrand();
}

sequence fasta_entry::getSequence2() const
{
    return m_sequence;
}

header fasta_entry::getHeader2() const
{
    return m_header;
}


bool fasta_entry::operator==(const fasta_entry& entry) const
{
    return m_sequence == entry.getSequence2() && m_header == entry.getHeader2();
}

void fasta::subsetFromInts(AOEbed &fileI) {
    m_content = std::vector <fasta_entry>();
    if(m_read == read_line) {
        while(!isEOF()) {
            fasta_entry entry = readFastaLine();
            for(const auto &entryB: fileI.getBedByID(entry.getChrom())) {
                fasta_entry entryT = entry.getSubset(entryB);
                if(!(entryT == fasta_entry())) {
                    m_content.push_back(entryT);
                } else {
                    std::cout << entryB.to_string() << std::endl;
                    std::cout << entry.getHeader() << std::endl;
                }
            }
        }
    } else {
        throw std::logic_error("function is for read_line type fasta");
    }
}

int fasta::size() const {
    return m_content.size();
}

fasta_entry fasta::getFastaByIndex(int index) const
{
    return m_content[index];
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
            if(ints[i].getStart() - ints[i-1].getStop() > 0) {
                results.push_back(bed_entry(m_header.getID(), ints[i-1].getStop(), ints[i].getStart()));
            }
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
