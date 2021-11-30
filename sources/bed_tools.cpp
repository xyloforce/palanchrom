#include "bed_tools.h"
#include <string>
#include <fstream>
#include <regex>
#include <iostream>

bed_entry::bed_entry ( std::string chrom, int start, int stop, std::string name, int score, char strand)
{
    m_chrom = chrom;
    m_start = start;
    m_stop = stop;
    m_name = name;
    m_score = score;
    m_strand = strand;
}

bed_entry::bed_entry ( std::string chrom, int start, int stop)
{
    m_chrom = chrom;
    m_start = start;
    m_stop = stop;
    m_name = ".";
    m_score = 0;
    m_strand = '.';
}

bed_entry::bed_entry()
{
    m_chrom = ".";
    m_start = 0;
    m_stop = 0;
    m_name = ".";
    m_score = 0;
    m_strand = '.';
}


int bed_entry::isInside(int pos, int size = 1) const
{
    int start(pos);
    int stop(pos + size);
// BEWARE : we don't count last base as being in the int so a seq of pos & size 1 is only pos
    if(start < m_start && stop <= m_start && start < m_stop && stop < m_stop) {
    // case -<->-|-|-----
        return 0;
    } else if(start < m_start && stop > m_start && start < m_stop && stop <= m_stop) {
    // case -<---|>|-----
        return 1;
    } else if(start >= m_start && stop > m_start && start < m_stop && stop <= m_stop) {
    // case -|---<>|-----
        return 2;
    } else if(start >= m_start && stop > m_start && start < m_stop && stop > m_stop) {
    // case -|---<|>-----
        return 3;
    } else if(start > m_start && stop > m_start && start >= m_stop && stop > m_stop) {
    // case -|---|<>-----
        return 4;
    } else if(start <= m_start && stop > m_start && start < m_stop && stop >= m_stop) {
    // case -<---||>-----
        return 5;
    } else {
        std::cout << "Int is " << m_start << ":" << m_stop << " and pos is " << pos << ":" << pos + size << std::endl;
        throw std::logic_error("Impossible combination of values");
    }
}

int bed_entry::isInside(bed_entry entry) const {
    int pos(entry.getStart());
    int size(entry.getStop() - entry.getStart());
    return isInside(pos, size);
}

int bed_entry::getStart() const
{
    return m_start;
}

int bed_entry::getStop() const
{
    return m_stop;
}

char bed_entry::getStrand() const
{
    return m_strand;
}

bool bed_entry::operator == (const bed_entry& entry) const
{
    return (m_start == entry.getStart() && m_stop == entry.getStop());
}

bool bed_entry::operator > (const bed_entry& entry) const
{
    return (m_start > entry.getStart());
}

bool bed_entry::operator < (const bed_entry& entry) const
{
    return (m_start < entry.getStart());
}

bool bed_entry::operator >= (const bed_entry& entry) const
{
    return (m_start >= entry.getStart());
}

bool bed_entry::operator <= (const bed_entry& entry) const
{
    return (m_start <= entry.getStart());
}

std::string bed_entry::getIDFull() const {
    std::string ID = std::to_string(m_start) + " " + std::to_string(m_stop) + " " + m_strand;
    return ID;
}

std::string bed_entry::getID() const {
    return m_chrom;
}

std::string bed_entry::getStringEntry() const {
    return m_chrom + '\t' + std::to_string(m_start) + '\t' + std::to_string(m_stop) + '\t' + m_name + '\t' + std::to_string(m_score) + '\t' + m_strand;
}

void bed_entry::setName(std::string name) {
    m_name = name;
}

std::string bed_entry::getName() const {
    return m_name;
}

AOE_entry::AOE_entry(std::string chrom, int start, int stop, char type, int zero) {
    m_chrom = chrom;
    m_start = start;
    m_stop = stop;
    m_type = type;
    m_zero = zero;
}

AOE_entry::AOE_entry() {
    m_chrom = "";
    m_start = 0;
    m_stop = 0;
    m_type = '\0';
    m_zero = 0;
}

AOE_entry::AOE_entry(bed_entry entry, int zero) {
    m_chrom = entry.getID();
    m_start = entry.getStart();
    m_stop = entry.getStop();
    m_type = entry.getStrand();
    m_zero = zero;
}

int AOE_entry::getRelativePos(int pos) const {
    int relativePos = m_zero - pos;
    if(m_type == 'L') {
        // its a left int so reverse pos
        relativePos = relativePos * -1;
    }
    return relativePos;
}

int AOE_entry::getZero() const {
    return m_zero;
}

bed::bed() {}

bed::bed ( std::string filename, bool read )
{
    if(read) {
        m_input = std::ifstream(filename);
        while(!m_input.eof()) {
        // file is chrom start stop name strand
            m_content.push_back(readBedLine());
        }
    } else {
        m_output = std::ofstream(filename);
    }
}

bed_entry bed::readBedLine() {
    char tchar = '\0';
    std::string tstart = "";
    std::string tstop = "";
    std::string chrom = "";
    int col = 1;
    std::string name = "";
    std::string tscore = "";
    char strand = '\0';
    int start = 0;
    int stop = 0;
    int score = 0;

    while(tchar != '\n' && !m_input.eof()) {
        // file is chrom start stop name strand
        m_input.get(tchar);
        
        if(tchar != '\t' && tchar != '\n') {
            switch(col) {
                case 1:
                    chrom += tchar;
                    break;
                case 2:
                    tstart += tchar;
                    break;
                case 3:
                    tstop += tchar;
                    break;
                case 4:
                    name += tchar;
                    break;
                case 5:
                    tscore += tchar;
                    break;
                case 6:
                    strand = tchar;
                    col ++;
                    break;
                default:
                    std::cout << "Skipping chars" << std::endl;
                    break;
            }
        } else if(tchar == '\t') {
            col ++;
        }
    }
    if (col == 3) {
        start = stoi(tstart);
        stop = stoi(tstop);
            
        return bed_entry(chrom, start, stop, ".", 0, '+');
    } else if (col < 3) {
        std::cout << "incorrect line, returning empty entry" << std::endl;
        return bed_entry();
    } else {
        start = stoi(tstart);
        stop = stoi(tstop);
        if(tscore != ".") {
            score = stoi(tscore);
        } else {
            score = 0;
        }
        return bed_entry(chrom, start, stop, name, score, strand);
    }
}

void bed::writeBedLine(bed_entry entry) {
    m_output << entry.getStringEntry() << '\n';
}

std::map <bed_entry, std::vector<bed_entry>> sorted_bed::overlap ( std::vector <bed_entry> currentInts, std::vector <bed_entry> pos)
{
    std::map <bed_entry, std::vector<bed_entry>> matchs;
    int A(0), B(currentInts.size());
    bool found = false;

    for(const auto &entry : pos) {
        while(A <= B && !found) {
            unsigned int index((A+B)/2);
            bed_entry current = currentInts[index];
            int status = current.isInside(entry.getStart(), entry.getStop() - entry.getStart());
            if(status == 0) {
                B = index - 1;
            } else if(status == 4) {
                A = index + 1;
            } else {
                // got an overlap
                matchs[entry].push_back(current);
                std::cout << "Searching for more overlaps" << std::endl;
                unsigned int indexA(index);
                unsigned int indexB(index);
                while((int)indexA - 1 > 0) {
                    indexA --;
                    bed_entry current(currentInts[indexA]);
                    std::cout << current.getStringEntry() << std::endl;
                    status = current.isInside(entry.getStart(), entry.getStop() - entry.getStart());
                    if(status != 0 && status != 4) {
                        matchs[entry].push_back(current);
                    } else {
                        break;
                    }
                }
                while(indexB + 1 < currentInts.size()) {
                    indexB ++;
                    bed_entry current = currentInts[indexB];
                    std::cout << current.getStringEntry() << std::endl;
                    status = current.isInside(entry.getStart(), entry.getStop() - entry.getStart());
                    if(status != 0 && status != 4) {
                        matchs[entry].push_back(current);
                    } else {
                        break;
                    }
                }
                found = true;
            }
        }
    }

    return matchs;
}

std::map <bed_entry, std::vector<bed_entry>> sorted_bed::getOverlap ( std::string chrom, std::vector <bed_entry> pos)
{
    std::vector <bed_entry> currentInts = getBedByID(chrom);
    return overlap(currentInts, pos);
}

// TODO create sorted_bed::isInsideSorted which takes vector of pos, sort it
// and keep last index to speed up search

bool sorted_bed::isInside(bed_entry entry) {
// usable with pos only
    std::vector <bed_entry> currentInts = getBedByID(entry.getID());
    int A(0), B(currentInts.size() - 1);
    bool found = false;

    while(A <= B && !found) {
        int index ((A+B) / 2);
        bed_entry selected(currentInts[index]);
        int val(selected.isInside(entry));
        if(val == 2) {
            found = true;
        } else if(val == 4) {
            A = index + 1;
        } else if(val == 0) {
            B = index - 1;
        } else {
            std::cout << "Val is : " << val << std::endl;
            throw std::logic_error("Unexpected overlap");
        }
    }
    return found;
}

std::map<std::string, bed_entry> bed::getBedByID ( std::string id ) const
{
    std::map <std::string, bed_entry> bed;
    for(const auto &entry : m_content) {
        if (entry.getID() == id) {
            bed[entry.getIDFull()] = entry;
        }
    }
    return bed;
}

sorted_bed::sorted_bed() {}

sorted_bed::sorted_bed(std::string filename) {
    m_input = std::ifstream(filename);
    int index = 0;
    while(!m_input.eof()) {
        bed_entry entry = readBedLine();
        if(!(entry == bed_entry())) {
            std::array <int, 3> key;
            key[0] = entry.getStart();
            key[1] = entry.getStop() - entry.getStart();
            key[2] = entry.getStrand();
            m_indexes[entry.getID()][key] = index;
            m_content.push_back(entry);
            index ++;
        }

        if(index % 1000 == 0) {
            std::cout << index << "         \r";
        }
    }
    std::cout << std::endl;
}

std::vector <bed_entry> sorted_bed::getBedByID(std::string id) {
    std::vector <bed_entry> output;
    for (const auto &pair : m_indexes[id]) {
        output.push_back(m_content[pair.second]);
    }
    std::sort(output.begin(), output.end());
    return output;
}

AOEbed::AOEbed(std::string filename) {
    m_input = std::ifstream(filename);
    int index(0);
    while(!m_input.eof()) {
        AOE_entry entry = readAOEline();
        if(!(entry == AOE_entry())) {
            m_content.push_back(entry);
            std::array <int, 2> pos;
            pos[0] = entry.getStart();
            pos[1] = entry.getStop();
            m_indexes[entry.getID()][pos] = index;
            index ++;
        }
    }
}

AOE_entry AOEbed::readAOEline() {
    char tchar = '\0';
    std::string tstart = "";
    std::string tstop = "";
    std::string chrom = "";
    int col = 1;
    std::string name = "";
    std::string tscore = "";
    char strand = '\0';
    int start = 0;
    int stop = 0;
    std::string tZero = "";
    int zero = 0;

    while(tchar != '\n' && !m_input.eof()) {
        // file is chrom start stop name strand
        m_input.get(tchar);

        if(tchar != '\t' && tchar != '\n') {
            switch(col) {
                case 1:
                    chrom += tchar;
                    break;
                case 2:
                    tstart += tchar;
                    break;
                case 3:
                    tstop += tchar;
                    break;
                case 4:
                    name += tchar;
                    break;
                case 5:
                    tscore += tchar;
                    break;
                case 6:
                    strand = tchar;
                    break;
                case 7:
                    tZero += tchar;
                    break;
                default:
                    std::cout << "Skipping chars" << std::endl;
                    break;
            }
        } else if(tchar == '\t') {
            col ++;
        }
    }
    try {
        start = stoi(tstart);
        stop = stoi(tstop);
        zero = stoi(tZero);
    } catch (std::invalid_argument exception) {
        std::cout << "Catched exception, creating empty line" << std::endl;
    }
    return AOE_entry(chrom, start, stop, strand, zero);
}

std::vector <AOE_entry> AOEbed::getBedByID(std::string id) {
    std::vector <AOE_entry> output;
    for (const auto &pair : m_indexes[id]) {
        output.push_back(m_content[pair.second]);
    }
    std::sort(output.begin(), output.end());
    return output;
}


std::map <bed_entry, std::vector<AOE_entry>> AOEbed::getOverlap (std::string chrom, std::vector <bed_entry> pos) {
    //objective : for a list of pos, return for each pos all ints that overlaps with it
    std::vector <AOE_entry> currentInts = getBedByID(chrom);
    std::map <bed_entry, std::vector<AOE_entry>> matchs;
    // need to convert AOE entries to bed and back ??
    std::vector <bed_entry> currentConv;
    for(auto entry: currentInts) {
        entry.setName(std::to_string(entry.getZero()));
        currentConv.push_back(entry);
    }
    std::map <bed_entry, std::vector<bed_entry>> tmp_matchs = overlap(currentConv, pos);
    for(const auto entry: tmp_matchs) {
        for(const auto entry2: entry.second) {
            matchs[entry.first].push_back(AOE_entry(entry2, std::stoi(entry2.getName())));
        }
    }
    return matchs;
}


// future alg is : for pos in fasta seq :
// if were in interval
// if interval is negative :
// put in place the complement of the base
