#include "bed_tools.h"
#include <string>
#include <fstream>
#include <regex>
#include <iostream>

bed_entry::bed_entry ( std::string chrom, int start, int stop, std::string name = ".", int score = 0, char strand = '.' )
{
    m_chrom = chrom;
    m_start = start;
    m_stop = stop;
    m_name = name;
    m_score = score;
    m_strand = strand;
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


int bed_entry::isInside(int pos, int size = 0) const
{
    if(m_start >= pos + size) {
    // case -<->-|-|-----
        return 0;
    } else if(m_start > pos && m_stop >= (pos + size)) {
    // case -<---|>|-----
        return 1;
    } else if(m_start <= pos && m_stop > (pos + size)) {
    // case -|---<>|-----
        return 2;
    } else if(m_start <= pos && m_stop > pos && m_stop <= (pos + size)) {
    // case -|---<|>-----
        return 3;
    } else if(m_stop <= pos) {
    // case -|---|<>-----
        return 4;
    } else if(m_start > pos && m_stop < (pos + size)) {
    // case -<---||>-----
        return 5;
    } else {
        std::cout << "Int is " << m_start << ":" << m_stop << " and pos is " << pos << ":" << pos + size << std::endl;
        throw std::logic_error("Impossible combination of values");
    }
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

std::string bed_entry::getIDFull() const {
    std::string ID = std::to_string(m_start) + " " + std::to_string(m_stop) + " " + m_strand;
    return ID;
}

std::string bed_entry::getID() const {
    return m_chrom;
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

std::vector <bed_entry> sorted_bed::inInt ( std::string chrom, std::vector <std::array <int, 3>> pos, bool stranded = false )
{
// match ints against pos
// BEWARE : pos must be on the same chrom
// BEWARE : actually CONSUMES the bed
// fix it
// works FOR NON OVERLAPING INTs only
    std::map <std::array<int, 3>, bed_entry, compareInts> currentInts = getBedByID(chrom);
    bed_entry output;
    std::vector <bed_entry> output_array;
    
    for(int i(0); i < pos.size(); i++) {
        std::array<int, 3> currentPos = pos[i];
        if(!stranded) {
            for(const auto &pair : currentInts) {
                int start(pos[i][0]);
                int size(pos[i][1]);
                int result = pair.second.isInside(start, size);
                if(result == 0) {
                    break; // smaller than smaller element : stop search
                } else if (result == 1) {
                    output = bed_entry(chrom, pair.second.getStart(), start+size); // overlap is between start of the bed int and the stop of the seeked int
                    break;
                } else if (result == 2) {
                    output = bed_entry(chrom, start, start + size); // int is IN the current int
                    break;
                } else if (result == 3) {
                    output = bed_entry(chrom, start, pair.second.getStop()); // stop is bigger than the current int
                    break;
                } else if (result == 4) {
                // bigger than this int : next
                } else if (result == 5) {
                // current int is inside the int
                    output = bed_entry(chrom, pair.second.getStart(), pair.second.getStop());
                    break;
                }
            }
            output_array.push_back(output);
        }
    }
    return output_array;
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
        std::array <int, 3> key;
        key[0] = entry.getStart();
        key[1] = entry.getStop() - entry.getStart();
        key[2] = entry.getStrand();
        m_indexes[entry.getID()][key] = index;
        m_content.push_back(entry);
        index ++;

        if(index % 1000 == 0) {
            std::cout << index << "         \r";
        }
    }
    std::cout << std::endl;
}

std::map <std::array <int, 3>, bed_entry, compareInts> sorted_bed::getBedByID(std::string id) {
    std::map <std::array <int, 3>, bed_entry, compareInts> output;
    for (const auto &pair : m_indexes[id]) {
        output[pair.first] = m_content[pair.second];
    }
    return output;
}

minimal_sorted_bed::minimal_sorted_bed(std::string filename) {
    m_input = std::ifstream(filename);
    int index = 0;
    while(!m_input.eof()) {
        std::tuple <int, std::string, int, int, char> inputLine = readBedLine();
        m_content.push_back(std::get <0>(inputLine));
        std::array <int, 3> key;
        key[0] = std::get <2>(inputLine);
        key[1] = std::get <3>(inputLine);
        key[3] = std::get <4>(inputLine);
        m_indexes[std::get <1>(inputLine)][key] = index;
        index ++;

        if(index % 100 == 0) {
            std::cout << index << "         \r";
        }
    }
    std::cout << "Finished loading file" << std::endl;
}

std::tuple <int, std::string, int, int, char> minimal_sorted_bed::readBedLine() {
    std::tuple <int, std::string, int, int, char> output;
    std::get <0>(output) = m_input.tellg();
    int col = 1;
    char tchar = '\0';
    char strand;
    std::string chrom = "";
    std::string tstart(""), tstop("");

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
                    break;
                case 5:
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
        std::get<1>(output) = chrom;
        std::get<2>(output) = stoi(tstart);
        std::get<3>(output) = stoi(tstop);
        std::get<4>(output) = '+';
    } else if (col < 3) {
        std::cout << "Empty line, returning empty entry" << std::endl;
        std::get<1>(output) = "";
        std::get<2>(output) = 0;
        std::get<3>(output) = 0;
        std::get<4>(output) = '+';
    } else {
        std::get<1>(output) = chrom;
        std::get<2>(output) = stoi(tstart);
        std::get<3>(output) = stoi(tstop);
        std::get<4>(output) = strand;
    }
    return output;
}

std::map <std::array <int, 3>, bed_entry> minimal_sorted_bed::getBedByID(std::string id) {
    std::map <std::array <int, 3>, bed_entry> output;
    for (const auto &pair : m_indexes[id]) {
        int posLine = m_content[pair.second];
        m_input.seekg(posLine, std::ios::beg);
        bed_entry line = bed::readBedLine();
        output[pair.first] = line;
    }
    return output;
}

// future alg is : for pos in fasta seq :
// if were in interval
// if interval is negative :
// put in place the complement of the base
