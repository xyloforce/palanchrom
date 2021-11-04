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
        return 0;
    } else if(m_start > pos && m_stop >= (pos + size)) {
        return 1;
    } else if(m_start <= pos && m_stop > (pos + size)) {
        return 2;
    } else if(m_start <= pos && m_stop > pos && m_stop <= (pos + size)) {
        return 3;
    } else if(m_stop <= pos) {
        return 4;
    } else if(m_start > pos && m_stop < (pos + size)) {
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

void bed::readBed ( std::string filename )
{
    std::ifstream bedFile(filename);
    std::string line;
    
    std::string chrom = "";
    std::string tstart = "";
    std::string tstop = "";
    int start = 0;
    int stop = 0;
    std::string name = "";
    std::string tscore = "";
    int score = 0;
    char strand;
    
    int col = 1;
    char tchar;
    
    while(!bedFile.eof()) {
        // file is chrom start stop name strand
        bedFile.get(tchar);
        
        if(tchar != '\n' && tchar != '\t') {
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
                default:
                    throw std::domain_error("incorrect number of cols");
                    break;
            }
        } else if(tchar == '\t') {
            col ++;
        } else if(tchar == '\n' && col == 6) {
            col = 1;
            start = stoi(tstart);
            stop = stoi(tstop);
            if(tscore != ".") {
                score = stoi(tscore);
            } else {
                score = 0;
            }
            bed_entry entry(chrom, start, stop, name, score, strand);
            m_content[chrom][start] = entry;
            
            chrom = "";
            tstart = "";
            tstop = "";
            name = "";
            tscore = "";
        } else if(tchar == '\n' && col == 3) {
            col = 1;
            start = stoi(tstart);
            stop = stoi(tstop);
            
            bed_entry entry(chrom, start, stop, ".", 0, '+');
            m_content[chrom][start] = entry;
            
            chrom = "";
            tstart = "";
            tstop = "";
        }
    }
}

bed::bed ( std::string filename, bool read )
{
    if(read) {
        readBed(filename);
    }
}

bed_entry bed::inInt ( std::string chrom, int pos, int size = 0 )
{
    int current = 5;
    bed_entry inInt;
    
    for(auto it = m_content[chrom].cbegin(); it != m_content[chrom].cend();) {
        current = it -> second.isInside(pos, size);
        if(current == 0) { // pos + size is smaller than start of int in a sorted array : pos is out of range
            break;
        } else if(current == 2) { // pos in the right int
            inInt = it -> second;
            break;
        } else if(current == 1 || current == 5) {
            // overlap start = need to adjust start to correct size => set start to start int
            bed_entry overlapCorrected(chrom, it -> second.getStart(), pos + size);
            inInt = overlapCorrected;
            break;
        } else if(current == 3) {
            // same but for stop
            bed_entry overlapCorrected(chrom, pos, it -> second.getStop(), ".", 0, it -> second.getStrand());
            inInt = overlapCorrected;
            break;
        } else if(current == 4) {
            it = m_content[chrom].erase(it); // since we're only incrementing no need to keep past int
        }
    }
    
    return inInt;
    // loop through key : if smaller than start : stop search
    // if bigger than stop : continue
    // if in it : stop
    
}

std::map<int, bed_entry> bed::getBedByID ( std::string id )
{
    return m_content[id];
}



// future alg is : for pos in fasta seq :
// if were in interval
// if interval is negative :
// put in place the complement of the base
