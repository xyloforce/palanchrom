#include "vcf_tools.h"
#include <fstream>
#include <regex>
#include <string>
#include <iostream>

vcf_entry::vcf_entry(std::string chrom, int pos, std::string id, std::string ref, char alt, int qual, std::string filter, std::string info)
{
    m_chrom = chrom;
    m_pos = pos;
    m_id = id;
    m_ref = ref;
    m_alt = alt;
    m_qual = qual;
    m_filter = filter;
    m_info = info;
}

vcf_entry::vcf_entry()
{
    m_chrom = ".";
    m_pos = 0;
    m_id = ".";
    m_ref = ".";
    m_alt = '.';
    m_qual = 0;
    m_filter = ".";
    m_info = ".";
}


char vcf_entry::get_ancestralBase() const
{
    return m_alt;
}

std::string vcf_entry::getAttributeString() const
{
    std::string result;
    result = m_chrom + std::to_string(m_pos) + m_id + m_ref + m_alt + std::to_string(m_qual) + m_filter + m_info;
    
    return result;
}

void vcf::vcf_read(std::string filename)
{
    std::ifstream vcfFile(filename);
    std::string line;
    char tchar;
    int col = 1;
    
    std::string chrom = "";
    std::string tpos = "";
    int pos = 0;
    std::string id = "";
    std::string ref = "";
    char alt;
    std::string tqual = "";
    int qual = 0;
    std::string filter = "";
    std::string info = "";
    
    while(!vcfFile.eof()) {
        // file is chrom pos id ref alt qual filter info
        vcfFile.get(tchar);
        
        if(tchar != '\n' && tchar != '\t' && tchar != '#') {
            switch(col) {
                case 1:
                    chrom += tchar;
                    break;
                case 2:
                    tpos += tchar;
                    break;
                case 3:
                    id += tchar;
                    break;
                case 4:
                    ref += tchar;
                    break;
                case 5:
                    alt = tchar;
                    break;
                case 6:
                    tqual += tchar;
                    break;
                case 7:
                    filter += tchar;
                    break;
                case 8:
                    info += tchar;
                    break;
                default:
                    std::cout << tchar << std::endl;
                    throw std::domain_error("incorrect number of cols");
                    break;
            }
        } else if(tchar == '\t') {
            col ++;
        } else if(tchar == '\n' && col == 8) {
            col = 1;
            pos = stoi(tpos);
            if(tqual != ".") {
                qual = stoi(tqual);
            } else {
                qual = 0;
            }
            vcf_entry entry(chrom, pos, id, ref, alt, qual, filter, info);
            m_content[chrom][pos] = entry;
        
            chrom = "";
            tpos = "";
            id = "";
            ref = "";
            tqual = "";
            filter = "";
            info = "";
        } else if(tchar == '#') {
            while(tchar != '\n') {
                vcfFile.get(tchar);
            }
        } else {
        }
    }
    // init m_content
}

vcf::vcf(std::string filename, bool read) {
    if(read) {
        vcf::vcf_read(filename);
    }
}

char vcf::isMuted ( std::string chrom, int pos, char base )
{
    if (m_content.find(chrom) != m_content.end()) {
        if (m_content[chrom].find(pos) != m_content[chrom].end()) {
            return m_content[chrom][pos].get_ancestralBase();
        } else {
            return base;
        }
    } else {
        return base;
    }
}

void vcf_writeline(std::string filename, vcf_entry entry_vcf) {
    std::ofstream outputFile(filename, std::ofstream::app);
    outputFile << entry_vcf.getAttributeString();
}
