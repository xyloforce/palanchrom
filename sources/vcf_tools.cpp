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

bool vcf_entry::operator == (const vcf_entry& entry) const
{
    return (m_pos == entry.getPos() && m_ref == entry.getRef() && m_alt == entry.getAlternate());
}

char vcf_entry::getAlternate() const
{
    return m_alt;
}

std::string vcf_entry::getRef() const
{
    return m_ref;
}

std::string vcf_entry::getChrom() const {
    return m_chrom;
}

std::string vcf_entry::getAttributeString() const
{
    std::string result;
    result = m_chrom + "\t" + std::to_string(m_pos) + "\t" + m_id + "\t" + m_ref + "\t" + m_alt + "\t" + std::to_string(m_qual) + "\t" + m_filter + "\t" + m_info;
    
    return result;
}

void vcf_entry::vcf_writeline(std::ofstream& output) const
{
    output << getAttributeString() << '\n';
}

int vcf_entry::getPos() const
{
    return m_pos;
}

int vcf_entry::getQual() const {
    return m_qual;
}


vcf_entry vcf::readVCFLine()
{
    std::string line;
    char tchar = '\0';
    int col = 1;
    
    std::string chrom = "";
    std::string tpos = "";
    long int pos = 0;
    std::string id = "";
    std::string ref = "";
    char alt;
    std::string tqual = "";
    int qual = 0;
    std::string filter = "";
    std::string info = "";
    
    while(tchar != '\n' && !m_input.eof()) {
        // file is chrom pos id ref alt qual filter info
        m_input.get(tchar);
        
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
        } else if(tchar == '#') {
            while(tchar != '\n') {
                m_input.get(tchar);
            }
        }
    }
    if(col > 6) {
        pos = stol(tpos);
        if(tqual != ".") {
            qual = stoi(tqual);
        } else {
            qual = 0;
        }
        return vcf_entry(chrom, pos, id, ref, alt, qual, filter, info);
    } else {
        std::cout << "Incorrect nmber of cols : " << col << " , skipping empty line" << std::endl;
        return vcf_entry();
    }
}

vcf::vcf(std::string filename, bool read) {
    if(read) {
        int index = 0;
        m_input = std::ifstream(filename);
        while(!m_input.eof()) {
            vcf_entry entry = readVCFLine();
            if(!(entry == vcf_entry())) {
                m_content.push_back(entry);
                std::tuple <int, std::string, char> description(entry.getPos(), entry.getRef(), entry.getAlternate());
                m_indexes[entry.getChrom()][description] = index;
                index ++;
                if(index % 10000 == 0) {
                    std::cout << index << "           \r";
                }
            }
        }
        std::cout << std::endl;
    } else {
        m_output = std::ofstream(filename);
        m_output <<"##fileformat=VCFv4.2\n";
        m_output <<"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    }
}
 
// std::string vcf::isMuted ( std::string chrom, int pos, std::string ref_bases )
// {
//     if (m_content.find(chrom) != m_content.end()) {
//         if (m_content[chrom].find(pos) != m_content[chrom].end()) {
//             if(m_content[chrom][pos].getRef() == toUpper(ref_bases)) {
//                 std::string tstring = "";
//                 tstring += m_content[chrom][pos].getAlternate();
//                 return tstring;
//             } else {
//                 std::cout<<"Pos : "<< pos << " base is "<< ref_bases << " and expected is " << m_content[chrom][pos].getRef() << std::endl;
//                 throw std::logic_error("ref doesnt match current base");
//             }
//         } else {
//             return ref_bases;
//         }
//     } else {
//         return ref_bases;
//     }
// }

void vcf::vcf_writeline(vcf_entry entry_vcf) {
    entry_vcf.vcf_writeline(m_output);
}

std::vector<vcf_entry> vcf::getVCFByID(std::string id)
{
    std::vector <vcf_entry> results;
    for (const auto &pair: m_indexes[id]) {
        results.push_back(m_content[pair.second]);
    }
    return results;
}

std::vector <vcf_entry> vcf::getVCFEntries() const {
    return m_content;
}

bool vcf::isEOF() const {
    return m_input.eof();
}
