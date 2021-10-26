#include "vcf_tools.h"
#include <fstream>
#include <regex>

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

std::string vcf_entry::get_ancestralBase() const
{
    return m_ref;
}

void vcf::vcf_read(std::string filename)
{
    std::ifstream vcfFile(filename);
    std::string line;
    
    while(getline(vcfFile, line)) {
        std::regex pattern ("(\\w+)\\t(\\d+)\\t(.+)\\t(\\w+)\\t(\\w)\\t(.+)\\t(\\w+)\\t(\\w+)");
        std::smatch results;
        if(std::regex_search(line, results, pattern)) {
            vcf_entry entry;
            // cmatch is now chrom pos id ref alt qual filter info
            if(results[6].str() == ".") {
                vcf_entry entry(results[1].str(), stoi(results[2].str()), results[3].str(), results[4].str(), results[5].str()[0], stoi(results[6].str()), results[7].str(), results[8].str());
            } else {
                vcf_entry entry(results[1].str(), stoi(results[2].str()), results[3].str(), results[4].str(), results[5].str()[0], 0, results[7].str(), results[8].str());
            }
            m_content[results[1].str()][stoi(results[2].str())] = entry;
        }
    }
    // init m_content
}



vcf::vcf(std::string filename, bool read) {
    if(read) {
        vcf::vcf_read(filename);
    }
}

std::string vcf::isMuted(std::string chrom, int pos, char base)
{
    if (m_content.find(chrom) != m_content.end()) {
        if (m_content[chrom].find(pos) != m_content[chrom].end()) {
            return m_content[chrom][pos].get_ancestralBase();
        } else {
            std::string tstring = "";
            tstring += base;
            return tstring;
        }
    } else {
        std::string tstring = "";
        tstring += base;
        return tstring;
    }
}


