#include "vcf_tools.h"
#include "bed_tools.h"

vcf_entry::vcf_entry(std::string chrom, int pos, std::string id, std::string ref, std::vector <std::string> alt, int qual, std::string filter, std::string info)
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
    m_alt = std::vector <std::string>();
    m_qual = 0;
    m_filter = ".";
    m_info = ".";
}

vcf_entry::vcf_entry(bed_entry entry) {
    m_chrom = entry.getChrom();
    m_pos = entry.getStart() + 1;
    m_id = ".";
    std::string name = entry.getName();
    bool isRef(true);
    int item(0);
    std::string ref;
    std::vector <std::string> alt(1);
    for(int i(0); i < name.size(); i++) {
        if(name[i] == ':') {
            isRef = false;
        }
        if(isRef) {
            ref += name[i];
        } else {
            if(name[i] == ',') {
                alt.push_back(std::string());
                item ++;
            } else if(name[i] != ':') {
                alt[item] += name[i];
            }
        }
    }

    m_ref = ref;
    m_alt = alt;
    m_qual = entry.getScore();
    m_info = ".";
}

bool vcf_entry::operator == (const vcf_entry& entry) const
{
    return (m_pos == entry.getPos() && m_ref == entry.getRef() && m_alt == entry.getAlternate());
}

bool vcf_entry::operator < (const vcf_entry& entry) const
{
    if(m_chrom != entry.getChrom()) {
        return m_chrom < entry.getChrom();
    } else {
        return m_pos < entry.getPos();
    }
}

std::vector <std::string> vcf_entry::getAlternate() const
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

std::string vcf_entry::to_string() const
{
    std::string result;
    std::string alt;
    for(int i(0); i < m_alt.size() ; i ++) {
        alt += m_alt[i] + ",";
    }
    alt.pop_back();
    result = m_chrom + "\t" + std::to_string(m_pos) + "\t" + m_id + "\t" + m_ref + "\t" + alt + "\t" + std::to_string(m_qual) + "\t" + m_filter + "\t" + m_info;
    
    return result;
}

void vcf_entry::vcf_writeline(std::ofstream& output) const
{
    output << to_string() << '\n';
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
    long int pos(0);
    std::string id = "";
    std::string ref = "";
    std::vector <std::string> alt(1);
    int item(0);
    std::string tqual = "";
    int qual(0);
    std::string filter = "";
    std::string info = "";
    
    bool warn = false;
    
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
                    if(tchar == ',') {
                        item ++;
                        alt.push_back(std::string());
                    } else {
                        alt[item] += tchar;
                    }
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
                    warn = true;
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
        if(warn) {
            std::cout << "VCF has more columns than default : " << col << std::endl;
        }
        return vcf_entry(chrom, pos, id, ref, alt, qual, filter, info);
    } else {
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
                m_indexes[entry.getChrom()].push_back(index);
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

void vcf::vcf_writeline(vcf_entry entry_vcf) {
    entry_vcf.vcf_writeline(m_output);
}

std::vector<vcf_entry> vcf::getVCFByChrom(std::string chrom)
{
    std::vector <vcf_entry> results;
    for (const int &index: m_indexes[chrom]) {
        results.push_back(m_content[index]);
    }
    return results;
}

std::vector <vcf_entry> vcf::getVCFEntries() const {
    return m_content;
}

vcf_entry vcf::getVCFEntry(int index) {
    return m_content[index];
}

bool vcf::isEOF() const {
    return m_input.eof();
}

std::vector <bed_entry> vcf::convertToBed(std::vector <vcf_entry> entries) {
    std::vector <bed_entry> converted;
    for(const auto &entry: entries) {
        converted.push_back(bed_entry(entry));
    }
    return converted;
}

std::vector <bed_entry> vcf::convertToBed() {
    std::vector <bed_entry> converted;
    for(const auto &entry: m_content) {
        converted.push_back(bed_entry(entry));
    }
    return converted;
}

std::vector <std::string> vcf::getChroms() const {
    std::vector <std::string> results;
    for(const auto &pair: m_indexes) {
        results.push_back(pair.first);
    }
    return results;
}
