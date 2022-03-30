#include "vcf_tools.h"
#include "bed_tools.h"


vcf_entry::vcf_entry(std::string chrom, int pos, std::string id, std::string ref, std::vector <std::string> alt, int qual, std::string filter, std::map <std::string, std::string> info)
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
    m_info = std::map <std::string, std::string> ();
}

vcf_entry::vcf_entry(bed_entry entry) {
//     std::cout << entry.getChrom() << std::endl;
    m_chrom = entry.getChrom();
    m_pos = entry.getStart() + 1;
//     std::cout << entry.getStart() + 1 << std::endl;
    m_id = ".";
    std::string name = entry.getName();
//     std::cout << entry.getName() << std::endl;
    bool isRef(true);
    int item(0);
    std::string ref("");
    std::vector <std::string> alt(1, std::string());
    for(int i(0); i < name.size(); i++) {
        if(name[i] == ':') {
            isRef = false;
        }
        if(isRef) {
            ref += name[i];
        } else {
            if(name[i] == ',') {
                alt.push_back(std::string(""));
                item ++;
            } else if(name[i] != ':') {
                alt[item] += name[i];
            }
        }
    }

    m_ref = ref;
    m_alt = alt;
    m_qual = entry.getScore();
    m_filter = ".";
    m_info = std::map <std::string, std::string>();
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
    std::string info;
    for(int i(0); i < m_alt.size() ; i ++) {
        alt += m_alt[i] + ",";
    }
    alt.pop_back();
    if(!m_info.empty()) {
        for(const auto &pair: m_info) {
            info += pair.first;
            info += "=";
            info += pair.second;
            info += ";";
        }
        info.pop_back();
    } else {
        info = ".";
    }
    result = m_chrom + "\t" + std::to_string(m_pos) + "\t" + m_id + "\t" + m_ref + "\t" + alt + "\t" + std::to_string(m_qual) + "\t" + m_filter + "\t" + info;
    
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

std::string vcf_entry::getInfoValue(std::string key) const {
     std::map <std::string, std::string>::const_iterator it = m_info.find(key);
    if (it != m_info.end()) {
        return it->second;
    } else {
        return "";
    }
}

vcf_entry vcf::readVCFLine(bool &warned)
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
    std::map <std::string, std::string> info;
    bool warn = false;
    bool isValue = false;
    std::string key = "";
    std::string value = "";

    bool parse_info = false;
    
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
                    switch(tchar) {
                        case ';':
                            isValue = false;
                            info[key] = value;
                            key = "";
                            value = "";
                            parse_info = false;
                            break;
                        case '=':
                            isValue = true;
                            parse_info = true;
                            break;
                        default:
                            if(isValue) {
                                value += tchar;   
                            } else {
                                key += tchar;
                            }
                            break;
                    }
                    break;
                default:
                    warn = true;
                    break;
            }
        } else if(tchar == '\t') {
            col ++;
            if(parse_info) {
                isValue = false;
                info[key] = value;
                key = "";
                value = "";
                parse_info = false;
            }
        } else if(tchar == '#') {
            while(tchar != '\n') { // read whole line until end then return empty entry
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
        if(warn && !warned) {
            std::cout << "VCF has more columns than default : " << col << std::endl;
            warned = true;
        }
        return vcf_entry(chrom, pos, id, ref, alt, qual, filter, info);
    } else {
        return vcf_entry();
    }
}

vcf::vcf(std::string filename, openType type)
{
    if(type == openType::read) {
        bool warned = false;
        int index = 0;
        m_input = std::ifstream(filename);
        while(!m_input.eof()) {
            vcf_entry entry = readVCFLine(warned);
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
    } else if(type == openType::read_line) {
        m_input = std::ifstream(filename);
    } else {
        m_output = std::ofstream(filename);
        m_output <<"##fileformat=VCFv4.2\n";
        m_output <<"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    }
}

vcf::vcf(std::vector<vcf_entry> values, std::map<std::string, std::vector <int>> indexes)
{
    m_content = values;
    m_indexes = indexes;
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

std::vector <vcf_entry> vcf::readVCFByChrom(std::string chrom) {
    bool warned = false;
    vcf_entry entry(readVCFLine(warned));
    std::vector <vcf_entry> results;
    int current(0);
    int count(0);
    if(entry.getChrom() == chrom) {
        while(entry.getChrom() == chrom) {
            results.push_back(entry);
            current = m_input.tellg();
            vcf_entry entry(readVCFLine(warned));
            count ++;
            if(count % 10000 == 0) {
                std::cout << count << "\r";
            }
            std::cout << std::endl;
        }
        if(entry.getChrom() != chrom) {
            m_input.seekg(current, std::ios::beg);
        }
    } else {
        throw std::logic_error("ID found doesn't match current - either badly sorted or fasta entry without any ref");
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
void vcf::delEntry(vcf_entry entry, bool updateIndexB)
{
    std::vector <int> indexes = m_indexes[entry.getChrom()];
    for(int i(0); i < indexes.size(); i++) {
        if(m_content[indexes[i]] == entry) {
            m_content.erase(m_content.begin() + indexes[i]);
            if(updateIndexB) {
                updateIndex();
            }
            break;
        }
    }
}

void vcf::updateIndex() {
    m_indexes.clear();
    for(int i(0); i < m_content.size(); i++) {
        m_indexes[m_content[i].getChrom()].push_back(i);
    }
}
