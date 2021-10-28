#ifndef DEF_PERS
#define DEF_PERS

#include <string>
#include <map>

class vcf_entry {
public:
    vcf_entry(std::string chrom, int pos, std::string id, std::string ref, char alt, int qual = 0, std::string filter = ".", std::string info = ".");
    vcf_entry();
    char get_ancestralBase() const;
    std::string getAttributeString() const;
private:
    std::string m_chrom;
    int m_pos;
    std::string m_id;
    std::string m_ref;
    char m_alt;
    int m_qual;
    std::string m_filter;
    std::string m_info;
};

class vcf {
public:
    vcf(std::string filename, bool read);
    void vcf_writeline(std::string filename);
//     void vcf_writelines(std::string filename);
    void vcf_read(std::string filename);
    char isMuted(std::string chrom, int pos, char base);
private:
    std::map <std::string, std::map<int, vcf_entry>> m_content;
    bool m_isInit;
};

#endif
