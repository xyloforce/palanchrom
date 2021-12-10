#ifndef DEF_VCF
#define DEF_VCF
#include <string>
#include <map>
#include <fstream>
#include <vector>

std::string toUpper(std::string lower);

class bed;
class bed_entry;

class vcf_entry {
public:
    vcf_entry(std::string chrom, int pos, std::string id, std::string ref, char alt, int qual = 0, std::string filter = ".", std::string info = ".");
    vcf_entry();
    vcf_entry(bed_entry entry);
    char getAlternate() const;
    int getQual() const;
    std::string getRef() const;
    std::string getAttributeString() const;
    std::string getChrom() const;
    int getPos() const;
    void vcf_writeline(std::ofstream& output) const;
    bool operator == (const vcf_entry& entry) const;
    bool operator < (const vcf_entry& entry) const;

private:
    std::string m_chrom;
    long m_pos;
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
    void vcf_writeline(vcf_entry entry_vcf);
//     void vcf_writelines(std::string filename);
    vcf_entry readVCFLine();
    // std::string isMuted(std::string chrom, int pos, std::string ref_bases);
    std::vector <vcf_entry> getVCFByID(std::string id);
    std::vector <vcf_entry> getVCFEntries() const;
    std::vector <std::string> getChroms() const;
    vcf_entry getVCFEntry(int index);
    bool isEOF() const;
    std::vector <bed_entry> convertToBed(std::vector <vcf_entry> entries);
    std::vector <bed_entry> convertToBed();
private:
    std::vector<vcf_entry> m_content;
    std::map <std::string, std::map <std::tuple <int, std::string, char>, int>> m_indexes;
    bool m_isInit;
    std::ifstream m_input;
    std::ofstream m_output;
};

#endif
