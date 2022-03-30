#ifndef DEF_VCF
#define DEF_VCF
#include <string>
#include <map>
#include <fstream>
#include <vector>
#include "openType.h"

std::string toUpper(std::string lower);

class bed;
class bed_entry;

class vcf_entry {
public:
    vcf_entry(std::string chrom, int pos, std::string id, std::string ref, std::vector <std::string> alt, int qual = 0, std::string filter = ".", std::map <std::string, std::string> info = std::map<std::string, std::string>());
    vcf_entry();
    vcf_entry(bed_entry entry);
    std::vector <std::string> getAlternate() const;
    int getQual() const;
    std::string getRef() const;
    std::string getInfoValue(std::string key) const;
    std::string to_string() const;
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
    std::vector <std::string> m_alt;
    int m_qual;
    std::string m_filter;
    std::map <std::string, std::string> m_info;
};

class vcf {
public:
    vcf(std::string filename, openType type);
    vcf(std::vector <vcf_entry> values, std::map <std::string, std::vector <int>> indexes);
    void vcf_writeline(vcf_entry entry_vcf);
//     void vcf_writelines(std::string filename);
    vcf_entry readVCFLine(bool &warned);
    // std::string isMuted(std::string chrom, int pos, std::string ref_bases);
    std::vector <vcf_entry> getVCFByChrom(std::string chrom);
    std::vector <vcf_entry> readVCFByChrom(std::string chrom);
    std::vector <vcf_entry> getVCFEntries() const;
    std::vector <std::string> getChroms() const;
    vcf_entry getVCFEntry(int index);
    void delEntry(vcf_entry entry, bool updateIndexB = true);
    bool isEOF() const;
    std::vector <bed_entry> convertToBed(std::vector <vcf_entry> entries);
    std::vector <bed_entry> convertToBed();
    void updateIndex();
private:
    std::vector<vcf_entry> m_content;
    std::map <std::string, std::vector<int>> m_indexes;
    std::ifstream m_input;
    std::ofstream m_output;
};

#endif
