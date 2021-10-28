#ifndef DEF_BED
#define DEF_BED
#include <string>
#include <map>

class bed_entry {
public:
    bed_entry(std::string chrom, int start, int stop, std::string name, int score, char strand);
    int isInside(int pos, int size) const;
    bed_entry();
    int getStart() const;
    int getStop() const;
    bool operator == (const bed_entry& entry) const;
    char getStrand() const;
private:
    std::string m_chrom;
    int m_start;
    int m_stop;
    std::string m_name;
    int m_score;
    char m_strand;
};

class bed {
public:
    bed(std::string filename, bool read);
    void readBed(std::string filename);
    bed_entry inInt(std::string chrom, int pos, int size);
private:
    std::map <std::string, std::map<int, bed_entry>> m_content;
    bool m_isInit;
};

#endif
