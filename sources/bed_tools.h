#ifndef DEF_BED
#define DEF_BED
#include <string>
#include <map>
#include <fstream>
#include <vector>

class bed_entry {
public:
    bed_entry(std::string chrom, int start, int stop, std::string name, int score, char strand);
    int isInside(int pos, int size) const;
    bed_entry();
    int getStart() const;
    int getStop() const;
    bool operator == (const bed_entry& entry) const;
    char getStrand() const;
    std::string getIDFull() const;
    std::string getID() const;
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
    bed();
    bed(std::string filename, bool read);
    void readBed();
    bed_entry redBedLine();
    bed_entry inInt(std::string chrom, int pos, int size);
    std::map <std::string, bed_entry> getBedByID(std::string id);
protected:
    std::vector <bed_entry> m_content;
    bool m_isInit;
    std::ifstream m_input;
    std::ofstream m_output;
};

class sorted_bed: public bed {
public:
    sorted_bed(std::string filename);
    void readBed();
    std::map <std::array <int, 3>, bed_entry> getBedByID(std::string id);
private:
    std::map <std::string, std::map <std::array <int, 3>, bed_entry>> m_indexes;
};

#endif
