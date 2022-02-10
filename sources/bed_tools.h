#ifndef DEF_BED
#define DEF_BED
#include "vcf_tools.h"
#include <string>
#include <fstream>
#include <regex>
#include <iostream>

class bed_entry {
public:
    bed_entry(std::string chrom, int start, int stop, std::string name, int score, char strand);
    bed_entry ( std::string chrom, int start, int stop);
    bed_entry();
    bed_entry(vcf_entry entry);

  // functions
    int isInside(int pos, int size) const;
    int isInside(bed_entry entry) const;

  // getters
    int getStart() const;
    int getStop() const;
    std::string getName() const;
    std::string getIDFull() const;
    std::string getChrom() const;
    std::string getStringEntry() const;
    int getScore() const;
    char getStrand() const;

  // setters
    void setName(std::string name);
    void setStrand(char strand);

  // operators
    bool operator == (const bed_entry& entry) const;
    bool operator > (const bed_entry& entry) const;
    bool operator < (const bed_entry& entry) const;
    bool operator >= (const bed_entry& entry) const;
    bool operator <= (const bed_entry& entry) const;

protected:
    std::string m_chrom;
    int m_start;
    int m_stop;
    std::string m_name;
    int m_score;
    char m_strand;
};

class AOE_entry: public bed_entry {
public:
    AOE_entry(std::string chrom, int start, int stop, char type, int zero);
    AOE_entry();
    AOE_entry(bed_entry entry, int zero);
    int getRelativePos(int pos) const;
    int getZero() const;
    char getType() const;

    bool operator == (const AOE_entry& entry) const;
private:
    int m_zero;
    char m_type;
};

class bed {
public:
    bed();
    bed(std::string filename, bool read);
    bed_entry readBedLine();
    std::map <std::string, bed_entry> getBedByID(std::string id) const;
    void writeBedLine(bed_entry entry);
    bed_entry getBedEntry(int index);
    std::vector <bed_entry> getEntries() const;
protected:
    std::vector <bed_entry> m_content;
    bool m_isInit;
    std::ifstream m_input;
    std::ofstream m_output;
};

class sorted_bed: public bed {
public:
    sorted_bed();
    sorted_bed(std::string filename);
    sorted_bed(std::vector <bed_entry> content);

  // functions
    void readBed();
    bool isInside(bed_entry entry);
    std::vector <bool> areInside (std::vector <bed_entry> entries);
    std::map <bed_entry, std::vector<bed_entry>> getOverlap ( std::string chrom, std::vector <bed_entry> pos);
    std::map <bed_entry, std::vector <bed_entry>> getOverlap (sorted_bed& toOverlap);

  // getters
    std::vector <bed_entry> getBedByID(std::string id);
    std::vector <std::string> getChroms();

protected:
    std::map <std::string, std::map <std::array <int, 3>, int>> m_indexes;

  // internal functions
    std::map <bed_entry, std::vector<bed_entry>> overlap ( std::vector <bed_entry> currentInts, std::vector <bed_entry> pos);
    std::vector <bed_entry> intersect(std::vector <bed_entry> source, std::vector <bed_entry> toIntersect, bool fullS = false, bool fullT = false);
};

class AOEbed: public sorted_bed {
public:
    AOEbed(std::string filename);
    AOEbed(std::vector <AOE_entry> content);
    AOE_entry readAOEline();
    std::map <bed_entry, std::vector<AOE_entry>> getOverlap (sorted_bed& entries);
    std::map <bed_entry, std::vector<AOE_entry>> getOverlap (vcf& entries);
    std::vector <AOE_entry> getBedByID(std::string id);
    AOE_entry getEntryByIndex(int index) const;
    std::vector <AOE_entry> getIntersects(sorted_bed& inputFile, bool fullI = false, bool fullF = false);
    std::vector <bed_entry> convertToBed(std::vector <AOE_entry> source) const;
    std::vector <bed_entry> convertToBed() const;
    std::vector <AOE_entry> convertBack(std::vector <bed_entry> source);
private:
    std::vector <AOE_entry> m_content;
    std::map <std::string, std::map <std::array <int, 2>, int>> m_indexes;
};

#endif
