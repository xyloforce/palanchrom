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
    std::string getStringEntry(char sep = '\t') const;
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
    std::string to_string() const;

    bool operator == (const AOE_entry& entry) const;
private:
    int m_zero;
    char m_type;
};

class bed {
public:
    bed();
    bed(std::string filename, openType type);
    bed_entry readBedLine(int count);
    std::map <std::string, bed_entry> getBedByID(std::string id) const;
    void writeBedLine(bed_entry entry);
    void writeFullBed();
    void setupAndWriteBed(std::string filename);
    bed_entry getBedEntry(int index);
    std::vector <bed_entry> getEntries() const;
    bool isEOF() const;
    void updateTags();
    bed_entry getBedByTag(std::string tag);
protected:
    std::map <std::string, int> m_tags;
    std::vector <bed_entry> m_content;
    bool m_isInit;
    std::ifstream m_input;
    std::ofstream m_output;
    openType m_type;
};

class sorted_bed: public bed {
public:
    sorted_bed();
    sorted_bed(std::string filename, int skip = 0);
    sorted_bed(std::vector <bed_entry> content);

  // functions
    void readBed();
    bool isInside(bed_entry entry);
    std::vector <bool> areInside (std::vector <bed_entry> entries);
    std::map <bed_entry, std::vector<bed_entry>> getOverlap ( std::string chrom, std::vector <bed_entry> pos);
    std::map <bed_entry, std::vector <bed_entry>> getOverlap (sorted_bed& toOverlap);

  // getters
    std::vector <bed_entry> getBedByID(std::string id) const;
    std::vector <std::string> getChroms() const;

protected:
    std::map <std::string, std::map <std::array <int, 3>, int>> m_indexes;

  // internal functions
    std::map <bed_entry, std::vector< bed_entry >> overlap(const std::vector< bed_entry > intsA, const std::vector< bed_entry > intsB);
    std::vector <bed_entry> intersect(std::vector <bed_entry> source, std::vector <bed_entry> toIntersect, bool fullS = false, bool fullT = false);
};

class AOEbed: public sorted_bed {
public:
    AOEbed(std::string filename, openType oType = read);
    AOEbed(std::vector <AOE_entry> content);
    AOE_entry readAOEline();
    int size() const;
    std::vector <AOE_entry> getContent() const;
    std::map <bed_entry, std::vector<AOE_entry>> getOverlap (sorted_bed& entries);
    std::map <bed_entry, std::vector<AOE_entry>> getOverlap (vcf& entries);
    std::map <bed_entry, std::vector<AOE_entry>> getOverlapLowMem (vcf& entries);
    std::vector <AOE_entry> getBedByID(std::string id);
    AOE_entry getEntryByIndex(int index) const;
    std::vector <AOE_entry> getIntersects(sorted_bed& inputFile, bool fullI = false, bool fullF = false);
    std::vector <AOE_entry> getIntersects(std::vector <bed_entry> input, std::vector <AOE_entry> frame, bool fullI = false, bool fullF = false);
    std::vector <AOE_entry> getIntersects(bed& inputFile, bool fullI = false, bool fullF = false);
    std::map <int, int> countVector(std::vector <AOE_entry> input);
//     std::map <std::string, std::map <int, int>> getCounts(bed &inputFile, bool fullI = false, bool fullF = false, std::string filename = "");
    void cutToMask(bed &mask, bool fullI = false, bool fullF = false);
    void cutToMask(sorted_bed &mask, bool fullI = false, bool fullF = false);
    std::vector <bed_entry> convertToBed(std::vector <AOE_entry> source) const;
    std::vector <bed_entry> convertToBed() const;
    std::vector <AOE_entry> convertBack(std::vector <bed_entry> source);
    void writeToFile(std::string filename, int limit = 0);
    void dumpAOE(int limit);
    void dumpAOE();
    void loadBlock(int size);
    std::vector <AOE_entry> getEntries();
    std::vector <AOE_entry> getSortedEntries();
private:
    std::vector <AOE_entry> m_content;
    std::map <std::string, std::map <std::array <int, 2>, int>> m_indexes;
};

void dump(std::vector <AOE_entry> &data, int limit);

#endif
