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
    std::string getID() const;
    std::string getStringEntry() const;
    char getStrand() const;

  // setters
    void setName(std::string name);

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
    void readBed();
    std::vector <bed_entry> getBedByID(std::string id);
    bool isInside(bed_entry entry);
    std::map <bed_entry, std::vector<bed_entry>> overlap ( std::vector <bed_entry> currentInts, std::vector <bed_entry> pos);
    std::map <bed_entry, std::vector<bed_entry>> getOverlap ( std::string chrom, std::vector <bed_entry> pos);
    std::vector <bool> areInside (std::vector <bed_entry> entries);
protected:
    std::map <std::string, std::map <std::array <int, 3>, int>> m_indexes;
};

class AOEbed: public sorted_bed {
public:
	AOEbed(std::string filename);
	AOE_entry readAOEline();
  std::map <bed_entry, std::vector<AOE_entry>> getOverlap(std::string chrom, std::vector <bed_entry> pos);
  std::vector <AOE_entry> getBedByID(std::string id);
private:
	std::vector <AOE_entry> m_content;
	std::map <std::string, std::map <std::array <int, 2>, int>> m_indexes;
};

#endif
