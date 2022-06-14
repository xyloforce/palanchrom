#include "bed_tools.h"
#include "vcf_tools.h"

bed_entry::bed_entry ( std::string chrom, int start, int stop, std::string name, int score, char strand)
{
    if(start > stop) {
        std::cout << start << "  " << stop << std::endl;
        std::cout << "Error : start > stop" << std::endl;
    }
    m_chrom = chrom;
    m_start = start;
    m_stop = stop;
    m_name = name;
    m_score = score;
    m_strand = strand;
}

bed_entry::bed_entry ( std::string chrom, int start, int stop)
{
    if(start > stop) {
        std::cout << start << "  " << stop << std::endl;
        std::cout << "Error : start > stop" << std::endl;
    }
    m_chrom = chrom;
    m_start = start;
    m_stop = stop;
    m_name = ".";
    m_score = 0;
    m_strand = '.';
}

bed_entry::bed_entry()
{
    m_chrom = ".";
    m_start = 0;
    m_stop = 0;
    m_name = ".";
    m_score = 0;
    m_strand = '.';
}

bed_entry::bed_entry(vcf_entry entry) {
    m_chrom = entry.getChrom();
    m_start = entry.getPos() - 1; // one-based vcf
    m_stop = entry.getPos(); // half open int
    std::string alt = "";
    for (int i(0); i < entry.getAlternate().size(); i ++) {
        alt += entry.getAlternate()[i];
        alt += ',';
    }
    alt.pop_back();
    m_name = entry.getRef() + ':' + alt + ':' + entry.getID();
    m_score = entry.getQual();
    m_strand = '+';
}


int bed_entry::isInside(int pos, int size = 1) const
{
    int start(pos);
    int stop(pos + size);
// BEWARE : we don't count last base as being in the int so a seq of pos & size 1 is only pos
    // <-> is start/stop
    // |---| is m_start / m_stop
    if(stop <= m_start) {
    // case -<->-|-|-----
        return 0;
    } else if(start >= m_stop) {
    // case -|---|<>-----
        return 4;
    } else if(start >= m_start && stop <= m_stop) {
    // case -|---<>|-----
        return 2;
    } else if(start >= m_start && start < m_stop && stop > m_stop) {
    // case -|---<|>-----
        return 3;
    } else if(start < m_start && stop > m_start && stop <= m_stop) {
    // case -<---|>|-----
        return 1;
    } else if(start <= m_start && stop > m_start && start < m_stop && stop >= m_stop) {
    // case -<---||>-----
        return 5;
    } else {
        std::cout << "Int is " << m_start << ":" << m_stop << " and pos is " << pos << ":" << pos + size << std::endl;
        throw std::logic_error("Impossible combination of values");
    }
}

int bed_entry::isInside(bed_entry entry) const {
    int pos(entry.getStart());
    int size(entry.getStop() - entry.getStart());
    return isInside(pos, size);
}

int bed_entry::getStart() const
{
    return m_start;
}

int bed_entry::getStop() const
{
    return m_stop;
}

char bed_entry::getStrand() const
{
    return m_strand;
}

int bed_entry::getScore() const {
    return m_score;
}

bool bed_entry::operator == (const bed_entry& entry) const
{
    return (m_start == entry.getStart() && m_stop == entry.getStop() && m_chrom == entry.getChrom());
}

bool bed_entry::operator > (const bed_entry& entry) const
{
    if(m_chrom == entry.getChrom()) {
        return (m_start > entry.getStart());
    } else {
        return (m_chrom > entry.getChrom());
    }
}

bool bed_entry::operator < (const bed_entry& entry) const
{
    if(m_chrom == entry.getChrom()) {
        if(m_start == entry.getStart()) {
            if(m_stop == entry.getStop()) {
                return (m_strand < entry.getStrand()) ;
            } else {
                return (m_stop < entry.getStop());
            }
        } else {
            return (m_start < entry.getStart());
        }
    } else {
        return (m_chrom < entry.getChrom());
    }
}

bool bed_entry::operator >= (const bed_entry& entry) const
{
    if(m_chrom == entry.getChrom()) {
        return (m_start >= entry.getStart());
    } else {
        return (m_chrom >= entry.getChrom());
    }
}

bool bed_entry::operator <= (const bed_entry& entry) const
{
    if(m_chrom == entry.getChrom()) {
        return (m_start <= entry.getStart());
    } else {
        return (m_chrom <= entry.getChrom());
    }
}

std::string bed_entry::getIDFull() const {
    std::string ID = std::to_string(m_start) + " " + std::to_string(m_stop) + " " + m_strand;
    return ID;
}

std::string bed_entry::getChrom() const {
    return m_chrom;
}

std::string bed_entry::getStringEntry() const {
    return m_chrom + '\t' + std::to_string(m_start) + '\t' + std::to_string(m_stop) + '\t' + m_name + '\t' + std::to_string(m_score) + '\t' + m_strand;
}

void bed_entry::setName(std::string name) {
    m_name = name;
}

void bed_entry::setStrand(char strand) {
    m_strand = strand;
}

std::string bed_entry::getName() const {
    return m_name;
}

AOE_entry::AOE_entry(std::string chrom, int start, int stop, char type, int zero) {
    if(start > stop) {
        std::cout << start << "  " << stop << std::endl;
        std::cout << "Error : start > stop" << std::endl;
    }
    m_chrom = chrom;
    m_start = start;
    m_stop = stop;
    m_type = type;
    m_zero = zero;
    if(m_type == 'R') {
        m_strand = '-';
    } else {
        m_strand = '+';
    }
}

AOE_entry::AOE_entry() {
    m_chrom = "";
    m_start = 0;
    m_stop = 0;
    m_type = '\0';
    m_strand = '\0';
    m_zero = 0;
}

AOE_entry::AOE_entry(bed_entry entry, int zero) {
    m_chrom = entry.getChrom();
    m_start = entry.getStart();
    m_stop = entry.getStop();
    m_strand = entry.getStrand();
    if(m_strand == '-') {
        m_type = 'R';
    } else {
        m_type = 'L';
    }
    m_zero = zero;
}

int AOE_entry::getRelativePos(int pos) const {
    int relativePos = pos - m_zero;
    if(m_type == 'R') {
        // ------[--{--]---|---[--}--]-----
        // {--]---| is left int : before ] it's "-"
        // |--[---} is right int : before [ it's "+"
        // its a right int so reverse pos
        relativePos = relativePos * -1;
    }
    return relativePos;
}

int AOE_entry::getZero() const {
    return m_zero;
}

char AOE_entry::getType() const {
    return m_type;
}


bed::bed() {}

bed::bed ( std::string filename, openType type )
{
    int count(0);
    if(type == openType::read) {
        m_input = std::ifstream(filename);
        while(!m_input.eof()) {
        // file is chrom start stop name strand
            m_content.push_back(readBedLine(count));
            count ++;
        }
    } else if(type == openType::read_line) {
        m_input = std::ifstream(filename);
    } else if(type == openType::write) {
        m_output = std::ofstream(filename);
    }
}

bed_entry bed::readBedLine(int count) {
    char tchar = '\0';
    std::string tstart = "";
    std::string tstop = "";
    std::string chrom = "";
    int col = 1;
    std::string name = ".";
    std::string tscore = "";
    char strand = '\0';
    int start = 0;
    int stop = 0;
    int score = 0;

    while(tchar != '\n' && !m_input.eof()) {
        // file is chrom start stop name strand
        m_input.get(tchar);
        
        if(tchar != '\t' && tchar != '\n') {
            switch(col) {
                case 1:
                    chrom += tchar;
                    break;
                case 2:
                    tstart += tchar;
                    break;
                case 3:
                    tstop += tchar;
                    break;
                case 4:
                    name += tchar;
                    break;
                case 5:
                    tscore += tchar;
                    break;
                case 6:
                    strand = tchar;
                    col ++;
                    break;
                default:
                    std::cout << "Skipping chars" << std::endl;
            }
        } else if(tchar == '\t') {
            col ++;
        }
    }
    if (col == 3) {
        start = stoi(tstart);
        stop = stoi(tstop);
        //name = ".";
        score = 0;
        strand = '+';
    } else if (col < 3) {
        std::cout << "incorrect line, returning empty entry" << std::endl;
        return bed_entry();
    } else {
        start = stoi(tstart);
        stop = stoi(tstop);
        if(tscore != "" && tscore != ".") {
            score = stoi(tscore);
        } else {
            score = 0;
        }
    }
    if(name == ".") {
        name = std::to_string(count);
    }
    return bed_entry(chrom, start, stop, name, score, strand);
}

void bed::writeBedLine(bed_entry entry) {
    m_output << entry.getStringEntry() << '\n';
}

bed_entry bed::getBedEntry(int index) {
    return m_content[index];
}

bool bed::isEOF() const {
    return m_input.eof();
}

std::map <bed_entry, std::vector<bed_entry>> sorted_bed::overlap ( std::vector <bed_entry> intsA, std::vector <bed_entry> intsB)
{
// BEWARE : output map <bedFromB>:vector<bedFromA> !!
    std::sort(intsB.begin(), intsB.end());
    std::map <bed_entry, std::vector<bed_entry>> matchs;
    int A(0);

    for(const auto &entry : intsB) {
        bool found = false;
        int B(intsA.size() - 1);
        while(A <= B && !found) {
            unsigned int index((A+B)/2);
            int status = intsA[index].isInside(entry.getStart(), entry.getStop() - entry.getStart());
            switch(status) {
                case 0:
                    B = index - 1;
                    break;
                case 4:
                    A = index + 1;
                    break;
                default:
                    // got an overlap
                    if(matchs.find(entry) != matchs.end()) {
                        std::cout << matchs.find(entry) -> first.getStringEntry() << std::endl;
                        std::cout << entry.getStringEntry() << std::endl;
                        throw std::logic_error("id is not unique : already listed in map");
                    }
                    matchs[entry].push_back(intsA[index]);
                    unsigned int indexA(index);
                    unsigned int indexB(index);
                    while((int)indexA - 1 > 0) {
                        indexA --;
                        status = intsA[indexA].isInside(entry.getStart(), entry.getStop() - entry.getStart());
                        if(status != 0 && status != 4) {
                            matchs[entry].push_back(intsA[indexA]);
                        } else {
                            break;
                        }
                    }
                    while(indexB + 1 < intsA.size()) {
                        indexB ++;
                        status = intsA[indexB].isInside(entry.getStart(), entry.getStop() - entry.getStart());
                        if(status != 0 && status != 4) {
                            matchs[entry].push_back(intsA[indexB]);
                        } else {
                            break;
                        }
                    }
                    found = true;
            }
        }
    }

    return matchs;
}

std::map <bed_entry, std::vector<bed_entry>> sorted_bed::getOverlap ( std::string chrom, std::vector <bed_entry> pos)
{
    std::vector <bed_entry> currentInts = getBedByID(chrom);
    return overlap(currentInts, pos);
}

std::vector <std::string> sorted_bed::getChroms() {
    std::vector <std::string> chroms;
    for(const auto &pair: m_indexes) {
        chroms.push_back(pair.first);
    }
    return chroms;
}

std::map <bed_entry, std::vector <bed_entry>> sorted_bed::getOverlap (sorted_bed& toOverlap) {
    std::map <bed_entry, std::vector <bed_entry>> result;
    for(const auto &chrom: getChroms()) {
        std::vector <bed_entry> currentChrom = getBedByID(chrom);
        std::cout << "Overlapping on chrom " << chrom << " with " << currentChrom.size() << " overlaps." << std::endl;
        std::map <bed_entry, std::vector <bed_entry>> tmp = toOverlap.getOverlap(chrom, currentChrom);
        std::cout << "Adding results" << std::endl;
        result.insert(tmp.begin(), tmp.end());
    }
    return result; // keys are bed entries from source, values are from intersected
}

std::vector <bool> sorted_bed::areInside (std::vector <bed_entry> entries) {
    std::vector <bool> matchs;
    std::sort(entries.begin(), entries.end());
    std::string lastChrom = "";
    std::vector <bed_entry> currentInts;
    int A(0);

    for(const auto &entry: entries) {
        if(lastChrom != entry.getChrom()) {
            currentInts = getBedByID(entry.getChrom());
            A = 0;
            lastChrom = entry.getChrom();
        }
        int B(currentInts.size() - 1); // impossible to access the size element : 0-based
        bool found = false;
        while(A <= B && !found) {
            int index ((A+B) / 2);
            int val(currentInts[index].isInside(entry));
            if(val == 2) {
                found = true;
            } else if(val == 4) {
                A = index + 1;
            } else if(val == 0) {
                B = index - 1;
            } else {
                std::cout << "Val is : " << val << std::endl;
                throw std::logic_error("Unexpected overlap");
            }
        }
        matchs.push_back(found);
    }
    return matchs;
}

bool sorted_bed::isInside(bed_entry entry) {
// usable with pos only
    std::vector <bed_entry> currentInts = getBedByID(entry.getChrom());
    int A(0), B(currentInts.size() - 1);
    bool found = false;

    while(A <= B && !found) {
        int index ((A+B) / 2);
        bed_entry selected(currentInts[index]);
        int val(selected.isInside(entry));
        if(val == 2) {
            found = true;
        } else if(val == 4) {
            A = index + 1;
        } else if(val == 0) {
            B = index - 1;
        } else {
            std::cout << "Val is : " << val << std::endl;
            throw std::logic_error("Unexpected overlap");
        }
    }
    return found;
}

std::map<std::string, bed_entry> bed::getBedByID ( std::string id ) const
{
    std::map <std::string, bed_entry> bed;
    for(const auto &entry : m_content) {
        if (entry.getChrom() == id) {
            bed[entry.getIDFull()] = entry;
        }
    }
    return bed;
}

std::vector <bed_entry> bed::getEntries() const {
    return m_content;
}

sorted_bed::sorted_bed() {}

sorted_bed::sorted_bed(std::string filename) {
    m_input = std::ifstream(filename);
    int index(0), count(0);
    while(!m_input.eof()) {
        bed_entry entry = readBedLine(count);
        count ++;
        if(!(entry == bed_entry())) {
            std::array <int, 3> key;
            key[0] = entry.getStart();
            key[1] = entry.getStop() - entry.getStart();
            key[2] = entry.getStrand();
            m_indexes[entry.getChrom()][key] = index;
            m_content.push_back(entry);
            index ++;
        }

        if(index % 10000 == 0) {
            std::cout << index << "         \r";
        }
    }
    std::cout << std::endl;
}

sorted_bed::sorted_bed(std::vector <bed_entry> content) {
    int index(0);
    m_isInit = false;
    for(const auto &entry: content) {
        m_content.push_back(entry);
        std::array <int, 3> key;
        key[0] = entry.getStart();
        key[1] = entry.getStop() - entry.getStart();
        key[2] = entry.getStrand();
        m_indexes[entry.getChrom()][key] = index;
        index++;
    }
}

std::vector <bed_entry> sorted_bed::getBedByID(std::string id) {
    std::vector <bed_entry> output;
    for (const auto &pair : m_indexes[id]) {
        output.push_back(m_content[pair.second]);
    }
    std::sort(output.begin(), output.end());
    return output;
}

AOEbed::AOEbed(std::string filename,openType oType) {
    m_input = std::ifstream(filename);
    m_content.clear();
    m_type = oType;
    int index(0);
    if(oType == read) {
        while(!m_input.eof()) {
            AOE_entry entry = readAOEline();

            if(!(entry == AOE_entry())) {
                m_content.push_back(entry);
                std::array <int, 2> pos = {entry.getStart(),  entry.getStop()};
                m_indexes[entry.getChrom()][pos] = index;
                index ++;
            }
        }
    } else if(oType == read_line) {
        std::cout << "opening in read_line mode" << std::endl;
    } else {
        throw std::logic_error("non-implemented mode");
    }
}

void AOEbed::loadBlock(int size)
{
    m_content.clear();
    m_indexes.clear();
    int count(0);
    if(m_type != read_line) {
        throw std::logic_error("type of handler is incorrect");
    }
    while(!isEOF() && count < size) {
        AOE_entry entry = readAOEline();
        if(!(entry == AOE_entry())) {
            std::array <int, 2> pos = {entry.getStart(), entry.getStop()};
            m_content.push_back(entry);
            m_indexes[entry.getChrom()][pos] = count;
            count ++;
        }
    }
}


AOEbed::AOEbed(std::vector <AOE_entry> content) {
    int index(0);
    m_content = std::vector <AOE_entry>();
    for(const auto &entry: content) {
        m_content.push_back(entry);
        std::array <int, 2> pos;
        pos[0] = entry.getStart();
        pos[1] = entry.getStop();
        m_indexes[entry.getChrom()][pos] = index;
        index ++;
    }
}

AOE_entry AOEbed::readAOEline() {
    char tchar = '\0';
    std::string tstart = "";
    std::string tstop = "";
    std::string chrom = "";
    int col = 1;
    std::string name = "";
    std::string tscore = "";
    char strand = '\0';
    int start = 0;
    int stop = 0;
    std::string tZero = "";
    int zero = 0;

    while(tchar != '\n' && !m_input.eof()) {
        // file is chrom start stop name strand
        m_input.get(tchar);

        if(tchar != '\t' && tchar != '\n') {
            switch(col) {
                case 1:
                    chrom += tchar;
                    break;
                case 2:
                    tstart += tchar;
                    break;
                case 3:
                    tstop += tchar;
                    break;
                case  4:
                    name += tchar;
                    break;
                case 5:
                    tscore += tchar;
                    break;
                case 6:
                    if(tchar == 'R' || tchar == 'L') {
                        strand = tchar;
                    } else {
                        std::cout << "Expected R or L, got " << tchar << std::endl;
                        throw std::domain_error("malformed file, incorrect type definition");
                    }
                    break;
                case 7:
                    tZero += tchar;
                    break;
                default:
                    std::cout << "Skipping chars" << std::endl;
            }
        } else if(tchar == '\t') {
            col ++;
        }
    }
    try {
        start = stoi(tstart);
        stop = stoi(tstop);
        zero = stoi(tZero);
        return AOE_entry(chrom, start, stop, strand, zero);
    } catch (std::invalid_argument exception) {
        return AOE_entry();
    }
}

std::vector <AOE_entry> AOEbed::getBedByID(std::string id) {
    std::vector <AOE_entry> output;
    for (const auto &pair : m_indexes[id]) {
        output.push_back(m_content[pair.second]);
    }
    std::sort(output.begin(), output.end());
    return output;
}


std::map <bed_entry, std::vector<AOE_entry>> AOEbed::getOverlap (sorted_bed& entries) {
    std::string lastChrom = "";
    std::map <bed_entry, std::vector<AOE_entry>> matchs;
    for(const std::string &chrom: entries.getChroms()) {
        std::vector <bed_entry> intsB = entries.getBedByID(chrom);
        std::vector <AOE_entry> intsA = getBedByID(chrom);
        // need to convert AOE entries to bed and back ??
        std::vector <bed_entry> convertedA(convertToBed(intsA));
        std::map <bed_entry, std::vector<bed_entry>> tmp_matchs = overlap(convertedA, intsB);
        for(const auto &entry: tmp_matchs) {
            matchs[entry.first] = convertBack(entry.second);
        }
    }
    return matchs;
}

std::map <bed_entry, std::vector<AOE_entry>> AOEbed::getOverlap (vcf& entries) {
//     std::string lastChrom = "";
    std::map <bed_entry, std::vector<AOE_entry>> matchs;
    for(const std::string &chrom: entries.getChroms()) {
        std::vector <bed_entry> intsB = entries.convertToBed(entries.getVCFByChrom(chrom));
        std::vector <AOE_entry> intsA = getBedByID(chrom);
        std::vector <bed_entry> convertedA(convertToBed(intsA));
        std::map <bed_entry, std::vector<bed_entry>> tmp_matchs = overlap(convertedA, intsB);
        for(const auto &entry: tmp_matchs) {
            matchs[entry.first] = convertBack(entry.second);
        }
    }
    return matchs;
}

std::map <bed_entry, std::vector<AOE_entry>> AOEbed::getOverlapLowMem (vcf& entries) {
    int count(1);
    bool warned = false;
    vcf_entry entryV = entries.readVCFLine(warned);
    std::string lastChrom(entryV.getChrom());
    std::vector <bed_entry> intsB;
    std::map <bed_entry, std::vector<AOE_entry>> matchs;

    while(!entries.isEOF()) {
        std::cout << "reading entries" << std::endl;
        while(count % 100000 != 0 && !entries.isEOF() && lastChrom == entryV.getChrom()) {
            intsB.push_back(bed_entry(entryV));
            entryV = entries.readVCFLine(warned);
            count ++;
        }
        std::cout << "preparing comparison for chrom " << lastChrom << std::endl;
        std::vector <AOE_entry> intsA = getBedByID(lastChrom);
        lastChrom = entryV.getChrom();
        std::vector <bed_entry> convertedA(convertToBed(intsA));
        std::cout << "overlapping" << std::endl;
        std::map <bed_entry, std::vector<bed_entry>> tmp_matchs = overlap(convertedA, intsB);
        std::cout << "saving matches" << std::endl;
        for(const auto &entryB: tmp_matchs) {
            if(matchs.find(entryB.first) != matchs.end()) { // value is in map
                std::vector<AOE_entry> convertedBack = convertBack(entryB.second);
                matchs[entryB.first].insert(matchs[entryB.first].end(), convertedBack.begin(), convertedBack.end());
            } else {
                matchs[entryB.first] = convertBack(entryB.second);
            }
        }
        intsB.clear();
        count ++;
//         std::cout << count << "\r";
    }
    intsB.push_back(bed_entry(entryV));
    std::vector <AOE_entry> intsA = getBedByID(entryV.getChrom());
    std::vector <bed_entry> convertedA(convertToBed(intsA));
    std::map <bed_entry, std::vector<bed_entry>> tmp_matchs = overlap(convertedA, intsB);
    for(const auto &entryB: tmp_matchs) {
        if(matchs.find(entryB.first) != matchs.end()) { // value is in map
            std::vector<AOE_entry> convertedBack = convertBack(entryB.second);
            matchs[entryB.first].insert(matchs[entryB.first].end(), convertedBack.begin(), convertedBack.end());
        } else {
            matchs[entryB.first] = convertBack(entryB.second);
        }
    }
    std::cout << std::endl;
    return matchs;
}

std::vector <bed_entry> sorted_bed::intersect (std::vector <bed_entry> source, std::vector <bed_entry> toIntersect, bool fullS, bool fullT) {
    std::vector <bed_entry> results;
    std::map <bed_entry, std::vector <bed_entry>> matches = overlap(source, toIntersect); // return map bed"toIntersect" : [bed"source"]
    for(const auto &entryToVector: matches) { // entryToVector is from toIntersect
        for(const auto &entry: entryToVector.second) { // entry is from source bed
            int status(entryToVector.first.isInside(entry));
            int start(0);
            int stop(0);
            switch(status) {
                case 1: // stop of entry is inside int && start of entry is outside
                    if(!fullT && !fullS) {
                        start = entryToVector.first.getStart();
                        stop = entry.getStop();
                    }
                    break;
                case 2: // entry is IN entryToVector.first
                    if(!fullT) {
                        start = entry.getStart();
                        stop = entry.getStop();
                    }
                    break;
                case 3:
                    if(!fullT && !fullS) {
                        start = entry.getStart();
                        stop = entryToVector.first.getStop();
                    }
                    break;
                case 5:
                    if(!fullS) {
                        start = entryToVector.first.getStart();
                        stop = entryToVector.first.getStop();
                    }
                    break;
                default:
                    throw std::logic_error("Cant find non-overlapping ints");
            }
            if(start != 0 || stop != 0) {
                results.push_back(bed_entry(entry.getChrom(), start, stop, entryToVector.first.getName(), entry.getScore(), entryToVector.first.getStrand()));
            }
        }
    }
    return results;
}

std::vector <AOE_entry> AOEbed::getIntersects(sorted_bed& inputFile, bool fullI, bool fullF) {
    std::vector <AOE_entry> results;
    for(const auto &chrom: inputFile.getChroms()) {
        std::vector <bed_entry> input = inputFile.getBedByID(chrom);
        std::vector <bed_entry> converted(convertToBed(getBedByID(chrom)));
        std::vector <bed_entry> tmp_intersect = intersect(input, converted, fullI, fullF); // if fullI set : return only
        std::vector <AOE_entry> intersect = convertBack(tmp_intersect);
        results.insert(results.end(), intersect.begin(), intersect.end());
    }
    return results;
}

std::vector <AOE_entry> AOEbed::getIntersects(bed& inputFile, bool fullI, bool fullF) {
    std::vector <AOE_entry> results;
    std::vector <bed_entry> converted;
    std::vector <bed_entry> input;
    int count(0);
    bed_entry entry = inputFile.readBedLine(count);
    count ++;
    std::string last_chrom(entry.getChrom());

    while(!inputFile.isEOF()) {
        while(entry.getChrom() == last_chrom && count % 100000 != 0) {
            input.push_back(entry);
            entry = inputFile.readBedLine(count);
            count ++;
        }
        converted = convertToBed(getBedByID(last_chrom));
        std::vector <bed_entry> tmp_intersect = intersect(input, converted, fullI, fullF);
        std::vector <AOE_entry> intersected = convertBack(tmp_intersect);
        results.insert(results.end(), intersected.begin(), intersected.end());
        input.clear();
        last_chrom = entry.getChrom();
        std::cout << count << " \r";
        count ++;
    }
    converted = convertToBed(getBedByID(last_chrom));
    std::vector <AOE_entry> intersected = convertBack(intersect(std::vector <bed_entry>(1, entry), converted, fullI, fullF));
    results.insert(results.end(), intersected.begin(), intersected.end());

    std::cout << std::endl;
    return results;
}

std::vector <bed_entry> AOEbed::convertToBed(std::vector <AOE_entry> source) const {
    std::vector <bed_entry> currentConv;
    for(auto entry: source) {
        entry.setName(std::to_string(entry.getZero()));
        // entry.setStrand(entry.getType()); no longer needed as AOE entries does now have a strand
        currentConv.push_back(entry);
    }
    return currentConv;
}

// works ONLY with bed created by convertToBed
std::vector <AOE_entry> AOEbed::convertBack(std::vector <bed_entry> source) {
    std::vector <AOE_entry> convertBack;
    for(auto entry: source) {
        convertBack.push_back(AOE_entry(entry, stoi(entry.getName())));
    }
    return convertBack;
}

std::vector<bed_entry> AOEbed::convertToBed() const
{
    return convertToBed(m_content);
}

AOE_entry AOEbed::getEntryByIndex(int index) const
{
    return m_content[index];
}

bool AOE_entry::operator==(const AOE_entry& entry) const
{
    return m_start == entry.getStart() && m_stop == entry.getStop() && m_chrom == entry.getChrom() && m_zero == entry.getZero() && entry.getType() == m_type;
}

std::string AOE_entry::to_string() const {
    return m_chrom + "\t" + std::to_string(m_start) + "\t" +  std::to_string(m_stop)  + "\t" + "." + "\t" + "." + "\t" + m_type + "\t" + std::to_string(m_zero);
}

// void AOEbed::writeToFile(std::string filename)
void AOEbed::writeToFile(std::string filename, int limit)
{
    std::ofstream fileH(filename);
    int count(0);
    for(const auto &entry: m_content) {
        fileH << entry.to_string() << "\n";
        count ++;
        if(count == limit) {
            break;
        }
    }
}

std::vector <AOE_entry> AOEbed::getEntries() {
    return m_content;
}

void AOEbed::dumpAOE(int limit) {
    std::cout << "Dumping..." << std::endl;
    writeToFile("dump.AOE", limit);
    m_content.erase(m_content.begin(), m_content.begin() + limit);
    std::cout << "Updating index..." << std::endl;
    m_indexes = std::map <std::string, std::map <std::array <int, 2>, int>>();
    for(int i(0); i < m_content.size(); i++) {
        std::array <int, 2> currInt = {m_content[i].getStart(), m_content[i].getStop()};
        m_indexes[m_content[i].getChrom()][currInt] = i;
    }
}

int AOEbed::size() const
{
    return m_content.size();
}

std::vector <AOE_entry> AOEbed::getContent() const
{
    return m_content;
}

void AOEbed::cutToMask(bed &mask, bool fullI, bool fullF)
{
    m_content = getIntersects(mask, fullI, fullF);
    std::cout << "Updating index..." << std::endl;
    m_indexes = std::map <std::string, std::map <std::array <int, 2>, int>>();
    for(int i(0); i < m_content.size(); i++) {
        std::array <int, 2> currInt = {m_content[i].getStart(), m_content[i].getStop()};
        m_indexes[m_content[i].getChrom()][currInt] = i;
    }
}

void AOEbed::cutToMask(sorted_bed &mask, bool fullI, bool fullF) {
    m_content = getIntersects(mask, fullI, fullF);
    m_indexes = std::map <std::string, std::map <std::array <int, 2>, int>>();
    for(int i(0); i < m_content.size(); i++) {
        std::array <int, 2> currInt = {m_content[i].getStart(), m_content[i].getStop()};
        m_indexes[m_content[i].getChrom()][currInt] = i;
    }
}

void dump(std::vector <AOE_entry> &data, int limit) {
    std::ofstream dumpH("dump.AOE");
    int count (0);
    for(const auto &entry: data) {
        dumpH << entry.to_string() << "\n";
        count ++;
        if(count == limit) {
            break;
        }
    }
    data.erase(data.begin(), data.begin()+limit);
}
