#include "bio_tools.h"
#include <iostream>

char reverseComp(char base) {
    base = toupper(base);
    switch(base) {
        case 'A':
            return 'T';
            break;
        case 'C':
            return 'G';
            break;
        case 'T':
            return 'A';
            break;
        case 'G':
            return 'C';
            break;
        case 'N':
            return 'N';
            break;
        default:
            std::cout << std::endl;
            std::cout << base << std::endl;
            throw std::domain_error("unexpected base value");
    }
}

std::string toUpper(std::string lower) {
    std::string result = "";
    for(unsigned int i(0); i<lower.size(); i++) {
        result += toupper(lower[i]);
    }
    return result;
}

std::array <int, 5> countBasesInSequence(fasta_entry entry) {
    std::array <int, 5> results;
    std::string seq = entry.getSequence();
    for(const auto &chara : seq) {
        switch(toupper(chara)) {
            case 'A':
                results[0] ++;
                break;
            case 'C':
                results[1] ++;
                break;
            case 'G':
                results[2] ++;
                break;
            case 'T':
                results[3] ++;
                break;
            case 'N':
                results[4] ++;
                break;
            default:
                std::cout << chara << std::endl;
                throw std::logic_error("Unexpected base");
        }
    }
    return results;
}

std::vector <bed_entry> matchPattern(std::string pattern, fasta_entry entry) {
    std::vector <bed_entry> values;
    int size(pattern.size());
    int wait(0);
    int endPos(size);
    for(int i(0); i<entry.getSize(); i++) {
        if(wait == 0) {
            if(entry.subsetEntry(i, i + size).getSequence() == pattern) {
                endPos = i + size -1;
                for(int j(1); j*size<entry.getSize(); j++) {
                    if(entry.subsetEntry(i + (size*(j-1)), i + (size*j)).getSequence() != pattern) {
                        break;
                    } else {
                        std::cout << j << std::endl;
                        endPos = i+size*j;
                    }
                }
                values.push_back(bed_entry(entry.getChrom(), i, endPos));
                i = endPos - size;
            } else {
                endPos = i + size;
            }
            wait = entry.searchChar(pattern[0], endPos - size + 1); // -2 bc pattern begin
        }
        wait --;
        // obvious solution : take substring for each pos and check
        // but safer : take substring -> check -> then search first letter of pattern -> skip search until finding it again

    }
    return values;
}
