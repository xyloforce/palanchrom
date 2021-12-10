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
