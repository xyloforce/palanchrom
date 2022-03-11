#include "bio_tools.h"
#include <iostream>
#include <set>

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

// std::vector <std::string> addN(std::string toAdd) {
//     std::vector <std::string> addedN;
//     for(int i(0); i < toAdd.size(); i++) {
//         std::string tmp(toAdd);
//         if(tmp[i] != 'N') {
//             tmp[i] = 'N';
//             addedN.push_back(tmp);
//         }
//     }
//     return addedN;
// }

// std::vector <std::string> addNEachPos(std::string toAdd) {
//     std::vector <std::string> results;
//     std::vector <std::string> tmp (1, toAdd);
//     std::vector <std::string> tmp2;
//     std::string allN;
//     for(int i(0); i < toAdd.size(); i++) {
//         allN += 'N';
//     }
//     for(int j(0); j < tmp.size(); j++) {
//         if(tmp[j] == allN) {
//             break;
//         }
//         tmp2 = addN(tmp[j]);
//         results.insert(results.end(), tmp2.begin(), tmp2.end());
//         tmp.insert(tmp.end(), tmp2.begin(), tmp2.end());
//     }
//     results.resize(std::distance(results.begin(), std::unique(results.begin(), results.end())));
//     results.push_back(toAdd);
//     return results;
// }

// std::string constructRegex(std::vector <std::string> patterns) {
//     std::string finalRegex = "";
//     for(int i(0); i < patterns.size(); i++) {
//         finalRegex += patterns[i];
//         if(i + 1 != patterns.size()) {
//             finalRegex += "|";
//         }
//     }
//     return finalRegex;
// }

std::string constructRegex(std::string pattern, bool addN) {
    // assumption : size is equal for all patterns
    std::string result;
    std::vector <std::set <char>> pos;
    int sizeMotif(0);
    for(int i(0); i < pattern.size(); i++) {
        std::cout << pattern[i] << std::endl;
        if(pattern[i] == ',') {
            sizeMotif = 0;
        } else if(pattern[i] == '-') {
            for(const auto &set: pos) {
                result += "[";
                for(const auto &character: set) {
                    result += character;
                }
                result += "]";
            }
            result += "|";
            pos.clear();
            sizeMotif = 0;
        } else {
            if(pos.size() <= sizeMotif) {
                pos.push_back(std::set <char> ());
                if(addN) {
                    pos[pos.size() - 1].insert('N');
                }
            }
            pos[sizeMotif].insert(pattern[i]);
            sizeMotif ++;
        }
    }
    for(const auto &set: pos) {
        result += "[";
        for(const auto &character: set) {
            result += character;
        }
        result += "]";
    }
    return result;
}
