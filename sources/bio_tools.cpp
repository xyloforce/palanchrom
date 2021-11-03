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
        default:
            std::cout << std::endl;
            std::cout << base << std::endl;
            throw std::domain_error("unexpected base value");
    }
}

std::string toUpper(std::string lower) {
    std::string result = "";
    for(int i(0); i<lower.size(); i++) {
        result += toupper(lower[i]);
    }
    return result;
}
