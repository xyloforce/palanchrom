#ifndef BIOTOOLS_INCLUDED
#define BIOTOOLS_INCLUDED

#include <string>
#include "fasta_tools.h"
#include "bed_tools.h"

char reverseComp(char base);
std::string reverseComp(std::string sequence);
std::string toUpper(std::string lower);
std::array <int, 5> countBasesInSequence(fasta_entry entry);
// std::vector <std::string> addNEachPos(std::string toAdd);
// std::vector <std::string> addN(std::string toAdd);
std::string constructRegex(std::string pattern, bool addN = false);
std::map <char, std::string> getArgs(std::vector<std::string> args);

#endif
