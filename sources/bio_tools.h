#ifndef BIOTOOLS_INCLUDED
#define BIOTOOLS_INCLUDED

#include <string>
#include "fasta_tools.h"

char reverseComp(char base);
std::string toUpper(std::string lower);
std::array <int, 5> countBasesInSequence(fasta_entry entry);

#endif
