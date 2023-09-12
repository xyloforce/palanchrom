#include <iostream>
#include "fasta_tools.hpp"
#include "bed_tools.hpp"
#include "tools.hpp"

int main(int argc, char* argv[]) {
    std::map <char, std::string> args = getArgs(std::vector <std::string> (argv, argv + argc));
    std::string fasta_filename, matchs_filename, non_matchs_filename, pattern;
    std::cout << "Starting" << std::endl;
    try {
        fasta_filename = args.at('f');
        matchs_filename = args.at('1');
        pattern = args.at('p');
    } catch (std::out_of_range) {
        std::cout << "Missing obligatory parameters. Parameters are :\n";
        std::cout << "\t+f fasta filename\n";
        std::cout << "\t+1 matchs filename\n";
        std::cout << "\t+p pattern\n";
        std::cout << "Optionnal arg :\n";
        std::cout << "\t+2 non-matchs filename\n";
        throw std::out_of_range("Missing arguments");
    }
    std::cout << "Reading file\n";
    fasta_file sequences(fasta_filename, read, standard);
    sequences.readWholeFile();
    bed_file matchs(matchs_filename, write);

    std::cout << "Looking for patterns\n";
    std::map<int, std::vector<std::shared_ptr <bio_entry>>> results(sequences.matchPatterns(pattern));

    std::cout << "Saving patterns\n";
    matchs.appendVector(results[0]);
    matchs.writeToFile();

    try {
        if(pattern.size() != 2) {
            std::cout << "Non-matchs only implemented for simple patterns of size 2" << std::endl;
        } else {
            non_matchs_filename = args.at('2');
            bed_file non_matchs(non_matchs_filename, write);
            non_matchs.appendVector(results[1]);
            non_matchs.writeToFile();
        }
    } catch(std::out_of_range) {
        std::cout << "Skipping non-matchs" << std::endl;
    }
    
    return 0;
}