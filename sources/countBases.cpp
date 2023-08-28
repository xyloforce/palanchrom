#include "fasta_tools.hpp"
#include "bed_tools.hpp"
#include "tools.hpp"
#include <iostream>

int main(int argc, char* argv[]) {
    std::map <char, std::string> args = getArgs(std::vector <std::string> (argv, argv + argc));
    std::string fasta_filename, AOE_filename, bed_filename, tsv_filename;
    std::cout << "Starting" << std::endl;
    try {
        fasta_filename = args.at('f');
        AOE_filename = args.at('a');
        bed_filename = args.at('b');
        tsv_filename = args.at('o');
    } catch (std::out_of_range) {
        std::cout << "Missing obligatory parameters. Parameters are :\n";
        std::cout << "\t+f fasta filename\n";
        std::cout << "\t+b bed filename\n";
        std::cout << "\t+o output filename\n";
        std::cout << "\t+a aoe filename\n";
        throw std::out_of_range("Missing arguments");
    }
    fasta_file ancestral(fasta_filename, read, standard);
    ancestral.readWholeFile();
    bed_file mask(bed_filename, read);
    mask.readWholeFile();
    AOE_file aoe(AOE_filename, read);
    aoe.readWholeFile();

    std::cout << "Intersecting" << std::endl;
    std::map <std::string, std::map <char, std::map <int, int>>> summed_values;
    for(const auto& entry: aoe.intersect(mask)) {
        std::string seq = dynamic_cast <fasta_entry*> (ancestral.getEntriesByChr(entry.result.getChr())[0]) -> subset(entry.result);
        for(int i(entry.result.getStart()); i < entry.result.getEnd(); i++) {
            summed_values[entry.result.getChr()][seq[i - entry.result.getStart()]][dynamic_cast <AOE_entry*>(entry.source) -> getRelativePos(i)] ++;
        }
    }
    std::cout << "Writing results" << std::endl;
    std::ofstream output_file(tsv_filename);
    for(const auto& chrToChar: summed_values) {
        for(const auto& charToPos: chrToChar.second) {
            for(const auto& posToCount: charToPos.second) {
                output_file << chrToChar.first << "\t" << charToPos.first  << "\t" << posToCount.first  << "\t" << posToCount.second << "\n";
            }
        }
    }
    return 0;
}