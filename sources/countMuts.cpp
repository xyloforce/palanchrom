#include "vcf_tools.hpp"
#include "bed_tools.hpp"
#include "tools.hpp"
#include <iostream>

int main(int argc, char* argv[]) {
    std::map <char, std::string> args = getArgs(std::vector <std::string> (argv, argv + argc));
    std::string vcf_filename, AOE_filename, bed_filename, tsv_filename;
    std::cout << "Starting" << std::endl;
    try {
        vcf_filename = args.at('v');
        AOE_filename = args.at('a');
        bed_filename = args.at('b');
        tsv_filename = args.at('o');
    } catch (std::out_of_range) {
        std::cout << "Missing obligatory parameters. Parameters are :\n";
        std::cout << "\t+v vcf filename\n";
        std::cout << "\t+b bed filename\n";
        std::cout << "\t+o output filename\n";
        std::cout << "\t+a aoe filename\n";
        throw std::out_of_range("Missing arguments");
    }
    vcf_file mutations(vcf_filename, read);
    bed_file mask(bed_filename, read);
    mask.readWholeFile();
    AOE_file aoe(AOE_filename, read);
    aoe.readWholeFile();
    
    aoe.apply_intersect(mask);
    std::map <std::string, std::map <std::string, std::map <int, int>>> summed_values;
    while(mutations.remainToRead()) {
        mutations.eraseAndLoad(); // load exactly one entry
        std::vector <intersect_results> results = mutations.intersect(aoe, false);
        if(results.size() > 1) {
            std::cout << "error : \n";
            for(const auto& entry: results) {
                std::cout << entry.hit -> getString() << '\n';
            }
            throw std::logic_error("more than one overlap");
        } else if(results.size() == 1) {
            std::string current_mutation = dynamic_cast<vcf_entry*> (results[0].source) -> getRef() + dynamic_cast<vcf_entry*> (results[0].source) -> getAlt();
            if(current_mutation[1] != 'N') {
                summed_values[mutations.getEntry(0) -> getChr()][current_mutation][dynamic_cast <const AOE_entry*>(results[0].hit) -> getRelativePos(results[0].result.getStart())] ++;
            }
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