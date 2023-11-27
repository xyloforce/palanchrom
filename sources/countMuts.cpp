#include "vcf_tools.hpp"
#include "bed_tools.hpp"
#include "tools.hpp"
#include <iostream>

std::string reverseString(const std::string& sequence) {
    std::string sequence_final;
    for(const auto nuc: sequence) {
        switch(nuc) {
            case 'A':
                sequence_final = sequence_final + 'T';
                break;
            case 'C':
                sequence_final = sequence_final + 'G';
                break;
            case 'G':
                sequence_final = sequence_final + 'C';
                break;
            case 'T':
                sequence_final = sequence_final + 'A';
                break;
            default:
                sequence_final = sequence_final + nuc;
        }
    }
    return sequence_final;
}

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
        std::cout << "Optionnal :" << std::endl;
        std::cout << "\t+s check strand of ints" << std::endl;
        std::cout << "\t+i amount of lines to load at a time" << std::endl;
        throw std::out_of_range("Missing arguments");
    }

    bool stranded(false);
    try {
        args.at('s');
        stranded = true;
    } catch(std::out_of_range) {
        std::cout << "Ignoring strands" << std::endl;
    }

    int nb_blocks(100000);
    try {
        stranded = std::stoi(args.at('i'));
    } catch(std::out_of_range) {
        std::cout << "Number of blocks not specified, default 100000" << std::endl;
    }

    vcf_file mutations(vcf_filename, read);
    bed_file mask(bed_filename, read);
    mask.readWholeFile();
    AOE_file aoe(AOE_filename, read);
    aoe.readWholeFile();
    
    std::cout << "Intersect with mask..." << std::endl;
    aoe.apply_intersect(mask, stranded);
    std::cout << "Intersect and count mutations..." << std::endl;
    std::map <int, std::map <std::string, std::map <char, int>>> summed_values;
    while(mutations.remainToRead()) {
        mutations.eraseAndLoadBlock(nb_blocks);
        std::vector <intersect_results> results = mutations.intersect(aoe, false);
        for(int i(0); i < results.size(); i++) {
            std::string current_mutation = dynamic_cast<vcf_entry*> (results[i].source) -> getAlt() + "_" + dynamic_cast<vcf_entry*> (results[i].source) -> getRef();
            if(current_mutation[1] != 'N') {
                if(stranded && results[i].hit -> getStrand() == '-') {
                    current_mutation = reverseString(current_mutation);
                }
                summed_values[dynamic_cast <const AOE_entry*>(results[i].hit) -> getRelativePos(results[i].result.getStart())][current_mutation][results[i].hit -> getStrand()] ++;
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