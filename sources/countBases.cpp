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
        tsv_filename = args.at('o');
    } catch (std::out_of_range) {
        std::cout << "Missing obligatory parameters. Parameters are :\n";
        std::cout << "\t+f fasta filename\n";
        std::cout << "\t+o output filename\n";
        std::cout << "\t+a aoe filename\n";
        std::cout << "Optionnal :" << std::endl;
        std::cout << "\t+b bed filename for masking ints\n";
        std::cout << "\t+s check strand of ints" << std::endl;
        std::cout << "\t+i count id by id instead of merging everything" << std::endl;
        throw std::out_of_range("Missing arguments");
    }
    fasta_file ancestral(fasta_filename, read, standard);
    ancestral.readWholeFile();

    bool nomask = false;
    try {
        bed_filename = args.at('b');
    } catch(std::out_of_range) {
        nomask = true;
    }
    
    AOE_file aoe(AOE_filename, read);
    aoe.readWholeFile();

    bool stranded(false);
    try {
        args.at('s');
        stranded = true;
    } catch(std::out_of_range) {
        std::cout << "Ignoring strands" << std::endl;
    }

    bool keep_ids(false);
    try {
        args.at('i');
        keep_ids = true;
    } catch(std::out_of_range) {
        std::cout << "Ignoring ids" << std::endl;
    }

    std::map <std::string, std::map <int, std::map <char, std::map <char, int>>>> summed_values;
    if(!nomask) {
        std::cout << "Intersecting" << std::endl;
        bed_file mask(bed_filename, read);
        mask.readWholeFile();
        for(const auto& entry: aoe.intersect(mask, stranded)) {
            std::string seq = dynamic_cast <fasta_entry*> (ancestral.getEntriesByChr(entry.result.getChr())[0]) -> subset(entry.result);
            std::string id = "none";
            for(int i(entry.result.getStart()); i < entry.result.getEnd(); i++) {
                int pos_seq(i - entry.result.getStart());
                int rel_pos(dynamic_cast <AOE_entry*>(entry.source) -> getRelativePos(i));
                if(entry.source -> getStrand() == '-') {
                    pos_seq = (seq.size() - 1) - pos_seq; // revert bc negative seq are counted backwards
                }
                if(keep_ids) {
                    id = entry.source -> getID();
                    rel_pos = 0;
                }
                summed_values[id][rel_pos][seq.at(pos_seq)][entry.source -> getStrand()] ++;
            }
        }
    } else {
        std::cout << "Counting..." << std::endl;
        for(const auto& entry: aoe.getEntries()) {
            std::string seq = dynamic_cast <fasta_entry*> (ancestral.getEntriesByChr(entry -> getChr())[0]) -> subset(*entry);
            std::string id = "none";
            for(int i(entry -> getStart()); i < entry -> getEnd(); i++) {
                int pos_seq(i - entry -> getStart());
                int rel_pos(dynamic_cast <AOE_entry*>(entry) -> getRelativePos(i));
                if(entry -> getStrand() == '-') {
                    pos_seq = (seq.size() - 1) - pos_seq; // revert bc negative seq are counted backwards
                }
                if(keep_ids) {
                    id = entry -> getID();
                    rel_pos = 0;
                }
                summed_values[id][rel_pos][seq.at(pos_seq)][entry -> getStrand()] ++;
            }
        }
    }
    
    std::cout << "Writing results" << std::endl;
    // aoe.apply_intersect(mask);
    // aoe.typeToWrite("dump2.aoe");
    // aoe.writeToFile();
    std::ofstream output_file(tsv_filename);
    for(const auto& idToChr: summed_values) {
        for(const auto& chrToChar: idToChr.second) {
            for(const auto& charToPos: chrToChar.second) {
                for(const auto& posToCount: charToPos.second) {
                    std::string tmp = "";
                    if(keep_ids) {
                        tmp += idToChr.first;
                    }
                    tmp += std::to_string(chrToChar.first) + "\t" + charToPos.first  + "\t" + posToCount.first  + "\t" + std::to_string(posToCount.second) + "\n";
                    output_file << tmp;
                }
            }
        }
    }
    
    return 0;
}