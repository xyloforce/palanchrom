#include "fasta_tools.hpp"
#include "vcf_tools.hpp"
#include "tools.hpp"
#include <iostream>

bool all(const std::vector <fasta_file>& files, int index, int pos) {
    char tchar(files[0].getEntry(index) -> getBase(pos));
    for(int i(1); i < files.size(); i++) {
        if(files[i].getEntry(index) -> getBase(pos) != tchar) {
            return false;
        }
    }
    return true;
}

int main(int argc, char* argv[]) {
    std::map <char, std::string> args = getArgs(std::vector<std::string>(argv, argv + argc));
    std::string fasta_base_filename, fasta_pair_filename, vcf_filename;
    std::vector <std::string> fasta_outgroups_filenames;
    std::cout << "Starting" << std::endl;
    try {
        fasta_base_filename = args.at('1');
        fasta_pair_filename = args.at('2');
        fasta_outgroups_filenames.push_back(args.at('3')); // need at last one outgroup
        vcf_filename = args.at('v');
    } catch (std::out_of_range) {
        std::cout << "Missing obligatory parameters. Parameters are :\n";
        std::cout << "\t+1 base fasta file\n";
        std::cout << "\t+2 pair fasta file\n";
        std::cout << "\t+3 outgroup fasta file\n";
        std::cout << "\t+v output file\n";
        std::cout << "Optionnal :\n";
        std::cout << "\t+4, ... , 10 : additionnal outgroups fasta files\n";
        throw std::out_of_range("Missing arguments");
    }
    for(int i(4); i < 10; i++) {
        try {
            fasta_outgroups_filenames.push_back(args.at(std::to_string(i)[0]));
        } catch(std::out_of_range) {
            break;
        }
    }
    std::cout << "Loading files" << std::endl;
    fasta_file base_fasta(fasta_base_filename, read, bedtools_stranded), pair_fasta(fasta_pair_filename, read, bedtools_stranded);
    base_fasta.readWholeFile();
    pair_fasta.readWholeFile();
    std::vector <fasta_file> outgroups;
    for(const std::string& filename: fasta_outgroups_filenames) {
        outgroups.push_back(fasta_file(filename, read, bedtools_stranded));
        outgroups[outgroups.size() -1].readWholeFile();
    }
    std::string chr, id("."), ancestral_base, filter(".");
    std::map <std::string, std::string> infos;
    std::vector <std::string> current_base;
    int qual(0);
    long start, start_entry(0);
    bool write_entry(false);
    std::cout << "Analyzing mutations" << std::endl;
    vcf_file output(vcf_filename, write);
    for(int i(0); i < base_fasta.getSize(); i++) {
        if(i % 1000 == 0) {
            std::cout << "treating entry " << i << " out of " << base_fasta.getSize() << "\r";
        }
        chr = base_fasta.getEntry(i) -> getChr();
        start_entry = base_fasta.getEntry(i) -> getStart();
        for(int j(0); j < base_fasta.getEntry(i) -> getSize(); j++) {
            infos.clear();
            start = start_entry + j;
            current_base = std::vector<std::string>(1, std::string(1, base_fasta.getEntry(i) -> getBase(j)));
            if(current_base[0][0] != 'N') { // no need to check for mutations on non-defined reference bases
                if(!all(outgroups, i, j)) { // difference in outgroups
                    ancestral_base = "N";
                    write_entry = true;
                } else if(pair_fasta.getEntry(i) -> getBase(j) != outgroups[0].getEntry(i) -> getBase(j)) {
                    ancestral_base = "N";
                    infos["reason"] = "mismatch";
                    write_entry = true;
                } else if(pair_fasta.getEntry(i) -> getBase(j) != current_base[0][0]) {
                    ancestral_base = std::string(1, pair_fasta.getEntry(i) -> getBase(j));
                    write_entry = true;
                }
                if(write_entry) {
                    write_entry = false;
                    infos["source"] = chr + "_" + std::to_string(start_entry) + "_" + std::to_string(base_fasta.getEntry(i) -> getEnd());
                    output.writeBioLine(vcf_entry(chr, start, start + 1, id, ancestral_base, current_base, qual, filter, infos));
                }
            }
        }
    }
    return(0);
}