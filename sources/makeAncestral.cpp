#include "vcf_tools.hpp"
#include "bed_tools.hpp"
#include "fasta_tools.cpp"
#include "tools.hpp"

int main(int argc, char* argv[]) {
    std::map <char, std::string> args = getArgs(std::vector <std::string> (argv, argv + argc));
    std::string vcf_filename, bed_filename, fasta_input_filename, fasta_output_filename;
    std::cout << "Starting" << std::endl;
    try {
        vcf_filename = args.at('v');
        bed_filename = args.at('b');
        fasta_input_filename = args.at('i');
        fasta_output_filename = args.at('o');
    } catch (std::out_of_range) {
        std::cout << "Missing obligatory parameters. Parameters are :\n";
        std::cout << "\t+v vcf filename\n";
        std::cout << "\t+b bed filename\n";
        std::cout << "\t+i input fasta\n";
        std::cout << "\t+o output fasta\n";
        throw std::out_of_range("Missing arguments");
    }

    std::cout << "Loading files" << std::endl;

    vcf_file mutations(vcf_filename, read);
    bed_file intervals(bed_filename, read);
    fasta_file input_fasta(fasta_input_filename, read, standard);
    // , output_fasta(fasta_output_filename, write, standard)

    intervals.readWholeFile(); // load the entries
    input_fasta.readWholeFile();

    std::cout << "Masking the entries" << std::endl;
    input_fasta.typeToWrite(fasta_output_filename);
    input_fasta.apply_mask(intervals);

    std::cout << "Applying mutations" << std::endl;
    bool ok(true);
    vcf_entry* vcf_converted(0);
    std::unique_ptr <bio_entry> entry;
    while(mutations.remainToRead()) {
        ok = true;
        try {
            entry = std::move(mutations.readLine());
            vcf_converted = dynamic_cast <vcf_entry*>(entry.get());
            if(!vcf_converted) {
                throw std::bad_cast();
            }
        } catch(std::out_of_range) {
            ok = false;
            std::cout << "skipped a empty entry, probably whiteline at the eof" << std::endl;
        } 
        if(ok) {
            fasta_entry* fasta_converted = dynamic_cast <fasta_entry*> (input_fasta.getEntriesByChr(vcf_converted -> getChr())[0]);
            if(fasta_converted -> getBase(vcf_converted -> getStart()) == vcf_converted -> getRef()[0]) { // if the fasta entry matching the vcf entry chr has the same base as the vcf entry reference field
                fasta_converted -> mutate(vcf_converted -> getAlt()[0], vcf_converted -> getStart());
            } else {
                std::cout << "Ref base is " << fasta_converted -> getBase(vcf_converted -> getStart()) << "at pos " << vcf_converted -> getStart() << " in fasta and " << vcf_converted -> getRef()[0] << " in vcf" << std::endl;
                std::cout << vcf_converted -> getString() << std::endl;
                throw std::invalid_argument("Reference base is invalid");
            }
        }
    }
    std::cout << "Writing results" << std::endl;
    input_fasta.writeToFile();
}