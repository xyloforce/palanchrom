#include "fasta_tools.hpp"
#include "bed_tools.hpp"
#include "tools.hpp"
#include <iostream>
#include <thread>
#include <mutex>
#include <array>
#include <utility>
#include <list>

std::mutex mutex_map;

void count_entry(fasta_file const& genome, bio_entry* const& entry, bool keep_ids, bool count_gc, std::map<std::string, std::map<int, std::map<char, std::map<char, int>>>>& summed_values)
{
    std::string seq = dynamic_cast <fasta_entry*> (genome.getEntriesByChr(entry->getChr())[0])->subset(*entry);
    std::string id = "none";
    for (int i(entry->getStart()); i < entry->getEnd(); i++) {
        int pos_seq(i - entry->getStart());
        int rel_pos(dynamic_cast <AOE_entry*>(entry)->getRelativePos(i));
        if (entry->getStrand() == '-') {
            pos_seq = (seq.size() - 1) - pos_seq; // revert bc negative seq are counted backwards
        }
        if (keep_ids) {
            id = entry->getID();
            // rel_pos = 0;
        }
        if (count_gc) {
            if (seq.at(pos_seq) == 'G' || seq.at(pos_seq) == 'C') {
				mutex_map.lock();
                summed_values[id][rel_pos]['S'][entry->getStrand()]++;
				mutex_map.unlock();
            }
            else {
                mutex_map.lock();
                summed_values[id][rel_pos]['W'][entry->getStrand()]++;
                mutex_map.unlock();
            }
        }
        else {
            mutex_map.lock();
            summed_values[id][rel_pos][seq.at(pos_seq)][entry->getStrand()]++;
            mutex_map.unlock();
        }
    }
}

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
        std::cout << "Optionnal :\n";
        std::cout << "\t+b bed filename for masking ints\n";
        std::cout << "\t+s check strand of ints" << std::endl;
        std::cout << "\t+i count id by id instead of merging everything" << std::endl;
        std::cout << "\t+m id : keep both, source or hit" << std::endl;
		std::cout << "\t+g count GC instead of bases" << std::endl;
        std::cout << "\t+t number of threads" << std::endl;
        throw std::out_of_range("Missing arguments");
    }
    fasta_file genome(fasta_filename, read, standard);
    genome.readWholeFile();

    bool nomask = false;
    try {
        bed_filename = args.at('b');
    } catch(std::out_of_range) {
        nomask = true;
    }

    bool count_gc = false;
    try  {
        args.at('g');
        count_gc = true;
    } catch(std::out_of_range) {
        std::cout << "Counting bases" << std::endl;
	}

    unsigned long long count_threads = 1;
    try {
        count_threads = std::stol(args.at('t')) * 2; // you can fit a lot more threads bc they are pretty fast
    } catch (std::out_of_range) {
        std::cout << "using one thread" << std::endl;
    }
    
    AOE_file aoe(AOE_filename, read);
    aoe.readWholeFile();
    
    id_status status = both;
    bool keep_ids(false);
    try {
        args.at('i');
        keep_ids = true;
        try {
            std::string statusS(args.at('m'));
            if(statusS == "source") {
                status = source;
            } else if(statusS == "hit") {
                status = hit;
            } else {
                status = both;
            }
        } catch (std::out_of_range) {
            status = both;
        }
    } catch(std::out_of_range) {
        std::cout << "Ignoring ids" << std::endl;
    }

    std::map <std::string, std::map <int, std::map <char, std::map <char, int>>>> summed_values;
    if(!nomask) {
        std::cout << "Intersecting" << std::endl;
        bed_file mask(bed_filename, read);
        mask.readWholeFile();
        aoe.apply_intersect(mask, false, status);
    }

    std::cout << "Counting..." << std::endl;
    
	std::list <std::thread> all_threads;
    int running_threads(0);
	for (const auto& entry : aoe.getEntries()) { // generate threads for each entry
        if (running_threads > count_threads) {
			// FIXME last job is trying to remove a non-existing thread
            // wait for the thread that you will replace to finish
            all_threads.front().join();
			all_threads.pop_front(); // remove the thread that finished
        }
		all_threads.emplace_back(std::thread(count_entry, std::ref(genome), entry, keep_ids, count_gc, std::ref(summed_values))); // replace the thread with a new one
		running_threads ++;
    }

    for(auto& thread: all_threads) { // wait for the remaining threads to finish
        if(thread.joinable()) {
            thread.join();
        }
	}

    std::cout << "Writing results" << std::endl;
    std::ofstream output_file(tsv_filename);
    for(const auto& idToChr: summed_values) {
        for(const auto& chrToChar: idToChr.second) {
            for(const auto& charToPos: chrToChar.second) {
                for(const auto& posToCount: charToPos.second) {
                    std::string tmp = "";
                    if(keep_ids) {
                        tmp += idToChr.first + "\t";
                    }
                    tmp += std::to_string(chrToChar.first) + "\t" + charToPos.first  + "\t" + posToCount.first  + "\t" + std::to_string(posToCount.second) + "\n";
                    output_file << tmp;
                }
            }
        }
    }
    
    return 0;
}
