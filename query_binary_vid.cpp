#include "util/file_piece.hh"

#include "util/file.hh"
#include "util/scoped.hh"
#include "util/string_piece.hh"
#include "util/tokenize_piece.hh"
#include "util/murmur_hash.hh"
#include "util/probing_hash_table.hh"
#include "util/usage.hh"

#include "helpers/vocabid.hh"
#include "helpers/quering.hh"
#include "helpers/hash.hh"

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h> //For finding size of file
#include <boost/functional/hash.hpp>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <ctime> //for timing.
#include <chrono>

int main(int argc, char* argv[]) {
	if (argc != 6) {
		// Tell the user how to run the program
		std::cerr << "Usage: " << argv[0] << " path_to_hashtable path_to_data_bin path_to_vocabid tablesize max_entry_size" << std::endl;
		return 1;
	}

	QueryEngine queries(argv[1], argv[2], argv[3], argv[4], argv[5]);

	//Interactive search
	std::cout << "Please enter a string to be searched, or exit to exit." << std::endl;
	while (true){
		bool found;
		std::string cinstr = "";
		getline(std::cin, cinstr);
		if (cinstr == "exit"){
			break;
		}else{
			//Time lookup
			std::clock_t c_start = std::clock();
			auto t_start = std::chrono::high_resolution_clock::now();

			std::pair<bool, std::vector<target_text>> query_result;

			//Actual lookup
			query_result = queries.query(StringPiece(cinstr));

			if (query_result.first) {
				std::cout << "Key found! TODO. Print result" << std::endl;
			} else {
				std::cout << "Key not found!" << std::endl;
			}

			//End timing
			std::clock_t c_end = std::clock();
			auto t_end = std::chrono::high_resolution_clock::now();

			//Print timing results
			std::cout << "CPU time used: "<< 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC<< " ms\n";
			std::cout << "Real time passed: "<< std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count()<< " ms\n";
		}
	}
	//clean up
	std::cout << "CLEANUP NOT WORKING, FIX!" << std::endl;
	//delete queries;

	util::PrintUsage(std::cout);

	return 0;
}