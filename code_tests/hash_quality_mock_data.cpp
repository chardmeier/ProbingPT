#include "util/file_piece.hh"

#include "util/file.hh"
#include "util/scoped.hh"
#include "util/string_piece.hh"
#include "util/tokenize_piece.hh"
#include "util/murmur_hash.hh"
#include "util/probing_hash_table.hh"
#include "util/usage.hh"

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <fstream>
#include <iostream>
#include <set> //For counting

//Struct for holding processed line
struct line_text {
	int value;
	StringPiece text;
};

//Ask if it's better to have it receive a pointer to a line_text struct
line_text splitLine(StringPiece textin);

uint64_t getHash(StringPiece text);

uint64_t getHash(StringPiece text) {
	std::size_t len = text.size();
	uint64_t key = util::MurmurHashNative(text.data(), len);
	return key;
}

line_text splitLine(StringPiece textin) {
	const char delim[] = " ||| ";
	StringPiece str; //Temp string container
	std::string temp; //Temp string for conversion
	int num;
	line_text output;

	//Tokenize
	util::TokenIter<util::MultiCharacter> it(textin, util::MultiCharacter(delim));
	//Get string
	str = *it;
	output.text = str;
	it++;

	//Get num
	str = *it;

	//Convert to int
	temp = str.as_string();
	num = atoi(temp.c_str());
	output.value = num;

	return output;
}

int main(int argc, char* argv[]){
	if (argc != 2) {
		// Tell the user how to run the program
		std::cerr << "Provided " << argc - 1 << " arguments, needed 1." << std::endl;
		std::cerr << "Usage: " << argv[0] << " path_to_phrasetable" << std::endl;
		return 1;
	}
	util::FilePiece filein(argv[1]);
	std::set<uint64_t> hashes;
	unsigned long setsize = hashes.size();
	unsigned int clashes = 0;
	unsigned long lines = 0;

	while(true){
		try{
			line_text newline = splitLine(filein.ReadLine());
			uint64_t hash = getHash(newline.text); //Get hash
			hashes.insert(hash); //Insert it in the set
			unsigned long newsize = hashes.size();
			if (newsize == setsize){
				clashes++; //We have a clash
			}else{
				setsize = newsize;
			}
			lines++;
			if (lines%10000 == 0){
				std::cout << "On line: " << lines << std::endl;
			}
		} catch (util::EndOfFileException e){
			std::cout << "End of file" << std::endl;
			break;
		}
	}

	std::cout << "Entries: " << lines << " Collisions: " << clashes << std::endl;
}