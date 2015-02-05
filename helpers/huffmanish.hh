#pragma once

//Huffman encodes a line and also produces the vocabulary ids
#include "hash.hh"
#include "line_splitter.hh"
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

//Sorting for the second
struct sort_pair {
	bool operator()(const std::pair<std::string, unsigned int> &left, const std::pair<std::string, unsigned int> &right) {
		return left.second > right.second; //This puts biggest numbers first.
	}
};

struct sort_pair_vec {
	bool operator()(const std::pair<std::vector<unsigned char>, unsigned int> &left, const std::pair<std::vector<unsigned char>, unsigned int> &right) {
		return left.second > right.second; //This puts biggest numbers first.
	}
};

class Huffman {
	unsigned long uniq_lines; //Unique lines in the file.

	//Containers used when counting the occurence of a given phrase
	std::map<std::string, unsigned int> target_phrase_words;
	std::map<std::vector<unsigned char>, unsigned int> word_all1;

	//Same containers as vectors, for sorting
	std::vector<std::pair<std::string, unsigned int> > target_phrase_words_counts;
	std::vector<std::pair<std::vector<unsigned char>, unsigned int> > word_all1_counts;

	//Huffman maps
	std::map<std::string, unsigned int> target_phrase_huffman;
	std::map<std::vector<unsigned char>, unsigned int> word_all1_huffman;

	//inverted maps
	std::map<unsigned int, std::string> lookup_target_phrase;
	std::map<unsigned int, std::vector<unsigned char> > lookup_word_all1;

	public:
		Huffman (const char *);
		void count_elements (const line_text &line);
		void assign_values();
		void serialize_maps(const char * dirname);
		void produce_lookups();

		std::vector<unsigned int> encode_line(const line_text &line);

		//encode line + variable byte ontop
		std::vector<unsigned char> full_encode_line(const line_text &line);

		//Getters
		const std::map<unsigned int, std::string> &get_target_lookup_map() const{
			return lookup_target_phrase;
		}
		const std::map<unsigned int, std::vector<unsigned char> > &get_word_all1_lookup_map() const{
			return lookup_word_all1;
		}

		unsigned long getUniqLines() {
			return uniq_lines;
		}
};

class HuffmanDecoder {
	std::map<unsigned int, std::string> lookup_target_phrase;
	std::map<unsigned int, std::vector<unsigned char> > lookup_word_all1;

public:
	HuffmanDecoder (const char *);
	HuffmanDecoder (std::map<unsigned int, std::string> *, std::map<unsigned int, std::vector<unsigned char> > *);

	//Getters
	const std::map<unsigned int, std::string> &get_target_lookup_map() const{
		return lookup_target_phrase;
	}
	const std::map<unsigned int, std::vector<unsigned char> > &get_word_all1_lookup_map() const{
		return lookup_word_all1;
	}

	inline std::string getTargetWordFromID(unsigned int id);

	std::string getTargetWordsFromIDs(std::vector<unsigned int> ids);

	target_text decode_line (std::vector<unsigned int> input, int num_scores);

	//Variable byte decodes a all target phrases contained here and then passes them to decode_line
	template<class Iterator>
	std::vector<target_text> full_decode_line (Iterator line_begin, Iterator line_end, int num_scores);
};

std::string getTargetWordsFromIDs(std::vector<unsigned int> ids, std::map<unsigned int, std::string> * lookup_target_phrase);

inline std::string getTargetWordFromID(unsigned int id, std::map<unsigned int, std::string> * lookup_target_phrase);

inline unsigned int reinterpret_float(float * num);

inline float reinterpret_uint(unsigned int * num);

std::vector<unsigned char> vbyte_encode_line(const std::vector<unsigned int> &line);
inline std::vector<unsigned char> vbyte_encode(unsigned int num);
inline unsigned int bytes_to_int(const std::vector<unsigned char> &number);

template<class Iterator>
std::vector<unsigned int> vbyte_decode_line(Iterator line_begin, Iterator line_end);

inline unsigned int bytes_to_int(const std::vector<unsigned char> &number){
	unsigned int retvalue = 0;
	std::vector<unsigned char>::const_iterator it = number.begin();
	unsigned char shift = 0; //By how many bits to shift

	while (it != number.end()) {
		retvalue |= (*it & 0x7f) << shift;
		shift += 7;
		it++;
	}

	return retvalue;
}

template<class Iterator>
std::vector<target_text> HuffmanDecoder::full_decode_line (Iterator lines_begin, Iterator lines_end, int num_scores){
	std::vector<target_text> retvector; //All target phrases
	std::vector<unsigned int> decoded_lines = vbyte_decode_line(lines_begin, lines_end); //All decoded lines
	std::vector<unsigned int>::iterator it = decoded_lines.begin(); //Iterator for them
	std::vector<unsigned int> current_target_phrase; //Current target phrase decoded

	short zero_count = 0; //Count home many zeroes we have met. so far. Every 3 zeroes mean a new target phrase.

	while(it != decoded_lines.end()){
		if (zero_count == 1) {
			//We are extracting scores. we know how many scores there are so we can push them
			//to the vector. This is done in case any of the scores is 0, because it would mess
			//up the state machine.
			for (int i = 0; i < num_scores; i++){
				current_target_phrase.push_back(*it);
				it++;
			}
		}

		if (zero_count == 3) {
			//We have finished with this entry, decode it, and add it to the retvector.
			retvector.push_back(decode_line(current_target_phrase, num_scores));
			current_target_phrase.clear(); //Clear the current target phrase and the zero_count
			zero_count = 0; //So that we can reuse them for the next target phrase
		}
		//Add to the next target_phrase, number by number.
		current_target_phrase.push_back(*it);
		if (*it == 0) {
			zero_count++;
		}
		it++; //Go to the next word/symbol
	}
	//Don't forget the last remaining line!
	if (zero_count == 3) {
		//We have finished with this entry, decode it, and add it to the retvector.
		retvector.push_back(decode_line(current_target_phrase, num_scores));
		current_target_phrase.clear(); //Clear the current target phrase and the zero_count
		zero_count = 0; //So that we can reuse them for the next target phrase
	}

	return retvector;

}

template<class Iterator>
std::vector<unsigned int> vbyte_decode_line(Iterator line_begin, Iterator line_end){
	std::vector<unsigned int> huffman_line;
	std::vector<unsigned char> current_num;

	for (Iterator it = line_begin; it != line_end; it++){
		current_num.push_back(*it);
		if ((*it >> 7) != 1) {
			//We don't have continuation in the next bit
			huffman_line.push_back(bytes_to_int(current_num));
			current_num.clear();
		}
	}
	return huffman_line;
}

