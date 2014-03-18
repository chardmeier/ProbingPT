#include "probing_hash_utils.hh"
#include "huffmanish.hh"
#include "hash.hh" //Includes line splitter
#include <sys/stat.h> //For finding size of file


char * read_binary_file(char * filename);

class QueryEngine {
	unsigned char * binary_mmaped; //The binari phrase table file
	std::map<unsigned int, std::string> vocabids;

	Table table;
	char *mem; //Memory for the table, necessary so that we can correctly destroy the object

	HuffmanDecoder decoder;

	size_t binary_filesize;
	size_t table_filesize;
	public:
		QueryEngine (const char *);
		~QueryEngine();
		std::pair<bool, std::vector<target_text> > query(StringPiece source_phrase);
		std::pair<bool, std::vector<target_text> > query(std::vector<uint64_t> source_phrase);
		void printTargetInfo(std::vector<target_text> target_phrases);
  	    const std::map<unsigned int, std::string> getVocab() const
  	    { return decoder.get_target_lookup_map(); }

};


