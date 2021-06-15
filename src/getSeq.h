#pragma once
#include <string>
#include <vector>
#include <unordered_map>

/*
#ifndef K_COUNT
typedef struct {
	int arr[];
}kcount;
#endif // !K_COUNT
*/


// remove the carriages from the fasta string and return ID and the whole seq
int seqFromFa_stream(std::string&, std::string&, const std::string&); 

int seqFromFa(std::string&, std::string &, const std::string&);


// get the entier sequene transfered to bit to save mem
typedef struct
{
	std::vector<uint64_t> seqBits;
	unsigned short seqLen;
} bitSequence;

// from a string to the bit structure
int seqBit(const std::string&, bitSequence&);
bitSequence seqBit(const std::string&);

// from the bit structure to sequence
int seqBack(const bitSequence&, std::string&);
std::string seqBack(const bitSequence&);

//check whether a sequence contain non base characters, such as "N"
int checkSeq(const std::string&);

//check line number of a file
short countlines(FILE* const);

//read sequences from fasta file and store in the kmerset
int storeKmersFa(std::unordered_map<uint64_t,unsigned short*>&, const std::string&, const short&, const short&, const short&);

//read sequences form fastq file and store in the kmerset
int storeKmersFq(std::unordered_map<uint64_t, unsigned short*>&, const std::string&, const short&, const short&, const short&);
