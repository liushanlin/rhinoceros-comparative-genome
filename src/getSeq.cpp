#include "getSeq.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include "kmer.h"

#ifndef BIT_MAX
	#define BIT_MAX 32
#endif

// remove the carriages from the fasta string and return ID and the whole seq using getline
int seqFromFa_stream(std::string& seqName, std::string& Seq, const std::string& str)
{
	std::string tmp;
	std::stringstream ss(str);
	std::getline(ss, seqName);
	while (std::getline(ss,tmp))
	{
		Seq += tmp;
	}
	if (0 == seqName.size()) {
		std::cerr << "no seqName obtained\n";
		return 0;
	}
	if (0 == Seq.size())
	{
		std::cerr << "no Seq obtained\n";
		return 0;
	}
	return 1;
}
// remove the carriages from the fasta string and return ID and the whole seq by for check
int seqFromFa(std::string& seqName, std::string& Seq, const std::string& str)
{
	int pos;
	int num = 0;
	Seq.reserve(str.size());
	for (size_t i = 0;i < str.size(); i++) {
		if ('\n' == str.at(i)) {
			num++;
			if (1 == num)
				pos = i;
		}
		else if (num >= 1)
			Seq.push_back(str.at(i));
	}
	Seq.shrink_to_fit();
	if (0 == num) {
		std::cerr << "do not find carriage return in the string file\n";
		return 0;
	}
	seqName = str.substr(0, pos+1);
	if (0 == Seq.size())
	{
		std::cerr << "no Seq obtained\n";
		return 0;
	}
	return 1;
}


int seqBit(const std::string & str, bitSequence& seq)
{
	int out = checkSeq(str);
	if (out == 0)
	{
		std::cerr << "the input sequence contain more than ATCG, please check for N or other character before\n";
		exit(0);
	}
	seq.seqLen = str.size();
	if (seq.seqLen > 0xffff)
	{
		std::cerr << "sequence length shoul be <=" << 0xffff << std::endl;
		exit(0);
	}
	for (size_t i = 0; i < seq.seqLen; i += BIT_MAX) {
		string sub = str.substr(i, BIT_MAX);
		uint64_t bit = seq2bit(sub);
		seq.seqBits.push_back(bit);
	}
	seq.seqBits.shrink_to_fit();
	return 1;
}


bitSequence seqBit(const std::string& str)
{
	bitSequence seq;
	int out = checkSeq(str);
	if (out == 0)
	{
		std::cerr << "the input sequence contain more than ATCG, please check for N or other character before\n";
		exit(0);
	}
	seq.seqLen = str.size();
	if (seq.seqLen > 0xffff)
	{
		std::cerr << "sequence length shoul be <=" << 0xffff << std::endl;
		exit(0);
	}
	for (size_t i = 0; i < seq.seqLen; i += BIT_MAX) {
		string sub = str.substr(i, BIT_MAX);
		uint64_t bit = seq2bit(sub);
		seq.seqBits.push_back(bit);
	}
	seq.seqBits.shrink_to_fit();
	return seq;
}


int seqBack(const bitSequence& bit, std::string& seq)
{
	for (size_t i = 0; i < bit.seqBits.size() - 1; i++) {
		seq += bit2seq(bit.seqBits[i], BIT_MAX);
	}
	int reminant = bit.seqLen % BIT_MAX;
	std::cout << "seqlenth is\t" << bit.seqLen << endl << "reminant is\t" << reminant << endl;
	if (0 == reminant)
	{
		seq += bit2seq(bit.seqBits.back(), BIT_MAX);
	}
	else
		seq += bit2seq(bit.seqBits.back(), reminant);
	return 1;
};

std::string seqBack(const bitSequence& bit)
{
	std::string seq;
	for (size_t i = 0; i < bit.seqBits.size() - 1; i++) {
		seq += bit2seq(bit.seqBits[i], BIT_MAX);
	}
	int reminant = bit.seqLen % BIT_MAX;
	std::cout << "seqlenth is\t" << bit.seqLen << endl << "reminant is\t" << reminant << endl;
	if (0 == reminant)
	{
		seq += bit2seq(bit.seqBits.back(), BIT_MAX);
	}
	else
		seq += bit2seq(bit.seqBits.back(), reminant);
	return seq;
};

//check whether a sequence contain non base characters, such as "N"
extern char alphabet[128];
int checkSeq(const string& seq)
{
	int check = 1;
	for (size_t i = 0; i < seq.size(); i++) {
		if (alphabet[(int)seq.at(i)] == 4) {
			check = 0;
			break;
		}
	}
	return check;
}

short countlines(FILE* const in)
{
	short num = 0;
	fseek(in, 0, SEEK_SET);
	char readline[1000];
	while (fgets(readline, 1000, in) != NULL)
	{
		if (readline[0] == '\n' || readline[0] == '\0')
			continue;
		num++;
	}
	return num;
}

int storeKmersFa(std::unordered_map<uint64_t,unsigned short*> &kset, const std::string& subfile, const short& ksize, const short& sampleN, const short& order)
{
	std::ifstream input;
	input.open(subfile);
	std::string line;
	while (std::getline(input, line, '>'))
	{
		if (line.size() == 0)
			continue;
		std::string id, seq;
		if (!seqFromFa(id, seq, line))
			return 0;
		size_t seqKmerEnd = seq.size() - ksize;
		for (size_t i = 0; i <= seqKmerEnd; i++)
		{
			std::string kseq = seq.substr(i, ksize);
			uint64_t kbit = seq2bit(kseq);
			unordered_map<uint64_t, unsigned short*>::iterator iter = kset.find(kbit);
			if (iter == kset.end())
			{
				unsigned short* arr = new unsigned short(sampleN);
				for (int i = 0; i < sampleN; i++) {
					arr[i] = 0;
				}
				arr[order]=1;
				kset[kbit] = arr;
			}
			else
			{
				if (iter->second[order] < 0xffff)
				iter->second[order]++;
			}

		}
	}
	input.close();
	/*
	FILE* input;
	if ((input = fopen(subfile.c_str(), "r")) == NULL)
	{
		perror("fopen");
		return 0;
	}
	fclose(input);
*/
	return 1;
}

int storeKmersFq(std::unordered_map<uint64_t, unsigned short*>& kset, const std::string& subfile, const short& ksize, const short& sampleN, const short& order)
{
	std::ifstream input;
	input.open(subfile);
	std::string line;
	while (std::getline(input, line))
	{
		std::string seq, ori, qua;
		std::getline(input, seq);
		std::getline(input, ori);
		std::getline(input, qua);
		if (line.size() == 0 || seq.size() == 0 || qua.size() != seq.size() || line[0] != '@')
			return 0;
		size_t seqKmerEnd = seq.size() - ksize;
		for (size_t i = 0; i <= seqKmerEnd; i++)
		{
			std::string kseq = seq.substr(i, ksize);
			uint64_t kbit = seq2bit(kseq);
			unordered_map<uint64_t, unsigned short*>::iterator iter = kset.find(kbit);
			if (iter == kset.end())
			{
				unsigned short* arr = new unsigned short(sampleN);
				for (int i = 0; i < sampleN; i++) {
					arr[i] = 0;
				}
				arr[order] = 1;
				kset[kbit] = arr;
			}
			else
			{
				if (iter->second[order] < 0xffff)
					iter->second[order]++;
			}

		}
	}
	input.close();
	return 1;
}

