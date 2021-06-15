#include "kmer.h"


// three arrays that convert bases to int and vus verse
char alphabet[128] = {
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 0, 4,
  4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 0, 4,
  4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
};

char bases[5] = {
  'A', 'C', 'G', 'T', 'N'
};

char c_bases[5] = {
  'T', 'G', 'C', 'A', 'N'
};

uint64_t seq2bit(string& kseq)
{
	uint64_t kbit = 0;

	for (size_t i = 0; i < kseq.size(); i++) {
		kbit = (kbit << 2) | alphabet[(int)kseq.at(i)];
	}

	return kbit;
}

//convert kmer-bit to kmer-seq
string bit2seq(uint64_t kbit, int kmerSize)
{
	string kseq;

	for (int i = 0; i < kmerSize; i++) {
		kseq.push_back(bases[(kbit >> (kmerSize - 1 - i) * 2) & 0x3]);
	}

	return kseq;
}

void reverse_complement(string& in_str, string& out_str)
{
	for (int i = in_str.size() - 1; i >= 0; i--) {
		out_str.push_back(c_bases[(int)alphabet[(int)in_str.at(i)]]);
	}
}

// hash the kbit
void Jenkin_hash(uint64_t& key)
{
	key += ~(key << 32);
	key ^= (key >> 22);
	key += ~(key << 13);
	key ^= (key >> 8);
	key += (key << 3);
	key ^= (key >> 15);
	key += ~(key << 27);
	key ^= (key >> 31);
//	return key;
}
