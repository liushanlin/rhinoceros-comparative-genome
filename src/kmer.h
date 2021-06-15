#pragma once

/*
quite some those codes get from the COPE software programed by Wei Fan, Binghang Liu, Jianying Yuan from BGI
You can obtain more info from its github linkage: 
*/

#include <iostream>
#include <cmath>

using namespace std;

//convert kmer-seq to kmer-bit
uint64_t seq2bit(string& kseq);

//convert kmer-bit to kmer-seq
string bit2seq(uint64_t kbit, int kmerSize);

// for sequence reverse complement
void reverse_complement(string& in_str, string& out_str);

// hash the key value
void Jenkin_hash(uint64_t&);
