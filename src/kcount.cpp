#include <cstdio>
#include <stdlib.h>
#include <iostream>
#include <unistd.h>
#include <string>
#include <unordered_map>
//#include <bits/stdc++.h>
#include "getSeq.h"
#include "kmer.h"

#define File_Max_length 1000

using namespace std;

void usage(void)
{
	string kmer_usage = "Usage: kcluter\n\
-i	input file list, one file per line\n\
-n	sample number \n\
-k	kmer size | 17 by default\n\
-o	basename of outputs | kclusterTest by default\n\
-t	input file type | fa or fq \n\
\n\
contact shanlin Liu for bug correction: shanlin1115@gmail.com\n\
";
	cerr << kmer_usage << endl;
	exit(-1);
}

int main(int argc, char* argv[])
{
	if (argc < 2)
	{
		usage();
	}
	int opt;
	unsigned short sample_num;
	string infile, basename, filetype;
	short ksize = 17;
	while ((opt = getopt(argc, argv, "i:o:n:k:t:")) != -1)
	{
		switch (opt)
		{
		case 'i':
			infile = optarg;
			break;
		case 'o':
			basename = optarg;
			break;
		case 'n':
			if (!isdigit(optarg[0])) {
				cerr << "-n option should be int\n";
				exit(-1);
			}
			sample_num = atoi(optarg);
			break;
		case 'k':
			ksize = atoi(optarg);
			break;
		case 't':
			filetype = optarg;
			break;
		}
	}
	if (0 == infile.length()) {
		usage();
	}
	if (0 == basename.length()) {
		basename = "kcluterTest";
	}
	if (sample_num < 2) {
		cerr << "try jellyfish for kmer count\n";
		exit(-1);
	}
	if (ksize > 31)
	{
		cerr << "kmer size should be lower than 31\n";
		exit(-1);
	}
	if (filetype.size() == 0 || (filetype != "fa" && filetype != "fq"))
	{
		cerr << "the input filetype is " << filetype << " must be set as fa or fq instead\n";
		exit(-1);
	}

	unordered_map<uint64_t, unsigned short*> kset;

	FILE* input;
	FILE* output;
	   	  
	if ((input = fopen(infile.c_str(), "r")) == NULL)
	{
		perror("fopen");
		exit(-1);
	}
	short lineNum = countlines(input);
	if (lineNum != sample_num)
	{
		cerr << "there are a total of " << lineNum << " files in the " << infile << " not equal to the sample number of " << sample_num << endl;
		exit(-1);
	}
	fseek(input, 0, SEEK_SET);
	short sampleOrder = 0;
	vector<string> fileList;
	char filelines[File_Max_length];
	while (fgets(filelines, File_Max_length,input) != NULL)
	{
		if (filelines[0] == '\n' || filelines[0] == '\0')
			continue;
		string individual_file(filelines);
		if (individual_file.back() == '\n')
			individual_file.pop_back();
		int check;
		if (filetype == "fa")
			check = storeKmersFa(kset, individual_file, ksize, sample_num, sampleOrder);
		else if (filetype == "fq")
			check = storeKmersFq(kset, individual_file, ksize, sample_num, sampleOrder);
		if (check == 1) {
			cout << "read " << individual_file << " in kmer set\t" << "correspond to array id of " << sampleOrder << endl;
			sampleOrder++;
			fileList.push_back(individual_file);
		}
		else 
		{	
			unordered_map<uint64_t, unsigned short*>::iterator iter = kset.begin();
			while (iter != kset.end())
			{
				delete iter->second;
				iter->second = NULL;
				kset.erase(iter++);
			}
			kset.clear();
			return 0;
		}
	}
	fclose(input);
	fileList.shrink_to_fit();
	
	string outfile = basename + ".count";
	if ((output = fopen(outfile.c_str(), "w")) == NULL)
	{
		perror("fopen");
		exit(-1);
	}
	fputs("kmerID\t", output);
	for (size_t i = 0; i < fileList.size(); i++) {
		fputs(fileList[i].c_str(), output);
		fputc('\t', output);
	}
	fputc('\n', output);

	for (auto iter = kset.begin(); iter != kset.end();)
	{
		string kmer = bit2seq(iter->first, ksize);
		fputs(kmer.c_str(), output);
		fputc('\t', output);
		for (size_t i = 0; i < sample_num;i++)
		{
			fprintf(output, "%d\t", iter->second[i]);
		}
		fputc('\n', output);
		delete iter->second;
		iter->second = NULL;
		kset.erase(iter++);
	}
	fclose(output);
	return 1;
}
