headers= getSeq.h kmer.h
objects= getSeq.o kmer.o

G++_OPTION= -Wall

all : kcount

kcount : $(objects) kcount.o
        g++ -o kcount $(objects) kcount.o $(G++_OPTION)

kcount.o  :   $(headers) kcount.cpp
        g++ -c kcount.cpp -o kcount.o $(G++_OPTION)

getSeq.o    :   $(headers) getSeq.cpp
        g++ -c getSeq.cpp -o getSeq.o $(G++_OPTION)

kmer.o  :   $(headers)  kmer.cpp
        g++ -c kmer.cpp -o kmer.o $(G++_OPTION)

.PHONY  :   clean
clean   :
        rm $(objects) kcount.o
