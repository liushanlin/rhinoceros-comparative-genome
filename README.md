# rhinoceros-comparative-genome

It contains all the important scripts used for a comparative genome study of the rhinceros species, including all the five extant rhinoceros species and three extinct rhinocerso species.

The scripts were utilized to extract genome alingment information from MAF files. For instance, it tries to add mapping-based consensus sequences (obtained from ANGSD -doFast function using a bam file) to the multiple sequence alignment (MSA) file that contains its references; it also tries to obtain CDS and amino acid sequences from MAF files according to the genome refence's annotation information (GFF file). 

Please check the tutorial for step by step guidiance.

## genomeAlignSynteny.sh
> 
It takes input of genome alignment generated by LAST (in format of MAF) and outputs alignment blocks in synteny. It takes advantage of programs provided by UCSC tools, such axtChain, netToAxt, et al.

## block2msa.pl
>
It reads whole genome alignments (file in MAF format generated using multiZ) and a species table which contians reference species and the query species. The query species are those that do not have de novo genome reference and have to obtain the genome sequences via concensus calling using functions like ANGSD -doFast

## exonerate2coordinate.pl
>
It takes inputs of log files generated from exonerate and the corresponding nucleodtie fasta files, and output the amino acid sequuences and cds identified by exonerate. In addition, it calculates the p-distance between identified exons and the reference exons and finds the best hits, which can be used for the following exon based sequence alignment.

## genelocationMaf2cds.pl
>
It reads gene annotation information (e.g., gene location and exon location information) and the genome alignment file (in format of MAF), and outputs amino acid sequences that located in the genome alignment region.
## exonerateFrameShiftPos.pl
>
It reads exonerate log file, the gene alignment file (maf file) and finds frameshift mutations shared between different species, and outputs their location information
