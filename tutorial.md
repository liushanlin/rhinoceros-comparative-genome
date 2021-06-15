## genome alignment

```
lastal -P10 -m100 -E0.05 horseGenomeDB sumatranPsmcRef.fasta | last-split -m1 > sumatran.maf
faToTwoBit sumatranPsmcRef.fasta sumatranPsmcRef.2bit
faSize -detailed sumatranPsmcRef.fasta > sumatranPsmcRef.size
sh genomeAlignSynteny.sh sumatran.maf Equu.size sumatranPsmcRef.size Equu.2bit sumatranPsmcRef.2bit Equu Sumatran
```


## consensus call for ancient samples

    angsd -dofasta 2 -explode 1 -basesPerLine 100 -i filtered.bam -out javanAngsd -doCounts 1 -trim 5 -minQ 20 -minMapQ 20 -setMinDepth 5
> Elasmotherim minDepth was set as 3 and minMapQ set as 10 due to lower coverage

## genome regions extraction via sliding windows of 100 kb from the genome alignment in maf format

    perl mafRegionExtract_v1.pl filtered.maf list.spInMaf 100 > blocks
    perl blocks2msa_v1.pl sp.table equu.chr.info blocks

## missing ratio check for each silding window

    perl missingStatistic.pl list.fasta sp.table > missingStatistic.txt

## filter out short and high missing alignments    

    perl filterMissing.pl missingStatistic.txt |awk '{print $1}' >list.qualified
    
## guided alignment 
    
    perl /groups/hologenomics/shanlin/rhinos/msa/guideAln_v1.pl chr10_28300000.fasta /groups/hologenomics/shanlin/rhinos/msa/sp.table
> show as an example

## check for effective length ratio without N

    ls chr*/chr*.aln >list.aln
    perl checkEffectiveRatio.pl list.aln >checkEffectiveRatio.txt
    awk '$3>=0.2{print $1}' checkEffectiveRatio.txt >list.ge2

## run raxml to infer gene trees and then obtain the species tree using Astral


# obtan orthologs from the genome alignment

## splite mafs according the gene location 

    mafsInRegion -outDir bedFile genesMaf tapirsFiltered.maf
    
## get fasta for all the targeted species and query species saperately

    perl maf2msa.pl sp.table list.maf

## exonerate 

    exonerate -m protein2genome -n 6 ../faaFromMaf/gene10000.faa ../msa/gene10000.target >gene10000.log 2>gene10000.err

## obtain cds and faa for each species asking exons that can be pairwise correlated

    perl exonerate2coordinateRhino_v1.pl ../faaFromMaf/gene10002.faa ../sp.table ../exonerate/gene10002.log ../msa/gene10002.target ../msa/gene10002.query gene10002 2>>err
> this script will get fna and faa sequences for each species, plus, it will generate a log file containing the p-distance information for each exon comparing to its corresponding best hit reference exon. The log file will be useful when applying the exon check in the next step

    perl preAlnExonSelect.pl gene10000
> this script will check the log file generated in the last step and output fna and faa sequences that share the same exon

## alignment for faa and fna

    perl guideAln_pro.pl gene10000_ec.faa sp.table
    
    grep ">" -c *.aln |awk -F\: '$2!=14{print $1}' |xargs rm
> some may contain errors and won't generate sequences for all the species and has to be removed

    pal2nal.pl ~/testRhino/cell/tapirsMultiz/faaAln/gene10000_ec.aln ~/testRhino/cell/tapirsMultiz/cdsFaa/gene10000_ec.fna  -output fasta > gene10000.fasta
> use pal2nal to align the cds sequences

## merged sequences and mask the conflicts from different reference genome

    perl refAlnCheck.pl gene10000.fasta gene10000.aln 1
> 1 means strict filter and will mask differences (inlcuding missing ones) as Ns

