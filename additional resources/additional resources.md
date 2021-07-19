**DNA sequence read of the historic/ancient samples**

The sequencing reads obtained for the four historic and ancient samples exhibited a length distribution as expected for ancient (a)DNA – concentrated around 50–60 bp (Poinar et al., 2006) (Figure 1), and all but the Siberian unicorn (~9✕) have genome-wide coverage > 10✕ (Figure 2). Furthermore, different libraries built for Javan rhinoceros generated reads of length on par with each other (Figure 3A), while the two DNA extracts and sequencing libraries of Merck’s rhinoceros produced different read length distribution patterns with the one sequenced using BGISeq platform containing more long reads (Figure 3B). However, we noted that the sequencing library for the BGISeq platform had an extremely high clonality level - reaching 97%, therefore, contributing little to the final genome coverage and we removed it from the following analyses. We further examined the consistency between the two sequence strategies (PE100 and PE50) for the woolly rhinoceros, and showed that both allow reliable SNP inference, and neither exhibit conflicts derived from systematic bias (Figure 4). Our analysis also showed that data with lower depth (library I) identified fewer variances as expected (Figure 4A), and ca. 1/3 of those incongruent variances that identified in library I were homozygous (Figure 4B) and those homozygous sites (incongruent between the two) could result from insufficient sequence depth (Figure 4C). Thus, we merged the two alignments for the following analyses.


![Figure 1] (https://github.com/liushanlin/rhinoceros-comparative-genome/blob/main/additional%20resources/readLengthDistributionAll.png)

**Figure 1 Read length distribution for the three extinct rhinoceros and the historic Javan rhinoceros sample.** Note: The current plot was estimated using one million randomly chosen reads from their raw data. Both the Merck’s rhinoceros and the Siberian unicorn data was generated using the Illumina sequencing platform, while the Javan and woolly rhinoceros data was sequenced using the BGISeq platform.

![Figure 2] ()

**Figure 2 Genome coverage for the historical and ancient samples.** Y-axes represent the proportion of the genome covered by any size fraction. Depth distribution was calculated using one million randomly-selected sites. Panel B shows a depth distribution of the Javan rhinoceros with the right plot representing a duplication removal strategy where all libraries were merged together before duplication removal. HQ includes only the reads with mapping quality ≥ 10.

![Figure 3] (https://github.com/liushanlin/rhinoceros-comparative-genome/blob/main/additional%20resources/lengthJavanStepPlain.png)

**Figure 3 Read length distribution of Javan (A) and Merck’s rhinoceros (B) sequences for which we built more than two libraries that were sequenced independently.** The middle panel shows a length distribution for which we binned the reads of length > 100 bp for the Illumina sequencing platform for a sake of fair comparison. We randomly chose two libraries out of the total 12 libraries built for Javan rhinoceros.

![Figure 4] (https://github.com/liushanlin/rhinoceros-comparative-genome/blob/main/additional%20resources/woollyDiffLib.png)

**Figure 4 Genomic variants identified using the two sequencing strategies (library I of PE100 and library II of PE50) for woolly rhinoceros.** A: Venn plot showing the proportion of variants obtained using the two different sequencing strategies; B: Proportion of the types (heterozygosity 0/1, labelled as heter, and ancient homozygosity 0/0 labelled as homoAncient) of the incongruent variances (17,217) generated in library I; C: Depth comparison for the different types of variances. nonOverlapHeter means the incongruent variances that were identified as heterozygous variances in library I, and nonOverlapHomo means those that were identifed as homozygous variances in library I.

**The *Stephanorhinus* and *Coelodonta* phylogenetic conundrum**

The phylogenetic relationship of *Stephanorhinus* and *Coelodonta* has been extensively studied over the last forty years (reviewed in Cappellini et al., 2019). Both morpho-anatomical and molecular data have been used to argue that *Stephanorhinus* represents a paraphyletic group that is closely related to Coelodonta. Mitochondrial genome-based phylogenetic studies (Kirillova et al., 2017; Margaryan et al., 2020) using a single representative per taxon have placed *Stephanorhinus kirchbergensis* as sister to *Coelodonta*, forming a clade with *Dicerorhinus sumatrensis* (Sumatran rhinoceros) as outgroup. Proteome sequencing data derived from an Early Pleistocene dental enamel from a *Stephanorhinus* tooth was used to argue that *Coelodonta* evolved from an early *Stephanorhinus* lineage, and that the *Stephanorhinus* genus is therefore paraphyletic, which was corroborated by Antoine et al ‘s study that utilised 278 morphological characters (Antoine et al., 2021). However, the close phylogenetic affinities between the *Stephanorhinus-Coelodonta* clade and *Dicerorhinus* recovered by all molecular analyses are at odds with morphology, as analyses of the latter consistently argue that *Stephanorhinus* and/or *Coelodonta* are constantly more closely related to the African rhinoceroses (e.g., Antoine et al., 2021). Thus, the cladogram differences between the morphological and molecular datasets could be attributed to either morphological homoplasy, or long-branch attraction, the latter of which could happen as only limited numbers of genes were analysed, and as many Miocene-Pliocene rhinocerotines cannot be included in molecular trees.


**Genome synteny and collinearity**

By conducting pairwise genome alignment between each de novo assembled rhinoceros genome and that of horse (Equus caballus), we did not detect any major within-scaffold segment reshuffling between the rhinoceros species, as well as for that between rhinoceroses and the horse. We then applied Circos (Krzywinski et al., 2009) to visualise the genome synteny and collinearity (Figure 5). The good within-scaffold collinearity between rhinoceros species to some extent guarantees the reliability and accuracy of the reference-based analyses for those historic and ancient samples.


![Figure 5](https://github.com/liushanlin/rhinoceros-comparative-genome/blob/main/additional%20resources/Equu_chr1.png)

**Figure 5. A demonstration example shows the pairwise genome alignments between the de novo-assembled genomes and horse using chromosome #1**. Species can be identified by the first letter of each scaffold ID: “w” represents white rhinoceros (red lines), “b” represents black rhinoceros (blue lines), “s” represents Sumatran rhinoceros (green lines) and “i” represents greater one-horned rhinoceros (orange lines). The figure only shows scaffolds of length > 5 Mb and alignment blocks of length > 100 kb.

> **Reference**
> 
> Krzywinski, M., Schein, J., Birol, I., Connors, J., Gascoyne, R., Horsman, D., Jones, S.J., and Marra, M.A. (2009). Circos: an information aesthetic for comparative genomics. Genome Res 19, 1639-1645.
