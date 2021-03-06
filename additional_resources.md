# **DNA sequence read of the historic/ancient samples**

The sequencing reads obtained for the four historic and ancient samples exhibited a length distribution as expected for ancient (a)DNA – concentrated around 50–60 bp (Poinar et al., 2006) (Figure 1), and all but the Siberian unicorn (~9✕) have genome-wide coverage > 10✕ (Figure 2). Furthermore, different libraries built for Javan rhinoceros generated reads of length on par with each other (Figure 3A), while the two DNA extracts and sequencing libraries of Merck’s rhinoceros produced different read length distribution patterns with the one sequenced using BGISeq platform containing more long reads (Figure 3B). However, we noted that the sequencing library for the BGISeq platform had an extremely high clonality level - reaching 97%, therefore, contributing little to the final genome coverage and we removed it from the following analyses. We further examined the consistency between the two sequence strategies (PE100 and PE50) for the woolly rhinoceros, and showed that both allow reliable SNP inference, and neither exhibit conflicts derived from systematic bias (Figure 4). Our analysis also showed that data with lower depth (library I) identified fewer variances as expected (Figure 4A), and ca. 1/3 of those incongruent variances that identified in library I were homozygous (Figure 4B) and those homozygous sites (incongruent between the two) could result from insufficient sequence depth (Figure 4C). Thus, we merged the two alignments for the following analyses.

<br/><br/>
![Figure 1](https://github.com/liushanlin/rhinoceros-comparative-genome/blob/main/additional%20resources/readLengthDistributionAll.png)
**Figure 1 Read length distribution for the three extinct rhinoceros and the historic Javan rhinoceros sample.** Note: The current plot was estimated using one million randomly chosen reads from their raw data. Both the Merck’s rhinoceros and the Siberian unicorn data was generated using the Illumina sequencing platform, while the Javan and woolly rhinoceros data was sequenced using the BGISeq platform.

<br/><br/>
![Figure 2](https://github.com/liushanlin/rhinoceros-comparative-genome/blob/main/additional%20resources/depth.ancient.png)
**Figure 2 Genome coverage for the historical and ancient samples.** Y-axes represent the proportion of the genome covered by any size fraction. Depth distribution was calculated using one million randomly-selected sites. Panel B shows a depth distribution of the Javan rhinoceros with the right plot representing a duplication removal strategy where all libraries were merged together before duplication removal. HQ includes only the reads with mapping quality ≥ 10.

<br/><br/>
![Figure 3](https://github.com/liushanlin/rhinoceros-comparative-genome/blob/main/additional%20resources/lengthJavanStepPlain.png)
**Figure 3 Read length distribution of Javan (A) and Merck’s rhinoceros (B) sequences for which we built more than two libraries that were sequenced independently.** The middle panel shows a length distribution for which we binned the reads of length > 100 bp for the Illumina sequencing platform for a sake of fair comparison. We randomly chose two libraries out of the total 12 libraries built for Javan rhinoceros.

<br/><br/>
![Figure 4](https://github.com/liushanlin/rhinoceros-comparative-genome/blob/main/additional%20resources/woollyDiffLib.png)
**Figure 4 Genomic variants identified using the two sequencing strategies (library I of PE100 and library II of PE50) for woolly rhinoceros.** A: Venn plot showing the proportion of variants obtained using the two different sequencing strategies; B: Proportion of the types (heterozygosity 0/1, labelled as heter, and ancient homozygosity 0/0 labelled as homoAncient) of the incongruent variances (17,217) generated in library I; C: Depth comparison for the different types of variances. nonOverlapHeter means the incongruent variances that were identified as heterozygous variances in library I, and nonOverlapHomo means those that were identifed as homozygous variances in library I.

---



# **Genome synteny and collinearity**

By conducting pairwise genome alignment between each de novo assembled rhinoceros genome and that of horse (Equus caballus), we did not detect any major within-scaffold segment reshuffling between the rhinoceros species, as well as for that between rhinoceroses and the horse. We then applied Circos (Krzywinski et al., 2009) to visualise the genome synteny and collinearity (Figure 5). The good within-scaffold collinearity between rhinoceros species to some extent guarantees the reliability and accuracy of the reference-based analyses for those historic and ancient samples.

<br/><br/>
![Figure 5](https://github.com/liushanlin/rhinoceros-comparative-genome/blob/main/additional%20resources/Equu_chr1.png)<br/><br/>
**Figure 5. A demonstration example shows the pairwise genome alignments between the de novo-assembled genomes and horse using chromosome #1**. Species can be identified by the first letter of each scaffold ID: “w” represents white rhinoceros (red lines), “b” represents black rhinoceros (blue lines), “s” represents Sumatran rhinoceros (green lines) and “i” represents greater one-horned rhinoceros (orange lines). The figure only shows scaffolds of length > 5 Mb and alignment blocks of length > 100 kb.

---

# **Demographic estimation using different parameters**
As a comparison, and due to the uncertainty of μ and g for the three extinct rhinoceros species, we also estimated changes in Ne using shorter generation times (i.e., 10) and different substitution rates (i.e., 1.2 x 10-8 (Goto et al., 2011) and 3.3 x 10-8 (Orlando et al., 2013)) (Figure 6).

<br/><br/>
![Figure 6](https://github.com/liushanlin/rhinoceros-comparative-genome/blob/main/additional%20resources/difPara.png)

**Figure 6. Past demographic trajectory for the three extinct rhinoceros species.** Each curve represents a different species and each panel represents different substitution rates (substitutions/site/generation, μ) and generation times (g).


---

# **Genes with frameshift mutations that could contribute to rhinoceros’ specific phenotypes**

Although it remains challenging to functionally validate the genes that showed signatures of loss-of-function, we uncovered frameshift mutations in several genes that are functionally well studied and could play vital roles in shaping rhinoceros morphological features. These include IFT43 (Intraflagellar Transport 43, Figure S7.1), that is involved in the formation and maintenance of cilia, which are important for the development and function of the light-sensitive tissue at the back of the eye (the retina) (Arts et al., 2011). Studies show that mutations in IFT43 cause non-syndromic recessive retinal degeneration (Biswas et al., 2017). We detected loss-of-function of MOGAT 3 (Monoacylglycerol O-acyltransferase 3, Figure 7) and multiple FGFs (Fibroblast Growth Factors) genes for all the rhinoceroses, both of which are reported to be evolutionarily critical to the development of hair follicles (Imamura, 2014; Lopesmarques et al., 2019); KDSR (3-ketodihydrosphingosine reductase, Figure 8.1) functions in the homeostasis of keratinisation, and variances in which were identified to be able to cause severe skin disorder characterised by neonatal onset of thick, scaly skin on the face and genitals (Boyden et al., 2017). Genes like DIPK1C (Divergent Protein Kinase Domain 1C), SLC7A13 (Solute Carrier Family 7 Member 13) that were reported to be related to abnormality of retinal pigmentation (Cody et al., 2009; Yahyaoui and Perezfrias, 2019), and mutations in ADAMTS (A Disintegrin and Metalloproteinase with Thrombospondin motifs) genes that have diverse roles in tissue morphogenesis were reported to cause skin and eye related diseases (Kelwick et al., 2015). 

<br/><br/>
![Figure 7](https://github.com/liushanlin/rhinoceros-comparative-genome/blob/main/additional%20resources/nonsense%20mutation.png)<br/><br/>
**Figure 7 Example loss-of-function inducing frameshift mutations shared across all rhino species**

# **Reference**

> Antoine, P.-O., Reyes, M.C., Amano, N., Claude, J., Bautista, A.P., Vos, J.d., and Ingicco, T. (2021). A new clade of rhinoceroses from the Pleistocene of Asia sheds light on mainland mammal dispersals to the Philippines. Zoological Journal of the Linnean Society, Advance online. https://doi.org/10.1093/zoolinnean/zlab009 <br/><br/>
> Arts, H.H., Bongers, E.M.H.F., Mans, D.A., Beersum, S.E.C.v., Oud, M.M., Bolat, E., Spruijt, L., Cornelissen, E.A.M., Schuurs-Hoeijmakers, J.H.M., Leeuw, N.d., et al. (2011). C14ORF179 encoding IFT43 is mutated in Sensenbrenner syndrome. Journal of Medical Genetics 48, 390-395.<br/><br/>
> Biswas, P., Duncan, J.L., Ali, M., Matsui, H., Naeem, M.A., Raghavendra, P.B., Frazer, K.A., Arts, H.H., Riazuddin, S., and Akram, J. (2017). A mutation in IFT43 causes non-syndromic recessive retinal degeneration. Human Molecular Genetics 26, 4741-4751.<br/><br/>
> Boyden, L.M., Vincent, N.G., Zhou, J., Hu, R., Craiglow, B.G., Bayliss, S.J., Rosman, I.S., Lucky, A.W., Diaz, L.A., and Goldsmith, L.A. (2017). Mutations in KDSR cause recessive progressive symmetric erythrokeratoderma. American Journal of Human Genetics 100, 978-984.<br/><br/>
> Cody, J.D., Heard, P., Crandall, A.C., Carter, E., Li, J., Hardies, L.J., Lancaster, J.L., Perry, B., Stratton, R.F., and Sebold, C. (2009). Narrowing critical regions and determining penetrance for selected 18q- phenotypes. American Journal of Medical Genetics Part A 149, 1421-1430.<br/><br/>
> Goto, H., Ryder, O.A., Fisher, A.R., Schultz, B., Kosakovsky Pond, S.L., Nekrutenko, A., and Makova, K.D. (2011). A massively parallel sequencing approach uncovers ancient origins and high genetic variability of endangered Przewalski's horses. Genome Biology and Evolution 3, 1096-1106.<br/><br/>
> Imamura, T. (2014). Physiological functions and underlying mechanisms of Fibroblast Growth Factor (FGF) family members: recent findings and implications for their pharmacological application. Biological & Pharmaceutical Bulletin 37, 1081-1089.<br/><br/>
> Kelwick, R., Desanlis, I., Wheeler, G.N., and Edwards, D.R. (2015). The ADAMTS (A Disintegrin and Metalloproteinase with Thrombospondin motifs) family. Genome Biology 16, 113-113.<br/><br/>
> Kirillova, I.V., Chernova, O.F., Der Made, J.V., Kukarskih, V.V., Shapiro, B., Der Plicht, J.V., Shidlovskiy, F.K., Heintzmann, P.D., Van Kolfschoten, T., and Zanina, O.G. (2017). Discovery of the skull of Stephanorhinus kirchbergensis (Jäger, 1839) above the Arctic Circle. Quaternary Research 88, 537-550.<br/><br/>
> Krzywinski, M., Schein, J., Birol, I., Connors, J., Gascoyne, R., Horsman, D., Jones, S.J., and Marra, M.A. (2009). Circos: an information aesthetic for comparative genomics. Genome Res 19, 1639-1645.<br/><br/>
> Margaryan, A., Sinding, M.S., Liu, S., Vieira, F.G., Chan, Y.L., Nathan, S.K.S.S., Moodley, Y., Bruford, M.W., and Gilbert, M.T.P. (2020). Recent mitochondrial lineage extinction in the critically endangered Javan rhinoceros. Zoological Journal of the Linnean Society. 190, 372–383.<br/><br/>
> Orlando, L., Ginolhac, A., Zhang, G., Froese, D., Albrechtsen, A., Stiller, M., Schubert, M., Cappellini, E., Petersen, B., and Moltke, I. (2013). Recalibrating Equus evolution using the genome sequence of an early Middle Pleistocene horse. Nature 499, 74-78.<br/><br/>
> Poinar, H.N., Schwarz, C., Qi, J., Shapiro, B., MacPhee, R.D.E., Buigues, B., Tikhonov, A., Huson, D.H., Tomsho, L.P., Auch, A., et al. (2006). Metagenomics to paleogenomics: large-scale sequencing of mammoth DNA. Science 311, 392-394.
