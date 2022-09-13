---
hide:
  - navigation
  - toc
---

![Alt text](logo/logo.png?raw=true "logo")

# TB/NTM species identification pipeline Version 1.0 (beta)

# Short Description
The pipeline is constructed using nextflow as workflow manager running in a docker container.
It is able to identify species of TB/nTM from positive MGIT cultures.
To do so it uses a hsp65 database for fast identification coupled with
a Metagenomic method using centrifuge to identify on genome level.
For TB it also is able to identify subspecies.
Results are presented in automated pdf and html reports.

# Authors code
```bash
J.P.M. Coolen
H. Severin
```

# Tool References
fastp
> Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu; fastp: an ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, Volume 34, Issue 17, 1 September 2018, Pages i884–i890, https://doi.org/10.1093/bioinformatics/bty560

centrifuge
> https://ccb.jhu.edu/software/centrifuge/manual.shtml

KMA
> Clausen, P.T.L.C., Aarestrup, F.M. & Lund, O. Rapid and precise alignment of raw reads against redundant databases with KMA. BMC Bioinformatics 19, 307 (2018). https://doi.org/10.1186/s12859-018-2336-6

SNP-IT
>Lipworth S, Jajou R, de Neeling A, et al. SNP-IT Tool for Identifying Subspecies and Associated Lineages of Mycobacterium tuberculosis Complex. Emerging Infectious Diseases. 2019;25(3):482-488. doi:10.3201/eid2503.180894.

gatk
> Van der Auwera GA & O'Connor BD. (2020). Genomics in the Cloud: Using Docker, GATK, and WDL in Terra (1st Edition). O'Reilly Media.

bwa-mem
> Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60. [PMID: 19451168]

delly
> Tobias Rausch, Thomas Zichner, Andreas Schlattl, Adrian M. Stuetz, Vladimir Benes, Jan O. Korbel. DELLY: structural variant discovery by integrated paired-end and split-read analysis. Bioinformatics. 2012 Sep 15;28(18):i333-i339. https://doi.org/10.1093/bioinformatics/bts378

# Citation
```bash
Currently in submission fase:

MGIT enriched shotgun metagenomics for routine identification of nontuberculous mycobacteria: a route to personalized healthcare

Jodie A. Schildkraut1*, Jordy P.M. Coolen1*, Heleen Severin1, Ellen Koenraad1, Nicole Aalders1, Willem J.G. Melchers1, Wouter Hoefsloot2, Heiman F.L. Wertheim1, Jakko van Ingen1
1)	Radboudumc Center for Infectious Diseases, Department of Medical Microbiology, Radboud University Medical Center, Nijmegen, the Netherlands
2)	Radboudumc Center for Infectious Diseases, Department of pulmonary diseases, Radboud University Medical Center, Nijmegen, the Netherlands

*Authors contributed equally
```
