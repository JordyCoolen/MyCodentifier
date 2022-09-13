![Alt text](docs/logo/logo.png?raw=true "logo")

# TB/NTM species identification pipeline Version 1.0 (beta)

# Short Description
The pipeline is constructed using nextflow as workflow manager running in a docker container.
It is able to identify species of TB/nTM from positive MGIT cultures.
To do so it uses a hsp65 database for fast identification coupled with
a Metagenomic method using centrifuge to identify on genome level.
For TB it also is able to identify subspecies.
Results are presented in automated pdf and html reports.

## installation notes:
### OPTION1 (build docker image from scratch)
```bash
install docker https://www.docker.com/products/docker-desktop
git clone https://github.com/JordyCoolen/MyCodentifier.git
cd MyCodentifier
docker build --rm -t jonovox/mycodentifier:1.0 ./
```

### or

### OPTION2 (pull docker image from docker hub)
```bash
install docker https://www.docker.com/products/docker-desktop
git clone https://github.com/JordyCoolen/MyCodentifier.git
cd MyCodentifier
docker pull jonovox/mycodentifier:1.0
```

### CONDA environments (prebuild)
download the newest release of the pipeline via:
https://surfdrive.surf.nl/files/index.php/s/azkU2t09zY1r9qB
and extract into the Mycodentifier folder

## Databases
- currently the databases are not provided
- we will provide in future versions the databases
- please contact us for more information

## create singularity .sif out of docker image
```bash
bash docker/create_singularity_image.sh <filename> <imagename>
```
example:
```bash
bash docker/create_singularity_image.sh mycodentifier_1_0 mycodentifier:1.0
```

## run TEST sample:
```bash
cd Mycodentifier
sh docker/run.sh mycodentifier jonovox/mycodentifier:1.0
cd /workflow
# test sample
nextflow run myco.nf --reads "/workflow/test/H37Rv_MB_R{1,2}.fastq.gz" --threads 8 --sampleName TEST \
--centrifuge_cont "/workflow/db/centrifuge_db/refv2" --centrifuge_WGS "/workflow/db/centrifuge_db/refv2" -resume
```

# Example usage:

## singularity use
```bash
singularity exec -B </path/to/newcode>:/workflow,</path/to/fastq.gz/files:/data mycodentifier_v1_72.sif \
nextflow run /workflow/myco.nf --threads 16 --outDir /workflow/output --snpEff_base /ifs/software/external/mmb/singularity_images/snpEff/ \
-resume --reads "/data/<samplename>_R{1,2}.fastq.gz" --sampleName <samplename>
```

```bash
singularity exec -B /ifs/data/mmb/nextseq/runs/200309_NB501730_0373_AHYC7YAFXY/demultiplex/merged_lanes:/data,/ifs/software/external/mmb/mycodentifier:/workflow \
mycodentifier_v1_72.sif nextflow run /workflow/myco.nf --reads '/data/20090667701-1nieuw_MB_R{1,2}.fastq.gz' --threads 16 --outDir /ifs/software/external/mmb/singularity_images \
--snpeff_base snpEff/ --sampleName 20090667701-1nieuw_MB
```

## Example of docker run:
```bash
docker run \
  -it \
  --rm \
  --name mycoprofiler2 \
  --mount type=bind,source=/location/to/centrifuge/contaminants/database,target=/workflow/db/centrifuge_contaminant,readonly \
  mycodentifier:latest /bin/bash
```

## Output files (TEST example)
```bash
output/TEST
├── GATK
│             └── raw.snps.indels_TEST.final_rg_MD.vcf
├── NC_000962.3
│             └── NC_000962.3_filt_genomic.dict
├── SNPFasta
│             └── TEST_altref.fasta.gz
├── SNPIT
│             ├── TEST_20220502142245_qc.json
│             └── TEST_20220502144828_qc.json
├── WGS_typing
│             ├── TEST_20220502142245.json
│             ├── TEST_20220502144828.json
│             └── TEST_report.file
├── annot
│             ├── TEST_snps_indels_filt_annot.vcf
│             └── snpEff_summary.html
├── bam_qc
│             ├── BAM_QC.txt
│             ├── TEST_20220502142245_qc.json
│             ├── TEST_20220502144828_qc.json
│             └── avg_dp.txt
├── bamfiles
│             ├── TEST.bam
│             ├── TEST.bam.bai
│             ├── TEST.final.bam
│             ├── TEST.final.bam.bai
│             ├── TEST.final_rg.bam -> /workflow/work/9e/5cb0668a42153d1a5d191b9ba9c917/TEST.final_rg.bam
│             ├── TEST.final_rg.bam.bai -> /workflow/work/9e/5cb0668a42153d1a5d191b9ba9c917/TEST.final_rg.bam.bai
│             ├── TEST.final_rg_MD.bam -> /workflow/work/07/61441a8a33aef5b0da3436460176df/TEST.final_rg_MD.bam
│             ├── TEST.final_rg_MD.bam.bai -> /workflow/work/07/61441a8a33aef5b0da3436460176df/TEST.final_rg_MD.bam.bai
│             ├── TEST_subsampled_2C.bam
│             └── TEST_subsampled_2C.bam.bai
├── centrifuge
│             ├── TEST_20220502142245.json
│             ├── TEST_20220502144828.json
│             ├── TEST_report.file
│             └── TEST_spreporter.out
├── contaminants
│             ├── TEST_20220502142245.json
│             ├── TEST_20220502144828.json
│             └── TEST_report.file
├── db
│             └── genomes
│                 └── NC_000962.3
│                     ├── NC_000962.3_filt_genomic.fa
│                     ├── NC_000962.3_filt_genomic.fa.fai
│                     ├── NC_000962.3_genomic.gbff
│                     └── NC_000962.3_genomic.gff
├── delly
│             └── TEST_delly.vcf.gz
├── fastp
│             ├── H37Rv_MB_R1.fastq_fastp.fastq.gz
│             ├── H37Rv_MB_R2.fastq_fastp.fastq.gz
│             ├── TEST.fastp_json
│             ├── TEST_20220502142245.json
│             └── TEST_20220502144828.json
├── filtered
│             ├── TEST_delly_filt.vcf
│             ├── TEST_filtered_snps.vcf
│             ├── TEST_flagged_snps.vcf
│             ├── TEST_gatk_delly_merged.sorted.vcf
│             └── raw.snps.indels_TEST.final_rg_MD_flagged_indels.vcf
├── hsp65
│             ├── TEST_20220502142245.json
│             ├── TEST_20220502144828.json
│             ├── hsp65.aln
│             └── hsp65.res
├── json
│             └── TEST.fastp_json
├── merged
│             └── TEST_snps_indels_filt.vcf
├── mutations
│             ├── TEST_20220502142245_qc.json
│             └── TEST_20220502144828_qc.json
└── report
    ├── TEST.html
    ├── TEST.pdf
    ├── TEST_20220502142245_qc.json
    └── TEST_20220502144828_qc.json
```

## Flow diagram of current pipeline:

![Alt text](docs/logo/flowchart.png?raw=true "Flowdiagram")

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
