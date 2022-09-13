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

![Alt text](logo/flowchart.png?raw=true "Flowdiagram")
