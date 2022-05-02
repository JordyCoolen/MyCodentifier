# test sample
nextflow run myco.nf --reads "/workflow/test/H37Rv_MB_R{1,2}.fastq.gz" --threads 8 --sampleName TEST \
--centrifuge_cont "/workflow/db/centrifuge_db/refv2" --centrifuge_WGS "/workflow/db/centrifuge_db/refv2" -resume