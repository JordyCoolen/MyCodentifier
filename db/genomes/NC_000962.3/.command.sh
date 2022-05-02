#!/bin/bash -ue
genome_dir="/workflow/db/genomes/NC_000962.3"
# get the genome files from NCBI if not already present in the target directory
fetch_refs.py --accession NC_000962.3 --table /workflow/db/centrifuge_db/conversion_table.txt --format all --outdir ${genome_dir}

# unzip all the genome files in the current directory (in parallel)
parallel 'zcat {} > {/.}' ::: ${genome_dir}/*.gz

# replace type-material INFO (2020-02-18, recently added by NCBI causes problems with snpEff)
sed -i 's/type-material/typematerial/g' NC_000962.3_genomic.gff
sed -i 's/collection-date/collectiondate/g' NC_000962.3_genomic.gff
sed -i 's/lat-lon/latlon/g' NC_000962.3_genomic.gff
sed -i 's/isolation-source/isolationsource/g' NC_000962.3_genomic.gff
sed -i 's/isolation-source/isolationsource/g' NC_000962.3_genomic.gff
sed -i 's/culture-collection/culturecollection/g' NC_000962.3_genomic.gff

# filter accession code chromosomal DNA (exlude Plasmidic DNA)
samtools faidx NC_000962.3_genomic.fna NC_000962.3 > NC_000962.3_filt_genomic.fa
samtools faidx NC_000962.3_filt_genomic.fa
