#!/usr/bin/env Rscript
list.of.packages <- c("jsonlite", "vcfR", "argparse")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(vcfR)
library(argparse)

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-f", "--fasta", type="character",
                    help="fasta file path")
parser$add_argument("-g", "--gff", type="character",
                    help="gff file path")
parser$add_argument("-v", "--vcf", type="character",
                    help="vcf path")

args <- parser$parse_args()

## Arg1 exception
if(is.null(args$fasta)) {
  cat("please provide a fasta file\n")
}
## Arg1 exception
if(is.null(args$gff)) {
  cat("please provide a gff file\n")
}
## Arg1 exception
if(is.null(args$vcf)) {
  cat("please provide a vcf file\n")
}


# Input the files.
vcf <- read.vcfR(args$vcf, verbose = FALSE)
fasta <- ape::read.dna(args$fasta, format = "fasta")
gff <- read.table(args$gff, sep="\t", quote="")

# Create a chromR object.
chrom <- create.chromR(name="Supercontig", vcf=vcf, seq=fasta, ann=gff, verbose=TRUE)

chrom <- masker(chrom, min_QUAL=0, min_DP=25, max_DP=650, min_MQ=59.5, max_MQ=60.5)
chrom <- proc.chromR(chrom, verbose = TRUE)

# name the plot jpg
# jpeg('vcf_plots.jpg')

##commented out
#chromoqc(chrom, dp.alpha = 22)

# reset plots
# dev.off()

# name the plot jpg
# jpeg('chrom_plots.jpg')

## commented out
#plot(chrom)

## commented out
#dev.off()

# print head from chrom object
# head(chrom)

## print vcf meta data
vcf

vcf_out <- capture.output(vcf)

basename = basename(tools::file_path_sans_ext(args$vcf))
fileout = paste(basename, "_summary_results.txt", sep="")

write(vcf_out, file = fileout,
      ncolumns = if(is.character(vcf)) 1 else 5,
      append = FALSE, sep = " ")

