#!/bin/env bash

# Download TGP
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20120131_omni_genotypes_and_intensities/Omni25_genotypes_2141_samples.b37.vcf.gz

# Convert to PLINK format, remove related individuals, MAF filter
vcftools --gzvcf Omni25_genotypes_2141_samples.b37.vcf.gz --keep TGP_unrel.txt --maf 0.01 --max-missing 0.95 --max-alleles 2 --min-alleles 2 --chr 1 --chr 2 --chr 3 --chr 4 --chr 5 --chr 6 --chr 7 --chr 8 --chr 9 --chr 10 --chr 11 --chr 12 --chr 13 --chr 14 --chr 15 --chr 16 --chr 17 --chr 18 --chr 19 --chr 20 --chr 21 --chr 22 --plink-tped --out TGP_plink

# Convert to PLINK BED/BIM/FAM
plink --tfile TGP_plink --make-bed --out TGP_1718

