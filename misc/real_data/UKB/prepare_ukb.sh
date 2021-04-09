#!/bin/env bash

ukb_prefix=""

# MAF Filter
plink --bfile $ukb_prefix --maf 0.01 --make-bed --out maf_01

# LD Pruning
plink --bfile maf_01 --indep-pairwise 50 kb 80 0.1
plink --bfile maf_01 --extract plink.prune.in --make-bed --out pruned_ukb

# Long-range LD removal
plink --bfile pruned_ukb --exclude long_range_ld.filter --make-bed --out ukb_no_long_ld

