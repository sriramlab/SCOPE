#!/bin/env bash

# Calculate frequencies and get FST
plink --bfile TGP_1718 --freq --out tgp_freqs
plink --bfile TGP_1718 --fst --within tgp_pops.txt --out tgp_pops
plink --bfile TGP_1718 --fst --within tgp_super_pops.txt --out tgp_super_pops

# Generate parameter files for simulation scripts
awk '{print $5}' tgp_super_pops.fst | paste tgp_freqs.frq - > tgp_superpop_param.txt
awk '{print $5}' tgp_pops.fst | paste tgp_freqs.frq - > tgp_pop_param.txt 
