#!/bin/env bash

# Download HO
wget https://reich.hms.harvard.edu/sites/reich.hms.harvard.edu/files/inline-files/EuropeFullyPublic.tar.gz
tar -xvzf EuropeFullyPublic.tar.gz
mv EuropeFullyPublic/data.* .
rm -r EuropeFullyPublic

# Convert to PLINK BED/BIM/FAM using EIGENSOFT
touch ho_plink.convertf_par

echo "genotypename: data.geno" >> ho_plink.convertf_par
echo "snpname: data.snp" >> ho_plink.convertf_par
echo "indivname: data.ind" >> ho_plink.convertf_par
echo "outputformat: PACKEDPED" >> ho_plink.convertf_par
echo "genooutfilename: ho.bed" >> ho_plink.convertf_par
echo "snpoutfilename: ho.bim" >> ho_plink.convertf_par
echo "indoutfilename: ho.fam" >> ho_plink.convertf_par

convertf -p ho_plink.convertf_par
awk '{if ($1 == 90) print "c"$0; else print $0}' ho.bim > ho2.bim
rm ho.bim
mv ho2.bim ho.bim

# Create final HO BED/BIM/FAM with filters
plink --bfile ho --remove ../ho_remove.samps --geno 0.01 --maf 0.05 --allow-extra-chr --make-bed --out ho_final --chr 1-26

