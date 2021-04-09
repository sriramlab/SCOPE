
# Get HGDP
wget https://hagsc.org/hgdp/data/hgdp.zip
unzip hgdp.zip
mv hgdp/* .
rm -rf __MACOSX
rmdir hgdp

# Get Sample Information
wget http://rosenberglab.stanford.edu/data/rosenberg2006ahg/SampleInformation.txt

# Remove column sums
head -n -1 SampleInformation.txt > SampleInformation2.txt
mv SampleInformation2.txt SampleInformation.txt

# Convert to MAP/PED
python HGDP_text_to_tped.py

# Convert to PLINK BED/BIM/FAM
plink --tfile HGDP_940 --make-bed --out HGDP_940

# Calculate HGDP frequencies
plink --bfile HGDP_940 --freq --out HGDP_940



