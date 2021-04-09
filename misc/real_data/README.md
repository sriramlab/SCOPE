## Real Datasets

We have included several scripts to regenerate the real datasets and parameters used in our manuscript.

### HGDP

This subdirectory contains a script to generate the Human Genome Diversity Project (HGDP) dataset. We also include `hgdp_param.txt.gz`, which is the parameter file used in our manuscript that can be given as input to our simulation scripts found in `misc/simulations` from the root of the repository.

### HO

This subdirectory contains a script to generate the Human Origins (HO) dataset. We have included `ho_remove.samps`, which contains manually curated samples to be removed (_e.g._ non-human samples). We also include a binary,`convertrf`, from the [EIGENSOFT suite](https://www.hsph.harvard.edu/alkes-price/software/) for converting the EIGENSTRAT format.

### TGP

This subdirectory contains a script to generate the 1000 Genomes Project (TGP) dataset. The script, `get_tgp.sh`, will download and convert the data to a PLINK binary fileset. The script, `gen_tgp_plink_metrics.sh`, will generate frequencies and F<sub>ST</sub> values based on population and superpopulation labels found in `tgp_pops.txt` and `tgp_super_pops.txt`, respectively. It will also generate the parameter files used by our simulations scripts found in `misc/simulations`. `TGP_unrel.txt` contains a manually curated list of unrelated individuals curated by the authors of [TeraStructure](https://github.com/StoreyLab/terastructure).

### UKB

This subdirectory contains a script to process the UK Biobank (UKB) dataset. The UKB dataset is not publicly available, so users will need to edit the script to specify the location of their UKB PLINK binary files for the array data. We also include `long_range_ld.filter`, which lists several SNPs involved in long-range linkage disequilibrium.

### Acknowledgements

We would like to acknowledge the authors of TeraStructure for providing [starter scripts](https://github.com/StoreyLab/terastructure/tree/master/scripts/data) for us in preparing these datasets.

