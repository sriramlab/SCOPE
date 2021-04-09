## Simulation Scripts

We provide our simulation scripts used in our manuscript. We require the following packages:

+ numpy
+ [pyplink](https://lemieuxl.github.io/pyplink/)
+ scipy
+ [parmapper](https://github.com/Jwink3101/parmapper)

The script `install_packages.sh` will install `pyplink` and `parmapper`. The script `generate_simulations_python.sh` will generate several of the simulations utilized in our manuscript. The parameter files and how to create them can be found in `misc/real_data` from the root directory. Do **not** blindly run `generate_simulations_python.sh` as several of the generated datasets are very large and can easily exhaust computational resources.

### Simulation Scripts

We provide the simulation scripts for both scenarios used in our manuscript. The script `simulateA.py` generates simulation under the Pritchard-Stephens-Donneley (PSD) model. See `python simulateA.py --help` for additional information. Please note that the parameter file (`freq_fst_file` argument) is non-standard. It is the last column of the output from `plink --fst` appended to the end of `plink --freq`. Please see the scripts and files in `misc/real_data` to generate the TGP and HGDP parameter files we used in our manuscript.

The other script `simulateB.py` generated simulations under a spatial model as described in [Ochoa and Storey 2020](https://doi.org/10.1371/journal.pgen.1009241). Please see the original manuscript or our manuscript for more details. This script requires similar arugments as our other simulation script, but has a few differences, which can be found by running `python simulateB.py --help`.

### PLINK Conversion Script

SCOPE can take the output from `plink --freq` as input to perform supervised analysis. The script `generate_plink_frq.py` will convert output from previous SCOPE runs to match the output from `plink --freq`. It requires the path to the PLINK BIM file, the outputted frequencies, and the output filename as input.

