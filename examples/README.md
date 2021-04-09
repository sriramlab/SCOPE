## Example Files

We provide example source files in `source_files`. The script `run_scope_examples.sh` will run SCOPE both unsupervised and supervised. This script may need to be modified to correctly point to the path of your SCOPE binary. Expected output from this script can be found in `results` under their respective folders.

### Source Files

The example files in `source_files` were generated using the `simulateA.py` script found in `misc/simulations` from the root of the repository. The PLINK formatted frequency file was generated using the `generate_plink_frq.py` script found in `misc/simulations` from the root of the repository.

#### Contents

+ Example PLINK binary files: `example_1k.bed/bim/fam`
+ True allele frequencies: `example_1k_allele_frequencies.txt`
+ True proportions: `example_1k_proportions.txt`
+ True frequencies in PLINK format: `example_1k.plink.freq`

### Results

The `results` subdirectory contains the expected output from running `run_scope_example.sh`.

+ `*_Phat.txt` are the estimated minor allele frequencies
+ `*_V.txt` is the estimated latent subspace
+ `*_Qhat.txt` is the estimated admixture proportions


