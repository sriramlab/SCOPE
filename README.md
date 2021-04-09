# SCOPE - (SCalable pOPulation structure inferencE)

SCOPE is a method for performing scalable population structure inference on biobank-scale genomic data. SCOPE utilizes a likelihood-free framework that involves estimation of the individual allele frequency (IAF) matrix through a modified version of principal component analysis (PCA) known as latent subspace estimation (LSE) followed by alternating least squares (ALS) to transform the estimated IAF matrix into ancestral allele frequencies and admixture proportions. SCOPE utilizes two major optimizations to enable scalable inference. Firstly, SCOPE uses randomized eigendecomposition to efficiently estimate the latent subspace. Second, SCOPE uses the Mailman algorithm for fast matrix-vector multiplication involving the genotype matrix. 

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

This code is based on contributions from the following sources:
* [Eigen](http://eigen.tuxfamily.org/) - C++ template library for linear algebra
* [Spectra](https://spectralib.org/) - C++ Library For Large Scale Eigenvalue Problems
* [ProPCA](https://github.com/sriramlab/ProPCA) - Scalable probabilistic PCA for large-scale genetic variation data 

## Prerequisites

The following packages are required on a Linux machine to compile and use SCOPE.

```
g++
cmake
make
```

## Installing

To install SCOPE, run the following commands:

```
git clone https://github.com/sriramlab/SCOPE.git
cd SCOPE
mkdir build
cd build
cmake ..
make
```

An example script can be found in the `examples` subdirectory to test SCOPE.

## Documentation for SCOPE

### Parameters

SCOPE can be run the from the command line using the following options. At minimum, SCOPE require the path to the PLINK binary prefix.

```
* genotype (-g) : Path to PLINK binary prefix
* frequencies (-freq) : Path to PLINK MAF file for supervision (default: none)
* num_evec (-k) : Number of latent populations (default: 5)
* max_iterations (-m) : Maximum number of iterations for ALS (default: 1000)
* convergence_limit (-cl) : Convergence threshold for LSE and ALS (default: 0.00001)
* output_path (-o) : Output prefix (default: scope_)
* nthreads (-nt): Number of threads to use (default: 1)
* seed (-seed): Seed to use (default: system time)
```

To perform supervised population structure inference, provide the `-freq` parameter. The file needed for this parameter can be generated using `plink --maf`. If no frequency file provided, SCOPE will perform unsupervised population structure inference.

### Output

SCOPE will output the following files:

* `scope_V.txt`: the estimated latent subspace from LSE
* `scope_Phat.txt`: the estimated minor allele frequencies for the latent populations
* `scope_Qhat.txt`: the estimated admixture proportions for each individual

Each column of `Phat.txt` corresponds to a row of `Qhat.txt`. If `Qhat.txt` is transposed, the columns will correspond to the columns of `Phat.txt`. If running SCOPE in supervised mode, the order of the colums in `Phat.txt` corresponds to the order displayed in the PLINK MAF file.


