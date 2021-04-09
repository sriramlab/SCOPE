# SCOPE - (SCalable pOPulation structure inferencE)

SCOPE is a method for performing scalable population structure inference on biobank-scale genomic data. SCOPE is able to perform 

### Prerequisites

The following packages are required on a linux machine to compile and use the software package.

```
g++
cmake
make
```

### Installing

Installing SCOPE is fairly simple. Just issue the following commands on a linux machine

```
git clone https://github.com/sriramlab/SCOPE.git
cd SCOPE
mkdir build
cd build
cmake ..
make
```

## Documentation for ProPCA

After compiling the executable propca is present in the build directory.
Running the propca is fairly simple and can be done in two different ways

* ``./propca -p <parameter_file>``
* ``./propca <various_command_line arguments>``

### Parameters

The values in the brackets are the command line flags for running the code without the parameter file.

```
* genotype (-g) : The path of the genotype file or plink bed file prefix
* num_evec (-k) : The number of eigen vectors to output (default: 5)
* l (-l) : The extra calculation to be performed so that k_effective  = k + l (default: num_evec)
* max_iterations (-m) : The maximum number of iterations to run the EM for (default: num_evec + 2)
* debug (-v) : Enabling debug mode to output various debug informations (default: false). Need to build with DEBUG=1 as described above for this flag to work.
* accuracy (-a) : Output the likelihood computation as a function of iterations (default: false)
* convergence_limit (-cl) : The value of the threshold telling the algorithm that it has converged (default: -1, meaning no auto termination condition )
* output_path (-o) : The output prefix along with the path where the results will be stored
* accelerated_em (-aem) : The flag stating whether to use accelerated EM or not (default: 0).
* var_normalize (-vn) : The flag stating whether to perform varinance normalization or not (default: false).
* fast_mode (-nfm) : The flag whether to use a fast mode for the EM algorithm(default: true). Note: Setting the -nfm (NOT fast_mode) at command line will use a slow version of EM.
* missing (-miss) : This flag states whether there is any missing data present in the genotype matrix or not. 
* text_version (-txt) : This flag makes the input genotype file to be in the text format as described below. If not used, plink format will be used. (default: false)
* memory_efficient (-mem) : The flag states whether to use a memory effecient version for the EM algorithm or not. The memory efficient version is a little slow than the not efficient version (default: false)
* nthreads (-nt): Number of threads to use (default: 1)
* seed (-seed): Seed to use (default: system time)

```

An example parameter file is provided in the examples directory.

You can run the code using the command:

```
../build/propca -p par.txt
``` 

The equivalent command to issue for running the same code from the examples directory is:

```
../build/propca -g example.geno -k 5 -l 2 -m 20 -a -cl 0.001 -o example_ -aem 1 -vn -nfm -txt
```

ProPCA wil generate three files containing the eigenvectors/principal components, projections, and eigenvalues.

#### Second:

The inout can be in the plink binary format, as descibed at [Plink BED](https://www.cog-genomics.org/plink/1.9/input#bed)

Make sure to set the text_version to false in the parameter file, or don't use the -txt command line flag, when running. 

## S

* [Eigen](http://eigen.tuxfamily.org/) - The Linear algebra library for C++

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
