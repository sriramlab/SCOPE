#ifndef PROPCA_PROPCA_H_
#define PROPCA_PROPCA_H_

#include <string>

#include "genotype.h"
#include "matmult.h"

//struct timespec t0;

struct options {
	std::string GENOTYPE_FILE_PATH;
	std::string OUTPUT_PATH;
	int max_iterations;
	int num_of_evec;
	bool getaccuracy;
	bool debugmode;
	bool var_normalize;
	int accelerated_em;
	int l;  // What is this?
	double convergence_limit;
	bool memory_efficient;
	bool fast_mode;
	bool missing;
	bool text_version;
	int nthreads;
	int seed;
	bool given_seed;
};


class ProPCA {
 public:
 	options command_line_opts;
 	
 	Genotype g;
	MatrixXdr geno_matrix; //(p,n)
	MatMult mm;

	int MAX_ITER;
	int k, p, n;
	int k_orig;

	MatrixXdr c; //(p,k)
	MatrixXdr x; //(k,n)
	MatrixXdr v; //(p,k)
	MatrixXdr means; //(p,1)
	MatrixXdr stds; //(p,1)

	clock_t total_begin; //= clock();

	//options command_line_opts;
	bool debug = false;
	bool check_accuracy = false;
	bool var_normalize = false;
	int accelerated_em = 0;
	double convergence_limit;
	bool memory_efficient = false;
	bool missing = false;
	bool fast_mode = true;
	bool text_version = false;
	int nthreads = 1;
	std::string output_path;

	ProPCA() {}
	ProPCA(int argc, char const *argv[]);

	void printCorrectUsage(void);

	int run();

	/* Run one iteration of EM when genotypes are not missing
	 * c_orig : p X k matrix
	 * Output: c_new : p X k matrix 
	 */
	MatrixXdr run_EM_not_missing(MatrixXdr &c_orig);
	MatrixXdr run_EM_missing(MatrixXdr &c_orig);
	MatrixXdr run_EM(MatrixXdr &c_orig);

	std::pair<double, double> get_error_norm(MatrixXdr &c);

	void print_vals();
};

#endif  // PROPCA_PROPCA_H_
