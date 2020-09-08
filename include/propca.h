#ifndef PROPCA_PROPCA_H_
#define PROPCA_PROPCA_H_

#include "genotype.h"
#include "matmult.h"

//struct timespec t0;

class ProPCA {
 public:
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

	ProPCA() {};

	std::pair<double, double> get_error_norm(MatrixXdr &c);


	/* Run one iteration of EM when genotypes are not missing
	 * c_orig : p X k matrix
	 * Output: c_new : p X k matrix 
	 */
	MatrixXdr run_EM_not_missing(MatrixXdr &c_orig);
	MatrixXdr run_EM_missing(MatrixXdr &c_orig);
	MatrixXdr run_EM(MatrixXdr &c_orig);

	void print_vals();

	int run(int argc, char const *argv[]);
};

#endif  // PROPCA_PROPCA_H_
