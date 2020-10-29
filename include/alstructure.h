#ifndef PROPCA_ALSTRUCTURE_H_
#define PROPCA_ALSTRUCTURE_H_

#include <random>
#include <string>
#include <fstream>

#include "genotype.h"
#include "matmult.h"


struct options {
	std::string GENOTYPE_FILE_PATH;
	std::string ROWSPACE_FILE_PATH;
	std::string INITIAL_FILE_PATH;
	std::string OUTPUT_PATH;
	int max_iterations;
	int num_of_evec;
	bool getaccuracy;
	bool debugmode;
	double convergence_limit;
	bool memory_efficient;
	bool fast_mode;
	bool missing;
	bool text_version;
	bool fhat_version;
	bool fhattrunc_version;
	int nthreads;
	unsigned int seed;
	bool given_seed;
};


class ALStructure {
 public:
 	options command_line_opts;
 	
 	Genotype g;
	MatrixXdr geno_matrix;  //(p,n)
	MatMult mm;

	Eigen::VectorXd D;      //(n,1) for LSE
	unsigned int nops;      // number ops for LSE

	MatrixXdr V;            // (n,k) for truncated ALS
	MatrixXdr Fhat;         // (p,n) for truncated ALS
	MatrixXdr Phat;         // (p,k) for truncated ALS
	MatrixXdr Qhat;         // (k,n) for truncated ALS
	MatrixXdr Qhat_old;     // (k,n) for truncated ALS
	MatrixXdr diff;         // (k,n) for truncated ALS

	double rmse;
	int MAX_ITER;
	int k, p, n;
	long long int niter;
	int seed;

	//std::ofstream fp;

	clock_t total_begin;  //= clock();

	bool debug;  // = false;
	double convergence_limit;
	bool memory_efficient;  // = false;
	bool missing;  // = false;
	bool fast_mode;  // = true;
	bool text_version;  // = false;
	bool fhat_version;  // = false;
	bool fhattrunc_version;  // = false;
	int nthreads;  // = 1;
	std::string output_path;

	ALStructure() {}
	ALStructure(int argc, char const *argv[]);

	void printCorrectUsage(void);

	void solve_for_Qhat(void);
	void solve_for_Phat(void);

	void initialize(std::default_random_engine &prng_eng);
	void truncated_alternating_least_squares(void);

	int run();

	void write_matrix(MatrixXdr &mat, const std::string file_name);
	void write_vector(Eigen::VectorXd &vec, const std::string file_name);

	unsigned int rows();  // rows for Spectra
	unsigned int cols();  // cols for Spectra
	void perform_op(const double* x_in, double* y_out);  // operation for Spectra

};

#endif  // PROPCA_ALSTRUCTURE_H_
