#ifndef PROPCA_PROPCA_H_
#define PROPCA_PROPCA_H_


#include "genotype.h"

//struct timespec t0;

class ProPCA {
 public:
 	Genotype g;
	MatrixXdr geno_matrix; //(p,n)

	int MAX_ITER;
	int k, p, n;
	int k_orig;

	MatrixXdr c; //(p,k)
	MatrixXdr x; //(k,n)
	MatrixXdr v; //(p,k)
	MatrixXdr means; //(p,1)
	MatrixXdr stds; //(p,1)

	// How to batch columns:
	int blocksize;
	double **partialsums;
	double *sum_op;

	// Intermediate computations in E-step.
	// Size = 3^(log_3(n)) * k
	double **yint_e;
	//  n X k
	double ***y_e;

	// Intermediate computations in M-step.
	// Size = nthreads X 3^(log_3(n)) * k
	double **yint_m;
	//  nthreads X log_3(n) X k
	double ***y_m;

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

	void multiply_y_pre_fast_thread(int begin, int end, MatrixXdr &op, int Ncol_op, double *yint_m, double **y_m, double *partialsums, MatrixXdr &res);
	void multiply_y_post_fast_thread(int begin, int end, MatrixXdr &op, int Ncol_op, double *yint_e, double **y_e, double *partialsums);

	/*
	 * M-step: Compute C = Y E 
	 * Y : p X n genotype matrix
	 * E : n K k matrix: X^{T} (XX^{T})^{-1}
	 * C = p X k matrix
	 *
	 * op : E
	 * Ncol_op : k
	 * res : C
	 * subtract_means :
	 */
	void multiply_y_pre_fast(MatrixXdr &op, int Ncol_op, MatrixXdr &res, bool subtract_means);

	/*
	 * E-step: Compute X = D Y 
	 * Y : p X n genotype matrix
	 * D : k X p matrix: (C^T C)^{-1} C^{T}
	 * X : k X n matrix
	 *
	 * op_orig : D
	 * Nrows_op : k
	 * res : X
	 * subtract_means :
	 */
	void multiply_y_post_fast(MatrixXdr &op_orig, int Nrows_op, MatrixXdr &res, bool subtract_means);

	void multiply_y_pre_naive_mem(MatrixXdr &op, int Ncol_op, MatrixXdr &res);
	void multiply_y_post_naive_mem(MatrixXdr &op, int Nrows_op, MatrixXdr &res);

	void multiply_y_pre_naive(MatrixXdr &op, int Ncol_op, MatrixXdr &res);

	void multiply_y_post_naive(MatrixXdr &op, int Nrows_op, MatrixXdr &res);

	void multiply_y_post(MatrixXdr &op, int Nrows_op, MatrixXdr &res, bool subtract_means);

	void multiply_y_pre(MatrixXdr &op, int Ncol_op, MatrixXdr &res, bool subtract_means);

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
