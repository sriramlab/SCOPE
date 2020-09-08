#ifndef PROPCA_MATMULT_H_
#define PROPCA_MATMULT_H_


#include "genotype.h"


class MatMult {
 public:
 	Genotype g;
	MatrixXdr geno_matrix; //(p,n)

	// How to batch columns:
	int blocksize;
	int hsegsize;  // = log_3(n)
	int hsize;
	int vsegsize;  // = log_3(p)
	int vsize;

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

	MatMult() {};
	MatMult(Genotype &gg, int k);

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

	void multiply_y_post(MatrixXdr &op, int Nrows_op, MatrixXdr &res, bool subtract_means);

	void multiply_y_pre(MatrixXdr &op, int Ncol_op, MatrixXdr &res, bool subtract_means);

	void clean_up();
};


#endif  // PROPCA_MATMULT_H_
