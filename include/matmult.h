#ifndef PROPCA_MATMULT_H_
#define PROPCA_MATMULT_H_


#include "genotype.h"


class MatMult {
 public:
 	Genotype g;
	MatrixXdr geno_matrix;  // (p,n)

	bool debug = false;
	bool var_normalize = false;
	bool memory_efficient = false;
	bool missing = false;
	bool fast_mode = true;
	int nthreads = 1;

	// How to batch columns:
	int blocksize;  // k
	int hsegsize;  // = log_3(n)
	int hsize;
	int vsegsize;  // = log_3(p)
	int vsize;

	double **partialsums;
	double *sum_op;

	// Intermediate computations in E-step.
	double **yint_e;  // Size = 3^(log_3(n)) * k
	double ***y_e;    // n X k

	// Intermediate computations in M-step.
	double **yint_m;  // Size = nthreads X 3^(log_3(n)) * k
	double ***y_m;    // nthreads X log_3(n) X k

	MatMult() {}

	MatMult(Genotype &xg,
			bool xdebug,
			bool xvar_normalize,
			bool xmemory_efficient,
			bool xmissing,
			bool xfast_mode,
			int xnthreads,
			int xk);

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
