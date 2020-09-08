#include "matmult.h"
#include "mailman.h"

#include <thread>

#include <Eigen/Dense>
#include <Eigen/Core>

MatMult::MatMult(Genotype &gg, int k) {
	g = gg;

	blocksize = k;

	hsegsize = g.segment_size_hori;  // = log_3(n)
	hsize = pow(3, hsegsize);
	vsegsize = g.segment_size_ver;  // = log_3(p)
	vsize = pow(3, vsegsize);

	partialsums = new double*[nthreads];
	yint_m = new double*[nthreads];
	for (int t = 0; t < nthreads; t++) {
		partialsums[t] = new double[blocksize];
		yint_m[t] = new double[hsize*blocksize];
		memset (yint_m[t], 0, hsize*blocksize * sizeof(double));
	}

	sum_op = new double[blocksize];

	yint_e = new double* [nthreads];
	for (int t = 0; t < nthreads; t++) {
		yint_e[t] = new double[hsize*blocksize];
		memset (yint_e[t], 0, hsize*blocksize * sizeof(double));
	}

	y_e  = new double**[nthreads];
	for (int t = 0 ; t < nthreads ; t++) {
		y_e[t]  = new double*[g.Nindv];
		for (int i = 0 ; i < g.Nindv ; i++) {
			y_e[t][i] = new double[blocksize];
			memset (y_e[t][i], 0, blocksize * sizeof(double));
		}
	}

	y_m = new double**[nthreads];
	for (int t = 0; t < nthreads; t++) {
		y_m[t] = new double*[hsegsize];
		for (int i = 0; i < hsegsize; i++) {
			y_m[t][i] = new double[blocksize];
		}
	}
}

void MatMult::multiply_y_pre_fast_thread(int begin, int end, MatrixXdr &op, int Ncol_op, double *yint_m, double **y_m, double *partialsums, MatrixXdr &res) {
	for (int seg_iter = begin; seg_iter < end; seg_iter++) {
		mailman::fastmultiply(g.segment_size_hori, g.Nindv, Ncol_op, g.p[seg_iter], op, yint_m, partialsums, y_m);
		int p_base = seg_iter * g.segment_size_hori;
		for (int p_iter = p_base; (p_iter < p_base + g.segment_size_hori) && (p_iter < g.Nsnp); p_iter++) {
			for (int k_iter = 0; k_iter < Ncol_op; k_iter++) {
				res(p_iter, k_iter) = y_m[p_iter - p_base][k_iter];
			}
		}
	}
}

void MatMult::multiply_y_post_fast_thread(int begin, int end, MatrixXdr &op, int Ncol_op, double *yint_e, double **y_e, double *partialsums) {
	for (int i = 0; i < g.Nindv; i++) {
		memset (y_e[i], 0, blocksize * sizeof(double));
	}

	for (int seg_iter = begin; seg_iter < end; seg_iter++) {
		mailman::fastmultiply_pre(g.segment_size_hori, g.Nindv, Ncol_op, seg_iter * g.segment_size_hori, g.p[seg_iter], op, yint_e, partialsums, y_e);
	}
}


void MatMult::multiply_y_pre_fast(MatrixXdr &op, int Ncol_op, MatrixXdr &res, bool subtract_means) {
	for (int k_iter = 0; k_iter < Ncol_op; k_iter++) {
		sum_op[k_iter] = op.col(k_iter).sum();
	}

	#if DEBUG == 1
		if (debug) {
			print_time();
			std::cout << "Starting mailman on premultiply" << std::endl;
			std::cout << "Nops = " << Ncol_op << "\t" << g.Nsegments_hori << std::endl;
			std::cout << "Segment size = " << g.segment_size_hori << std::endl;
			std::cout << "Matrix size = " << g.segment_size_hori << "\t" << g.Nindv << std::endl;
			std::cout << "op = " <<  op.rows() << "\t" << op.cols() << std::endl;
		}
	#endif

	//TODO: Memory Effecient SSE FastMultipy

	nthreads = (nthreads > g.Nsegments_hori) ? g.Nsegments_hori: nthreads;

	std::thread th[nthreads];
	int perthread = g.Nsegments_hori / nthreads;
	// std::cout << g.Nsegments_hori << "\t" << nthreads << "\t" << perthread << std::endl;
	int t = 0;
	for (; t < nthreads - 1; t++) {
	// std::cout << "Launching thread " << t << std::endl;
		th[t] = std::thread(&MatMult::multiply_y_pre_fast_thread, this, t * perthread , (t+1)*perthread, std::ref(op), Ncol_op, yint_m[t], y_m[t], partialsums[t], std::ref(res));
	}

	th[t] = std::thread(&MatMult::multiply_y_pre_fast_thread, this, t * perthread , g.Nsegments_hori  - 1, std::ref(op), Ncol_op, yint_m[t], y_m[t], partialsums[t], std::ref(res));

	for (int t = 0; t < nthreads; t++) {
		th[t].join();
	}

	// for(int seg_iter = 0; seg_iter < g.Nsegments_hori - 1; seg_iter++){
	//	mailman::fastmultiply ( g.segment_size_hori, g.Nindv, Ncol_op, g.p[seg_iter], op, yint_m, partialsums, y_m);
	//	int p_base = seg_iter * g.segment_size_hori;
	//	for(int p_iter=p_base; (p_iter < p_base + g.segment_size_hori) && (p_iter < g.Nsnp) ; p_iter++ ){
	//		for(int k_iter = 0; k_iter < Ncol_op; k_iter++)
	//			res(p_iter, k_iter) = y_m [p_iter - p_base][k_iter];
	//	}
	//}

	int last_seg_size = (g.Nsnp % g.segment_size_hori != 0) ? g.Nsnp % g.segment_size_hori : g.segment_size_hori;
	mailman::fastmultiply(last_seg_size, g.Nindv, Ncol_op, g.p[g.Nsegments_hori-1], op, yint_m[0], partialsums[0], y_m[0]);
	int p_base = (g.Nsegments_hori - 1) * g.segment_size_hori;
	for (int p_iter = p_base; (p_iter < p_base + g.segment_size_hori) && (p_iter < g.Nsnp); p_iter++) {
		for (int k_iter = 0; k_iter < Ncol_op; k_iter++) {
			res(p_iter, k_iter) = y_m[0][p_iter - p_base][k_iter];
		}
	}

	#if DEBUG == 1
		if (debug) {
			print_time();
			std::cout << "Ending mailman on premultiply" << std::endl;
		}
	#endif

	if (!subtract_means) {
		return;
	}

	for (int p_iter = 0; p_iter < g.Nsnp; p_iter++) {
		for (int k_iter = 0; k_iter < Ncol_op; k_iter++) {
			res(p_iter, k_iter) = res(p_iter, k_iter) - (g.get_col_mean(p_iter) * sum_op[k_iter]);
			if (var_normalize) {
				res(p_iter, k_iter) = res(p_iter, k_iter) / (g.get_col_std(p_iter));
			}
		}
	}
}


void MatMult::multiply_y_post_fast(MatrixXdr &op_orig, int Nrows_op, MatrixXdr &res, bool subtract_means) {
	MatrixXdr op;
	op = op_orig.transpose();

	if (var_normalize && subtract_means) {
		for (int p_iter = 0; p_iter < g.Nsnp; p_iter++) {
			for (int k_iter = 0; k_iter < Nrows_op; k_iter++) {
				op(p_iter, k_iter) = op(p_iter, k_iter) / (g.get_col_std(p_iter));
			}
		}
	}

	#if DEBUG == 1
		if (debug) {
			print_time();
			std::cout << "Starting mailman on postmultiply" << std::endl;
		}
	#endif

	int Ncol_op = Nrows_op;

	nthreads = (nthreads > g.Nsegments_hori) ? g.Nsegments_hori: nthreads;

	std::thread th[nthreads];
	int perthread = g.Nsegments_hori / nthreads;
	// std::cout << "post: " << g.segment_size_hori << "\t" << g.Nsegments_hori << "\t" << nthreads << "\t" << perthread << std::endl;
	int t = 0;
	for (; t < nthreads - 1; t++) {
	// std::cout << "Launching " << t << std::endl;
		th[t] = std::thread(&MatMult::multiply_y_post_fast_thread, this, t * perthread, (t+1) * perthread, std::ref(op), Ncol_op, yint_e[t], y_e[t], partialsums[t]);
	}
	// std::cout << "Launching " << t << std::endl;
	th[t] = std::thread(&MatMult::multiply_y_post_fast_thread, this, t * perthread, g.Nsegments_hori - 1, std::ref(op), Ncol_op, yint_e[t], y_e[t], partialsums[t]);
	for (int t = 0; t < nthreads; t++) {
		th[t].join();
	}
	// std::cout << "Joined "<< std::endl;

	// int seg_iter;
	// for(seg_iter = 0; seg_iter < g.Nsegments_hori-1; seg_iter++){
	// 	mailman::fastmultiply_pre (g.segment_size_hori, g.Nindv, Ncol_op, seg_iter * g.segment_size_hori, g.p[seg_iter], op, yint_e, partialsums[0], y_e);
	// }

	for (int t = 1; t < nthreads; t++) {
		for (int n_iter = 0; n_iter < g.Nindv; n_iter++) {
			for (int k_iter = 0; k_iter < Ncol_op; k_iter++) {
				y_e[0][n_iter][k_iter] += y_e[t][n_iter][k_iter];
			}
		}
	}

	int last_seg_size = (g.Nsnp % g.segment_size_hori != 0) ? g.Nsnp % g.segment_size_hori : g.segment_size_hori;
	mailman::fastmultiply_pre(last_seg_size, g.Nindv, Ncol_op, (g.Nsegments_hori-1) * g.segment_size_hori,
							  g.p[g.Nsegments_hori-1], op, yint_e[0], partialsums[0], y_e[0]);

	for (int n_iter = 0; n_iter < g.Nindv; n_iter++)  {
		for (int k_iter = 0; k_iter < Ncol_op; k_iter++) {
			res(k_iter, n_iter) = y_e[0][n_iter][k_iter];
			y_e[0][n_iter][k_iter] = 0;
		}
	}

	#if DEBUG == 1
		if (debug) {
			print_time();
			std::cout << "Ending mailman on postmultiply" << std::endl;
		}
	#endif


	if (!subtract_means) {
		return;
	}

	double *sums_elements = new double[Ncol_op];
	memset (sums_elements, 0, Nrows_op * sizeof(int));

	for (int k_iter = 0; k_iter < Ncol_op; k_iter++) {
		double sum_to_calc = 0.0;
		for (int p_iter = 0; p_iter < g.Nsnp; p_iter++) {
			sum_to_calc += g.get_col_mean(p_iter) * op(p_iter, k_iter);
		}
		sums_elements[k_iter] = sum_to_calc;
	}
	for (int k_iter = 0; k_iter < Ncol_op; k_iter++) {
		for (int n_iter = 0; n_iter < g.Nindv; n_iter++) {
			res(k_iter, n_iter) = res(k_iter, n_iter) - sums_elements[k_iter];
		}
	}
}

void MatMult::multiply_y_pre_naive_mem(MatrixXdr &op, int Ncol_op, MatrixXdr &res) {
	for (int p_iter = 0; p_iter < g.Nsnp; p_iter++) {
		for (int k_iter = 0; k_iter < Ncol_op; k_iter++) {
			double temp = 0;
			for (int n_iter = 0; n_iter < g.Nindv; n_iter++) {
				temp+= g.get_geno(p_iter, n_iter, var_normalize) * op(n_iter, k_iter);
			}
			res(p_iter, k_iter) = temp;
		}
	}
}

void MatMult::multiply_y_post_naive_mem(MatrixXdr &op, int Nrows_op, MatrixXdr &res) {
	for (int n_iter = 0; n_iter < g.Nindv; n_iter++) {
		for (int k_iter = 0; k_iter < Nrows_op; k_iter++) {
			double temp = 0;
			for(int p_iter = 0; p_iter < g.Nsnp; p_iter++)
				temp += op(k_iter, p_iter) * (g.get_geno(p_iter, n_iter, var_normalize));
			res(k_iter, n_iter) = temp;
		}
	}
}

void MatMult::multiply_y_post(MatrixXdr &op, int Nrows_op, MatrixXdr &res, bool subtract_means) {
    if (fast_mode) {
        multiply_y_post_fast(op, Nrows_op, res, subtract_means);
    } else {
		if(memory_efficient)
			multiply_y_post_naive_mem(op, Nrows_op, res);
		else
			res = op * geno_matrix;
	}
}

void MatMult::multiply_y_pre(MatrixXdr &op, int Ncol_op, MatrixXdr &res, bool subtract_means) {
    if (fast_mode) {
        multiply_y_pre_fast(op, Ncol_op, res, subtract_means);
    } else {
		if (memory_efficient) {
			multiply_y_pre_naive_mem(op, Ncol_op, res);
		} else {
			res = geno_matrix * op;
		}
	}
}

void MatMult::clean_up() {
	delete[] sum_op;

	for (int t = 0; t < nthreads; t++) {
		delete[] yint_e[t];
	}
	delete[] yint_e;

	for (int t = 0; t < nthreads; t++) {
		delete[] yint_m[t];
		delete[] partialsums[t];
	}
	delete[] yint_m;
	delete[] partialsums;

	for (int t = 0; t < nthreads; t++) {
		for (int i  = 0; i < hsegsize; i++) {
			delete[] y_m[t][i];
		}
		delete[] y_m[t];
	}
	delete[] y_m;

	for (int t = 0; t < nthreads; t++) {
		for (int i  = 0; i < g.Nindv; i++) {
			delete[] y_e[t][i];
		}
		delete[] y_e[t];
	}
	delete[] y_e;	
}
