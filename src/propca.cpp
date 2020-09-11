/** 
 All of this code is written by Aman Agrawal 
 (Indian Institute of Technology, Delhi)
*/

#include "propca.h"

#include <iostream>
#include <fstream>
#include <iomanip>

#include <chrono>
#include "time.h"

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/SVD>
#include <Eigen/QR>

#include "config.h"
#include "storage.h"


struct timespec t0;


ProPCA::ProPCA(int argc, char const *argv[]) {
	// Set default values
	command_line_opts.num_of_evec = 5;
	command_line_opts.max_iterations = command_line_opts.num_of_evec + 2;
	command_line_opts.getaccuracy = false;
	command_line_opts.debugmode = false;
	command_line_opts.OUTPUT_PATH = "propca_";
	bool got_genotype_file = false;
	command_line_opts.var_normalize = false;
	command_line_opts.l = command_line_opts.num_of_evec;
	command_line_opts.accelerated_em = 0;
	command_line_opts.convergence_limit = -1.0;
	command_line_opts.memory_efficient = false;
	command_line_opts.fast_mode = true;
	command_line_opts.missing = false;
	command_line_opts.text_version = false;
	command_line_opts.nthreads = 1;
	command_line_opts.seed = -1;
	command_line_opts.given_seed = false;

	if (argc < 3) {
		std::cout << "Correct Usage is " << argv[0] << " -p <parameter file>" << std::endl;
		exit(-1);
	}

	if (strcmp(argv[1], "-p") == 0) {
		// Read arguments from configuration file
		std::string cfg_filename = std::string(argv[2]);
		ConfigFile cfg(cfg_filename);
		got_genotype_file = cfg.keyExists("genotype");
		command_line_opts.num_of_evec = cfg.getValueOfKey<int>("num_evec", 5);
		command_line_opts.max_iterations = cfg.getValueOfKey<int>("max_iterations", command_line_opts.num_of_evec + 2);
		command_line_opts.getaccuracy = cfg.getValueOfKey<bool>("accuracy", false);
		command_line_opts.debugmode = cfg.getValueOfKey<bool>("debug", false);
		command_line_opts.l = cfg.getValueOfKey<int>("l", command_line_opts.num_of_evec);
		command_line_opts.OUTPUT_PATH = cfg.getValueOfKey<std::string>("output_path", std::string("fastppca_"));
		command_line_opts.GENOTYPE_FILE_PATH = cfg.getValueOfKey<std::string>("genotype", std::string(""));
		command_line_opts.convergence_limit = cfg.getValueOfKey<double>("convergence_limit", -1.0);
		command_line_opts.var_normalize = cfg.getValueOfKey<bool>("var_normalize", false);
		command_line_opts.accelerated_em = cfg.getValueOfKey<int>("accelerated_em", 0);
		command_line_opts.memory_efficient = cfg.getValueOfKey<bool>("memory_efficient", false);
		command_line_opts.fast_mode = cfg.getValueOfKey<bool>("fast_mode", true);
		command_line_opts.missing = cfg.getValueOfKey<bool>("missing", false);
		command_line_opts.text_version = cfg.getValueOfKey<bool>("text_version", false);
		command_line_opts.nthreads = cfg.getValueOfKey<int>("nthreads", 1);
		command_line_opts.seed = cfg.getValueOfKey<int>("seed", -1);
		command_line_opts.given_seed = command_line_opts.seed >= 0 ? true: false;
	} else {
		// Read arguments from standard input
		bool got_max_iter = false;
		for (int i = 1; i < argc; i++) {
			if (i + 1 != argc) {
				if (strcmp(argv[i], "-g") == 0) {
					command_line_opts.GENOTYPE_FILE_PATH = std::string(argv[i+1]);
					got_genotype_file = true;
					i++;
				} else if (strcmp(argv[i], "-o") == 0) {
					command_line_opts.OUTPUT_PATH = std::string(argv[i+1]);
					i++;
				} else if (strcmp(argv[i], "-k") == 0) {
					command_line_opts.num_of_evec = atoi(argv[i+1]);
					i++;
				} else if (strcmp(argv[i], "-m") == 0) {
					command_line_opts.max_iterations = atoi(argv[i+1]);
					got_max_iter = true;
					i++;
				} else if (strcmp(argv[i], "-nt") == 0) {
					command_line_opts.nthreads = atoi(argv[i+1]);
					i++;
				} else if (strcmp(argv[i], "-seed") == 0) {
					command_line_opts.seed = atoi(argv[i+1]);
					command_line_opts.given_seed = command_line_opts.seed >= 0 ? true: false;
					i++;
				} else if (strcmp(argv[i], "-l") == 0) {
					command_line_opts.l = atoi(argv[i+1]);
					i++;
				} else if (strcmp(argv[i], "-cl") == 0) {
					command_line_opts.convergence_limit = atof(argv[i+1]);
					i++;
				} else if (strcmp(argv[i], "-aem") == 0) {
					command_line_opts.accelerated_em = atof(argv[i+1]);
					i++;
				} else if (strcmp(argv[i], "-v") == 0) {
					command_line_opts.debugmode = true;
				} else if (strcmp(argv[i], "-vn") == 0) {
					command_line_opts.var_normalize = true;
				} else if (strcmp(argv[i], "-a") == 0) {
					command_line_opts.getaccuracy = true;
				} else if (strcmp(argv[i], "-mem") == 0) {
					command_line_opts.memory_efficient = true;
				} else if (strcmp(argv[i], "-miss") == 0) {
					command_line_opts.missing = true;
				} else if (strcmp(argv[i], "-nfm") == 0) {
					command_line_opts.fast_mode = false;
				} else if (strcmp(argv[i], "-txt") == 0) {
					command_line_opts.text_version = true;
				} else {
					std::cout << "Not Enough or Invalid arguments" << std::endl;
					printCorrectUsage();
					exit(-1);
				}
			} else if (strcmp(argv[i], "-v") == 0) {
				command_line_opts.debugmode = true;
			} else if (strcmp(argv[i], "-a") == 0) {
				command_line_opts.getaccuracy = true;
			} else if (strcmp(argv[i], "-vn") == 0) {
				command_line_opts.var_normalize = true;
			} else if (strcmp(argv[i], "-mem") == 0) {
				command_line_opts.memory_efficient = true;
			} else if (strcmp(argv[i], "-nfm") == 0) {
				command_line_opts.fast_mode = false;
			} else if (strcmp(argv[i], "-miss") == 0) {
				command_line_opts.missing = true;
			} else if (strcmp(argv[i], "-txt") == 0) {
				command_line_opts.text_version = true;
			}
		}
		if (!got_max_iter)
			command_line_opts.max_iterations = command_line_opts.num_of_evec + 2;
	}

	if (got_genotype_file == false) {
		std::cout << "Genotype file missing" << std::endl;
		printCorrectUsage();
		exit(-1);
	}
}


void ProPCA::printCorrectUsage(void) {
	std::cout << "Correct Usage: "
			  << "run_propca \\\n"
			  << "    -g <genotype file> \\\n"
			  << "    -k <number of eigenvectors> \\\n"
			  << "    -m <maximum number of iterations> \\\n"
			  << "    -v (for debug mode) \\\n"
			  << "    -a (for getting accuracy)\n"
			  << std::endl;
}


std::pair<double, double> ProPCA::get_error_norm(MatrixXdr &c) {
	Eigen::HouseholderQR<MatrixXdr> qr(c);
	MatrixXdr Q;
	Q = qr.householderQ() * MatrixXdr::Identity(p, k);
	MatrixXdr q_t(k, p);
	q_t = Q.transpose();
	MatrixXdr b(k, n);
	// Need this for subtracting the correct mean in case of missing data
	if (missing) {
		mm.multiply_y_post(q_t, k, b, false);
		// Just calculating b from seen data
		MatrixXdr M_temp(k, 1);
		M_temp = q_t * means;
		for (int j = 0; j < n; j++) {
			MatrixXdr M_to_remove(k, 1);
			M_to_remove = MatrixXdr::Zero(k, 1);
			for (int i = 0; i < g.not_O_j[j].size(); i++) {
				int idx = g.not_O_j[j][i];
				M_to_remove = M_to_remove + (Q.row(idx).transpose() * g.get_col_mean(idx));
			}
			b.col(j) -= (M_temp - M_to_remove);
		}
	} else {
		mm.multiply_y_post(q_t, k, b, true);
	}

	Eigen::JacobiSVD<MatrixXdr> b_svd(b, Eigen::ComputeThinU | Eigen::ComputeThinV);
	MatrixXdr u_l, d_l, v_l;
	if (fast_mode) {
        u_l = b_svd.matrixU();
	} else {
        u_l = Q * b_svd.matrixU();
    }
	v_l = b_svd.matrixV();
	d_l = MatrixXdr::Zero(k, k);
	for (int kk = 0; kk < k; kk++) {
		d_l(kk, kk) = (b_svd.singularValues())(kk);
	}

	MatrixXdr u_k, v_k, d_k;
	u_k = u_l.leftCols(k_orig);
	v_k = v_l.leftCols(k_orig);
	d_k = MatrixXdr::Zero(k_orig, k_orig);
	for(int kk = 0; kk < k_orig; kk++)
		d_k(kk, kk)  =(b_svd.singularValues())(kk);

	MatrixXdr b_l, b_k;
    b_l = u_l * d_l * (v_l.transpose());
    b_k = u_k * d_k * (v_k.transpose());

    if (fast_mode) {
        double temp_k = b_k.cwiseProduct(b).sum();
        double temp_l = b_l.cwiseProduct(b).sum();
        double b_knorm = b_k.norm();
        double b_lnorm = b_l.norm();
        double norm_k = (b_knorm * b_knorm) - (2 * temp_k);
        double norm_l = (b_lnorm * b_lnorm) - (2 * temp_l);
        return std::make_pair(norm_k, norm_l);
    } else {
        MatrixXdr e_l(p, n);
        MatrixXdr e_k(p, n);
        for (int p_iter = 0; p_iter < p; p_iter++) {
            for (int n_iter = 0; n_iter < n; n_iter++) {
                e_l(p_iter, n_iter) = g.get_geno(p_iter, n_iter, var_normalize) - b_l(p_iter, n_iter);
                e_k(p_iter, n_iter) = g.get_geno(p_iter, n_iter, var_normalize) - b_k(p_iter, n_iter);
            }
        }

        double ek_norm = e_k.norm();
        double el_norm = e_l.norm();
        return std::make_pair(ek_norm, el_norm);
    }
}


/* Run one iteration of EM when genotypes are not missing
 * c_orig : p X k matrix
 * Output: c_new : p X k matrix 
 */
MatrixXdr ProPCA::run_EM_not_missing(MatrixXdr &c_orig) {
	#if DEBUG == 1
		if (debug) {
			print_time();
			std::cout << "Enter: run_EM_not_missing" << std::endl;
		}
	#endif

	// c_temp : k X p matrix: (C^T C)^{-1} C^{T}
	MatrixXdr c_temp(k, p);
	MatrixXdr c_new(p, k);
	c_temp = ((c_orig.transpose()*c_orig).inverse()) * (c_orig.transpose());

	#if DEBUG == 1
		if (debug) {
			print_timenl();
		}
	#endif

	/* E-step: Compute X = D Y 
	* Y : p X n genotype matrix
 	* D : k X p matrix: (C^T C)^{-1} C^{T}
 	* X : k X n matrix
 	*  x_fn: X
 	*  c_temp: D 
 	*/
	MatrixXdr x_fn(k, n);
	mm.multiply_y_post(c_temp, k, x_fn, true);

	#if DEBUG == 1
		if (debug) {
			print_timenl();
		}
	#endif

	// x_temp : n X k matrix X^{T} (XX^{T})^{-1}
	MatrixXdr x_temp(n, k);
	x_temp = (x_fn.transpose()) * ((x_fn * (x_fn.transpose())).inverse());

	/* M-step: C = Y E
	 * Y : p X n genotype matrix
	 * E : n K k matrix: X^{T} (XX^{T})^{-1}
	 * C = p X k matrix
	 * c_new : C 
	 * x_temp : E 
	 */
	mm.multiply_y_pre(x_temp, k, c_new, true);

	#if DEBUG == 1
		if (debug) {
			print_time();
			std::cout << "Exiting: run_EM_not_missing" << std::endl;
		}
	#endif

	return c_new;
}


MatrixXdr ProPCA::run_EM_missing(MatrixXdr &c_orig) {
	MatrixXdr c_new(p, k);
	MatrixXdr mu(k, n);

	// E step
	MatrixXdr c_temp(k, k);
	c_temp = c_orig.transpose() * c_orig;

	MatrixXdr T(k, n);
	MatrixXdr c_fn;
	c_fn = c_orig.transpose();
	mm.multiply_y_post(c_fn, k, T, false);

	MatrixXdr M_temp(k, 1);
	M_temp = c_orig.transpose() * means;

	for (int j = 0; j < n; j++) {
		MatrixXdr D(k, k);
		MatrixXdr M_to_remove(k, 1);
		D = MatrixXdr::Zero(k, k);
		M_to_remove = MatrixXdr::Zero(k, 1);
		for (int i = 0; i < g.not_O_j[j].size(); i++) {
			int idx = g.not_O_j[j][i];
			D = D + (c_orig.row(idx).transpose() * c_orig.row(idx));
			M_to_remove = M_to_remove + (c_orig.row(idx).transpose() * g.get_col_mean(idx));
		}
		mu.col(j) = (c_temp-D).inverse() * (T.col(j) - M_temp + M_to_remove);
	}

	#if DEBUG == 1
		if (debug) {
			std::ofstream x_file;
			x_file.open((output_path + std::string("x_in_fn_vals.txt")).c_str());
			x_file << std::setprecision(15) << mu << std::endl;
			x_file.close();
		}
	#endif


	// M step

	MatrixXdr mu_temp(k, k);
	mu_temp = mu * mu.transpose();
	MatrixXdr T1(p, k);
	MatrixXdr mu_fn;
	mu_fn = mu.transpose();
	mm.multiply_y_pre(mu_fn, k, T1, false);
	MatrixXdr mu_sum(k, 1);
	mu_sum = MatrixXdr::Zero(k, 1);
	mu_sum = mu.rowwise().sum();

	for (int i = 0; i < p; i++) {
		MatrixXdr D(k, k);
		MatrixXdr mu_to_remove(k, 1);
		D = MatrixXdr::Zero(k, k);
		mu_to_remove = MatrixXdr::Zero(k, 1);
		for (int j = 0; j < g.not_O_i[i].size(); j++) {
			int idx = g.not_O_i[i][j];
			D = D + (mu.col(idx) * mu.col(idx).transpose());
			mu_to_remove = mu_to_remove + (mu.col(idx));
		}
		c_new.row(i) = (((mu_temp-D).inverse()) * (T1.row(i).transpose() - (g.get_col_mean(i) * (mu_sum-mu_to_remove)))).transpose();
		double mean;
		mean = g.get_col_sum(i);
		mean = mean -  (c_orig.row(i) * (mu_sum-mu_to_remove))(0, 0);
		mean = mean * 1.0 / (n - g.not_O_i[i].size());
		g.update_col_mean(i, mean);
	}

	// IMPORTANT: Update the value of means variable present locally, so that for next iteration, updated value of means is used.
	for (int i = 0; i < p; i++) {
		means(i, 0) = g.get_col_mean(i);
		// Also updating std, just for consistency, though, it is not used presently.
		stds(i, 0) = g.get_col_std(i);
	}

	return c_new;
}

MatrixXdr ProPCA::run_EM(MatrixXdr &c_orig) {
	if (missing) {
		return run_EM_missing(c_orig);
	} else {
		return run_EM_not_missing(c_orig);
	}
}

void ProPCA::print_vals() {
	Eigen::HouseholderQR<MatrixXdr> qr(c);
	MatrixXdr Q;
	Q = qr.householderQ() * MatrixXdr::Identity(p, k);
	MatrixXdr q_t(k, p);
	q_t = Q.transpose();
	MatrixXdr b(k, n);

	// Need this for subtracting the correct mean in case of missing data
	if (missing) {
		mm.multiply_y_post(q_t, k, b, false);
		// Just calculating b from seen data
		MatrixXdr M_temp(k, 1);
		M_temp = q_t * means;
		for (int j = 0; j < n; j++) {
			MatrixXdr M_to_remove(k, 1);
			M_to_remove = MatrixXdr::Zero(k, 1);
			for (int i = 0; i < g.not_O_j[j].size(); i++) {
				int idx = g.not_O_j[j][i];
				M_to_remove = M_to_remove + (Q.row(idx).transpose()*g.get_col_mean(idx));
			}
			b.col(j) -= (M_temp - M_to_remove);
		}
	} else {
		mm.multiply_y_post(q_t, k, b, true);
	}

	Eigen::JacobiSVD<MatrixXdr> b_svd(b, Eigen::ComputeThinU | Eigen::ComputeThinV);
	MatrixXdr u_l, v_l, u_k, v_k, d_k;
	u_l = b_svd.matrixU();
	v_l = b_svd.matrixV();
	u_k = u_l.leftCols(k_orig);
	v_k = v_l.leftCols(k_orig);

	std::ofstream evec_file;
	evec_file.open((output_path + std::string("evecs.txt")).c_str());
	evec_file << std::setprecision(15) << Q*u_k << std::endl;
	evec_file.close();
	std::ofstream eval_file;
	eval_file.open((output_path + std::string("evals.txt")).c_str());
	for(int kk = 0; kk < k_orig; kk++)
		eval_file << std::setprecision(15) << (b_svd.singularValues())(kk)*(b_svd.singularValues())(kk) / g.Nsnp << std::endl;
	eval_file.close();

	std::ofstream proj_file;
	proj_file.open((output_path + std::string("projections.txt")).c_str());
	proj_file << std::setprecision(15) << v_k << std::endl;
	proj_file.close();
	if (debug) {
		std::ofstream c_file;
		c_file.open((output_path + std::string("cvals.txt")).c_str());
		c_file << std::setprecision(15) << c << std::endl;
		c_file.close();

		std::ofstream means_file;
		means_file.open((output_path + std::string("means.txt")).c_str());
		means_file << std::setprecision(15) << means << std::endl;
		means_file.close();

		d_k = MatrixXdr::Zero(k_orig, k_orig);
		for(int kk = 0; kk < k_orig; kk++)
			d_k(kk, kk)  = (b_svd.singularValues())(kk);
		MatrixXdr x_k;
		x_k = d_k * (v_k.transpose());
		std::ofstream x_file;
		x_file.open((output_path + std::string("xvals.txt")).c_str());
		x_file << std::setprecision(15) << x_k.transpose() << std::endl;
		x_file.close();
	}
}


int ProPCA::run() {
	auto start = std::chrono::system_clock::now();

	clock_t io_begin = clock();
    clock_gettime(CLOCK_REALTIME, &t0);
    //clock_gettime(CLOCK_REALTIME, &timespec);

	std::pair<double, double> prev_error = std::make_pair(0.0, 0.0);
	double prevnll = 0.0;

	// TODO: Memory Effecient Version of Mailman

	memory_efficient = command_line_opts.memory_efficient;
	text_version = command_line_opts.text_version;
    fast_mode = command_line_opts.fast_mode;
	missing = command_line_opts.missing;
    MAX_ITER =  command_line_opts.max_iterations;
	k_orig = command_line_opts.num_of_evec;
	debug = command_line_opts.debugmode;
	check_accuracy = command_line_opts.getaccuracy;
	var_normalize = command_line_opts.var_normalize;
	accelerated_em = command_line_opts.accelerated_em;
	nthreads = command_line_opts.nthreads;
	output_path = std::string(command_line_opts.OUTPUT_PATH);

	if (text_version) {
		if (fast_mode) {
			g.read_txt_mailman(command_line_opts.GENOTYPE_FILE_PATH, missing);
		} else {
			g.read_txt_naive(command_line_opts.GENOTYPE_FILE_PATH, missing);
		}
	} else {
		g.read_plink(command_line_opts.GENOTYPE_FILE_PATH, missing, fast_mode);
	}

	// TODO: Implement these codes.
	if (missing && !fast_mode) {
		std::cout << "Missing version works only with mailman i.e. fast mode"
				  << std::endl
				  <<" EXITING..."
				  << std::endl;
		exit(-1);
	}
	if (fast_mode && memory_efficient) {
		std::cout << "Memory effecient version for mailman EM not yet implemented" << std::endl;
		std::cout << "Ignoring Memory effecient Flag" << std::endl;
	}
	if (missing && var_normalize) {
		std::cout << "Missing version works only without variance normalization"
				  << std::endl
				  << "EXITING..."
				  << std::endl;
		exit(-1);
	}

	k = k_orig + command_line_opts.l;
	k = static_cast<int>(ceil(k / 10.0) * 10);
	command_line_opts.l = k - k_orig;
	p = g.Nsnp;
	n = g.Nindv;
	convergence_limit = command_line_opts.convergence_limit;

	bool toStop = false;
	if(convergence_limit != -1) {
		toStop = true;
	}

	if (command_line_opts.given_seed) {
		srand(command_line_opts.seed);
	} else {
		srand((unsigned int) time(0));
	}

	c.resize(p, k);
	x.resize(k, n);
	v.resize(p, k);
	means.resize(p, 1);
	stds.resize(p, 1);

	if (!fast_mode && !memory_efficient) {
		geno_matrix.resize(p, n);
		g.generate_eigen_geno(geno_matrix, true, var_normalize);
	}

	clock_t io_end = clock();

	// TODO: Initialization of c with gaussian distribution
	c = MatrixXdr::Random(p, k);

	mm = MatMult(g, geno_matrix, debug, var_normalize, memory_efficient,
				 missing, fast_mode, nthreads, k);

	for (int i = 0; i < p; i++) {
		means(i, 0) = g.get_col_mean(i);
		stds(i, 0) = g.get_col_std(i);
	}

	std::ofstream c_file;
	if (debug) {
		c_file.open((output_path + std::string("cvals_orig.txt")).c_str());
		c_file << std::setprecision(15) << c << std::endl;
		c_file.close();
		std::cout << "Read Matrix" << std::endl;
	}

	std::cout << "Running on Dataset of "
			  << g.Nsnp << " SNPs and "
			  << g.Nindv << " Individuals"
			  << std::endl;

	#if SSE_SUPPORT == 1
		if (fast_mode)
			std::cout << "Using Optimized SSE FastMultiply" << std::endl;
	#endif

	clock_t it_begin = clock();
	for (int i = 0; i < MAX_ITER; i++) {
		MatrixXdr c1, c2, cint, r, v;
		double a, nll;
		if (debug) {
			print_time();
			std::cout << "*********** Begin epoch " << i << "***********" << std::endl;
		}
		if (accelerated_em != 0) {
			#if DEBUG == 1
				if (debug) {
					print_time();
					std::cout << "Before EM" << std::endl;
				}
			#endif
			c1 = run_EM(c);
			c2 = run_EM(c1);
			#if DEBUG == 1
				if (debug) {
					print_time();
					std::cout << "After EM but before acceleration" << std::endl;
				}
			#endif
			r = c1 - c;
			v = (c2 - c1) - r;
			a = -1.0 * r.norm() / (v.norm());
			if (accelerated_em == 1) {
				if (a > -1) {
					a = -1;
					cint = c2;
				} else {
					cint = c - 2*a*r + a*a*v;
					nll = get_error_norm(cint).second;
					if (i > 0) {
						while (nll > prevnll && a < -1) {
							a = 0.5 * (a-1);
							cint = c - 2*a*r +(a*a*v);
							nll = get_error_norm(cint).second;
						}
					}
				}
				c = cint;
			} else if (accelerated_em == 2) {
				cint = c - 2*a*r + a*a*v;
				c = cint;
				// c = run_EM(cint);
			}
		} else {
			c = run_EM(c);
		}

		if (accelerated_em == 1 || check_accuracy || toStop) {
			std::pair<double, double> e = get_error_norm(c);
			prevnll = e.second;
			if (check_accuracy) {
				std::cout << "Iteration " << i + 1 << "  " << std::setprecision(15) << e.first << "  " << e.second << std::endl;
			}

			if (std::abs(e.first-prev_error.first) <= convergence_limit) {
				std::cout << "Breaking after " << i + 1 << " iterations" << std::endl;
				break;
			}
			prev_error = e;
		}
		if (debug) {
			print_time();
			std::cout << "*********** End epoch " << i << "***********" << std::endl;
		}
	}
	clock_t it_end = clock();

    print_vals();

    mm.clean_up();

	clock_t total_end = clock();
	double io_time = static_cast<double>(io_end - io_begin) / CLOCKS_PER_SEC;
	double avg_it_time = static_cast<double>(it_end - it_begin) / (MAX_ITER * 1.0 * CLOCKS_PER_SEC);
	double total_time = static_cast<double>(total_end - total_begin) / CLOCKS_PER_SEC;
	std::cout << "IO Time:  " << io_time << "\nAVG Iteration Time:  " << avg_it_time << "\nTotal runtime:   " << total_time << std::endl;

	std::chrono::duration<double> wctduration = std::chrono::system_clock::now() - start;
	std::cout << "Wall clock time = " <<  wctduration.count() << std::endl;

	return 0;
}
