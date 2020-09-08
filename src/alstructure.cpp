//#include <bits/stdc++.h>

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>
#include "time.h"

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/SVD>
#include <Eigen/QR>

#include "alstructure.h"
#include "config.h"
#include "helper.h"
#include "storage.h"


struct timespec t0;


// https://stackoverflow.com/questions/34247057/how-to-read-csv-file-and-assign-to-eigen-matrix
template<typename M>
M load_tsv (const std::string & path) {
	std::ifstream indata;
	indata.open(path);
	std::string line;
	std::vector<double> values;
	int rows = 0;
	while (std::getline(indata, line)) {
		std::stringstream lineStream(line);
		std::string cell;
		while (std::getline(lineStream, cell, '\t')) {
			values.push_back(std::stod(cell));
		}
		++rows;
	}
	return Eigen::Map<const Eigen::Matrix<typename M::Scalar, M::RowsAtCompileTime, M::ColsAtCompileTime, Eigen::RowMajor>>(values.data(), rows, values.size()/rows);
}


double fix_interval(double x) {
	/* Random function in Eigen generates numbers from a 
	 * uniform distribution on the interval [-1, 1].
	 * We need to map these values to the interval [0, 1].
	 * To map x in [a, b] to [c, d]
	 * f(x) = c + ((d-c)/(b-a)) * (x-a)
	 */
	return (0.5 * (x + 1.0));
}


double divide_by_two(double x) {
	return x / 2.0;
}


double truncate_with_epsilon(double x) {
	double epsilon = 0.0000000001;
	if (x <= 0.0)
		return epsilon;
	if (x >= 1.0)
		return 1.0 - epsilon;
	return x;
}


double truncate_xxx(double x) {
	if (x < 0.0)
		return 0.0;
	if (x > 1.0)
		return 1.0;
	return x;
}


void project_onto_simplex(std::vector<double> &data) {
	std::vector<size_t> inds(data.size());

	std::iota(inds.begin(), inds.end(), 0);

	//std::stable_sort(inds.begin(), inds.end(), [&data](size_t i, size_t j) { return data[i] > data[j]; });

	double tmpsum = 0;
	double tmax;
	bool bget = false;

	for(int i = 1; i < data.size(); i++){
		tmpsum = tmpsum + data[inds[i - 1]];
		tmax = (tmpsum - 1.0) / i;
		if (tmax >= data[inds[i]]) {
			bget = true;
			break;
		}
	}

	if (!(bget)) {
		tmpsum = tmpsum + data[inds[data.size() - 1]];
		tmax = (tmpsum - 1) / data.size();
	}

	for (int i = 0; i < data.size(); i++) {
		data[i] = data[i] - tmax;
		if (data[i] < 0.0) {
			data[i] = 0.0;
		}
	}
}


ALStructure::ALStructure(int argc, char const *argv[]) {
	// Set default values
	command_line_opts.num_of_evec = 5;
	command_line_opts.debugmode = false;
	command_line_opts.OUTPUT_PATH = "alstructure_";
	bool got_genotype_file = false;
	bool got_rowspace_file = false;
	command_line_opts.convergence_limit = 0.00001;  // Used by R code
	command_line_opts.max_iterations = 1000;        // Used by R code
	command_line_opts.memory_efficient = false;
	command_line_opts.fast_mode = true;
	command_line_opts.missing = false;
	command_line_opts.text_version = false;
	command_line_opts.fhat_version = false;
	command_line_opts.fhattrunc_version = false;
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
		command_line_opts.debugmode = cfg.getValueOfKey<bool>("debug", false);
		command_line_opts.OUTPUT_PATH = cfg.getValueOfKey<std::string>("output_path", std::string("fastppca_"));
		command_line_opts.GENOTYPE_FILE_PATH = cfg.getValueOfKey<std::string>("genotype", std::string(""));
		command_line_opts.convergence_limit = cfg.getValueOfKey<double>("convergence_limit", -1.0);
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
				} else if (strcmp(argv[i], "-r") == 0) {
					command_line_opts.ROWSPACE_FILE_PATH = std::string(argv[i+1]);
					got_rowspace_file = true;
					i++;
				} else if (strcmp(argv[i], "-i") == 0) {
					command_line_opts.INITIAL_FILE_PATH = std::string(argv[i+1]);
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
				} else if (strcmp(argv[i], "-cl") == 0) {
					command_line_opts.convergence_limit = atof(argv[i+1]);
					i++;
				} else if (strcmp(argv[i], "-v") == 0) {
					command_line_opts.debugmode = true;
				} else if (strcmp(argv[i], "-mem") == 0) {
					command_line_opts.memory_efficient = true;
				} else if (strcmp(argv[i], "-miss") == 0) {
					command_line_opts.missing = true;
				} else if (strcmp(argv[i], "-nfm") == 0) {
					command_line_opts.fast_mode = false;
				} else if (strcmp(argv[i], "-txt") == 0) {
					command_line_opts.text_version = true;
				} else if (strcmp(argv[i], "-fhat") == 0) {
					command_line_opts.fhat_version = true;
				} else if (strcmp(argv[i], "-fhattrunc") == 0) {
					command_line_opts.fhattrunc_version = true;
				} else {
					std::cout << "Not Enough or Invalid arguments" << std::endl;
					printCorrectUsage();
					exit(-1);
				}
			} else if (strcmp(argv[i], "-v") == 0) {
				command_line_opts.debugmode = true;
			} else if (strcmp(argv[i], "-mem") == 0) {
				command_line_opts.memory_efficient = true;
			} else if (strcmp(argv[i], "-nfm") == 0) {
				command_line_opts.fast_mode = false;
			} else if (strcmp(argv[i], "-miss") == 0) {
				command_line_opts.missing = true;
			} else if (strcmp(argv[i], "-txt") == 0) {
				command_line_opts.text_version = true;
			} else if (strcmp(argv[i], "-fhat") == 0) {
				command_line_opts.fhat_version = true;
			} else if (strcmp(argv[i], "-fhattrunc") == 0) {
				command_line_opts.fhattrunc_version = true;
			}
		}
	}

	if (got_genotype_file == false) {
		std::cout << "Genotype file missing" << std::endl;
		printCorrectUsage();
		exit(-1);
	}
}


void ALStructure::printCorrectUsage(void) {
	std::cout << "Correct Usage: "
			  << "run_alstructure \\\n"
			  << "    -g <genotype file> \\\n"
			  << "    -r <rowspace file> \\\n"
			  << "    -m <maximum number of iterations> \\\n"
			  << "    -v (for debug mode) \\\n"
			  << std::endl;
}


void ALStructure::solve_for_Qhat() {
	if (fhat_version || fhattrunc_version) {
		Qhat = (((Phat.transpose() * Phat).inverse() * Phat.transpose()) * Fhat);
	}

	MatrixXdr temp_kxp(k, p);
	temp_kxp = ((Phat.transpose() * Phat).inverse()) * Phat.transpose();

	// temp_kxn = temp_kxp * X, where X is the p X n genotype matrix
	MatrixXdr temp_kxn(k, n);
	mm.multiply_y_post(temp_kxp, k, temp_kxn, false);

	temp_kxn = temp_kxn.unaryExpr(&divide_by_two);

	Qhat = (temp_kxn * V) * V.transpose();
}


void ALStructure::solve_for_Phat() {
	if (fhat_version || fhattrunc_version) {
		Phat = (Fhat * (Qhat.transpose() * (( Qhat * Qhat.transpose() ).inverse()) ));
	}

	MatrixXdr temp_nxk(n, k);
	temp_nxk = V * (V.transpose() * ((Qhat.transpose() * (Qhat * Qhat.transpose()).inverse())));

	mm.multiply_y_pre(temp_nxk, k, Phat, false);

	Phat = Phat.unaryExpr(&divide_by_two);
}


void ALStructure::initialize() {
	std::ofstream fp;

	if (command_line_opts.given_seed) {
		srand(seed);
	} else {
		seed = (unsigned int) time(0);
		srand(seed);
	}

	std::cout << "Initializing Phat using seed " << seed << std::endl;

	Phat = MatrixXdr::Random(p, k);
	Phat = Phat.unaryExpr(&fix_interval);
	if (debug) {
		fp.open((command_line_opts.OUTPUT_PATH + "Phat_0.txt").c_str());
		fp << std::setprecision(15) << Phat << std::endl;
		fp.close();
	}
}


void ALStructure::truncated_alternating_least_squares() {
	std::ofstream fp;

	solve_for_Qhat();
	if (debug) {
		fp.open((command_line_opts.OUTPUT_PATH + "Qhat_0.txt").c_str());
		fp << std::setprecision(15) << Qhat << std::endl;
		fp.close();
	}

	solve_for_Phat();
	if (debug) {
		fp.open((command_line_opts.OUTPUT_PATH + "Phat_1.txt").c_str());
		fp << std::setprecision(15) << Phat << std::endl;
		fp.close();
	}

	for (niter = 1; niter < MAX_ITER; niter++) {
		Qhat_old = Qhat;

		solve_for_Qhat();
		if (debug) {
			fp.open((command_line_opts.OUTPUT_PATH + "Qhat_" + std::to_string(niter) + ".txt").c_str());
			fp << std::setprecision(15) << Qhat << std::endl;
			fp.close();
		}

		std::vector<double> col;
		col.resize(k);
		for (int c_iter = 0; c_iter < n; c_iter++) {
 			//VectorXd::Map(&col[0], d) = Qhat.col(c_iter);
			for (int r_iter = 0; r_iter < k; r_iter++) {
				col[r_iter] = Qhat(r_iter, c_iter);
			}

			project_onto_simplex(col);

			for (int r_iter = 0; r_iter < k; r_iter++) {
				Qhat(r_iter, c_iter) = col[r_iter];
			}
		}

		if (debug) {
			fp.open((command_line_opts.OUTPUT_PATH + "Qhat_" + std::to_string(niter) + "_w_constraints.txt").c_str());
			fp << std::setprecision(15) << Qhat << std::endl;
			fp.close();
		}

		solve_for_Phat();
		if (debug) {
			fp.open((command_line_opts.OUTPUT_PATH + "Phat_" + std::to_string(niter+1) + ".txt").c_str());
			fp << std::setprecision(15) << Phat << std::endl;
			fp.close();
		}

		Phat = Phat.unaryExpr(&truncate_with_epsilon);
		if (debug) {
			fp.open((command_line_opts.OUTPUT_PATH + "Phat_" + std::to_string(niter+1) + "_w_constraints.txt").c_str());
			fp << std::setprecision(15) << Phat << std::endl;
			fp.close();
		}

		diff = Qhat - Qhat_old;
 		rmse = diff.norm() / sqrt(n * k);
		std::cout << "Iteration " << niter+1 << "  -- RMSE " << std::setprecision(15) << rmse << std::endl;
		if ((rmse <= convergence_limit) || std::isnan(rmse)) {
			std::cout << "Breaking after " << niter+1 << " iterations" << std::endl;
			break;
		}
	}
}


int ALStructure::run() {
	std::ofstream fp;

	memory_efficient = command_line_opts.memory_efficient;
	text_version = command_line_opts.text_version;
	fhat_version = command_line_opts.fhat_version;
	fhattrunc_version = command_line_opts.fhattrunc_version;
	fast_mode = command_line_opts.fast_mode;
	missing = command_line_opts.missing;
	MAX_ITER =  command_line_opts.max_iterations;
	convergence_limit = command_line_opts.convergence_limit;
	debug = command_line_opts.debugmode;
	nthreads = command_line_opts.nthreads;
	output_path = std::string(command_line_opts.OUTPUT_PATH);
	seed = command_line_opts.seed;

	auto start = std::chrono::system_clock::now();

	clock_t io_begin = clock();
	clock_gettime(CLOCK_REALTIME, &t0);

	// Read genotype matrix X
	if (text_version) {
		if (fast_mode) {
			g.read_txt_mailman(command_line_opts.GENOTYPE_FILE_PATH, missing);
		} else {
			g.read_txt_naive(command_line_opts.GENOTYPE_FILE_PATH, missing);
		}
	} else {
		g.read_plink(command_line_opts.GENOTYPE_FILE_PATH, missing, fast_mode);
	}

	p = g.Nsnp;
	n = g.Nindv;

	// TODO: Implement these codes.
	if (missing) {
		std::cout << "Missing version not yet implemented!" << std::endl;
		exit(-1);
	}

	if (fast_mode && memory_efficient) {
		std::cout << "Memory effecient version for mailman EM not yet implemented" << std::endl;
		std::cout << "Ignoring Memory effecient Flag" << std::endl;
	}

	if (!fast_mode && !memory_efficient) {
		geno_matrix.resize(p, n);
		g.generate_eigen_geno(geno_matrix, false, false);
		if (debug) {
			fp.open((command_line_opts.OUTPUT_PATH + "X.txt").c_str());
			fp << std::setprecision(15) << geno_matrix << std::endl;
			fp.close();
		}
	}

	// Read eigenvectors of the n x n matrix: G = (1/m) * (X^T X - D)
	V = load_tsv<MatrixXdr>(command_line_opts.ROWSPACE_FILE_PATH);
	k = V.cols();

	if (V.rows() != n) {
		std::cout << "Dimensions of genotype matrix and rowspace matrix do not agree!" << std::endl;
		exit(-1);
	}

	clock_t io_end = clock();

	mm = MatMult(g, geno_matrix, debug, false, memory_efficient, missing, fast_mode, nthreads, k);

	// Compute Fhat = (1/2) * (X V V^T)
	if (fhat_version || fhattrunc_version) {
		std::cout << "Explicitly computing Fhat";
		MatrixXdr temp_pxk(p, k);
		mm.multiply_y_pre(V, k, temp_pxk, false);

		if (debug) {
			fp.open((command_line_opts.OUTPUT_PATH + "XV.txt").c_str());
			fp << std::setprecision(15) << temp_pxk << std::endl;
			fp.close();
		}

		Fhat = temp_pxk * V.transpose();
		Fhat = Fhat.unaryExpr(&divide_by_two);
		if (debug) {
			fp.open((command_line_opts.OUTPUT_PATH + "Fhat.txt").c_str());
			fp << std::setprecision(15) << Fhat << std::endl;
			fp.close();
		}

		if (command_line_opts.fhattrunc_version) {
			std::cout << " -- truncated";
			Fhat = Fhat.unaryExpr(&truncate_xxx);
			if (debug) {
				fp.open((command_line_opts.OUTPUT_PATH + "Fhat_truncated.txt").c_str());
				fp << std::setprecision(15) << Fhat << std::endl;
				fp.close();
			}
		}
		std::cout << std::endl;
	}

	std::cout << "Running on Dataset of " << p << " SNPs and " << n << " Individuals" << std::endl;

	#if SSE_SUPPORT == 1
		if (fast_mode)
			std::cout << "Using Optimized SSE FastMultiply" << std::endl;
	#endif

	clock_t it_begin = clock();

	if (std::string(command_line_opts.INITIAL_FILE_PATH) != "") {
		std::cout << "Using initial Phat provided" << std::endl; 
		Phat = load_tsv<MatrixXdr>(command_line_opts.INITIAL_FILE_PATH);
	} else {
		initialize();
	}

	truncated_alternating_least_squares();

	// Restart!
	if ((niter < 10) && (std::isnan(rmse))) {
		if (std::string(command_line_opts.INITIAL_FILE_PATH) == "") {
			command_line_opts.given_seed = false;
			initialize();
			truncated_alternating_least_squares();
		}
	}

	clock_t it_end = clock();

	fp.open((command_line_opts.OUTPUT_PATH + "Phat.txt").c_str());
	fp << std::setprecision(15) << Phat << std::endl;
	fp.close();

	fp.open((command_line_opts.OUTPUT_PATH + "Qhat.txt").c_str());
	fp << std::setprecision(15) << Qhat << std::endl;
	fp.close();

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
