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

#include <Spectra/SymEigsSolver.h>

#include "alstructure.h"
#include "config.h"
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

template<typename M>
M read_plink_freq_file (const std::string &path, int k) {
	/*
	Reads from plink.frq.strat file to initialize P matrix
	Returns rows/k x k matrix, where rows in number of rows
		in file
	*/
        std::ifstream indata;
        indata.open(path);
        std::string line;
        std::vector<double> values;
        int rows = 1;
        std::getline(indata, line);
        while (std::getline(indata, line)) {
                std::stringstream lineStream(line);
                std::string cell;
                std::vector<std::string> seglist;
                while (lineStream >> cell){
                        seglist.push_back(cell);
                }
                values.push_back(std::stod(seglist[5]));
                ++rows;
        }
        return Eigen::Map<const Eigen::Matrix<typename M::Scalar, M::RowsAtCompileTime, M::ColsAtCompileTime, Eigen::RowMajor>>(values.data(), rows/k, k);
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

	std::stable_sort(inds.begin(), inds.end(), [&data](size_t i, size_t j) { return data[i] > data[j]; });

	double tmpsum = 0;
	double tmax;
	bool bget = false;

	for (int i = 1; i < data.size(); i++) {
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
	bool got_freq_file = false;
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
	nops = 1;


	if (argc < 3) {
		printCorrectUsage();
		//std::cout << "Correct Usage is " << argv[0] << " -p <parameter file>" << std::endl;
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
				} else if (strcmp(argv[i], "-freq") == 0) {
					command_line_opts.FREQ_FILE_PATH = std::string(argv[i+1]);
					got_freq_file = true;
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
			  << "    -k <latent dimension> \\\n"
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


void ALStructure::initialize(std::default_random_engine &prng_eng) {
	if (!command_line_opts.given_seed) {
		seed = static_cast<unsigned int>(time(NULL));
	}

	std::cout << "Initializing Phat using seed " << seed << std::endl;

	// srand(seed);
	prng_eng.seed(seed);

	// Phat = MatrixXdr::Random(p, k);
	// Phat = Phat.unaryExpr(&fix_interval);
	std::uniform_real_distribution<double> dis(0, 1);
	Phat = MatrixXdr::Zero(p, k).unaryExpr([&](float dummy){return dis(prng_eng);});

	if (debug) write_matrix(Phat, std::string("Phat_0_") + std::to_string(seed) +  std::string(".txt"));
}


void ALStructure::truncated_alternating_least_squares(bool projection_mode) {

	if (projection_mode){
		std::cout << "Solving for Q using provided frequencies" << std::endl;
		solve_for_Qhat();
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
		return;
	}

	solve_for_Qhat();
	if (debug) write_matrix(Qhat, std::string("Qhat_0.txt"));

	solve_for_Phat();
	if (debug) write_matrix(Phat, std::string("Phat_1.txt"));

	for (niter = 1; niter < MAX_ITER; niter++) {
		Qhat_old = Qhat;

		solve_for_Qhat();
		if (debug) write_matrix(Qhat, std::string("Qhat_" + std::to_string(niter) + ".txt"));

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
		if (debug) write_matrix(Qhat, std::string("Qhat_" + std::to_string(niter) + "_w_constraints.txt"));

		solve_for_Phat();
		if (debug) write_matrix(Phat, std::string("Phat_" + std::to_string(niter+1) + ".txt"));

		Phat = Phat.unaryExpr(&truncate_with_epsilon);
		if (debug) write_matrix(Phat, std::string("Phat_" + std::to_string(niter+1) + "_w_constraints.txt"));

		diff = Qhat - Qhat_old;
 		rmse = diff.norm() / sqrt(n * k);
		std::cout << "Iteration " << niter+1 << "  -- RMSE " << std::setprecision(15) << rmse << std::endl;
		if ((rmse <= convergence_limit) || std::isnan(rmse)) {
			std::cout << "Breaking after " << niter+1 << " iterations" << std::endl;
			break;
		}
	}
}


void ALStructure::write_matrix(MatrixXdr &mat, const std::string file_name) {
	std::ofstream fp;
	fp.open((command_line_opts.OUTPUT_PATH + file_name).c_str());
	fp << std::setprecision(15) << mat << std::endl;
	fp.close();
}

void ALStructure::write_vector(Eigen::VectorXd &vec, const std::string file_name) {
	std::ofstream fp;
	fp.open((command_line_opts.OUTPUT_PATH + file_name).c_str());
	fp << std::setprecision(15) << vec << std::endl;
	fp.close();
}

int ALStructure::run() {

	total_begin = clock();

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
	k = command_line_opts.num_of_evec;

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
		if (debug) write_matrix(geno_matrix, "X.txt");
	}

	clock_t io_end = clock();

	std::cout << "Running on Dataset of " << p << " SNPs and " << n << " Individuals" << std::endl;

	#if SSE_SUPPORT == 1
		if (fast_mode)
			std::cout << "Using Optimized SSE FastMultiply" << std::endl;
	#endif

	clock_t it_begin = clock();

	mm = MatMult(g, geno_matrix, debug, false, memory_efficient, missing, fast_mode, nthreads, k);

	if (std::string(command_line_opts.ROWSPACE_FILE_PATH) != "") {
		std::cout << "Using provided V" << std::endl;

	// Read eigenvectors of the n x n matrix: G = (1/m) * (X^T X - D)
		V = load_tsv<MatrixXdr>(command_line_opts.ROWSPACE_FILE_PATH);
		if (k != V.cols()) {
			k = V.cols();
			std::cout << "Mismatch between column number of provided V and provided k!" << std::endl;
			std::cout << "Changing k to number of columns in V" << std::endl;
		}

		if (V.rows() != n) {
			std::cout << "Dimensions of genotype matrix and rowspace matrix do not agree!" << std::endl;
			exit(-1);
		}
	} else {
		std::cout << "Performing latent subspace estimation" << std::endl;

		// Calculate D matrix
		D.resize(g.Nindv);
		for (int i = 0; i < g.Nindv; ++i) {
			D[i] =  2 * g.rowsum[i] - g.rowsqsum[i];
		}
		if (debug) write_vector(D, "D.txt");

		// Calculate V
		Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE, ALStructure> eigs(this, k, k * 2 + 1);
		eigs.init();
		eigs.compute(MAX_ITER, convergence_limit);

		if (eigs.info() == Spectra::SUCCESSFUL) {
			V = eigs.eigenvectors();
			write_matrix(V, "V.txt");
			if (debug) {
				Eigen::VectorXd evals = eigs.eigenvalues().array() / (n-1);
				write_vector(evals, "evals.txt");
			}
			std::cout << "Latent subspace esimation completed after " << nops << " iterations" << std::endl;
		}
		else {
			throw new std::runtime_error(
				std::string("Spectra eigendecomposition unsucessful") + ", status" + std::to_string(eigs.info())
			);
		}
	}

	// Compute Fhat = (1/2) * (X V V^T)
	if (fhat_version || fhattrunc_version) {
		std::cout << "Explicitly computing Fhat";
		MatrixXdr temp_pxk(p, k);
		mm.multiply_y_pre(V, k, temp_pxk, false);

		// if (debug) write_matrix(temp_pxk, "XV.txt");

		Fhat = temp_pxk * V.transpose();
		Fhat = Fhat.unaryExpr(&divide_by_two);
		if (debug) write_matrix(Fhat, "Fhat.txt");

		if (command_line_opts.fhattrunc_version) {
			std::cout << " -- truncated";
			Fhat = Fhat.unaryExpr(&truncate_xxx);
			if (debug) write_matrix(Fhat, "Fhat_truncated.txt");
		}
		std::cout << std::endl;
	}

	// Create pseudo random number generator (PRNG) engine
	std::default_random_engine prng_eng{};

	if (std::string(command_line_opts.INITIAL_FILE_PATH) != "") {
		std::cout << "Using initial Phat provided" << std::endl; 
		Phat = load_tsv<MatrixXdr>(command_line_opts.INITIAL_FILE_PATH);
	} 
	else if (std::string(command_line_opts.FREQ_FILE_PATH) != ""){
		Phat = read_plink_freq_file<MatrixXdr>(command_line_opts.FREQ_FILE_PATH,k);
		MAX_ITER = 1;
	}
	else {
		initialize(prng_eng);
	}

	if (std::string(command_line_opts.FREQ_FILE_PATH) != ""){
		truncated_alternating_least_squares(true);
	}
	else {
		truncated_alternating_least_squares();

		// Try restarting with new seed!
		if (std::isnan(rmse) && (std::string(command_line_opts.INITIAL_FILE_PATH) == "")) {
			command_line_opts.given_seed = false;
			for (int xx = 0; xx < 5; xx++) {
				initialize(prng_eng);
				truncated_alternating_least_squares();
				if (!std::isnan(rmse)) {
					break;
				}
			}
		}
	}

	clock_t it_end = clock();

	write_matrix(Phat, "Phat.txt");
	write_matrix(Qhat, "Qhat.txt");

	mm.clean_up();

	clock_t total_end = clock();
	double io_time = static_cast<double>(io_end - io_begin) / CLOCKS_PER_SEC;
	//double avg_it_time = static_cast<double>(it_end - it_begin) / (MAX_ITER * 1.0 * CLOCKS_PER_SEC);
	double total_time = static_cast<double>(total_end - total_begin) / CLOCKS_PER_SEC;
	std::cout << "Completed!" << std::endl;
	std::cout << "IO Time:  " << io_time << std::endl;
	std::cout << "Total runtime:   " << total_time << std::endl;

	std::chrono::duration<double> wctduration = std::chrono::system_clock::now() - start;
	std::cout << "Wall clock time = " <<  wctduration.count() << std::endl;

	return 0;
}

unsigned int ALStructure::cols() {
	return n;
}

unsigned int ALStructure::rows() {
	return n;
}

void ALStructure::perform_op(const double* x_in, double* y_out) {
	// Performs ((Xv)^T X)^T - Dv
	MatrixXdr x = Eigen::Map<const Eigen::VectorXd> (x_in, n);
	Eigen::Map<Eigen::VectorXd> y(y_out, n);

	MatrixXdr temp_px1(p,1);
	mm.multiply_y_pre(x,1,temp_px1,false);
	temp_px1.transposeInPlace(); // 1xp

	MatrixXdr temp_1xn(1,n);
	mm.multiply_y_post(temp_px1,1,temp_1xn,false);

	y.noalias() = temp_1xn.transpose() - D.cwiseProduct(x);
	nops++;
}
