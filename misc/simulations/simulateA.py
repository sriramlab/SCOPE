#!/bin/env python3

import sys
import numpy as np
import pyplink
import argparse
import parmapper

def draw_allele_freqs(i,k,maf,fst):
	"""
	Draws allele frequencies using
	BD model

	Parameters
	i: index to use
	k: number of draws/populations
	maf: MAF vector from real data
	fst: FST level from real data

	returns: length-k vector of allele freqs
	"""
	freq = maf[i]
	fst_val = fst[i]
	if fst_val <= 0 or fst_val == 1 or freq == 0 or freq == 1:
		return(np.repeat(0,k))
	else:
		alpha = ((1-fst_val)/fst_val)*freq
		beta = ((1-fst_val)/fst_val)*(1-freq)
		return(np.random.beta(alpha,beta,k))

def seeded_binom(i,p):
	'''
	Seeded binomial function for parallelization

	paramters
	i: index
	p: probability matrix

	returns: binomial draws
	'''
	np.random.seed(i)
	return(np.random.binomial(2,p[i]))

if __name__=='__main__':

	parser = argparse.ArgumentParser(description="Genotype simulation under BN-PSD model (Simulation A)")

	#parser = argparse.ArgumentParser(description="Genotype simulation under BN-PSD model (Simulation A)",epilog="Performs simulations under population genetics models. Allele frequencies are drawn from the Balding-Nichols models using parameters provided in the freq_fst_file. The freq_fst_file is assumed to come from real data and generated from PLINK using the --fst flag. Population proportions for each individual are generated from the Pritchard-Stephens-Donnelly model using a two-stage Dirichlet draw. We assume S regions with K latent populations. Each individual is drawn from a region's population vector. This corresponds to Scenario A in the TeraStructure paper (2016).")

	parser.add_argument("n_ind",help="number of individuals",type=int)
	parser.add_argument("n_snps",help="number of SNPs",type=int)
	parser.add_argument("freq_fst_file",help="MAF and FST file (PLINK --fst, --freq) file",type=str)
	parser.add_argument("outprefix",help="prefix for output files")
	parser.add_argument("-S","--regions",help="number of regions (default: 50)",type=int,default=50)
	parser.add_argument("-K","--latent",help="number of latent populations (default: 6)",type=int,default=6)
	parser.add_argument("-a","--alpha",help="initial dirichlet parameter (default: 0.2)",type=float,default=0.2)
	parser.add_argument("-g","--gamma",help="second-level dirichlet scale (default: 50)",type=float,default=50)
	parser.add_argument("-c","--chunk",help="chunk size (default: 10000)",type=int,default=10000)
	parser.add_argument("-s","--seed",help="seed (default: 1)",type=int,default=1)
	parser.add_argument("-np","--processes",help="number of processes to spawn (default: 1)",type=int,default=1)
	args = parser.parse_args()

	S = args.regions
	K = args.latent
	alpha = args.alpha
	gamma = args.gamma
	n_ind = args.n_ind
	n_snps = args.n_snps
	freq_fst_file = args.freq_fst_file
	chunk_size = args.chunk
	outprefix = args.outprefix
	seed = args.seed
	n_cpu = args.processes

	print("Simulation A")
	print("-------------")
	print("Parameters")
	print("-------------")
	print("S: {}\nK: {}\nalpha: {}\ngamma: {}\nn_ind: {}\nn_snps: {}\nfreq_fst_file: {}\nchunk_size: {}\noutprefix: {}\nseed: {}\nprocesses: {}".format(S,K,alpha,gamma,n_ind,n_snps,freq_fst_file,chunk_size,outprefix,seed,n_cpu))
	print("-------------")
	
	np.random.seed(seed)
	
	print("Generating individual proportions...")
	n_region = int(n_ind/S)
	q = np.random.dirichlet(np.repeat(alpha,K),size=S)

	pop_thetas = np.apply_along_axis(np.random.dirichlet,1,gamma*q,n_region)
	pop_thetas = pop_thetas.reshape(pop_thetas.shape[0]*pop_thetas.shape[1],pop_thetas.shape[2])

	np.savetxt("{}_proportions.txt".format(outprefix),pop_thetas,delimiter='\t')	

	print("Successfully generated individual proportions")

	print("Drawing allele frequencies...")

	print("Using information from {}".format(freq_fst_file))
	maf, fst = np.genfromtxt(freq_fst_file,skip_header=1,missing_values='NA',filling_values=0,usecols=[4,6],unpack=True)

	snp_ind = np.random.choice(list(range(len(maf))),size=n_snps,replace=True)

	pop_allele_freq = np.array(list(map(lambda x: draw_allele_freqs(x,K,maf,fst),snp_ind)))
	
	del snp_ind
	del maf
	del fst

	np.savetxt("{}_allele_frequencies.txt".format(outprefix),pop_allele_freq,delimiter='\t')

	print("Successfully drew allele frequencies")

	print("Drawing genotypes...")	

	pop_thetas = np.transpose(pop_thetas)

	with pyplink.PyPlink(outprefix,mode="w") as bed:
		for i in range(int(np.ceil(n_snps/chunk_size))):	
			ind_freq = np.matmul(pop_allele_freq[i*chunk_size:(i+1)*chunk_size,:],pop_thetas)
			ind_genos = list(parmapper.parmap(lambda x: seeded_binom(x,ind_freq),range(ind_freq.shape[0]),N=n_cpu))
			#ind_genos = np.random.binomial(2,ind_freq)
			del ind_freq
			for j in range(len(ind_genos)):
				bed.write_genotypes(ind_genos[j])
				#bed.write_genotypes(ind_genos[j,:])
			del ind_genos
			print("{}/{} SNPs written".format(np.minimum((i+1)*chunk_size,n_snps),n_snps))

	print("Finished drawing genotypes")

	bim_file = np.stack((np.repeat(1,n_snps),list(range(1,n_snps+1)),np.repeat(0,n_snps),list(range(1,n_snps+1)),np.repeat(1,n_snps),np.repeat(2,n_snps)),axis=1)

	np.savetxt("{}.bim".format(outprefix),bim_file,fmt='%i',delimiter='\t')

	fam_file = np.stack((list(range(1,n_ind+1)),list(range(1,n_ind+1)),np.repeat(0,n_ind),np.repeat(0,n_ind),np.repeat(0,n_ind),np.repeat(0,n_ind)),axis=1)

	np.savetxt("{}.fam".format(outprefix),fam_file,fmt='%i',delimiter='\t')

	print("PLINK files written to {}.(bed,bim,fam)".format(outprefix))

	print("Population allele frequencies saved to {}_allele_freqs.txt".format(outprefix))

	print("Individual proportions saved to {}_proportions.txt".format(outprefix))

	print("Finished simulation")

	exit(0)

