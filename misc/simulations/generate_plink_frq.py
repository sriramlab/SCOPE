#!/bin/env python3

import sys

if __name__ == '__main__':

	if len(sys.argv) != 4:
		print("Usage: python3 generate_plink_frq.py plink_bim true_freqs outfile",file=sys.stderr)
		exit(1)
	
	outfile = open(sys.argv[3],'w+')
	bim_file = open(sys.argv[1],'r+')
	freqs_file = open(sys.argv[2],'r+')

	bim = bim_file.readlines()
	freqs = freqs_file.readlines()

	bim_file.close()
	freqs_file.close()

	outfile.write("CHR\tSNP\tCLST\tA1\tA2\tMAF\tMAC\tNCHROBS\n")

	for i in range(len(bim)):
		line = bim[i].rstrip().strip().split()
		CHR = line[0]
		SNP = line[1]
		A1 = line[4]
		A2 = line[5]
		snp_freq = freqs[i].rstrip().strip().split()
		for j in range(len(snp_freq)):
			outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(CHR,SNP,j+1,A1,A2,snp_freq[j],"N","N"))

	outfile.close()

	exit(0)
