#!/bin/env python3

echo Run \"bash install_packages.sh\" once before running simulations
echo ---------------------------

mkdir simulationA
mkdir simulationB

python3 simulateA.py -np 6 -c 100000 10000 1000000 tgp_superpop_param.txt simulationA/n10k_1m_snps 
python3 simulateA.py -np 6 -c 100000 100000 1000000 tgp_superpop_param.txt simulationA/n100k_1m_snps
python3 simulateA.py -np 6 -c 10000 1000000 1000000 tgp_superpop_param.txt simulationA/n1m_1m_snps

python3 simulateB.py -np 6 -c 100000 10000 100000 tgp_superpop_param.txt simulationB/n10k_100k_snps
python3 simulateB.py -np 6 -c 100000 10000 1000000 tgp_superpop_param.txt simulationB/n10k_1m_snps

python3 simulateA.py -np 6 -c 10000 1000000 1000000 tgp_superpop_param.txt simulationA/n1m_1m_snps

