#!/bin/env bash

# Path to SCOPE binary
scope_binary=../build/scope

# Run SCOPE unsupervised
$scope_binary -g source_files/example_1k -k 6 -seed 12345 -o unsupervised_example_1k_

# Run SCOPE supervised
$scope_binary -g source_files/example_1k -k 6 -seed 12345 -freq source_files/example_1k.plink.freq -o supervised_example_1k_
