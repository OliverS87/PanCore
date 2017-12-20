# Run parsnp on a collection of sequences (+ a reference)
# Generate sets with sequentially increasing size:
# Start with 2, add more to go to 5, 10, 20, 50, 100, 200, 300, 400

import sys
import os
from random import shuffle
from SimpleParSNP2 import simple_parsnp


if len(sys.argv[1:]) != 8:
    print("Usage: python3 {0} reference sample_dir distance_value cpu_count out_folder substract sets prefix".format(
        sys.argv[0]))
    exit(3)
ref_path = sys.argv[1]
dist_val = int(sys.argv[3])
cpu_count = int(sys.argv[4])
out_folder = sys.argv[5]
assemblies_dir = sys.argv[2]
substract = abs(int(sys.argv[6]))*-1
nr_of_sets = int(sys.argv[7])
prefix = sys.argv[8]
file_list = [os.path.join(assemblies_dir, file) for file in os.listdir(assemblies_dir)
                      if os.path.isfile(os.path.join(assemblies_dir, file))]

# Run a core analysis for all sample
first_simple_snp = simple_parsnp()
first_simple_snp.set_dist(dist_val)
first_simple_snp.set_size(21)
first_simple_snp.set_reference(ref_path)
first_simple_snp.set_threads(cpu_count)
first_simple_snp.set_prefix("{0}_{1}".format(prefix, len(file_list)))
first_simple_snp.add_files(file_list)
first_simple_snp.run_parsnp(out_folder, False, False)

for set in range(0, nr_of_sets):
    this_file_list = [item for item in file_list]
    shuffle(this_file_list)
    this_file_list = this_file_list[:substract]
    while this_file_list:
        set_simple_snp = simple_parsnp()
        set_simple_snp.set_dist(dist_val)
        set_simple_snp.set_size(21)
        set_simple_snp.set_reference(ref_path)
        set_simple_snp.set_threads(cpu_count)
        set_simple_snp.set_prefix("{0}_R{1}_{2}".format(prefix, set, len(this_file_list)))
        set_simple_snp.add_files(this_file_list)
        set_simple_snp.run_parsnp(out_folder, False, False)
        # Reduce file list
        shuffle(this_file_list)
        this_file_list = this_file_list[:substract]


