# Calculate clusters from SimpleParSNPs ic.stat file
# a) Cluster based on rearrangements
# b) Cluster based on length differences in IC regions
# c) Cluster based on mash all-vs-all run on unaligned sequences
from os import path
from os import mkdir as mkdir
import subprocess
import math

class Cluster:
    def __init__(self, ic_stat_p, min_nr_clstr):
        self.ic_stat_p = ic_stat_p
        self.min_nr_clstr = min_nr_clstr


    # Calculate pairwise ANI with mash
        # First, split up useqs in individual files, one for each si
    def cluster_ani(self, cpu):
        # Retrieve path to .unalign file from the ic_stat file (same folder, different filename)
        out_p, filename = path.split(self.ic_stat_p)
        # Modify the filename to access the .unalign file
        filename_unalign = filename.replace(".ic.csv", ".unalign")
        # Define file name to mash sketch and mash dist results
        mash_sketch_p = path.join(out_p, filename.replace(".ic.csv", ".sketch"))
        mash_dist_p = path.join(out_p, filename.replace(".ic.csv", ".msh.dist"))
        # Each round of clustering gets its unique prefix, retrieve it here
        prefix = filename.split(".")[0]
        # Unaligned sequence file is split into individual files, one for each si
        # Define path to the temporary folder holding these individual files
        useqs_split_p = path.join(out_p, "{0}_useqs".format(prefix))
        # Try to create the folder for the temporary .unalign split files
        try:
            mkdir(useqs_split_p)
        except FileExistsError:
            pass
        # File name to global .unalign file
        useq_p = path.join(out_p, filename_unalign)
        # Hold links to individual split files in a dict
        sep_files_dict = {}
        # To determine optimum parameters for mash ANI estimation, we need to know the number of nucleotides
        # in the largest useq file
        useq_lengths = {}
        with open(useq_p, "r") as useq_f:
            for line in useq_f:
                if line.startswith(">"):
                    data = line.split()
                    si = int(data[0].split(":")[0][1:])
                    try:
                        active_out_f = sep_files_dict[si]
                    except KeyError:
                        active_out_f = open(path.join(useqs_split_p, "useq_{0}.faa".format(si)), "w")
                        sep_files_dict[si] = active_out_f
                elif line.startswith(">") or line.startswith("="):
                    active_out_f = None
                    continue
                # All other cases: nt seqs
                try:
                    useq_lengths[si] += len(line.strip())
                except KeyError:
                    useq_lengths[si] = len(line.strip())
                active_out_f.write(line)
        # Close all output files
        [f.close() for f in sep_files_dict.values()]
        # Remove reference from useq files
        sep_files_dict.pop(1)
        useq_lengths.pop(1)
        # Find longest useq file:
        max_length = max(useq_lengths.values())
        # calculate k-mer length for mash distance estimation
        k = math.ceil(math.log10(99*max_length)/math.log10(4))
        # Run mash to sketch all useqs into one sketch
        run = subprocess.run("mash sketch -p {0} -k {1} -s {2} -o {3} {4}".
                             format(cpu, k, 5000, mash_sketch_p,
                                    " ".join([f.name for f in sep_files_dict.values()])), shell=True)
        if run.returncode != 0:
            return run.returncode
        # Run mash to evaluate the pairwise distance
        #  ./mash dist  -p 1 -t talia.msh talia.msh > talia.out
        run = subprocess.run("mash dist -p {0} -t {1}.msh {1}.msh > {2}".
                             format(cpu, mash_sketch_p, mash_dist_p), shell=True)
        if run.returncode != 0:
            return run.returncode
        # Modify mash output, reduce identifier to si (currently it is path)
        modified_out = []
        with open(mash_dist_p, "r") as mash_dist_f:
            # First correct header
            header = next(mash_dist_f)
            data = header.split()
            modified_out.append([cell.rsplit("_")[-1].replace(".faa", "") for cell in data])
            for line in mash_dist_f:
                data = line.split()
                data[0] = data[0].rsplit("_")[-1].replace(".faa", "")
                modified_out.append(data)
        with open(mash_dist_p, "w") as new_mash_dist_f:
            for line in modified_out:
                new_mash_dist_f.write("\t".join(line)+"\n")
        # Next, run an rscript to do the clustering
        cluster_p = path.join(out_p, filename.replace(".ic.", ".clstr."))
        png_p =  path.join(out_p, filename.replace(".ic.csv", ".png"))
        run = subprocess.run("Rscript {5} {0} {1} {2} {3} {4}".format(
            mash_dist_p, cluster_p, png_p, self.min_nr_clstr, prefix,
        path.join(out_p, "mash_ani_clustering.r")), shell=True)
        return run.returncode

    def cluster_length(self):
        p_ic_region = {}
        max_si = 0
        with open(self.ic_stat_p, "r") as ic_stat_f:
            # Skip header
            header = next(ic_stat_f)
            for line in ic_stat_f:
                data = line.strip().split(",")
                cb1 = int(data[1])
                cb2 = int(data[2])
                si = int(data[0])
                max_si = max(max_si, si)
                type = data[6]
                if type == "P":
                    try:
                        p_ic_region[(cb1, cb2)].append(line)
                    except KeyError:
                        p_ic_region[(cb1, cb2)] = [line]
        # Keep only those ic regions where all assemblies are present
        p_ic_region = {key: val for key, val in p_ic_region.items() if len(val) == max_si}
        # Define output path
        out_p, filename = path.split(self.ic_stat_p)
        filename = filename.replace(".ic.", ".len.")
        # Write filtered ic regions to file
        with open(path.join(out_p, filename), "w") as out_f:
            out_f.write(header)
            for val in p_ic_region.values():
                for item in val:
                    out_f.write(item)
        # Next, run an rscript to do the clustering
        cluster_filename = filename.replace(".len.", ".clstr.")
        png_filename = filename.replace(".csv", ".png")
        run = subprocess.run("Rscript {4} {0} {1} {2} {3}".format(
            path.join(out_p, filename), path.join(out_p, cluster_filename),
            path.join(out_p, png_filename), self.min_nr_clstr,
        path.join(out_p, "iclength_deviation_eucl_cluster.r")), shell=True)
        return run.returncode


    def cluster_rearrangement(self):
        # Identify CB Nr. signaling 5" and 3" end
        min_cb = None
        max_cb = None
        del_ic_regions = {}
        with open(self.ic_stat_p, "r") as ic_stat_f:
            # Skip header
            next(ic_stat_f)
            for line in ic_stat_f:
                data = line.strip().split(",")
                cb1 = int(data[1])
                cb2 = int(data[2])
                si = int(data[0])
                icb_id = data[1] + "." + data[2]
                type = data[6]
                try:
                    del_ic_region = del_ic_regions[(cb1, cb2)]
                    del_ic_region[0] += 1
                    del_ic_region[1].add(type)
                    del_ic_region[2].append((si, icb_id, type))
                except KeyError:
                    del_ic_regions[(cb1, cb2)] = [1, set(type), [(si, icb_id, type)]]
        # Find max SI
        max_si = max([item[0] for sublist in del_ic_regions.values() for item in sublist[2]])
        # Filter regions, keep only those where a D occurs
        del_ic_regions = {key: val for key, val in del_ic_regions.items() if "D" in val[1]}
        # Filte regions, keep only those where there is one occurence per SI
        del_ic_regions = {key: val[2] for key, val in del_ic_regions.items() if val[0] == max_si}
        # Sort regions by region and si
        del_ic_regions = sorted([val for sublist in del_ic_regions.values() for val in sublist],
                                key=lambda x: (x[1], x[0]))
        # Create the feature vector matrix for jaccard clustering of IC regions
        feature_vector = {(reg[0], reg[1]): "1" if reg[2] == "P" else "0" for reg in del_ic_regions}
        sorted_ic_reg_list = sorted(list(set([item[1] for item in del_ic_regions])))
        # Define output path
        out_p, filename = path.split(self.ic_stat_p)
        filename = filename.replace(".ic.", ".fdel.")
        with open(path.join(out_p, filename), "w") as feat_out_f:
            # Write header
            feat_out_f.write("si," + ",".join(sorted_ic_reg_list) + "\n")
            for si in range(1, max_si + 1):
                feat_out_f.write(str(si))
                for ic_reg in sorted_ic_reg_list:
                    feat_out_f.write("," + feature_vector[(si, ic_reg)])
                feat_out_f.write("\n")
        # Next, run an rscript to do the clustering
        cluster_filename = filename.replace(".fdel.",".clstr.")
        run = subprocess.run("Rscript {3} {0} {1} {2}".format(path.join(out_p, filename),
                                                            path.join(out_p, cluster_filename),
                                                              self.min_nr_clstr,
                                                              path.join(out_p, "rearrangement_jac_cluster.r")),
                             shell=True)
        return run.returncode

