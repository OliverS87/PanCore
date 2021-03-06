# Group assemblies into cluster
# Methods use ParCores intercore region stat table or the unaligned sequence xmfa file
# Cluster methods:
# a) Cluster based on rearrangements
# b) Cluster based on mash all-vs-all run on unaligned sequence
# c) Cluster based on length differences in IC regions
from os import path
from os import mkdir as mkdir
from os import remove
from os import rmdir
import subprocess
import math


# Each cluster method requires
# - Path to intercore region stat table (path to unaligned sequence file is derived from it)
# - Minimum number of multi-assembly cluster
# - Create a plot of each cluster result?
# - Debug mode? I.e. keep temporary files?
class Cluster:
    def __init__(self, ic_stat_p, min_nr_clstr, make_plots, debug):
        self.ic_stat_p = ic_stat_p
        self.min_nr_clstr = min_nr_clstr
        self.make_plots = make_plots
        self.debug = debug

    # a) Cluster based on sequence rearrangements:
    # Algorithm looks for intercore regions that are missing (with respect to the reference), i.e.
    # -> deleted intercore regions
    def cluster_rearrangement(self):
        # Find intercore regions that are not found in assemblies != reference
        del_ic_regions = {}
        # Open and parse the tabular intercore region stat file
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
        # Find max SI (we need this value for some iterations later)
        max_si = max([item[0] for sublist in del_ic_regions.values() for item in sublist[2]])
        # From the list of all intercore regions, keep only those where a deletion occurs for
        # at least ony assembly
        del_ic_regions = {key: val for key, val in del_ic_regions.items() if "D" in val[1]}
        # The cluster algorithm cannot deal well with 'NA's. Remove all intercore regions for whom
        # we don't have information for all assemblies
        del_ic_regions = {key: val[2] for key, val in del_ic_regions.items() if val[0] == max_si}
        # Sort intercore regions by ID and assembly ID
        del_ic_regions = sorted([val for sublist in del_ic_regions.values() for val in sublist],
                                key=lambda x: (x[1], x[0]))
        # Convert the dict into a presence/absence matrix:
        # For each assembly ID, we list if intercore region X is present(1) or absent (0)
        feature_vector = {(reg[0], reg[1]): "1" if reg[2] == "P" else "0" for reg in del_ic_regions}
        # Sort intercore regions by name (not necessary, but makes plot look nicer)
        sorted_intercore_ids = sorted(list(set([item[1] for item in del_ic_regions])))
        # Define output path for feature vector: Same as intercore tabular stat file
        out_p, filename = path.split(self.ic_stat_p)
        filename = filename.replace(".ic.", ".fdel.")
        prefix = filename.split(".")[0]
        # Write feature vector to file
        with open(path.join(out_p, filename), "w") as feat_out_f:
            # Write header
            feat_out_f.write("si," + ",".join(sorted_intercore_ids) + "\n")
            for si in range(1, max_si + 1):
                feat_out_f.write(str(si))
                for ic_reg in sorted_intercore_ids:
                    feat_out_f.write("," + feature_vector[(si, ic_reg)])
                feat_out_f.write("\n")
        # Next, run an rscript to do the clustering
        # First, define the output file name for the cluster results
        cluster_filename = filename.replace(".fdel.",".clstr.")
        # Should the Rscript make plots for the clustering?
        # If so, define the output path
        if self.make_plots:
            png_filepath = path.join(out_p, filename.replace(".fdel.csv",".png"))
        else:
            png_filepath = "NA"
        # Run Rscript
        run = subprocess.run("Rscript --vanilla {5} {0} {1} {2} {3} {4}".format(path.join(out_p, filename),
                                                            path.join(out_p, cluster_filename),
                                                                                png_filepath,
                                                                                self.min_nr_clstr,
                                                                                prefix,
                                                              path.join(out_p, "rearrangement_jac_cluster.r")),
                             stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
        # Remove temporary files
        if not self.debug:
            try:
                remove(path.join(out_p, filename))
            except FileNotFoundError:
                pass
        return run.returncode

    # Cluster based on length differences between intercore regions with the same label
    # -> Cluster based on insertions/deletions
    def cluster_length(self):
        # Only look at intercore regions that are present in every sequence
        # otherwise the length would be 'NA" and that is only upsetting R
        p_ic_region = {}
        # We need the number of sequences later for an iteration
        max_si = 0
        # Parse the intercore stat table file
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
        prefix = filename.split(".")[0]
        filename = filename.replace(".ic.", ".len.")
        # Write filtered ic regions to file
        with open(path.join(out_p, filename), "w") as out_f:
            out_f.write(header)
            for val in p_ic_region.values():
                for item in val:
                    out_f.write(item)
        # Next, run an rscript to do the clustering
        cluster_filename = filename.replace(".len.", ".clstr.")
        if self.make_plots:
            png_filepath = path.join(out_p, filename.replace(".len.csv", ".png"))
        else:
            png_filepath = "NA"
        run = subprocess.run("Rscript --vanilla {5} {0} {1} {2} {3} {4}".format(
            path.join(out_p, filename), path.join(out_p, cluster_filename),
            png_filepath, self.min_nr_clstr, prefix,
        path.join(out_p, "iclength_deviation_eucl_cluster.r")), stderr=subprocess.PIPE,
            stdout=subprocess.PIPE, shell=True)
        # Clean up?
        if not self.debug:
            try:
                remove(path.join(out_p, filename))
            except FileNotFoundError:
                pass
        return run.returncode


    # Calculate pairwise ANI with mash
    # First, split up useqs in individual files, one for each si
    def cluster_ani(self, cpu, all_seqs):
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
        # Define the path to the .unalign file
        useq_p = path.join(out_p, filename_unalign)
        # We next split the .unalign file into individual files, one for each SI
        # Hold links to individual split files in a dict
        sep_files_dict = {}
        # To determine optimum parameters for mash ANI estimation, we need to know the number of nucleotides
        # in the largest useq file
        useq_lengths = {}
        # use all sequences for ANI clustering?
        if not all_seqs:
            # Find lowest and highest core block ID
            core_block_ids = []
            with open(useq_p, "r") as useq_f:
                for line in useq_f:
                    if line.startswith(">1:"):
                        core_block_ids.extend([int(item) for item in line.split()[2][6:].split(".")])
            contig_end_ids = [min(core_block_ids), max(core_block_ids)]
        else:
            contig_end_ids = []
        with open(useq_p, "r") as useq_f:
            for line in useq_f:
                if line.startswith(">"):
                    data = line.split()
                    # If not including all seqs, check whether this seq is from between two core blocks
                    if not all_seqs:
                        cb1, cb2 = [int(item) for item in data[2][6:].split(".")]
                        if cb1 in contig_end_ids or cb2 in contig_end_ids:
                            active_out_f = None
                            continue
                    si = int(data[0].split(":")[0][1:])
                    try:
                        active_out_f = sep_files_dict[si]
                    except KeyError:
                        active_out_f = open(path.join(useqs_split_p, "useq_{0}.faa".format(si)), "w")
                        sep_files_dict[si] = active_out_f
                elif line.startswith("="):
                    active_out_f = None
                    continue
                # All other cases: nt seqs
                if active_out_f:
                    try:
                        useq_lengths[si] += len(line.strip())
                    except KeyError:
                        useq_lengths[si] = len(line.strip())
                    try:
                        active_out_f.write(line)
                    except AttributeError:
                        pass

        # Close all output files
        [f.close() for f in sep_files_dict.values()]
        # Remove reference from useq files
        ref_useq = sep_files_dict.pop(1)
        remove(ref_useq.name)
        useq_lengths.pop(1)
        # Find longest useq file:
        max_length = max(useq_lengths.values())
        # calculate k-mer length for mash distance estimation
        k = math.ceil(math.log10(99*max_length)/math.log10(4))
        # Run mash to sketch all useqs into one sketch
        # Pipe stdout/stderr to prevent printing to screen
        run = subprocess.run("mash sketch -p {0} -k {1} -s {2} -o {3} {4}".
                             format(cpu, k, 5000, mash_sketch_p,
                                    " ".join([f.name for f in sep_files_dict.values()])),
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        # Remove temporary useq files
        if not self.debug:
            [remove(file.name) for file in sep_files_dict.values()]
            rmdir(useqs_split_p)
        # Run mash to evaluate the pairwise distance
        #  ./mash dist  -p 1 -t talia.msh talia.msh > talia.out
        run = subprocess.run("mash dist -p {0} -t {1}.msh {1}.msh > {2}".
                             format(cpu, mash_sketch_p, mash_dist_p), shell=True)
        # If mash sketch or mash dist encountered an error, the return code of mash dist will be
        # != 0
        # Only proceed with clustering when mash ran successfull
        if run.returncode == 0:
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
            # Create a plot of this clustering?
            if self.make_plots:
                png_p = path.join(out_p, filename.replace(".ic.csv", ".png"))
            else:
                png_p = "NA"
            run = subprocess.run("Rscript --vanilla {5} {0} {1} {2} {3} {4}".format(
                mash_dist_p, cluster_p, png_p, self.min_nr_clstr, prefix,
            path.join(out_p, "mash_ani_clustering.r")), stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
        # Clean up
        if not self.debug:
            try:
                remove(mash_sketch_p+".msh")
            except FileNotFoundError:
                pass
            try:
                remove(mash_dist_p)
            except FileNotFoundError:
                pass
        return run.returncode




