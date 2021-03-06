# Takes the result of a clustering and runs separte SimpleParSNP runs for each cluster
import os
import shutil
import sys

from PanCoreLibs.cluster import Cluster
from PanCoreLibs.iclength_deviation_eucl_cluster import IclengthClusterRscript
from PanCoreLibs.mash_ani_clustering import MashAnoClusteringRscript
from PanCoreLibs.rearrangement_jac_cluster import RearrangementJacCluster
from ParCore import ParCore
from PanCoreLibs.reduce_input import ReduceInput
from PanCoreLibs.random_cluster import RandomCluster
from PanCoreLibs.build_xmfa import BuildXmfa

def clean_up(outpath, prefix, keep_all_core, debug):
    if debug:
        return
    # Remove parsnp tmp out if it still exists
    shutil.rmtree(os.path.join(outpath, prefix), ignore_errors=True)
    # Remove parsnp config file if it still exists
    try:
        os.remove(os.path.join(outpath, "{0}parsnp_config.ini".format(prefix)))
    except FileNotFoundError:
        pass
    if not keep_all_core:
        try:
            os.remove(os.path.join(outpath, "{0}.xmfa".format(prefix)))
        except FileNotFoundError:
            pass
    try:
        os.remove(os.path.join(outpath, "{0}.ic.csv".format(prefix)))
    except FileNotFoundError:
        pass
    try:
        os.remove(os.path.join(outpath, "{0}.fdel.csv".format(prefix)))
    except FileNotFoundError:
        pass
    try:
        os.remove(os.path.join(outpath, "{0}.clstr.csv".format(prefix)))
    except FileNotFoundError:
        pass
    try:
        os.remove(os.path.join(outpath, "{0}.len.csv".format(prefix)))
    except FileNotFoundError:
        pass
    try:
        os.remove(os.path.join(outpath, "{0}.unalign".format(prefix)))
    except FileNotFoundError:
        pass

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("reference", help="Full length reference sequence file")
    parser.add_argument("sample_folder", help="Folder where all sample files are located")
    parser.add_argument("-p", "--prefix", default="pansnp", help="Prefix for all output files")
    parser.add_argument("-o", "--outdir", default="pansnp", help="Output path")
    parser.add_argument("-c", "--cpu", default=1, type=int, help="Number of CPU cores")
    parser.add_argument("-s", "--size", default=21, type=int, help="Minimum core block size")
    parser.add_argument("-d", "--distance", default=30, type=int,
                        help="Maximum distance between two MUMS in a core block")
    #parser.add_argument("-u", "--unaligned", action="store_true", help="Output unaligned regions?")
    parser.add_argument("-a", "--all_core", action="store_true",
                        help="Output core block alignment for each sample subset")
    parser.add_argument("-b", "--debug", action="store_true",
                        help="Keep temporary files")
    parser.add_argument("-i", "--plot", action="store_true",
                        help="Plot cluster for each sample subset")
    parser.add_argument("-r", "--random", type=int,
                        help="For each clustering create x random clusterings with clusters of same size.")
    parser.add_argument("-e", "--reduce", action="store_true",
                        help="Reduce input after each round of clustering to the unclustered parts")
    parser.add_argument("-m", '--method', choices=["r", "sa", "sc", "l"],
                        help="Cluster by 'R'earrangements, 'S'imilarity 'A'll or between 'C'ore blocks "
                             "or 'L'ength", default="sc")
    parser.add_argument("-l", "--cluster", default=2, type=int, help="Max. number of multi-element cluster created during each cycle")
    args = parser.parse_args()

    ref_p = args.reference
    samples_folder = args.sample_folder
    out_p = args.outdir
    dist_param = args.distance
    cpu_count = args.cpu
    size_param = args.size
    cluster_method = args.method
    keep_all_core = args.all_core
    plot = args.plot
    reduce = args.reduce
    debug = args.debug
    if args.random:
        random = max(0, args.random)
    else:
        random = 0
    # Create at least two multi-cluster per iteration
    min_cluster = max(2, args.cluster)
    prefix = args.prefix
    # Create output folder
    if not os.path.isdir(out_p):
        os.makedirs(out_p, exist_ok=True)
    # All the consecutive alignments created by PanCore are
    # merged into one large xmfa file at the end
    # Overlaps between alignments are resolved, with that alignment with more
    # sequences "winning"
    build_xmfa = BuildXmfa(out_p, prefix)
    # Initialize the cluster queue
    # Each element is a list of files that are going to be core clustered by parsnp
    # each item comes in this format (id[file1,file2,file3,...])
    # Each parsnp run gives rise to more cluster, these are appended to the end of the queue
    # A cluster comes with a boolean subset? If subset? is FALSE, this cluster will not
    # be split into subcluster. This is used for randomized cluster.
    parsnp_queue = [(prefix,[], ref_p, True)]
    # Get all files in sample folder
    [parsnp_queue[0][1].append(os.path.join(samples_folder, item)) for item in os.listdir(samples_folder)
     if os.path.isfile(os.path.join(samples_folder, item))]
    # Clustering may involve somecomputations with R
    # R is called via Rscript
    # To avoid unnecesary problems with absolute and relative paths,
    # the R scripts are written into the defined output_folder and deleted again
    # when Pansnp terminates
    if cluster_method == "l":
        rscript = IclengthClusterRscript(out_p)
    elif cluster_method == "sa" or cluster_method == "sc":
        rscript = MashAnoClusteringRscript(out_p)
    else:
        rscript = RearrangementJacCluster(out_p)
    rscript.write_script()
    # Load Reduce input module
    reducer = ReduceInput(out_p)
    while parsnp_queue:
        print("Length of queue: {0}".format(len(parsnp_queue)))
        # 1. Get the first element of the queue and run parsnp on it
        prefix, file_list, this_ref, subsetting = parsnp_queue.pop(0)
        # if there is only one sequence to core-analyse, don't do it
        if len(file_list) <= 1:
            print("Skipping {0} because length <= 1".format(prefix))
            continue
        # Prepare simpleparsnp run
        sp = ParCore()
        sp.set_dist(dist_param)
        sp.set_reference(this_ref)
        sp.set_threads(cpu_count)
        sp.add_files(file_list)
        sp.set_prefix(prefix)
        # Run parsnp, if clustering by mash/ani, create the unaligned file
        # Else, create the intracore region stat file
        if cluster_method == "sa" or cluster_method == "sc":
            sp_rc = sp.run_parcore(out_p, True, False)
        else:
            sp_rc = sp.run_parcore(out_p, (False or reduce), True)
        if sp_rc != 0:
            print("ParCore for {0} failed.\n{1}".format(prefix, sp_rc))
            clean_up(out_p, prefix, keep_all_core, debug)
            continue
        # If the current cluster shall not be split into subcluster, stop here and continue with the next
        # cluster in the queue
        if not subsetting:
            continue
        # Is the reduce input flag set?
        if reduce:
            # Map the indices of this rounds .xmfa and .unalign back to the
            # coordinate system of the original input
            if reducer.sequences_added():
                reducer.map_back(file_list, os.path.join(out_p, "{0}.xmfa".format(prefix)),
                                                         os.path.join(out_p, "{0}.unalign".format(prefix)))
            # Create new input sample sequences from the unaligned part of this input
            reducer.useq_to_input(file_list, os.path.join(out_p, "{0}.unalign".format(prefix)))
            reducer.mask_reference(this_ref, os.path.join(out_p, "{0}.xmfa".format(prefix)), prefix)
            file_list = [os.path.join(out_p, "new_input", os.path.basename(item)) for item in file_list]
            next_ref = os.path.join(out_p, "new_input", "{0}_ref.faa".format(prefix))
        else:
            next_ref = this_ref
        # Add the xmfa alignment for this cluster to the global xmfa alignment:
        # The global alignment removes all sequence regions from the alignment that are
        # part of any alignment in any parent node of this node
        build_xmfa.add_xmfa(os.path.join(out_p, "{0}.xmfa".format(prefix)))
        # Split the assemblies associated with this node into subsets
        # The cluster method is by
        # a) Rearrangements: Missing intercore regions indicate sequence rearrangements
        # b) Average nucleotide identity: ANI between the sequences, either full length sequences or the
        #   region between core blocks
        # c) Insertions/deletions: Differences in intercore region lengths
        # Data used for the clustering is the tabular intercore stat file (length and rearrangement) or
        # the unaligned sequences file (ANI) produced by ParCore.
        # Every cluster method returns the cluster result in tabular format
        # If the cluster method fails, an return code != is returned
        cluster = Cluster(os.path.join(out_p, "{0}.ic.csv".format(prefix)), min_cluster, plot, debug)
        if cluster_method == "r":
            clustering = cluster.cluster_rearrangement()
        elif cluster_method == "sa" or cluster_method == "sc":
            clustering = cluster.cluster_ani(cpu_count, cluster_method == "sa")
        else:
            clustering = cluster.cluster_length()
        if clustering != 0:
            print("clustering for {0} failed.".format(prefix))
            clean_up(out_p, prefix, keep_all_core, debug)
            continue
        # Open and parse the cluster file
        cluster_db = {}
        with open(os.path.join(out_p, "{0}.clstr.csv".format(prefix)), "r") as cluster_f:
            # Skip header
            next(cluster_f)
            for line in cluster_f:
                data = line.split()
                try:
                    cluster_db[int(data[1])].append(file_list[int(data[0])-2])
                except KeyError:
                    cluster_db[int(data[1])] = [file_list[int(data[0])-2]]
        # Sort cluster by length
        cluster_list = sorted(cluster_db.values(), key=lambda x: -len(x))
        # If there is only one cluster, stop here
        if len(cluster_list) == 1:
            clean_up(out_p, prefix, keep_all_core, debug)
            continue
        # Create some random clusters?
        for r in range(0, random):
            rc = RandomCluster()
            random_cluster_list = rc.randomize(cluster_list)
            # Add random cluster to the queue
            for i, clstr in enumerate(random_cluster_list):
                parsnp_queue.append(("{0}_{1}_R{2}".format(prefix, i, r), clstr, next_ref, False))
        # Add cluster to the queue
        for i, clstr in enumerate(cluster_list):
            parsnp_queue.append(("{0}_{1}".format(prefix, i), clstr, next_ref, True))
        clean_up(out_p, prefix, keep_all_core, debug)
    # Remove rscript that was used for clustering
    rscript.remove_script()
    # Build the global xmfa alignment file
    build_xmfa.generate_combined_xmfa()
    # That's it, mate :) Nothing left to do
