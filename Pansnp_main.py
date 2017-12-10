# Takes the result of a clustering and runs separte SimpleParSNP runs for each cluster
import os
import shutil
import sys

from pansnp_libs.cluster import Cluster
from pansnp_libs.iclength_deviation_eucl_cluster import IclengthClusterRscript
from pansnp_libs.mash_ani_clustering import MashAnoClusteringRscript
from SimpleParSNP import SimpleParSNP


def clean_up(outpath, prefix):
    return
    # Remove parsnp tmp out if it still exists
    shutil.rmtree(os.path.join(outpath, prefix), ignore_errors=True)
    # Remove parsnp config file if it still exists
    try:
        os.remove(os.path.join(outpath, "{0}parsnp_config.ini".format(prefix)))
    except FileNotFoundError:
        pass
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


if __name__ == '__main__':
    args = sys.argv
    # usage: python3 run_cluster.py cluster.csv prev_xmfa ref_path sample_folder out_path
    if len(args[1:]) != 7:
        print('usage: python3 {0} sample_folder ref_path out_path dist cpu cluster[r|l|s] min_cluster'.format(
            args[0]))
        exit(1)
    ref_p = args[2]
    samples_folder = args[1]
    out_p = args[3]
    dist_param = int(args[4])
    cpu_count = int(args[5])
    if args[6].lower() == "r":
        cluster_method = "rec"
    elif args[6].lower() == "s":
        cluster_method = "sim"
    else:
        cluster_method = "len"
    min_cluster = int(args[7])


    # Initialize the cluster queue
    # Each element is a list of files that are going to be core clustered by parsnp
    # each item comes in this format (id[file1,file2,file3,...])
    # Each parsnp run gives rise to more cluster, these are appended to the end of the queue
    parsnp_queue = [("talia",[])]
    # Get all files in sample folder
    [parsnp_queue[0][1].append(os.path.join(samples_folder, item)) for item in os.listdir(samples_folder)
     if os.path.isfile(os.path.join(samples_folder, item))]
    # Clustering may involve somecomputations with R
    # R is called via Rscript
    # To avoid unnecesary problems with absolute and relative paths,
    # the R scripts are written into the defined output_folder and deleted again
    # when Pansnp terminates
    if cluster_method == "len"
        rscript = IclengthClusterRscript(out_p)
    elif cluster_method == "sim":
        rscript = MashAnoClusteringRscript(out_p)
    elif cluster_method == "rec":
        rscript = None
    rscript.write_script()
    while parsnp_queue:
        print("Length of queue: {0}".format(len(parsnp_queue)))
        # 1. Get the first element of the queue and run parsnp on it
        prefix, file_list = parsnp_queue.pop(0)
        # if there is only one sequence to core-analyse, don't do it
        if len(file_list) <= 1:
            print("Skipping {0} because length <= 1".format(prefix))
            continue
        # Prepare simpleparsnp run
        sp = SimpleParSNP()
        sp.set_dist(dist_param)
        sp.set_reference(ref_p)
        sp.set_threads(cpu_count)
        sp.add_files(file_list)
        sp.set_prefix(prefix)
        # Run parsnp, do not create unaligned file
        # But create stat file
        sp_rc = sp.run_parsnp(out_p, True, True)
        if sp_rc != 0:
            print("parsnp for {0} failed.\n{1}".format(prefix, sp_rc))
            clean_up(out_p, prefix)
            continue
        # Convert the IC stat file
        cluster = Cluster(os.path.join(out_p, "{0}.ic.csv".format(prefix)), min_cluster)
        if cluster_method == "rec":
            clustering = cluster.cluster_rearrangement()
        elif cluster_method == "sim":
            clustering = cluster.cluster_ani(cpu_count)
        else:
            clustering = cluster.cluster_length()
        if clustering != 0:
            print("clustering for {0} failed.".format(prefix))
            clean_up(out_p, prefix)
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
            clean_up(out_p, prefix)
            continue
        # Else, add cluster to the queue
        for i, clstr in enumerate(cluster_list):
            parsnp_queue.append(("{0}_{1}".format(prefix, i), clstr))
        clean_up(out_p, prefix)
    rscript.remove_script()