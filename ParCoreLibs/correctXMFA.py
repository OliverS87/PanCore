from os import path
from datetime import datetime as dt


class CorrectXMFA:
    # Required instance values:
    # - parsnp distance value
    # - path to xmfa raw output
    # - min SI and max SI
    # - Global contig start/stop for each contig for each SI
    def __init__(self, dist_param, xmfa_path, min_si, max_si, input_contigs_dict, write_log):
        # The sequence index offset depends on the contig index nr.
        #  and the distance parameter used by parsnp:
        # (contig_id-1)*(distance+10)
        # The (+10) is added here
        self.dist_param = int(dist_param) + 10
        self.xmfa_path = xmfa_path
        self.min_si = int(min_si)
        self.max_si = int(max_si)
        self.input_contigs_dict = input_contigs_dict
        self.write_log = write_log

    # Main function to call all steps necessary to correct the
    # raw XMFA file in out_path
    def correct_xmfa(self):
        # 1. Correct index error in header
        # Function returns a dict with the corrected headers
        # or an error code
        start = dt.now()
        # corrected_header contains only the sequence headers
        # the actual nucleotide sequences have not been parsed
        # Each sequence has a unique combination of assembly ID and cluster ID
        # With this pair of values it is possible to match the corrected header sequence
        # back to the corresponding nucleotide sequence
        # This prevents the need for keeping all nucleotide sequences in memory
        corrected_header = self.__correct_header__()
        time_elapsed = (dt.now()-start).total_seconds()
        # Check if an error code was returned
        try:
            int(corrected_header)
            self.write_log("Error while correcting xmfa header. Bummer :(\n")
            return corrected_header
        # if not, continue
        except TypeError:
            self.write_log("Corrected xmfa header in {0} seconds.\n".format(time_elapsed))
            pass
        # Search in header for clusters that are redundant, i.e. already contained in
        # a longer clustered sequence
        start = dt.now()
        redundant_cluster = self.__find_redundant_clusters__(corrected_header)
        # Filter header list, remove header that belong to a redundant cluster
        corrected_header = {key:val for key, val in corrected_header.items() if val[4] not in redundant_cluster}
        time_elapsed = (dt.now() - start).total_seconds()
        self.write_log("Removed {0} redundant cluster in {1} seconds.\n".format(len(redundant_cluster), time_elapsed))
        # Write new xmfa file
        start = dt.now()
        self.__write_new_xmfa__(corrected_header)
        time_elapsed = (dt.now() - start).total_seconds()
        self.write_log("Wrote cluster alignment file in {0} seconds.\n".format(time_elapsed))
        # Return the header information
        # Other functions can use this list instead of reading in the xmfa file again
        return corrected_header


    # Read in the header lines from the raw xmfa output
    # and apply correction formula.
    # Return dict with header sequences, sorted by si and cluster
    def __correct_header__(self):
        # Reformat the input_contigs_dict:
        # input_contigs[si]=[[si, contig_id, global_start, global_stop, contig_start, contig_stop, contig_length,
        # strand, left_coreblock, right_coreblock],...]
        # Into this format:
        # contig_bounds_corr[si] = assembly_contig_bounds_corr[(global_start,global_stop)] = ContigID
        contig_bounds_corr = {key: {(item[2], item[3]): item[1] for item in val} for key, val
                              in self.input_contigs_dict.items()}
        # Calculate contig length from end-pos minus start-pos +1
        # Final result is a list of contig length sorted by contig ID, starting with contig ID 1
        # These lists are stored in contig_lengths dict under their si
        contig_lengths = {}
        for si in range(self.min_si, self.max_si + 1):
            contig_lengths[si] = [abs(item[0][1] - item[0][0]) + 1 for item in sorted(contig_bounds_corr[si].items(),
                                                                                      key=lambda x: x[1])]
        # Correct xmfa file header
        # Store each corrected header in core_header dict
        # Format: core_header[(si, clusterID)] = header-line
        # Some core fragments are redundant, i.e. they appear a second time inside a larger fragment
        # Those header are removed at a later step
        core_header = {}
        # Iterate through raw XMFA file, look only at header lines
        with open(self.xmfa_path, "r") as xmfa_f:
            for line in xmfa_f:
                if line.startswith(">"):
                    # Get old start/stop indices
                    data = line.split()
                    cluster_id = int(data[2].replace("cluster", ""))
                    old_start = int(data[0][data[0].index(":") + 1:data[0].index("-")])
                    old_stop = int(data[0][data[0].index("-") + 1:])
                    # Sometimes there is no sx annotation for contig id
                    # Assume that it is contig 1
                    try:
                        contig_id = int(data[3][1:data[3].index(":")])
                    except ValueError:
                        contig_id = 1
                    seq_id = int(data[0][1:data[0].index(":")])
                    # Use the contig id and the offset constant to calculate the correct start/stop indices
                    new_start = old_start - (contig_id - 1) * self.dist_param
                    new_stop = old_stop - (contig_id - 1) * self.dist_param
                    cluster_length = abs(new_stop - new_start) + 1
                    cluster_id = int(data[2].replace("cluster", ""))
                    # The contig index may have changed with the new indices
                    assembly_contig_bounds = contig_bounds_corr[seq_id]
                    this_contig_bounds = [item for item in assembly_contig_bounds.keys() if
                                          item[0] <= new_start and item[1] >= new_stop - 1]
                    if len(this_contig_bounds)  != 1:
                        print("Failed to find: {0} {1} in\n{2}".format(new_start,
                                                                       new_stop, assembly_contig_bounds.keys()))
                        return -1
                    new_contig_id = assembly_contig_bounds[this_contig_bounds[0]]
                    # Calculate new start position on contig:
                    # This is the start position on merged, minus all previous contig lengths.
                    new_contig_start = new_start - sum(contig_lengths[seq_id][:new_contig_id - 1])
                    new_header = (seq_id, new_start, new_stop, data[1],cluster_id,new_contig_id, new_contig_start)
                    core_header[(seq_id, cluster_id)] = new_header
        return core_header

    # Some clusters are redundant, i.e. they are contained completely within a larger cluster
    # All cluster indices are sorted by start and stop values.
    def __find_redundant_clusters__(self, core_header):
        redundant_clusters = []
        core_fragment_indices = {si:[] for si in range(self.min_si, self.max_si + 1)}
        # Format for core_fragments: {si:[(start,stop,clstrID), (start,stop,clstrID), ...]}
        [core_fragment_indices[item[0]].append((item[1], item[2], item[4])) for item in core_header.values()]
        for si in range(self.min_si, self.max_si + 1):
            # Sort ascending by start and descending by stop, i.e. for clusters with the same start index
            # the largest cluster is in the first position
            core_fragments = sorted(core_fragment_indices[si], key=lambda x: (x[0], -x[1]))
            # Put all non-redundant clusters in core_fragments_no_red
            # The first cluster in the sorted list cannot be redundant, so add it right now
            core_fragments_no_red = [core_fragments[0]]
            # Compare each element in the sorted core_fragments list with the last element in the
            # non_redundant clusters list. If the first element in the core_fragments list starts before
            # the last non_redundant cluster ends, it is contained within that non_redundant cluster
            # and needs to be filtered out
            for frag in core_fragments[1:]:
                if frag[0] <= core_fragments_no_red[-1][1]:
                    redundant_clusters.append(frag[2])
                else:
                    core_fragments_no_red.append(frag)
        return redundant_clusters

    def __write_new_xmfa__(self, corr_core_header):
        # Define out path
        new_xmfa_p = self.xmfa_path+".corr"
        with open(self.xmfa_path, "r") as xmfa_f, open(new_xmfa_p, "w") as new_xmfa_f:
            remove_clstr = False
            for line in xmfa_f:
                if line.startswith(">"):
                    # Retrieve corrected header line
                    data = line.split()
                    seq_id = int(data[0][1:data[0].index(":")])
                    cluster_id = int(data[2].replace("cluster", ""))
                    try:
                        header_data = corr_core_header[(seq_id, cluster_id)]
                        new_header = ">{0}:{1}-{2} {3} cluster{4} s{5}:p{6}\n".format(*header_data)
                        new_xmfa_f.write(new_header)
                        remove_clstr = False
                    except KeyError:
                        remove_clstr = True
                else:
                    if not remove_clstr:
                        new_xmfa_f.write(line)
