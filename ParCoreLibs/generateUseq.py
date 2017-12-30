# Generate Useq from a parsnp XMFA file
from Bio import SeqIO
from os import path
from datetime import datetime as dt
class GenerateUseq:
    def __init__(self, seq_db, xmfa_header, input_contigs, input_contigs_names, out_path, prefix, write_log):
        # Seq_db matches SI to file path
        self.seq_db = seq_db
        # Define smallest and largest Sequence Index, needed later to iterate through each assembly
        self.min_si = min(self.seq_db.keys())
        self.max_si = max(self.seq_db.keys())
        # Header from the xmfa file, so we don't have to read in the whole xmfa again
        self.xmfa_header = xmfa_header
        # Two DB with information about all contigs used as input (incl. ref)
        self.input_contigs = input_contigs
        self.input_contigs_names = input_contigs_names
        # Define out folder
        self.out_path = out_path
        self.prefix = prefix
        # Add messages to log file
        self.write_log = write_log

    # Fragments the contigs in input_contigs by core regions found in the xmfa file
    def __fragment_contigs(self, clusters):
        for si in range(self.min_si, self.max_si + 1):
            # Dict cluster by contig and order them by global start
            this_si_clusters = {}
            for clstr in clusters[si]:
                try:
                    this_si_clusters[clstr[1]].append(clstr)
                except KeyError:
                    this_si_clusters[clstr[1]] = [clstr]
            # Sort clusters
            this_si_clusters = {key: sorted(val, key=lambda x: x[2]) for key, val in this_si_clusters.items()}
            # this_si_clusters = sorted(this_si_clusters, key = lambda x: (x[1],x[2]))
            # Retrieve all contigs for this si and store them by contig id
            this_si_contigs = {}
            for contig in self.input_contigs[si]:
                try:
                    this_si_contigs[contig[1]][0].append(contig)
                except KeyError:
                    # Append two lists to this contig SI/Contig ID combo:
                    # First one is for the unprocessed contigs
                    # Second one is for the already fragmented
                    # During fragmentation, each fragment is only looked at once,
                    # saving some computational time
                    this_si_contigs[contig[1]] = [[contig], []]
            # Iterate over all core regions for this si
            # Clusters are sorted by contig ID, start with first contig ID
            for ctg_id, ctg_clusters in this_si_clusters.items():
                this_si_contig = this_si_contigs[ctg_id]
                for cluster in ctg_clusters:
                    while this_si_contig[0] and not (this_si_contig[0][0][2] <= cluster[2]
                                                     and this_si_contig[0][0][3] >= cluster[3]):
                        this_si_contig[1].append(this_si_contig[0][0])
                        this_si_contig[0].pop(0)
                    # Break if there are no more fragments left for this contig id
                    if not this_si_contig[0]:
                        break
                    host_ctg_fragment = this_si_contig[0][0]
                    # Remove core block region from contig fragment, thereby splitting it into two fragments
                    # Fragments can be of size <1, these fragments are stored as size 0 fragments
                    # Format of contig fragments remains the same, as defined for input_contigs above.
                    # left_coreblock is the current core_block_seq ID for the new_right_ctg_fragment and
                    # for new_left_ctg_fragment it is right_coreblock. The other two ?_coreblock values
                    # remain unchanged (remain either contig-end or other core block)
                    new_left_ctg_fragment = (
                    host_ctg_fragment[0], host_ctg_fragment[1], host_ctg_fragment[2], cluster[2] - 1,
                    host_ctg_fragment[4], cluster[4] - 1, max(cluster[2] - host_ctg_fragment[2], 0),
                    host_ctg_fragment[7], host_ctg_fragment[8], cluster[8])
                    new_right_ctg_fragment = (
                    host_ctg_fragment[0], host_ctg_fragment[1], cluster[3] + 1, host_ctg_fragment[3],
                    cluster[5] + 1, host_ctg_fragment[5],
                    max(host_ctg_fragment[5] - cluster[5] + 2, 0),
                    host_ctg_fragment[7], cluster[8], host_ctg_fragment[9])
                    # Remove original ctg fragment from this_si_contigs list
                    this_si_contig[0].pop(0)
                    # Fragment to the right of the core region could become further fragmented by the next core region
                    # Fragment to the left is not, because clusters were sorted by start position
                    # The cut-out core block region is not part of the contig fragments anymore
                    this_si_contig[1].append(new_left_ctg_fragment)
                    this_si_contig[0].insert(0, new_right_ctg_fragment)
                # After all core regions for this si and contig id are processed,
                # return fragmented contig to input_contigs list
                this_si_contigs[ctg_id] = [this_si_contig[0] + this_si_contig[1]]
                # print(this_si_contigs[ctg_id])
            self.input_contigs[si] = [item for sublist in this_si_contigs.values() for item in sublist[0]]

    def __intracore_sorting__(self):
        # Sort each unaligned region by flanking blocks (core block ID or 5"/3" ctg end)
        useq_regions = {}
        for fragments in self.input_contigs.values():
            for fragment in fragments:
                key = tuple(sorted([fragment[8], fragment[9]]))
                try:
                    useq_regions[key].append(fragment)
                except KeyError:
                    useq_regions[key] = [fragment]
        return useq_regions

    def __generate_useq_file__(self, useq_regions, reference_intracore_blocks_order, non_reference_intracore_blocks):
        # Load all input contig sequences into memory. This fastens the creation of the unalign file.
        input_sequences = {}
        for si in range(self.min_si, self.max_si + 1):
            input_file_path = self.seq_db[si]
            input_sequences[si] = {rec.id: rec.seq for rec in SeqIO.parse(input_file_path, "fasta")}
        # Define output file name
        useq_p = path.join(self.out_path, "{0}.unalign".format(self.prefix))
        # Generate useq file - Write the intracluster regions appearing in the reference first
        with open(useq_p, "w") as useq_f:
            for ic_block_id in reference_intracore_blocks_order:
                # Sort occurences by si
                this_ic_occurences = sorted(useq_regions[ic_block_id], key=lambda x: x[0])
                # Write IC region for each si (except for those cases where SI does not have this region)
                for occ in this_ic_occurences:
                    fragment_seq =\
                        input_sequences[occ[0]][self.input_contigs_names[(occ[0], occ[1])]][occ[4]:occ[5] + 1]
                    # Do not write zero length fragments
                    if not fragment_seq:
                        continue
                    useq_f.write(">{0}:{1}-{2} + intCls{3}.{4} s{5}:p{6}\n".format(occ[0], occ[2], occ[3],
                                                                                   ic_block_id[0], ic_block_id[1],
                                                                                   occ[1], occ[4]))
                    useq_f.write(str(fragment_seq) + "\n")
                useq_f.write("=\n")
        # Next, add all the intracore blocks to the .unalign file that do not appear in the reference
        with open(useq_p, "a") as useq_f:
            for ic_block_id in non_reference_intracore_blocks:
                # Sort occurences by si
                this_ic_occurences = sorted(useq_regions[ic_block_id], key=lambda x: x[0])
                # Check if a sequence was added for this IC region
                # Sometimes seqs are of length zero and are not added
                added_seq = False
                for occ in this_ic_occurences:
                    fragment_seq = \
                        input_sequences[occ[0]][self.input_contigs_names[(occ[0], occ[1])]][occ[4]:occ[5] + 1]
                    # Do not write zero length fragments
                    if not fragment_seq:
                        continue
                    useq_f.write(">{0}:{1}-{2} + intCls{3}.{4} s{5}:p{6}\n".format(occ[0], occ[2], occ[3],
                                                                                   ic_block_id[0], ic_block_id[1],
                                                                                   occ[1], occ[4]))
                    useq_f.write(str(fragment_seq) + "\n")
                    added_seq = True
                if added_seq:
                    useq_f.write("=\n")

    def __generate_ic_file__(self, useq_regions, clusters_sizes, left_end_cluster, right_end_cluster):
        # Define output file name
        ic_p = path.join(self.out_path, "{0}.ic.csv".format(self.prefix))
        with open(ic_p, "w") as csv_f:
            csv_f.write("SI,CB.L,CB.R,LEN.CB.L,LEN.CB.R,LEN.IC,TYPE,SR\n")
            # for CBL,CBR in zip(range(left_end_cluster+1, right_end_cluster-1), range(left_end_cluster+2, right_end_cluster)):
            for CBL, CBR in useq_regions.keys():
                # Sort occurences by si
                try:
                    this_ic_occurences = sorted(useq_regions[(CBL, CBR)], key=lambda x: x[0])
                except KeyError:
                    this_ic_occurences = []
                # First, test whether this IC region is present on the reference
                if this_ic_occurences and this_ic_occurences[0][0] == 1:
                    this_ic_region_lenght_on_ref = this_ic_occurences[0][6]
                    csv_f.write("{0},{1},{2},{3},{4},{5},{6},{7}\n".format(1, CBL, CBR, clusters_sizes.get((1, CBL), 0),
                                                                           clusters_sizes.get((1, CBR), 0),
                                                                           this_ic_region_lenght_on_ref, "P", 1))
                    this_ic_occurences.pop(0)
                    # Now check for the presence of this IC for each other SI
                    # If it is present, write occurence to file and remove the occurence from the ic_occurences stack
                    for si in range(self.min_si + 1, self.max_si + 1):
                        if this_ic_occurences and this_ic_occurences[0][0] == si:
                            # Have a min ic region length of 1 to avoid division by zero
                            s_r_len_ratio = this_ic_occurences[0][6] / max(this_ic_region_lenght_on_ref, 1)
                            csv_f.write("{0},{1},{2},{3},{4},{5},{6},{7}\n".format(si, CBL, CBR,
                                                                                   clusters_sizes.get((si, CBL), 0),
                                                                                   clusters_sizes.get((si, CBR), 0),
                                                                                   this_ic_occurences[0][6], "P",
                                                                                   s_r_len_ratio))
                            this_ic_occurences.pop(0)
                        elif this_ic_occurences and this_ic_occurences[0][0] != si:
                            # The transition is present in the reference but not in this si
                            # Check if the transition is blocked because the sequence is fragmented into two
                            # contigs at this position
                            contig_end_transitions = []
                            if CBL != right_end_cluster:
                                try:
                                    contig_end_transitions.extend([useq_region_tuple for useq_region_tuple in
                                                                   useq_regions[(left_end_cluster, CBL)]
                                                                   if useq_region_tuple[0] == si])
                                except KeyError:
                                    pass
                            if CBR != right_end_cluster:
                                try:
                                    contig_end_transitions.extend([useq_region_tuple for useq_region_tuple in
                                                                   useq_regions[(left_end_cluster, CBR)]
                                                                   if useq_region_tuple[0] == si])
                                except KeyError:
                                    pass
                            if CBL != left_end_cluster:
                                try:
                                    contig_end_transitions.extend([useq_region_tuple for useq_region_tuple in
                                                                   useq_regions[(CBL, right_end_cluster)]
                                                                   if useq_region_tuple[0] == si])
                                except KeyError:
                                    pass
                            if CBR != left_end_cluster:
                                try:
                                    contig_end_transitions.extend([useq_region_tuple for useq_region_tuple in
                                                                   useq_regions[(CBR, right_end_cluster)]
                                                                   if useq_region_tuple[0] == si])
                                except KeyError:
                                    pass
                            # Remove duplicate elements from contig_end_transitions:
                            # contig_end_transitions = list(set(contig_end_transitions))
                            if len(contig_end_transitions) > 4:
                                pass
                            if contig_end_transitions:
                                for transition in contig_end_transitions:
                                    csv_f.write(
                                        "{0},{1},{2},{3},{4},{5},{6},{7}\n".format(transition[0], transition[-2],
                                                                                   transition[-1],
                                                                                   clusters_sizes.get(
                                                                                       (transition[0], transition[-2]),
                                                                                       0),
                                                                                   clusters_sizes.get(
                                                                                       (transition[0], transition[-1]),
                                                                                       0),
                                                                                   transition[-4], "LE" if
                                                                                   transition[
                                                                                       -2] == left_end_cluster else "RE",
                                                                                   0))
            # Even if one of the CB has a contig end as a neighbor, this does not mean that this transition
            # is because of a genome fragmentation at this point. Example:
            # Reference: [CB9]-[CB10]-[CB11]
            # Sample: [CB9]-[CB10]-[CB20] and [CB11]-[RE]
            # This is still a deletion for [CB10]-[CB11] because [CB10] does not have a contig end as neighbor
            # [[si, contig_id, global_start, global_stop, contig_start, contig_stop, contig_length,
            #                   strand, left_coreblock, right_coreblock]
                            contig_end_transitions = [[item[-1], item[-2]] for item in contig_end_transitions]
                            contig_end_transitions = [item for sublist in contig_end_transitions for item in sublist]
                            if not (CBL in contig_end_transitions and CBR in contig_end_transitions):
                                csv_f.write("{0},{1},{2},{3},{4},{5},{6},{7}\n".format(si, CBL, CBR,
                                                                                       clusters_sizes.get((si, CBL), 0),
                                                                                       clusters_sizes.get((si, CBR), 0),
                                                                                       0, "D", 0))

                else:
                    # This IC region is not present on the reference
                    this_ic_region_lenght_on_ref = 0
                    csv_f.write("{0},{1},{2},{3},{4},{5},{6},{7}\n".format(1, CBL, CBR, clusters_sizes.get((1, CBL), 0),
                                                                           clusters_sizes.get((1, CBR), 0),
                                                                           this_ic_region_lenght_on_ref, "NP", 0))
                    for si in range(self.min_si + 1, self.max_si + 1):
                        if this_ic_occurences and this_ic_occurences[0][0] == si:
                            csv_f.write("{0},{1},{2},{3},{4},{5},{6},{7}\n".format(si, CBL, CBR,
                                                                                   clusters_sizes.get((si, CBL), 0),
                                                                                   clusters_sizes.get((si, CBR), 0),
                                                                                   this_ic_occurences[0][6], "N", 0))
                            this_ic_occurences.pop(0)
                        else:
                            csv_f.write("{0},{1},{2},{3},{4},{5},{6},{7}\n".format(si, CBL, CBR,
                                                                                   clusters_sizes.get((si, CBL), 0),
                                                                                   clusters_sizes.get((si, CBR), 0), 0,
                                                                                   "NP", 0))
    def generate_useqs(self, generate_useq_file, generate_ic_file):
        # First, find the smallest and largest cluster nr.
        # min-1 and max+1 become "virtual" clusters signaling the end of a contig
        start = dt.now()
        cluster_ids = [item[4] for item in self.xmfa_header.values()]
        left_end_cluster = min(cluster_ids) - 1
        right_end_cluster = max(cluster_ids) + 1
        # Update this information in input_contigs
        for val in self.input_contigs.values():
            for contig in val:
                contig[-2] = left_end_cluster
                contig[-1] = right_end_cluster
        # Extract core locations from xmfa header, store by seqid
        # Store in this format:
        # clusters[si] = [[si, contig_id, global_start, global_stop, contig_start, contig_stop, fragment_length,
        #                   strand, ClusterID]]
        clusters = {si:[] for si in range(self.min_si, self.max_si+1)}
        for header in self.xmfa_header.values():
            global_start = header[1]
            global_stop = header[2]
            contig_id = header[5]
            contig_start = header[6]
            seq_id = header[0]
            cluster_id = header[4]
            strand = header[3]
            fragment_length = abs(global_stop - global_start) + 1
            contig_stop = contig_start + fragment_length - 1
            clusters[seq_id].append((seq_id, contig_id, global_start, global_stop,
                                     contig_start, contig_stop, fragment_length, strand, cluster_id))
        # Create another dict just for cluster sizes
        # Format: {(si, clstrID):cluster_size}
        clusters_sizes = {(key, item[-1]): item[-3] for key, val in clusters.items() for item in val}

        # Some contigs may be without core block. Find and report this unaligned contigs
        aligned_contigs = set([(ctg[0], ctg[1]) for val in clusters.values() for ctg in val])
        all_contigs = set([(ctg[0], ctg[1]) for val in self.input_contigs.values() for ctg in val])
        unaligned_contigs = all_contigs - aligned_contigs

        # Convert unaligned_contig list into a dict, si are the keys
        unaligned_contig_dict = {si: [] for si in range(self.min_si, self.max_si + 1)}
        [unaligned_contig_dict[ctg[0]].append(ctg[1]) for ctg in unaligned_contigs]
        time_elapsed = (dt.now() - start).total_seconds()
        self.write_log("Identified {0} unaligned contigs in {1} seconds.\n".format(len(unaligned_contigs),
                                                                                   time_elapsed))
        # Fragment input contigs into unaligned fragments with the previously identified core regions
        start = dt.now()
        self.__fragment_contigs(clusters)
        time_elapsed = (dt.now() - start).total_seconds()
        self.write_log("Calculated the indices of unaligned regions in {0} seconds.\n".format(time_elapsed))
        # Sort fragments in input_contigs by intracore ID
        start = dt.now()
        useq_regions = self.__intracore_sorting__()

        # Find "reference" order of unaligned blocks, as given by the reference genome
        # In many cases there is a linear range of core region ids, starting at 1 till max(core region id)
        # But sometimes it doesn't start at one and core regions in between are missing
        reference_intracore_blocks_order = [tuple(sorted((item[8], item[9]))) for item in sorted(self.input_contigs[1],
                                                                                                 key=lambda x: x[2])]

        # Some flanking orientations may not occur in the reference  but in some or all of the assemblies
        all_intracore_blocks = [tuple(sorted((frag[8], frag[9]))) for fragments in self.input_contigs.values() for frag in
                                fragments]
        non_reference_intracore_blocks = set(all_intracore_blocks) - set(reference_intracore_blocks_order)
        # sort non_reference_intracore_blocks by the number of their occurences
        non_reference_intracore_blocks = sorted(list(non_reference_intracore_blocks),
                                                key=lambda x: -len(useq_regions[x]))

        time_elapsed = (dt.now() - start).total_seconds()
        self.write_log("Sorted the unaligned regions by intra-core blocks in {0} seconds.\n".format(time_elapsed))
        # Create file with all unaligned sequences, sorted by intracore region?
        if generate_useq_file:
            start = dt.now()
            self.__generate_useq_file__(useq_regions, reference_intracore_blocks_order, non_reference_intracore_blocks)
            time_elapsed = (dt.now() - start).total_seconds()
            self.write_log("Wrote unaligned sequence output in {0} seconds.\n".format(time_elapsed))
        # Generate the IC region stat file?
        if generate_ic_file:
            start = dt.now()
            self.__generate_ic_file__(useq_regions, clusters_sizes, left_end_cluster, right_end_cluster)
            time_elapsed = (dt.now() - start).total_seconds()
            self.write_log("Wrote intra-core region stat file in {0} seconds.\n".format(time_elapsed))