# Pansnp module: Build an "overall" alignment out of all the individual alignments
# created at each subset core-analysis step.
# Pansnp employs a hierarchical approach on pan-genome computation. The same hierarchical structure
# is employed to stepwise build the global alignment file.
# Each sequence subset has n parent nodes. When adding a new alignment to the global alignment, we
# need to check the new alignment for overlaps with alignments in the parent node. Although the order
# in which we check is irrelevant, we do this along the path from the root to the new alignment node.
# If there is an overlap for a cluster in the new alignment node, we split those cluster into two fragments,
# one to the left and one to the right of the overlap. Sequence indices are adjusted after each split.
#
# This module has additional functionality to write down all the sequence fragments that are not part of
# any alignment.
from os import path
from os import makedirs
from Bio import SeqIO


# class ClusterObject represents a single cluster from a ParCore core alignment
# Each instance has one sequence and one set of indices per assembly. Assemblies
# are identified by the filename of their input file. Each assembly is represented
# exactly once. Multiple sequences from the same assembly are not possible (and anyway
# not supported by ParCore).
class ClusterObject:
    def __init__(self, cluster_id):
        self.cluster_id = cluster_id
        # Store sequence fragments in a dict
        self.sequences = {}
        # Store indices in a dict. Format:
        # (global_start, global_stop, strand, contigID, contig_start)
        self.indices = {}
        # Flag that indicates a zero-length cluster
        self.zero_length = False

    # Add a sequence fragment for one assembly. If there is already a
    # sequence stored for that assembly, the new sequence string will be appended
    # To replace a sequence, use replace_sequence()
    def add_sequence(self, filename, seq):
        # Append sequence if there is already a sequence string for 'filename"
        # If not, add new sequence as new element to the dict
        try:
            self.sequences[filename] += seq.strip()
        except KeyError:
            self.sequences[filename] = seq.strip()

    # Replace sequence string for 'filename'
    def replace_sequence(self, filename, seq):
        self.sequences[filename] = seq.strip()

    # Add/replace index for a 'filename'
    def add_index(self, filename, index):
        self.indices[filename] = index

    # Getter/setter
    def get_id(self):
        return self.cluster_id

    def set_id(self, new_id):
        self.cluster_id = new_id

    def get_indices(self):
        return self.indices

    def is_zero_length(self):
        return self.zero_length

    # Write this cluster to file
    # Order of sequences is defined in fn_order
    # File object is file_obj
    def write_alignment(self, fn_order, file_obj):
        for fn in fn_order:
            # Retrieve indices and sequence for this filename
            index = self.indices[fn]
            seq = self.sequences[fn]
            # Write header
            file_obj.write(">{0}:{1}-{2} Cluster{3} {4} s{5}:p{6}\n".format(fn, index[0], index[1],
                                                                            self.cluster_id, index[2],
                                                                            index[3], index[4]))
            # Write sequence in formatted style:
            # Linebreak every 80 nt
            seq_pointer = 0
            while seq_pointer < len(seq):
                file_obj.write(seq[seq_pointer:seq_pointer+80] + "\n")
                seq_pointer += 80

    # Split this cluster into two cluster based on the indices of a second cluster, parent_cluster,
    # that overlaps, at least some parts, of this cluster. The cluster will be split into a
    # fragment located downstream of parent_cluster (left cluster) and a fragment
    # located upstream of parent_cluster (right cluster). This cluster becomes the left fragment,
    # the function returns the right fragment as new cluster.
    # The left fragment may have size < 1 afterwards, i.e. global start > global stop
    # If right fragment has size <1, a None object is returned
    # The split is done based on the indices stored in this cluster and in parent_cluster for one
    # assembly, identified by file_name. For all other assemblies in this cluster the sequence string
    # is cut to the same length as for the sequence string associated with file_name. The indices for those
    # other sequences are adjusted based on the number of nucleotides removed at the 5' and 3' end,
    # not counting gaps.
    def split(self, parent_cluster, file_name):
        # First, retrieve indices for file_name for parent and child cluster
        index_parent = parent_cluster.get_indices()[file_name]
        index_child = self.indices[file_name]
        # Get a list of all the file names in this cluster, the 'child' cluster, that are unequal to file_name
        other_file_names = [key for key in self.indices.keys() if key != file_name]
        # Get global start, stop indices for parent and child cluster
        # Global refers to an index-system that consideres all contigs in the input file to be merged
        # together into one large sequence.
        parent_global_start, parent_global_stop = index_parent[:2]
        child_global_start, child_global_stop = index_child[:2]
        # Left and right fragment could disappear if the parent global start or stop is positioned beyond
        # the childs global start or stop
        # Define zero to be the minimum fragment length, otherwise we could ran into some trouble with
        # Pythons slicing scheme.
        # The left fragment is on the 5' side of the overlap region and the right fragment
        # is on the 3' side of the overlap region for strands in '+' orientation. For strands in '-'
        # orientation we need to swap this around, so that the left fragment comes from the 3' side and
        # vice versa.
        if index_child[2] == "+":
            left_fragment_length = max(0, parent_global_start-child_global_start)
            right_fragment_length = max(0, child_global_stop - parent_global_stop)
        else:
            left_fragment_length = max(0, child_global_stop - parent_global_stop)
            right_fragment_length = max(0, parent_global_start -  child_global_start)
        # Fragment length refers to the ungapped sequence
        # Calculate the actual position on the sequence where to make the cut
        # This is independent of the strand orientation.
        child_sequence = self.sequences[file_name]
        # Check if there are gaps in the child sequence, i.e. the actual sequence is longer than
        # indicated by the global start/stop of the child sequence
        gaps = (len(child_sequence) > child_global_stop-child_global_start+1)
        # The left sequence ends at left_fragment_length -1 , so we cut to the right of
        # left_sequence_cut_pos (Python has zero-centric coordinates for lists).
        left_sequence_cut_pos = left_fragment_length-1

        # If there are gaps, we have to move the sequence cut position to the right
        # until the length of the sequence (with gaps removed) matches the targeted length
        if gaps:
            while len(child_sequence[:left_sequence_cut_pos+1].replace("-","")) < left_fragment_length:
                left_sequence_cut_pos += 1
        # For the right sequence we cut to the left of right_sequence_cut_pos
        right_sequence_cut_pos = len(child_sequence)-right_fragment_length
        if gaps:
            while len(child_sequence[right_sequence_cut_pos:].replace("-","")) < right_fragment_length:
                right_sequence_cut_pos -= 1

        # Create a new cluster, this will be the part of this cluster to the right of the parental overlap
        # Get the index vector for the right cluster. Format:
        # (global_start, global_stop, strand, contigID, contig_start)
        # Get orientation of right fragment
        r_strand = index_child[2]
        # Define global start/stop of right fragment for '+' or '-' orientation
        if r_strand == "+":
            r_global_start = parent_global_stop + 1
            r_global_stop = child_global_stop
        else:
            r_global_start = child_global_start
            r_global_stop = parent_global_start - 1
        # Contig ID remains unchanged
        r_contigID = index_child[3]
        # On-contig start position is easily adjusted for the changed global start:
        # Just calculate the offset from the old and new global start position
        r_contig_start = r_global_start-index_child[0]+index_child[4]
        r_index = (r_global_start, r_global_stop, r_strand, r_contigID, r_contig_start)
        # Only create a cluster when it will be > 0 nt
        right_cluster = None
        if r_global_start <= r_global_stop:
            # Get the sequence for the right cluster
            r_sequence = child_sequence[right_sequence_cut_pos:]
            # Create a unqiue ID for this new cluster
            # The actual ID is not important as it is only used internally as a key
            # But it is important that this ID is unique to prevent accidental overwriting
            right_cluster = ClusterObject(str(self.cluster_id) + "R" + str(r_contigID) + str(r_global_start)
                                          + str(r_global_stop))
            right_cluster.add_index(file_name, r_index)
            right_cluster.replace_sequence(file_name, r_sequence)
            # So far, this cluster only includes the sequence for file_name. We now need to include
            # all the other file_names present in this child cluster: We know how long the sequence needs
            # to be, so we cut each sequence at right_sequence_cut_pos. We then need to adjust the
            # indices associated with this other_file_name (ofn). Global_stop remains the same
            # but we need to adjust global_start by what we cut away. Important: gap positions in the
            # removed sequence portion do not count for the index adjusting.
            for ofn in other_file_names:
                other_index = self.indices[ofn]
                other_sequence = self.sequences[ofn]
                other_r_strand = other_index[2]
                # Cut the sequence to the left of right_sequence_cut_pos
                other_r_sequence = other_sequence[right_sequence_cut_pos:]
                # Calculate the actual number of nucleotides cut away
                removed_nucleotides_count = len(other_sequence[:right_sequence_cut_pos].replace("-", ""))
                # Adjust the global start by the actual number of nucleotides included in this fragment
                # For '-' strands: adjust global stop
                if other_r_strand == '+':
                    other_r_global_start = other_index[0]+removed_nucleotides_count
                    # Global stop remains unaffected
                    other_r_global_stop = other_index[1]
                else:
                    other_r_global_stop = other_index[1] - removed_nucleotides_count
                    # Global start remains unaffected
                    other_r_global_start = other_index[0]
                # If we removed everything from the child_sequence, the global stop will
                # become smaller than the global_start
                # This zero-length cluster is then of no use, so we will discard it
                if other_r_global_stop < other_r_global_start:
                    right_cluster = None
                    break
                other_r_contigID = other_index[3]
                other_r_contig_start = other_index[4] + other_r_global_start - other_index[0]
                other_r_index = (other_r_global_start, other_r_global_stop, other_r_strand,
                                 other_r_contigID, other_r_contig_start)
                # Add modified index and sequence for this other_file_name to the cluster
                # Before continuing with the next _other_file_name
                right_cluster.add_index(ofn, other_r_index)
                right_cluster.replace_sequence(ofn, other_r_sequence)

        # Modify the indices of this cluster, which will become the left-side fragment
        # Again, the strand orientation is important
        l_strand = index_child[2]
        if l_strand == "+":
            # Global start positions remain unchanged, only the global stop changes
            # as the left fragment stops right before the start of the overlap region
            l_global_start = child_global_start
            l_global_stop = parent_global_start-1
        else:
            # Global stop remains unchanged, global start is adjusted as the left fragment starts
            # next to the overlap region
            l_global_start = parent_global_stop + 1
            l_global_stop = child_global_stop
        l_contigID = index_child[3]
        # Adjust on-contig start
        l_contig_start = l_global_start-index_child[0] + index_child[4]
        l_index = (l_global_start, l_global_stop, l_strand, l_contigID, l_contig_start)
        # if the left-side fragment becomes size <= 0, stop here as this cluster will be removed soon
        if l_global_start > l_global_stop:
            self.zero_length = True
        # Otherwise, adjust the sequence string and indices for the left side cluster for file_name
        else:
            # Get the sequence for the left cluster
            l_sequence = child_sequence[:left_sequence_cut_pos+1]
            self.add_index(file_name, l_index)
            self.replace_sequence(file_name, l_sequence)
            # As for the right_cluster, we need to reduce the sequence for each other file_name
            # to the length of l_sequence and subsequently adjust their indices
            for ofn in other_file_names:
                other_index = self.indices[ofn]
                other_sequence = self.sequences[ofn]
                other_l_strand = other_index[2]
                # Cut the sequence to the right of left_sequence_cut_pos
                other_l_sequence = other_sequence[:left_sequence_cut_pos+1]
                # Calculate the actual number of nucleotides cut away
                removed_nucleotides_count = len(other_sequence[left_sequence_cut_pos+1:].replace("-", ""))
                if other_l_strand == '+':
                    # Global start remains unaffected
                    other_l_global_start = other_index[0]
                    # Adjust global_stop by subtracting all removed nucleotides
                    other_l_global_stop = other_index[1] - removed_nucleotides_count
                else:
                    # Global stop remains unaffected
                    other_l_global_stop = other_index[1]
                    # Adjust global_start by subtracting all removed nucleotides
                    other_l_global_start = other_index[0] + removed_nucleotides_count
                # Check that there are some nucleotides left in this fragment
                if other_l_global_start > other_l_global_stop:
                    self.zero_length = True
                    break
                other_l_contigID = other_index[3]
                other_l_contig_start = other_l_global_start-other_index[0] + other_index[4]
                other_l_index = (other_l_global_start, other_l_global_stop, other_l_strand,
                                 other_l_contigID, other_l_contig_start)
                self.add_index(ofn, other_l_index)
                self.replace_sequence(ofn, other_l_sequence)
        # Finally, return the new right_cluster
        return right_cluster


# Represent an entire alignment (i.e. a node in the hierarchical PanCore graph)
# as a collection of ClusterObjects
class AlignmentNode:
    def __init__(self, name):
        self.name = name
        # Collect all clusters in a dict
        self.clusters = {}
        # Remember lowest and highest cluster ID
        self.min_clusterID = -1
        self.max_clusterID = -1

    # Getter/setter
    def add_clstr(self, clstr):
        self.clusters[clstr.get_id()] = clstr

    def get_clstr(self):
        return self.clusters

    def rename_cluster(self, old_name, new_name):
        renamed_cluster = self.clusters.pop(old_name)
        renamed_cluster.set_id(new_name)
        self.clusters[new_name] = renamed_cluster

    # Return all file-names represented in this node
    def get_filenames(self):
        # Pick any cluster since all clusters have the same file names
        return list(list(self.clusters.values())[0].get_indices().keys())

    # Remove all clusters from this node that are of zero length
    def remove_zero_length_clstr(self):
        self.clusters = {clstr.get_id(): clstr for clstr in self.clusters.values() if not clstr.is_zero_length()}

    # Renumerate clusters starting with next_clstr_id
    # Sort cluster by global start position in ref_fn file
    # Returns last clstr_id assigned
    def renumerate_clusters(self, next_clstr_id, ref_fn):
        # Sort all clusters by global start of the reference
        clstr_list = sorted(self.clusters.values(), key=lambda x: x.get_indices()[ref_fn][0])
        # Clear current dict of clusters (clusters are still saved in clstr_list)
        self.clusters.clear()
        self.min_clusterID = next_clstr_id
        # iterate through sorted cluster list, add clusters together with their new ID to cluster dict
        for clstr in clstr_list:
            clstr.set_id(next_clstr_id)
            self.clusters[next_clstr_id] = clstr
            next_clstr_id += 1
        self.max_clusterID = next_clstr_id - 1
        return next_clstr_id

    # Write alignments for all clusters in this node to file
    # The order of sequences in these alignment blocks is given by filename_order
    def write_alignments(self, filename_order, out_file_obj):
        # Filter filename_order to include only those file names actually present in this node
        filename_order = [fn for fn in filename_order if fn in self.get_filenames()]
        # Iterate through all clusters
        for clusterID in range(self.min_clusterID, self.max_clusterID + 1):
            self.clusters[clusterID].write_alignment(filename_order, out_file_obj)
            # Write a '=' to signal end of this cluster alignment
            out_file_obj.write("=\n")


# Main class for the construction of combined xmfa alignment files
# and combined unaligned sequences files
class BuildXmfa:
    def __init__(self, out_p, prefix):
        # Create output folder
        if not path.isdir(out_p):
            makedirs(out_p, exist_ok=True)
        self.out_p = out_p
        # Collect all nodes in a dict
        self.al_tree = {}
        # Each PanCore subset alignment output filename starts with a prefix
        self.prefix = prefix
        # Remember name of reference file
        self.ref_file = None

    # Once all subsets are processed by PanCore we can write the combined xmfa core alignment file
    def generate_combined_xmfa(self):
        # First, we define the order in which we will add the alignment nodes to the output file
        # We walk the alignment tree in breadth-first order
        # On each level, sort nodes numerically, 0_0 < 0_1 < 1_0...
        node_labels = sorted([[int(x) for x in item.split("_") if x] for item in self.al_tree.keys()],
                             key=lambda x: (len(x), x))
        # Reformat into string format
        node_labels = ["_".join([str(x) for x in id_list]) for id_list in node_labels]
        node_labels = ["_" + label if label else label for label in node_labels]
        # Rename clusters -> Assign each a new ID, starting with zero
        # renumerate cluster returns the next new_clstr_id, so we can simple reassign
        # that variable
        new_clstr_id = 0
        for node_label in node_labels:
            a_node = self.al_tree[node_label]
            new_clstr_id = a_node.renumerate_clusters(new_clstr_id, self.ref_file)
        # When writing the alignments to file, we need to have a defined order of file_names
        # Only the root node has all file_names
        # Take the fn from the root-node, sort alphabetically
        # but put the reference in first place
        root_node = self.al_tree[node_labels[0]]
        root_files = sorted(root_node.get_filenames())
        root_files.remove(self.ref_file)
        root_files.insert(0, self.ref_file)
        # Now walk along the tree, for each node trigger the writing of all alignments in that node to out_f
        with open(path.join(self.out_p, "{0}_combined.xmfa".format(self.prefix)), "w") as out_f:
            for node_label in node_labels:
                self.al_tree[node_label].write_alignments(root_files, out_f)


    # Extract all unaligned sequences for one input assembly file
    # Needs to be finished in a future version
    def generate_useq(self, input_f_path):
        # Extract filename from input_f_path
        filename = path.split(input_f_path)[-1]
        # First, find all nodes in the tree where this filename appears
        nodes = [node for node in self.al_tree.values() if filename in node.get_filenames()]
        # For each node, get all the indices for this filename
        # Since indices do not overlap we can simple collect all indices from all nodes and combine them
        clusters = [list(node.get_clstr().values()) for node in nodes]
        this_fn_indices = [item.get_indices()[filename] for clusterlist in clusters for item in clusterlist]
        print(clusters)
        # Parse the input file
        input_seqs_records = SeqIO.parse(input_f_path, "fasta")
        input_seqs_dict = {}
        # Store sequences in a dict, key is the contig nr, starting with one
        for id, rec in enumerate(input_seqs_records, 1):
            seq = rec.seq
            input_seqs_dict[id] = (len(seq), seq)

    # Add a new xmfa file to the alignment node tree
    # Each PanCore subset creates one alignment xmfa file
    # Here, this xmfa file is converted into a alignment node
    # Due to the hierarchical approach of PanCore, each node has zero or more
    # parent nodes. Nucleotides contained in a core block in any parent node cannot be
    # in a core block in this child node. The new alignments are therefore checked against all parental
    # alignments to remove regions covered in parental alignments.
    def add_xmfa(self, xmfa_p):
        # Create a new alignment node for this xmfa file
        # The label is the filename:
        # pancore_0_1_0 -> _0_1_0 is the label
        # For the root, the label is an empty string ''
        node_label = path.basename(xmfa_p)[len(self.prefix):].rsplit(".")[0]
        # Create an, so far empty, new alignment node
        new_al_node = AlignmentNode(node_label)
        # First, map the sequence index (SI) from the xmfa  to the file name
        # SI change with each new alignment, so the file name is
        # the only unique identifier we have
        si_fn_dict = {}
        si = None
        with open(xmfa_p, "r") as xmfa_f:
            for line in xmfa_f:
                if not line.startswith("#"):
                    break
                if line.startswith("##SequenceIndex"):
                    si = int(line.split()[1])
                if line.startswith('##SequenceFile'):
                    si_fn_dict[si] = line.split()[1]
        # The reference file used by PanCore always has SI == 1
        self.ref_file = si_fn_dict[1]
        # Convert xmfa file into cluster objects and add each
        # Cluster object to the new alignment node
        # Iterate only once through xmfa file. The start of a new cluster
        # is detected by a '='
        next_cluster = None
        # Keep track of the filename a sequence string belongs to
        # We need this information to add the string to the correct assembly in the cluster object
        current_fn = None
        with open(xmfa_p, "r") as xmfa_f:
            for line in xmfa_f:
                if line.startswith("#"):
                    continue
                if line.startswith("="):
                    new_al_node.add_clstr(next_cluster)
                    next_cluster = None
                    continue
                if line.startswith(">"):
                    data = line.split()
                    cluster_id = int(data[2].replace("cluster", ""))
                    global_start = int(data[0][data[0].index(":") + 1:data[0].index("-")])
                    global_stop = int(data[0][data[0].index("-") + 1:])
                    seq_id = int(data[0][1:data[0].index(":")])
                    strand = data[1]
                    contig_id, contig_start = [int(item[1:]) for item in data[3].split(":")]
                    if not next_cluster:
                        next_cluster = ClusterObject(cluster_id)
                    current_fn = si_fn_dict[seq_id]
                    # (global_start, global_stop, strand, contigID, contig_start)
                    index = (global_start, global_stop, strand, contig_id, contig_start)
                    next_cluster.add_index(current_fn, index)
                else:
                    if next_cluster and current_fn:
                        next_cluster.add_sequence(current_fn, line)
        # Add this cluster node to the tree
        self.al_tree[node_label] = new_al_node
        # If this is the xmfa file for the root node, i.e. the first alignment performed,
        # we are finished here
        if not node_label:
            return 0
        # Next, we need to substract everything from this alignment node that is already covered
        # by all parent nodes
        # For each cluster, we check whether it is partially or totally covered by any cluster in
        # any parent node. We walk along the tree towards the new alignment node, comparing only
        # against those nodes on our path. The path is described by the node label. For example:
        # pansnp_0_1_0_1.xmfa -> from root to child node 0, to child node 1, to child node zero. This
        # alignment node is child node one to the last child node.
        # Each sequence comes from an assembly, represented here as a file name
        # There are always less file names in a child node compared to the parent node
        # For each file name, we check whether its sequences in this node are already contained
        # in the parent node. In a perfect world, it should be sufficient to check only one file name,
        # e.g. the reference. This would only work if each nucleotide always aligns to the same nucleotie
        # in the reference. But due to the progressive nature of this hierarchical approach, we cannot
        # take this for granted. We therefore have to check every sequence from each file name for
        # presence in the parent. Ideally, this should not be much computational effort. For each file name,
        # if there is an overlap with the parent node, all sequences from all files are trimmed. After the first
        # round of trimming for the first file name, not so much ( or at best no work) should be left to do
        # for all other file names.
        # Get the path along the tree we walk from the root to the new nie
        tree_path = []
        # Make a copy of the node_label
        new_node_label = node_label
        while new_node_label:
            new_node_label = new_node_label.rsplit("_", maxsplit=1)[0]
            tree_path.append(new_node_label)
        tree_path.reverse()
        # Walk along the tree
        for parent_label in tree_path:
            try:
                parent_node = self.al_tree[parent_label]
            except KeyError:
                # If we can't find the node, some alignment files are missing
                # return error code
                return -1
            # Get the child_node back out of the tree
            child_node = self.al_tree[node_label]
            # Get all the clusters for this child alignment node
            # and parent node
            child_clstr = list(child_node.get_clstr().values())
            parent_clstr = list(parent_node.get_clstr().values())
            # Get a list of all assemblies (i.e. file names) involved in this child node
            # If there are no assemblies in this child node left, our work is done here
            # All of the sequences covered by the child alignment are seen in
            # a parent alignment before
            try:
                child_filenames = list(list(child_clstr)[0].get_indices().keys())
            except IndexError:
                return
            # For each filename, check for overlaps between the child and the parent clusters
            # If there is an overlap, split the sequences in the separate cluster
            # Ensure that the reference file name is in first position
            child_filenames.remove(self.ref_file)
            child_filenames.insert(0, self.ref_file)
            for fn in child_filenames:
                # Get the list of cluster in the child node
                # We already got this list above to extract the file names
                # But we need to retrieve this list at every iteration
                # because clusters may have been removed during the previous loop
                child_clstr = list(child_node.get_clstr().values())
                # Get the indices for each cluster
                # First, we need to do some reformatting
                # Format:
                # [(ClusterID0, [(filename1, index1), (filename2, index2),...]),
                # (ClusterID1, [(filename1, index1), (filename2, index2),...]),
                # (ClusterID2, [(filename1, index1), (filename2, index2),...])...]
                child_indices = [(clstr.get_id(), list(clstr.get_indices().items())) for clstr in child_clstr]
                parent_indices = [(clstr.get_id(), list(clstr.get_indices().items())) for clstr in parent_clstr]
                # Reformat:
                # [[ClusterID0, filename1, index1], ([ClusterID0, filename2, index2],[ClusterID1, filename1, index1]..],
                child_indices = [[item[0], *index] for item in child_indices for index in item[1]]
                parent_indices = [[item[0], *index] for item in parent_indices for index in item[1]]
                # Reformat indices into this format and leave only those indices matching current file name in list
                # [contigID, global_start, global_stop, cluster_id]
                child_indices = [(index[2][3], index[2][0], index[2][1], index[0]) for index in child_indices
                                 if index[1] == fn]
                parent_indices = [(index[2][3], index[2][0], index[2][1], index[0]) for index in parent_indices
                                  if index[1] == fn]

                # Finally, sort cluster indices by contigID and global_start
                child_indices = sorted(child_indices, key=lambda x: (x[0], x[1]))
                parent_indices = sorted(parent_indices, key=lambda x: (x[0], x[1]))

                # Now walk through the list of parent and child indices until we reach the end of
                # the parent or child indices
                while parent_indices and child_indices:
                    # Get the first element of parent and child list
                    child_index = child_indices[0]
                    parent_index = parent_indices[0]
                    # If the first child_index contig ID is lower than the first parent_index contig_id,
                    # move along the child indices chain
                    if child_index[0] < parent_index[0]:
                        child_indices.pop(0)
                    # Opposite case: first parent_index contig ID is lower than the first child_index contig_id
                    # move along the parent indices chain
                    elif parent_index[0] < child_index[0]:
                        parent_indices.pop(0)
                    # Last case: both child and parent index are from the same contig
                    else:
                        # Does the child fragment end before the parent fragment starts?
                        # If so, walk along the child chain
                        if child_index[2] < parent_index[1]:
                            child_indices.pop(0)
                        # Or does the parent fragment end before the child fragment starts?
                        # If so, walk along the parent chain
                        elif parent_index[2] < child_index[1]:
                            parent_indices.pop(0)
                        # In any other case, there is some overlap between this child and parent fragment
                        # This means that we have to split the child fragment into two fragments,
                        # One fragment to the left of the overlap, one to the right. Obviously, depending on the
                        # type of the overlap, one or both of the fragments may be of size zero
                        else:
                            # First, retrieve the overlapping cluster from the child and parent node
                            child_overlap_cluster = child_node.get_clstr()[child_index[3]]
                            parent_overlap_cluster = parent_node.get_clstr()[parent_index[3]]
                            # Split the child cluster
                            # The child cluster becomes the left split, the function returns the right split
                            # as new cluster
                            right_split = child_overlap_cluster.split(parent_overlap_cluster, fn)
                            # The child indices at pos 0 are not valid anymore due to the splitting of the
                            # associated cluster
                            child_indices.pop(0)
                            # For at maximum one of the fragments of the split of this child cluster
                            # there may be more overlap with other parent clusters
                            # It depends on the strand orientation whether it will be the left or the right
                            # fragment of the split
                            # Try to extract the left and right fragment start positions
                            if child_overlap_cluster.is_zero_length():
                                left_frag_start = None
                            else:
                                left_frag_start = child_overlap_cluster.get_indices()[fn][0]
                            try:
                                right_frag_start = right_split.get_indices()[fn][0]
                            except AttributeError:
                                right_frag_start = None
                            # Now determine whether we should add the left or the right fragment indices
                            # back to the child indices chain. We only need to add the fragment from the 3' side
                            # of the split as there cannot be another overlap between the 5' side fragment
                            # and any other parental cluster due to the sorted nature of the indices lists.
                            add_left_fragment = False
                            add_right_fragment = False
                            if left_frag_start and right_frag_start:
                                if left_frag_start > right_frag_start:
                                    add_left_fragment = True
                                else:
                                    add_right_fragment = True
                            elif left_frag_start:
                                add_left_fragment = True
                            elif right_frag_start:
                                add_right_fragment = True
                            remaining_fragment_index = None
                            if add_left_fragment:
                                left_frag_index = child_overlap_cluster.get_indices()[fn]
                                left_split_id = child_overlap_cluster.get_id()
                                remaining_fragment_index = (left_frag_index[3], left_frag_index[0],
                                                            left_frag_index[1], left_split_id)
                            if add_right_fragment:
                                right_split_index = right_split.get_indices()[fn]
                                right_split_id = right_split.get_id()
                                remaining_fragment_index = (right_split_index[3], right_split_index[0],
                                                            right_split_index[1], right_split_id)
                            # Add the 3' side fragment index to the very start of the child indices list
                            if remaining_fragment_index:
                                child_indices.insert(0, remaining_fragment_index)
                            # If the split created a right_split cluster of size > 0 nt,
                            # add it do this child_nodes cluster dict
                            if right_split:
                                child_node.add_clstr(right_split)
                # Finally, after scanning through all parent and child indices,
                # remove any cluster from this node that are of zero length
                child_node.remove_zero_length_clstr()





if __name__ == '__main__':
    import sys
    input_xmfas = sys.argv[1:]
    build_xmfa = BuildXmfa("./build_test", "cd50_rerun")
    for ix, input_xmfa in enumerate(input_xmfas):
        print("Adding {0}".format(ix))
        build_xmfa.add_xmfa(input_xmfa)
    build_xmfa.generate_combined_xmfa()
    #build_xmfa.generate_useq("Clostridium_reference_genome_noplasmid.fna")


