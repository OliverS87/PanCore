# Build an xmfa file by stepwise adding new xmfa alignments
# New alignments are checked for overlaps with present alignments
# Overlaps are trimmed so that only novel aligned regions are added to the xmfa file
# Each incoming alignment has a parent alignment. Overlaps are only checked against
# the immediate parent
from os import path
from os import makedirs
# Represent a single cluster in an xmfa file
class ClusterObject:
    def __init__(self, cluster_id):
        self.cluster_id = cluster_id
        self.sequences = {}
        # Indices:
        # (global_start, global_stop, strand, contigID, contig_start)
        self.indices = {}

    def add_sequence(self, filename, seq):
        try:
            self.sequences[filename] += seq.strip()
        except KeyError:
            self.sequences[filename] = seq.strip()

    def add_index(self, filename, index):
        self.indices[filename] = index

    def get_id(self):
        return self.cluster_id

    def get_indices(self):
        return self.indices

    # Split this cluster based on the indices of one file_name
    # into two cluster: this cluster becomes the left fragment, the function
    # returns the right fragment
    # The left fragment may have size < 1 afterwards, i.e. global start > global stop
    # If right fragment has size <1, a None object is returned
    # Each sequence in a ClusterObject has the same length
    # But: There may be gaps that are not reflected by the global start/stop values
    # When adjusting the global start/stop and contig start for removed sequence material,
    # remove gap positions in the sequence from the calculation
    def split(self, parent_cluster, file_name):
        # First, retrieve indices for file_name for parent and child cluster
        index_parent = parent_cluster.get_indices()[file_name]
        index_child = self.indices[file_name]
        # Get a list of all the file names in this cluster that are unequal to file_name
        other_file_names = [key for key in self.indices.keys() if key != file_name]
        # Calculate how long the sequence for the left fragment and right fragment will be
        parent_global_start, parent_global_stop = index_parent[:2]
        child_global_start, child_global_stop = index_child[:2]
        # Left and right fragment could disapear if the parent global start or stop is positioned beyond
        # the childs global start or stop
        # Define zero to be the minimum fragment length
        left_fragment_length = max(0, parent_global_start-child_global_start)
        right_fragment_length = max(0, child_global_stop - parent_global_stop)
        # Sequence length refers to the ungapped sequence
        # Calculate the actual position on the sequence where to make the cut
        child_sequence = self.sequences[file_name]
        # Check if there are gaps in the sequence
        gaps = (len(child_sequence) == abs(child_global_stop-child_global_start)+1)
        left_sequence_cut_pos = left_fragment_length-1
        if gaps:
            while len(child_sequence[:left_sequence_cut_pos+1].replace("-","")) < left_fragment_length:
                left_sequence_cut_pos += 1
        right_sequence_cut_pos = len(child_sequence)-right_fragment_length
        if gaps:
            while len(child_sequence[right_sequence_cut_pos:].replace("-","")) < right_fragment_length:
                right_sequence_cut_pos -= 1
        # Remember how many nucleotides and gaps are going to be cutoff at the left and right
        left_sequence_length = left_fragment_length
        right_sequence_length = len(child_sequence) - right_sequence_cut_pos + 1
        # Create a new cluster, this will be the part of this cluster to the right of the parental overlap
        # Get the index vector for the right cluster
        # (global_start, global_stop, strand, contigID, contig_start)
        r_global_start = parent_global_start + 1
        r_global_stop = child_global_stop
        r_strand = index_child[2]
        r_contigID = index_child[3]
        r_contig_start = r_global_start-index_child[0]+index_child[4]
        r_index = (r_global_start, r_global_stop, r_strand, r_contigID, r_contig_start)
        # Only create a cluster when it will be > 0 nt
        if r_global_start <= r_global_stop:
            # Get the sequence for the right cluster
            r_sequence = child_sequence[len(child_sequence) - right_sequence_length:]
            right_cluster = ClusterObject(self.cluster_id + "R")
            right_cluster.add_index(file_name, r_index)
            right_cluster.add_sequence(file_name, r_sequence)
            # Now add all the other sequences from the other file names to the cluster
            # If for any sequence a cluster size <1 results, stop and forget this cluster
            zero_length_cluster_R = False
            for ofn in other_file_names:
                # (global_start, global_stop, strand, contigID, contig_start)
                other_index = self.indices[ofn]
                other_sequence = self.sequences[ofn]
                # Reduce the sequence to the right fragment sequence length calculated above
                other_r_sequence = other_sequence[len(other_sequence)-right_sequence_length:]
                # Calculate the actual number of nucleotides in this sequence fragment
                other_r_nucleotide_count = len(other_r_sequence.replace("-", ""))
                other_r_global_start = other_index[1]-other_r_nucleotide_count+1
                other_r_global_stop = other_index[1]
                if other_r_global_stop < other_r_global_start:
                    zero_length_cluster_R = True
                    break
                other_r_strand = other_index[2]
                other_r_contigID = other_index[3]
                other_r_contig_start = other_index[4] + other_r_global_start - other_index[0]
                other_r_index = (other_r_global_start, other_r_global_stop, other_r_strand,
                                 other_r_contigID, other_r_contig_start)
                right_cluster.add_index(ofn, other_r_index)
                right_cluster.add_sequence(ofn, other_r_sequence)
        # Modify the indices of this cluster, which will become the left-handed fragment
        l_global_start = child_global_start
        l_global_stop = parent_global_start-1








# Represent an entire alignment (i.e. a node in the hierarchical clustering graph
# as a collection of ClusterObjects
class AlignmentNode:
    def __init__(self, name):
        self.name = name
        self.clusters = {}

    def add_clstr(self, clstr):
        self.clusters[clstr.get_id()] = clstr

    def get_clstr(self):
        return self.clusters

class BuildXmfa:
    def __init__(self, out_p, prefix):
        # Create output folder
        if not path.isdir(out_p):
            makedirs(out_p, exist_ok=True)
        self.out_f = open(path.join(out_p, "{0}_combined.xmfa".format(prefix)), "w")
        # Collect all nodes in a dict
        self.al_tree = {}
        self.prefix = prefix

    def close_file(self):
        self.out_f.close()

    def add_xmfa(self, xmfa_p):
        # Create a new alignment node for this xmfa file
        node_label = path.basename(xmfa_p)[len(self.prefix):].rsplit(".")[0]
        new_al_node = AlignmentNode(node_label)
        # First, map the SI to the file name
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
        # Convert xmfa file into cluster objects and add each
        # Cluster object to the new alignment node
        next_cluster = None
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
                    next_cluster.add_sequence(current_fn, line)
        # Add this cluster node to the tree
        self.al_tree[node_label] = new_al_node
        # If this is the xmfa file for the root node, i.e. the first alignment performed,
        # we are finished here
        if not node_label:
            return 0
        # Next, we need to substract everything from this alignment node that is already covered
        # by the parent node
        # For each cluster, we check whether it is partially or totally covered by a cluster in
        # the parent node. Each sequence comes from an assembly, represented here as a file name
        # There are always less file names in a child node compared to the parent node
        # For each file name, we check whether its sequences in this node are already contained
        # in the parent node. In a perfect world, it should be sufficient to check only one file name,
        # e.g. the reference. This would only work if each nucleotide always aligns to the same nucleotie
        # in the reference. But due to the progressive nature of thus hierarchcial approach, we cannot
        # take this for granted. We therefore have to check every sequence from each file name for
        # presence in the parent. Idealy, this should not be much computational effort. For each file name,
        # if there is an overlap with the reference, all sequences from all files are trimmed. After the first
        # round of trimming for the first file name, not so much ( or at best no work) should be left to do
        # for all other files.
        # First, get the indices for this child node and the parent node
        parent_label = node_label.rsplit("_",maxsplit=1)[0]
        try:
            parent_node = self.al_tree[parent_label]
        except KeyError:
            return 1
        child_node = self.al_tree[node_label]
        # Get all the clusters for this child alignment node
        # and parent node
        child_clstr = list(child_node.get_clstr().values())
        parent_clstr = list(parent_node.get_clstr().values())
        # Get a list of all assemblies (i.e. file names) involved in this child node
        child_filenames = list(child_clstr)[0].get_indices().keys()
        # For each filename, check for overlaps between the child and the parent clusters
        for fn in child_filenames:
            print('Checking {0}'.format(fn))
            # Get the updated list of child nodes
            child_clstr = list(child_node.get_clstr().values())
            # Get the indices for each cluster
            # Format:
            # [(ClusterID0, [(filename1, index1), (filename2, index2),...]),
            # (ClusterID1, [(filename1, index1), (filename2, index2),...]),
            # (ClusterID2, [(filename1, index1), (filename2, index2),...])...]
            child_indices = [(clstr.get_id(), list(clstr.get_indices().items())) for clstr in child_clstr]
            parent_indices = [(clstr.get_id(), list(clstr.get_indices().items())) for clstr in parent_clstr]
            print("A:{0}".format(child_indices[0]))
            # Reformat:
            # [[ClusterID0, filename1, index1], ([ClusterID0, filename2, index2],[ClusterID1, filename1, index1]...],
            child_indices = [[item[0], *index] for item in child_indices for index in item[1]]
            parent_indices = [[item[0], *index] for item in parent_indices for index in item[1]]
            print("B:{0}".format(child_indices[0]))
            # Reformat indices into this format and leave only those indices matching current file name in list
            # [contigID, global_start, global_stop, cluster_id]
            child_indices = [(index[2][3], index[2][0], index[2][1], index[0]) for index in child_indices
                             if index[1] == fn]
            parent_indices = [(index[2][3], index[2][0], index[2][1], index[0]) for index in parent_indices
                              if index[1] == fn]

            print("C:{0}".format(child_indices[0]))
            # And sort them by contigID and global_start
            child_indices = sorted(child_indices, key=lambda x: (x[0], x[1]))
            parent_indices = sorted(parent_indices, key=lambda x: (x[0], x[1]))
            print("D:{0}".format(child_indices[0]))
            # Now walk through the list of parent and child indices until we reach the end of the parent or child indices
            while parent_indices and child_indices:
                child_index = child_indices[0]
                parent_index = parent_indices[0]
                # If the first child_index contig ID is lower than the first parent_index contig_id,
                # move along the child indices chain
                if child_index[0] < parent_index[0]:
                    child_indices.pop(0)
                    print("a")
                # Opposite case: move along parent indices chain
                elif parent_index[0] < child_index[0]:
                    parent_indices.pop(0)
                    print("b")
                # Last case: both child and parent index are from the same contig
                else:
                    # Does the child fragment end before the parent fragment starts?
                    # If so, walk along the child chain
                    if child_index[2] < parent_index[1]:
                        child_indices.pop(0)
                        print("c")
                    # Or does the parent fragment end before the child fragment starts?
                    # If so, walk along the parent chain
                    elif parent_index[2] < child_index[1]:
                        parent_indices.pop(0)
                        print("d")
                    # In any other case, there is some overlap between this child and parent fragment
                    # This means that we have to split the child fragment into two fragments,
                    # One fragment to the left of the overlap, one to the right. Obviously, depending on the
                    # type of the overlap, one or both of the fragments may be of size zero
                    else:
                        print('Overlap detected')
                        # First, retrieve the overlapping cluster from the child and parent node
                        child_overlap_cluster = child_node.get_clstr()[child_index[3]]
                        parent_overlap_cluster = parent_node.get_clstr()[parent_index[3]]
                        print(child_index)
                        print(child_overlap_cluster.get_indices())
                        print(parent_index)
                        print(parent_overlap_cluster.get_indices())
                        child_indices.pop(0)
                        parent_indices.pop(0)




if __name__ == '__main__':
    import sys
    input_xmfa_1 = sys.argv[1]
    input_xmfa_2 = sys.argv[2]
    build_xmfa = BuildXmfa("./build_test", "talia")
    print("Adding 1")
    build_xmfa.add_xmfa(input_xmfa_1)
    print("Adding 2")
    build_xmfa.add_xmfa(input_xmfa_2)



