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
        left_sequence_cut_pos = left_fragment_length
        while len(child_sequence[:left_sequence_cut_pos].replace("-","")) < left_fragment_length:
            left_sequence_cut_pos += 1
        #right_sequence_cut_pos =






# Represent an entire alignment (i.e. a node in the hierarchical clustering graph
# as a collection of ClusterObjects
class AlignmentNode:
    def __init__(self, name):
        self.name = name
        self.clusters = {}

    def add_clstr(self, clstr):
        self.clusters[clstr.get_id()] = clstr

class BuildXmfa:
    def __init__(self, out_p, prefix):
        # Create output folder
        if not path.isdir(out_p):
            makedirs(out_p, exist_ok=True)
        self.out_f = open(path.join(out_p, "{0}_combined.xmfa".format(prefix)), "w")
        # Intialize progressive alignment tree
        # Name of root node is ""
        self.al_tree = {"": AlignmentNode("")}
        self.prefix = prefix

    def close_file(self):
        self.out_f.close()

    def add_xmfa(self, xmfa_p):
        # Create a new alignment node for this xmfa file
        node_label = path.basename(xmfa_p)[len(self.prefix)+1:].rsplit(".")[0]
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
        # Next, we need to substract everything from this alignment node that is already covered
        # by the parent node




if __name__ == '__main__':
    import sys
    input_xmfa = sys.argv[1]
    build_xmfa = BuildXmfa("./build_test", "talia")
    build_xmfa.add_xmfa(input_xmfa)





