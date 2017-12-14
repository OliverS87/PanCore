# Reduce input sample sequences to the regions that were not clustered
# in a previous parsnp run
# Provide functionality to map the results of the next round of parsnp
# on the reduced data to the original sequence indices
from os import path
from os import makedirs
from Bio import SeqIO

class ReduceInput:
    def __init__(self, original_input_sample_folder, pansnp_out_dir):
        self.index_mapper = {}
        self.original_input_sample_folder = original_input_sample_folder
        self.new_input_f = path.join(pansnp_out_dir, "new_input")
        makedirs(self.new_input_f, exist_ok=True)

    def sequences_added(self):
        return len(self.index_mapper) > 0

    # Two tasks:
    # 1. Split the useq into individual files, one for each si, these are used as input sample in the
    #   next pansnp round
    # 2. Store the header of the useqs in a dict, we need this header later to map the results of the
    #   next pansnp round back to the original reference
    def useq_to_input(self, file_list, useq_p):
        # Clear previous versions of new input files from the dict
        for file in file_list:
            filename = path.basename(file)
            try:
                del self.index_mapper[filename]
                print("Deleted {0}".format(filename))
            except KeyError:
                print("Could not delete {0}".format(filename))
        file_dict = {}
        active_out_f = None
        with open(useq_p, "r") as useq_f:
            for line in useq_f:
                if line.startswith(">"):
                    si = int(line.strip().split(":")[0][1:])

                    if si == 1:
                        active_out_f = None
                    else:
                        file_name = path.basename(file_list[si - 2])
                        try:
                            active_out_f = file_dict[si]
                        except KeyError:
                            active_out_f = open(path.join(self.new_input_f, file_name), "w")
                            file_dict[si] = active_out_f
                        index_data = line.split()
                        index_data[0] = [int(item) for item in index_data[0].split(":")[1].split("-")]
                        index_data[3] = [int(item[1:]) for item in index_data[3].split(":")]
                        # Format: global start, global stop, contigID, contig start
                        index_data = [*index_data[0], *index_data[3]]
                        try:
                            self.index_mapper[file_name].append(index_data)
                        except KeyError:
                            self.index_mapper[file_name] = [index_data]
                if line.startswith("="):
                    active_out_f = None
                if active_out_f:
                    active_out_f.write(line)
        # Close files
        [f.close() for f in file_dict.values()]

    def __map_index__(self, line, file_list):
        if line.startswith(">") and not line.startswith(">1:"):
            # Break down the header line:
            data = line.split()
            cluster_id = data[2]
            global_start, global_stop = [int(item) for item in data[0].split(":")[1].split("-")]
            fragment_length = abs(global_stop-global_start)+1
            si = int(data[0][1:data[0].index(":")])
            strand = data[1]
            contigID, contig_start = [int(item[1:]) for item in data[3].split(":")]
            # Pick file name
            file_name = path.basename(file_list[si-2])
            # Get index data for the input sequence
            # Format: global start, global stop, contigID, contig start
            org_index_data = self.index_mapper[file_name][contigID-1]
            mapped_global_start = org_index_data[0]+contig_start
            mapped_global_stop = org_index_data[0]+contig_start + fragment_length-1
            mapped_contigID = org_index_data[2]
            mapped_contig_start = org_index_data[3] + contig_start
            line = ">{0}:{1}-{2} {3} {4} s{5}:p{6}\n".format(si, mapped_global_start,
                                                             mapped_global_stop, strand,
                                                             cluster_id, mapped_contigID, mapped_contig_start)
        return line
    # Map the indices of the xmfa and useq file back to the original reference index
    def map_back(self, file_list, xmfa_file, useq_file):
        for file in (xmfa_file, useq_file):
            with open(file, "r") as f:
                f_content = f.readlines()
            # Map all indices back
            f_content_mapped = [self.__map_index__(line, file_list) for line in f_content]
            with open(file, "w") as f:
                f.writelines(f_content_mapped)


    # Modify the reference: Mask the regions involved in an alignment with 'E's
    # Replacement with 'N" or any other non-ACGT letter from the nucleotide alphabet
    # causes problems with parsnp: Indices are suffering from another offset, that cannot
    # be corrected by SimpleParSNPs usual correction methods
    def mask_reference(self, ref_file_p, xmfa_p, prefix):
        ref_seq = list(str(next(SeqIO.parse(ref_file_p, "fasta")).seq))
        with open(xmfa_p, "r") as xmfa_f:
            for line in xmfa_f:
                if line.startswith(">1:"):
                    data = line.split()
                    global_start, global_stop = [int(item) for item in data[0].split(":")[1].split("-")]
                    ref_seq[global_start:global_stop+1] = (global_stop+1-global_start)*'E'
        with open(path.join(self.new_input_f, "{0}_ref.faa".format(prefix)),"w") as new_ref:
            new_ref.write(">Ref:{0}\n".format(prefix))
            new_ref.write("".join(ref_seq))

if __name__ == '__main__':
    file_list = ["SRR1015292.fasta", "SRR1015290.fasta", "SRR1015291.fasta",
                 "SRR1015296.fasta", "SRR1015297.fasta", "SRR1015299.fasta", "SRR1015295.fasta"]
    file_list = [path.join("ab_assemblies_subset", item) for item in file_list]
    useq = "ab/parsnp.unalign"
    ri = ReduceInput("ab_assemblies_subset", "ab")
    #ri.useq_to_input(file_list, useq)
    #ri.map_back(file_list, "ab/cleo2.xmfa", "ab/cleo2.unalign")
    ri.mask_reference("NC_017162.gb.fna", "ab/parsnp.xmfa", "talia")