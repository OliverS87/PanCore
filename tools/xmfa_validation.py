# xmfa_validation:
# Check if the index positions found in xmfa file sequence headers match exactly with the original input file
# Based on these indices, a sequence substring is extracted from the input file and compared to the sequence
# in the xmfa file (raw string-vs-string comparison). Gaps are removed from the aligned sequence before comparison.
# The scripts output file validation.log counts the number of correctly matches indices (separate for global indices
# and contig start indices). It also states how many sequences span across contig boundaries.
# cmdline args: parSNP.xmfa  path_to_samples_folder reference.faa cpu_count out_path temporary_data_path
import sys
import os
import json
from Bio import SeqIO
if len(sys.argv[1:])!= 6:
    print("Usage: python3 {0} parSNP.xmfa path_to_samples_folder ref.faa cpu_count out_path "
          "temporary_data_path".format(sys.argv[0]))
    exit(3)
xmfa_file = sys.argv[1]
samples_path = sys.argv[2]
ref_file = sys.argv[3]
cpu_count = sys.argv[4]
# Create out folder
out_path = sys.argv[5]
if not os.path.isdir(out_path):
     os.mkdir(out_path)
# Create a folder where all merged assembly files can be stored
merged_path = sys.argv[6]
if not os.path.isdir(merged_path):
     os.mkdir(merged_path)


# Match sequence index (SI) with actual file name
seq_db = {}
with open(xmfa_file, "r") as xmfa:
    for line in xmfa:
        if line.startswith("##SequenceIndex"):
            si = int(line.split()[1])
            continue
        if line.startswith("##SequenceFile"):
            seq_db[si] = line.split()[1]
# Define the first and last SI, needed to iterate over all SI (i.e. all input assemblies)
min_si = min(seq_db.keys())
max_si = max(seq_db.keys())

# Convert each multifasta sample file into one merged fasta file
# Also store the contig boundaries in the merged file in a separate dict, i.e. contig start/stop in the merged assembly
# Dict keys are the filenames of the assembly files
# Check if a previous generated contig_bounds dict is present in the merged_path folder
if os.path.isfile(os.path.join(merged_path, "contig_bounds.dict")):
    with open(os.path.join(merged_path, "contig_bounds.dict"), "r") as dict_f:
        contig_bounds = json.load(dict_f)
        for filename, si_dict in contig_bounds.items():
            contig_bounds[filename] = {int(key):value for key,value in si_dict.items()}
    contig_bounds_dict_found = True
else:
    contig_bounds = {}
    contig_bounds_dict_found = False
for si in range(min_si, max_si+1):
    if si != 1:
        file_path = os.path.join(samples_path, seq_db[si])
    else:
        file_path = os.path.join(ref_file)
    # If a merged assembly file was created before and the contig_bounds dict was generated before, skip this step
    if os.path.isfile(os.path.join(merged_path, seq_db[si]+".merged.faa")) and contig_bounds_dict_found:
        print("Found merged file for {0}".format(seq_db[si]))
        continue
    with open(os.path.join(merged_path, seq_db[si]+".merged.faa"), "w") as single_faa_file:
            this_contig_bounds = {}
            next_contig_start = 0
            contig_index = 0
            single_faa_file.write(">{0} merged\n".format(seq_db[si]))
            for record in SeqIO.parse(file_path, "fasta"):
                contig_index +=1
                contig_length = len(record.seq)
                this_contig_bounds[contig_index]=(next_contig_start, next_contig_start+contig_length-1)
                next_contig_start = next_contig_start + contig_length
                single_faa_file.write("{0}\n".format(record.seq))
            contig_bounds[seq_db[si]]=this_contig_bounds

# Dump contig_bounds dict to file
with open(os.path.join(merged_path, "contig_bounds.dict"), "w") as dfile:
    json.dump(contig_bounds, dfile)



# Convert xmfa file to multiple multifasta files, one for each si
# Define output files, one for each si
query_files = {}
for index in range(min_si, max_si+1):
    query_files[index] = open(os.path.join(out_path, "seq{0}_query.faa".format(index)), "w")
# Define a reference for the active output file
active_output_file = None
# Iterate through xmfa file
with open(xmfa_file, "r") as xmfa:
    for line in xmfa:
        if line.startswith(">"):
            # Extract information from header
            header = line.strip()
            this_seq_si = int(header[1:header.index(":")])
            # Define active output file
            active_output_file = query_files[this_seq_si]
            # Rewrite header
            active_output_file.write(line)
            continue
        if line.startswith("#") or line.startswith("="):
            continue
        else:
            active_output_file.write(line)
# Close query files
[qfile.close() for qfile in query_files.values()]

# Define a reverse complement function
def rev_compl(sequence):
    try:
        return "".join([compl_nt(nt)for nt in sequence[::-1].upper()])
    except TypeError:
        print("Type error: {0}".format(sequence))
        return ""


def compl_nt(nt):
    if nt=="A":
        return "T"
    if nt=="T":
        return "A"
    if nt=="G":
        return "C"
    if nt=="C":
        return "G"
    else:
        return nt


# For each SI, we have an individual query file containing all the sequences for this particular SI found in the xmfa
# We next try to retrieve this sequence from the original assembly (based on contig ID and contig start value as stored
# in the header) and the merged assembly (based on start/stop indices as stored in the header)

# Start with merged assembly file
xmfa_position_matched_merged = 0
xmfa_position_not_matched_merged = 0
xmfa_within_ctg_bounds = []

with open(os.path.join(out_path, "xmfa_unmatched_merged.faa"), "w") as xmfa_unmatched_file:
    for index in range(min_si, max_si+1):
        # Load the assembly/reference sequence (reference for si == 1)
        orig_genome = next(SeqIO.parse(os.path.join(merged_path, seq_db[index] + ".merged.faa"), "fasta")).seq
        # For each sequence in the xmfa file, retrieve the sequence designated by its header coordinates from the
        # (merged) original assembly file and compare both sequences for identity
        for record in SeqIO.parse(os.path.join(out_path, "seq{0}_query.faa".format(index)), "fasta"):
            start_stop = [int(item) for item in record.description.split()[0].split(":")[1].split("-")]
            try:
                contig_id = int(record.description.split()[3].split(":")[0][1:])
            except ValueError:
                contig_id = -1
            # Test if sequence start/stop and contig ID start/stop are within this constraints:
            # contig_start <= seq_start <= seq_stop <= contig_stop
            try:
                xmfa_within_ctg_bounds.append(contig_bounds[seq_db[index]][contig_id][0] <= start_stop[0] <= start_stop[1] \
                                    <= contig_bounds[seq_db[index]][contig_id][1])
            except KeyError as ke:
                print(ke)
                xmfa_within_ctg_bounds.append(False)
            # Test if sequence was taken from the reverse strand
            if record.description.split()[1]=="-":
                sequence = rev_compl(record.seq)
            else:
                sequence = str(record.seq.upper())
            # Retrieve the sequence by index from the original input file
            try:
                orig_seq = orig_genome[start_stop[0]:start_stop[1]+1].upper()
            # If index wrong and sequence cannot be retrieved, value of sequence becomes "ERROR"
            except IndexError:
                orig_seq = "ERROR"
            # Sequences in XMFA may have gaps, condense sequences by removing gaps,
            # otherwise we can't compare them to the sequences from the ungapped input
            if sequence.replace("-","") == orig_seq:
                xmfa_position_matched_merged+=1
            else:
                xmfa_position_not_matched_merged+=1
                xmfa_unmatched_file.write(">{0}\n{1}\n".format(record.description, sequence))
                xmfa_unmatched_file.write(">ORG_seq_{0}_{1}_{2}\n{3}\n".format(index, start_stop[0],
                                                                               start_stop[1], orig_seq))




# Next check correctness of on-contig start and contig ID of each sequence stored in the xmfa, proceed from
# si = 1 till max(si)
xmfa_position_matched_contig = 0
xmfa_position_not_matched_contig = 0
with open(os.path.join(out_path, "xmfa_unmatched_contig.faa"), "w") as xmfa_unmatched_file:
    for index in range(min_si, max_si+1):
        # Load the assembly/reference sequence (reference for si == 1)
        if index != 1:
            file_path = os.path.join(samples_path, seq_db[index])
        else:
            file_path = os.path.join(ref_file)
        # Store contig sequences in a list, ordered by their order in the file
        orig_contigs = [rec.seq for rec in SeqIO.parse(file_path, "fasta")]
        # For each sequence in the xmfa file, retrieve the sequence designated by its header coordinates from the
        # contig original assembly file and compare both sequences for identity
        for record in SeqIO.parse(os.path.join(out_path, "seq{0}_query.faa".format(index)), "fasta"):
            start_stop = [int(item) for item in record.description.split()[0].split(":")[1].split("-")]
            seq_length = abs(start_stop[1]-start_stop[0]+1)
            try:
                contig_id = int(record.description.split()[3].split(":")[0][1:])
                contig_start = int(record.description.split()[3].split(":")[1][1:])
            except ValueError:
                contig_id = -1
                contig_start = 0
            # Test if sequence was taken from the reverse strand
            if record.description.split()[1] == "-":
                sequence = rev_compl(record.seq)
            else:
                sequence = str(record.seq.upper())
            # Retrieve the sequence by index from the original input file
            try:
                orig_seq = orig_contigs[contig_id-1][contig_start:contig_start+seq_length]
            # If index wrong and sequence cannot be retrieved, value of sequence becomes "ERROR"
            except IndexError:
                orig_seq = "ERROR"
            # Sequences in XMFA may have gaps, condense sequences by removing gaps,
            # otherwise we can't compare them to the sequences from the ungapped input
            if sequence.replace("-","") == orig_seq:
                xmfa_position_matched_contig += 1
            else:
                xmfa_position_not_matched_contig += 1
                xmfa_unmatched_file.write(">{0}\n{1}\n".format(record.description, sequence))
                xmfa_unmatched_file.write(">ORG_seq_{0}_{1}_{2}_{3}\n{4}\n".format(index, contig_id,
                                                                               contig_start, seq_length, orig_seq))






with open(os.path.join(out_path, "validation.log"), "w") as log_file:
    log_file.write("Matched {0} xmfa sequences ({1} unmatched (merged))\n".format(xmfa_position_matched_merged,
                                                                         xmfa_position_not_matched_merged))
    log_file.write("Xmfa seqs within contig bounds: {0} out of {1}\n".format(sum(xmfa_within_ctg_bounds),
                                                                             len(xmfa_within_ctg_bounds)))
    log_file.write("Matched {0} xmfa sequences ({1} unmatched (contig))\n".format(xmfa_position_matched_contig,
                                                                         xmfa_position_not_matched_contig))





# Clean up
for index in range(min_si, max_si+1):
    os.remove(os.path.join(out_path, "seq{0}_query.faa".format(index)))

