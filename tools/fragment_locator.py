# Fragment locator:
# Find the position of fragments found in parSNPs .xmfa file in the original sample input files.
# Sequence fragments in .xmfa are matched back to the original assembly file using blastn
# A blast DB is build from the assembly files, two DB for each file: One in which assembly contigs are merged into one
# sequence in the order of their appearance in the file and one DB with contigs stored as separate sequences
# The indices of the blast hits are checked against the indices reported in parSNPs output file.
# parsSNPs output has two types of indices: Global and Contig. Global indices see the input assembly file as one large
# sequence, with all contigs merged together in the order of their appearance. Contig indices refer to an actual
# position on that specific contig.
# These indices are checked:
# parsnp.xmfa: Global start position of a fragment and a fragments contig start position.
# The script creates two output files: contig_offsets.txt  and merged_offsets.txt for index global start index
# offsets and contig start index offsets. Columns of the tsv-formatted output files are:
# SequenceIndex ContigID Strand-orientation index-offset and sequence length


import sys
import os
import subprocess

# Check for correct number of cmdline arguments
# Expected format:
# parSNP.xmfa path_to_samples_folder reference.faa cpu_count out_path
if len(sys.argv[1:]) != 5:
    print("Usage: python3 {0} parSNP.xmfa path_to_samples_folder reference.faa "
          "cpu_count out_path".format(sys.argv[0]))
    exit(3)

# Parse paths to input arguments from sys.argv
xmfa_file = sys.argv[1]
samples_path = sys.argv[2]
ref_file = sys.argv[3]
cpu_count = sys.argv[4]

# Create the output folder if it does not exist
out_path = sys.argv[5]
if not os.path.isdir(out_path):
    os.makedirs(out_path, exist_ok=True)

# Check if blast is in $PATH
run = subprocess.run(["blastn -version"],shell=True, stdin=None,stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE,close_fds=True)
if run.returncode != 0:
    print("Can't find blast executables in $PATH :-(")
    exit(4)
else:
    print("Found {0}".format(run.stdout.decode()))

# Sequence headers in parSNPs output use a sequence index (SI) to identify assemblies
# Match sequence index (SI) with actual file name in seq_db dict
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

# Convert xmfa file to multiple multi-fasta files, one for each SI
# These files will be used as queries against the blastdb later.
# Define output files, one for each SI
query_files = {index:open(os.path.join(out_path, "seq{0}_query.faa".format(index)), "w")
               for index in range(min_si, max_si+1)}

# active_output_file points to the file handler currently in use
active_output_file = None

# Iterate once through the xmfa file
# Everytime we encounter a header line, the active output file is switched according to the SI in the header
# The header has to be rewritten so that it can be parsed by blastn
with open(xmfa_file, "r") as xmfa:
    for line in xmfa:
        if line.startswith(">"):
            # Extract information from header
            header = line.strip()
            this_seq_si = int(header[1:header.index(":")])
            contig_start = header.split(":p")[1]
            contig_index = header.split()[3].split(":")[0]
            strand_orientation = header.split()[1]
            overall_start = header.split()[0].split(":")[1].split("-")[0]
            overall_stop = header.split()[0].split(":")[1].split("-")[1]
            # Define active output file
            active_output_file = query_files[this_seq_si]
            # Rewrite header
            active_output_file.write(">lcl|{0}_{1}_{2}_{3}_{4}\n".format(overall_start,
                                                                 overall_stop, strand_orientation, contig_start,
                                                                 contig_index))
        # Write sequence data to the active output file
        elif not line.startswith("#") and not line.startswith("="):
            active_output_file.write(line)
# Close query files
[qfile.close() for qfile in query_files.values()]

# Create merged assembly files: All contigs are joined into one sequence in the order of their appearance in the
# assembly file
for si in seq_db.keys():
    # All non-reference assembly files
    if si != 1:
        with open(os.path.join(samples_path, seq_db[si]), "r") as multi_faa_file:
            with open(os.path.join(out_path, seq_db[si]+".merged.faa"), "w") as single_faa_file:
                single_faa_file.write(">{0} merged\n".format(seq_db[si]))
                for line in multi_faa_file:
                    if not line.startswith(">") and not len(line.strip()) == 0:
                        single_faa_file.write(line)
    # Reference file is in a different location
    if si == 1:
        with open(os.path.join(ref_file), "r") as multi_faa_file:
            with open(os.path.join(out_path, "ref_merged.faa"), "w") as single_faa_file:
                single_faa_file.write(">{0} merged\n".format(seq_db[si]))
                for line in multi_faa_file:
                    if not line.startswith(">") and not len(line.strip()) == 0:
                        single_faa_file.write(line)


# For each SI, run the fragments from the .xmfa file against the blast DBs created from the input files
# si == 1 points to the reference file which is stored in a different location then the other assembly files
for si in range(min_si, max_si+1):
    # For all non-reference assemblies
    if si != 1:
        # Create blast_db from one assembly file (keeping contigs separate)
        subprocess.run(["makeblastdb", "-in", os.path.join(samples_path, seq_db[si]), "-input_type", "fasta", "-dbtype",
                       "nucl", "-parse_seqids", "-out", os.path.join(out_path, "contig_blastdb_{0}".format(si))])
        # Run this SI query file against the contig-level blast db
        subprocess.run(["blastn", "-db", os.path.join(out_path, "contig_blastdb_{0}".format(si)), "-query",
                       os.path.join(out_path, "seq{0}_query.faa".format(si)), "-parse_deflines", "-outfmt",
                       "6 qseqid sstart send mismatch", "-num_threads", cpu_count, "-out",
                        os.path.join(out_path, "contig_blastres_{0}".format(si))])
        # Create blast_db from one assembly file (contigs merged into one sequence)
        subprocess.run(["makeblastdb", "-in", os.path.join(out_path, seq_db[si]+".merged.faa"), "-input_type", "fasta",
                       "-dbtype", "nucl", "-parse_seqids", "-out",
                       os.path.join(out_path, "merged_blastdb_{0}".format(si))])
        # Run this SI query file against the merged-level blast db
        subprocess.run(["blastn", "-db", os.path.join(out_path, "merged_blastdb_{0}".format(si)), "-query",
                       os.path.join(out_path, "seq{0}_query.faa".format(si)), "-parse_deflines", "-outfmt",
                       "6 qseqid sstart send mismatch", "-num_threads", cpu_count,
                        "-out", os.path.join(out_path, "merged_blastres_{0}".format(si))])
    # Only for the reference assembly
    # The reference has to be a continuous sequence so the nr. of contigs is one
    # We still perform a blast search for a "merged" reference sequence, which is identical to the
    # contig-level reference sequence. Blast results for both runs should be identical, so this is a good
    # check to ensure that there are no systematic errors in this approach
    if si == 1:
        # Create blast_db from one assembly file (keeping contigs separate)
        subprocess.run(["makeblastdb", "-in", ref_file, "-input_type", "fasta", "-dbtype",
                       "nucl", "-parse_seqids", "-out", os.path.join(out_path, "contig_blastdb_{0}".format(si))])
        # Run this SI query file against the contig-level blast db
        subprocess.run(["blastn", "-db", os.path.join(out_path, "contig_blastdb_{0}".format(si)), "-query",
                       os.path.join(out_path, "seq{0}_query.faa".format(si)), "-parse_deflines", "-outfmt",
                       "6 qseqid sstart send mismatch", "-num_threads", cpu_count, "-out",
                        os.path.join(out_path, "contig_blastres_{0}".format(si))])
        # Create blast_db from one assembly file (contigs merged into one sequence)
        subprocess.run(["makeblastdb", "-in", os.path.join(out_path, "ref_merged.faa"), "-input_type", "fasta",
                       "-dbtype", "nucl", "-parse_seqids", "-out",
                       os.path.join(out_path, "merged_blastdb_{0}".format(si))])
        # Run this SI query file against the merged-level blast db
        subprocess.run(["blastn", "-db", os.path.join(out_path, "merged_blastdb_{0}".format(si)), "-query",
                       os.path.join(out_path, "seq{0}_query.faa".format(si)), "-parse_deflines", "-outfmt",
                       "6 qseqid sstart send mismatch", "-num_threads", cpu_count,
                        "-out", os.path.join(out_path, "merged_blastres_{0}".format(si))])

# Compare the on-contig start position as reported in the .xmfa file with the actual on-contig start position
# returned by blastn. Blast runs were done indepentently for each SI. We calculate eventual differences between the
# reported and actual start position for each fragment and store them in dicts, sorted by SI and contig ID.
contig_offset = {}
for si in range(min_si, max_si+1):
    # Store all offsets for this SI
    si_contig_offset = {}
    with open(os.path.join(out_path, "contig_blastres_{0}".format(si)), "r") as contig_file:
        # Not all blast hits should be considered, some need to be filtered out. These are:
        # Erroneous hits: Hits were subject and query length do not match 100% and/or where mismatches occur
        filtered_file = {}
        for line in contig_file:
            # Split header of query file into separate columns
            # original format:
            # overall_start, overall_stop, strand_orientation, contig_start, contig_index
            data = line.replace("_", " ").split()
            # Each line is now: qseqid subjectstart subjectend mismatch
            # with qseqid split into overall_start, overall_stop, strand_orientation, contig_start, contig_index
            if len(data) != 8:
                continue
            # Check for misplaced blast hits:
            # 1. Length of aligned seq and blast hit doesnt't match 100%:
            if abs(int(data[0])-int(data[1])) != abs(int(data[5])-int(data[6])):
                continue
            # 2. Check for mismatches
            if int(data[7]) != 0:
                continue
            # Calculate difference between sequence start position on contig as reported by parsnp and actual start
            # position as reported by blast
            # For sequence fragments aligned in reverse complement, use stop position reported by blast
            if data[2] == "+":
                offset = abs(int(data[3])-int(data[5])+1)
            else:
                offset = abs(int(data[3])-int(data[6])+1)
            # For queries that have multiple hits, keep only that one with the lowest offset (abs)
            # Store complete lines, line-id is query_id (unsplit)
            line_id = line.split()[0]
            if line_id not in filtered_file:
                filtered_file[line_id] = (line, offset)
            # If offset stored in filtered_file for this query ID is larger than the offset calculated in this loop,
            # replace the stored line with the current line
            if filtered_file[line_id][1] > offset:
                filtered_file[line_id] = (line, offset)
        # Retrieve all lines that were not filtered out
        filtered_lines = [item[0] for item in filtered_file.values()]
        # For each contig_id, store all offsets in a list in this format:
        # [(offset, reported_start_on_contig, length_of_aligned_fragment), ...]
        # For sequence fragments aligned in rev. compl use the subject_stop_blast value
        for line in filtered_lines:
            data = line.replace("_", " ").split()
            contig_id = data[4]
            if not contig_id in si_contig_offset:
                # For each contig_id store offsets on the + and - oriented alignment in two separate lists
                si_contig_offset[contig_id] = ([], [])
            if data[2]=="+":
                si_contig_offset[contig_id][0].append((int(data[3])-int(data[5])+1,
                                                      int(data[3]), abs(int(data[0])-int(data[1]))+1))
            else:
                si_contig_offset[contig_id][1].append((int(data[3])-int(data[6])+1,
                                                      int(data[3]), abs(int(data[0])-int(data[1]))+1))
    # Add the offsets for this SI to the overall collection of offsets
    contig_offset[si] = si_contig_offset

# Write offsets to a tsv file
# Format:
# SI ContigID Strand Offset SeqSize
with open(os.path.join(out_path, "contig_offsets.txt"), "w") as out_file:
    for si in range(min_si, max_si+1):
        si_contig_offset = contig_offset[si]
        # Sort by contig id, with s1 < s2 < s3...
        offsets_sorted = sorted(si_contig_offset.items(), key=lambda x: x[0])
        # For each contig id, write the collected + and - offsets to file
        # Format: seq_id, contig_id, strand, offset, length of fragment
        for offset in offsets_sorted:
            contig_id = offset[0]
            for plusstrand_offset in offset[1][0]:
                out_file.write("{0}\t{1}\t+\t{2}\t{3}\n".
                               format(si, contig_id, plusstrand_offset[0],plusstrand_offset[2]))
            for minusstrand_offset in offset[1][1]:
                out_file.write("{0}\t{1}\t-\t{2}\t{3}\n".
                               format(si, contig_id, minusstrand_offset[0],minusstrand_offset[2]))

# Compare the global start/stop indices reported in the .xmfa output to the actual start/stop indices as reported by
# blastn.
merged_offset = {}
for si in range(min_si, max_si+1):
    si_merged_offset = {}
    with open(os.path.join(out_path, "merged_blastres_{0}".format(si)), "r") as merged_file:
        filtered_file = {}
        for line in merged_file:
            # Split header of query file into separate columns
            # original format:
            # overall_start, overall_stop, strand_orientation, contig_start, contig_index
            # Each line is now: qseqid subjectstart subjectend mismatch
            # with qseqid split into overall_start, overall_stop, strand_orientation, contig_start, contig_index
            data = line.replace("_", " ").split()
            if len(data) != 8:
                continue
            # Check for misplaced blast hits:
            # 1. Length of aligned seq and blast hit doesnt't match 100%:
            if abs(int(data[0])-int(data[1])) != abs(int(data[5])-int(data[6])):
                continue
            # 2. Check for mismatches
            if int(data[7]) != 0:
                continue
            # Calculate difference between sequence start position on contig as reported by parsnp and actual start
            # position as reported by blast
            # For sequence fragments aligned in reverse complement, use stop position reported by blast
            if data[2] == "+":
                offset = abs(int(data[0])-int(data[5])+1)
            else:
                offset = abs(int(data[0])-int(data[6])+1)
            line_id = line.split()[0]
            # For queries that have multiple hits, keep only that one with the lowest offset (abs)
            if line_id not in filtered_file:
                filtered_file[line_id] = (line, offset)
            # If offset stored in filtered_file for this query ID is larger than the offset calculated in this loop,
            # replace the stored line with the current line
            if filtered_file[line_id][1] > offset:
                filtered_file[line_id] = (line, offset)
        # Retrieve lines from filtered_file
        filtered_lines = [item[0] for item in filtered_file.values()]
        for line in filtered_lines:
            data = line.replace("_", " ").split()
            # Sort offsets by contig ID and strand
            contig_id = data[4]
            if not contig_id in si_merged_offset:
                si_merged_offset[contig_id] = ([], [])
            if data[2]=="+":
                si_merged_offset[contig_id][0].append((int(data[0])-int(data[5])+1,
                                                      int(data[0]), abs(int(data[0])-int(data[1]))+1))
            else:
                si_merged_offset[contig_id][1].append((int(data[0])-int(data[6])+1,
                                                      int(data[0]), abs(int(data[0])-int(data[1]))+1))
    # Add the offsets for this SI to the overall collection of offsets
    merged_offset[si] = si_merged_offset

# Write offsets to TSV file
# Format: SeqID ContigID Strand Offset SeqSize
with open(os.path.join(out_path, "merged_offsets.txt"), "w") as out_file:
    for si in range(min_si, max_si+1):
        si_merged_offset = merged_offset[si]
        # Sort by contig_id
        offsets_sorted = sorted(si_merged_offset.items(), key=lambda x: x[0])
        for offset in offsets_sorted:
            contig_id = offset[0]
            for plusstrand_offset in offset[1][0]:
                out_file.write("{0}\t{1}\t+\t{2}\t{3}\n".
                               format(si, contig_id, plusstrand_offset[0],plusstrand_offset[2]))
            for minusstrand_offset in offset[1][1]:
                out_file.write("{0}\t{1}\t-\t{2}\t{3}\n".
                               format(si, contig_id, minusstrand_offset[0],minusstrand_offset[2]))
