# Edit header lines in fasta files to remove special characters ':' and '-' not tolerated by mugsy
import sys
fasta_file_list = sys.argv[1:]
for fasta_p in fasta_file_list:
    file_content = []
    with open(fasta_p, "r") as fasta_f:
        for line in fasta_f:
            if line.startswith(">"):
                line = line.replace(":","")
                line = line.replace("-", "")
            file_content.append(line)
    with open(fasta_p, "w") as fasta_f:
        for line in file_content:
            fasta_f.write(line)

