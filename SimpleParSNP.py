# SimpleParSNP class
# Instances configure the .ini file required by parsnp
# and run parsnp
# Also provides functionality to correct parsnps xmfa index bug
# and to create a file containing all unaligned sequences
from os import path
import os
import sys
import subprocess
from simpleparsnp_libs import correctXMFA
from simpleparsnp_libs import generateUseq
from Bio import SeqIO
from datetime import datetime as dt


class SimpleParSNP:
    def __init__(self):
        # Set defaults for parsnp parameters:
        # Max. allowed distance between two seeds during LCB formation
        self.dist = 30
        # Min. cluster size
        self.min_cluster_size = 21
        # Nr. of threads
        self.threads = 1
        # Path to reference
        self.ref_p = ""
        # Path to input files
        self.files_p = []
        # Match SI to file path
        self.seq_db = {}
        self.min_si = 0
        self.max_si = 0
        # Prefix for all output file names
        self.prefix = ""
        # Store some information about each contig in the set of input files
        self.input_contigs_dict = {}
        self.input_contigs_org_names = {}
        # Global log file
        self.log_p = None

    def set_reference(self, new_ref_p):
        if path.isfile(new_ref_p):
            self.ref_p = new_ref_p

    def add_files(self, new_files_p):
        for new_file_p in new_files_p:
            if path.isfile(new_file_p):
                self.files_p.append(new_file_p)
        # After adding new files, rebuild contig DB
        self.update_seqdb()
        self.update_contigdb()

    def set_threads(self, thread_nr):
        try:
            int(thread_nr)
        except ValueError:
            return
        self.threads = thread_nr

    def set_dist(self, new_dist):
        try:
            int(new_dist)
        except ValueError:
            return
        self.dist = new_dist

    def set_prefix(self, new_pf):
        # Dots are not allowed in prefix, replace with underscore
        self.prefix = new_pf.replace(".","_")

    def generate_ini(self, out_dir):
        ini_p = path.join(out_dir, "{0}parsnp_config.ini".format(self.prefix))
        with open(ini_p, "w") as ini_f:
            # Write file paths
            ini_f.write(";Parsnp configuration File\n")
            ini_f.write(";\n")
            ini_f.write("[Reference]\n")
            ini_f.write("file={0}\n".format(self.ref_p))
            ini_f.write("reverse=0\n")
            ini_f.write("[Query]\n")
            for file_p in enumerate(self.files_p, 1):
                ini_f.write("file{0}={1}\n".format(file_p[0], file_p[1]))
                ini_f.write("reverse{0}=0\n".format(file_p[0]))
            # Write standard MUM settings
            ini_f.write("[MUM]\nanchors=1.0*(Log(S))\nanchorfile=\nanchorsonly=0\ncalcmumi=0\n"
                        "mums=1.1*(Log(S))\nmumfile=\nfilter=1\nfactor=2.0\nextendmums=0\n")
            # Write LCB settings
            ini_f.write("[LCB]\nrecombfilter=0\n")
            ini_f.write("cores={0}\n".format(self.threads))
            ini_f.write("diagdiff=0.12\ndoalign=2\n")
            ini_f.write("c={0}\n".format(self.min_cluster_size))
            ini_f.write("d={0}\n".format(self.dist))
            ini_f.write("q=30\np=15000000\nicr=0\nunaligned=0\n")
            # Write output settings
            ini_f.write("[Output]\n")
            ini_f.write("outdir={0}\n".format(path.join(out_dir, self.prefix)))
            ini_f.write("prefix={0}\nshowbps=0\n".format(self.prefix))

    # Update seq DB dict that matches each file name to an SI as returned later by parsnp
    def update_seqdb(self):
        # Reference
        self.seq_db[1] = self.ref_p
        # Sample files
        for file_p in enumerate(self.files_p, 2):
            self.seq_db[file_p[0]] = file_p[1]
        # Define smallest and largest Sequence Index, needed later to iterate through each assembly
        self.min_si = min(self.seq_db.keys())
        self.max_si = max(self.seq_db.keys())

    # Update the information about contigs in the input set
    # We don't have the information about the lowest and highest cluster nr, yet
    # Use None as placeholders for this information
    def update_contigdb(self):
        for si in range(self.min_si, self.max_si + 1):
            self.input_contigs_dict[si] = []
            input_file_path = self.seq_db[si]
            next_start_global = 0
            contig_count = 0
            for rec in SeqIO.parse(input_file_path, "fasta"):
                contig_count += 1
                self.input_contigs_org_names[(si, contig_count)] = rec.id
                seq_len = len(rec.seq)
                self.input_contigs_dict[si].append([si, contig_count,
                                          next_start_global, next_start_global + seq_len - 1, 0,
                                          seq_len - 1, seq_len, "+", None, None])
                next_start_global += seq_len

    # Create a log file
    def create_log(self, out_dir):
        self.log_p = path.join(out_dir, "{0}.log".format(self.prefix))
        # Create  log file
        with open(self.log_p, "w") as log_f:
            log_f.write("---simpleparsnp_libs V0.0.1---\n")

    # Write to log file
    def write_log(self, message):
        with open(self.log_p, "a") as log_f:
            log_f.write(message)

    def run_parsnp(self, out_dir, generate_useq, generate_icstats):
        start = dt.now()
        if not path.isdir(path.join(out_dir, self.prefix)):
            os.makedirs(path.join(out_dir, self.prefix), exist_ok=True)
        self.create_log(out_dir)
        self.generate_ini(out_dir)
        self.write_log("Running parsnp on {0} samples using {1} as reference.\n".format(len(self.files_p),
                                                                                        path.basename(self.ref_p)))
        # Write parsnp parameters to log file
        self.write_log("Parsnp parameters:\nMaxDist={0},MinClstrSize={1},MaxDiagDiff=0.12\n".
                       format(self.dist, self.min_cluster_size))
        # Write file names to log file
        for i, p in enumerate(self.files_p, 2):
            self.write_log("si{0}:{1}\n".format(i, path.basename(p)))
        run = subprocess.run("./parsnp.simple " + path.join(out_dir, "{0}parsnp_config.ini".format(self.prefix)),
                             shell=True, stdin=None,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             close_fds=True, executable="/bin/bash")
        # Return RC if failed
        if run.returncode != 0:
            self.write_log("Parsnp run failed. Bummer :(\n")
            return run.returncode
        # parsnp creates its own log file, join its content to "our" log file
        with open(os.path.join(out_dir, self.prefix, "parsnpAligner.log"), "r") as parsnp_log_f:
            for line in parsnp_log_f:
                if line.startswith("Mum anchor size"):
                    break
            for line in parsnp_log_f:
                self.write_log(line)

        # Rename raw xmfa output file
        #os.rename(path.join(out_dir, "parsnpAligner.xmfa"), path.join(out_dir, "{0}.xmfa".format(self.prefix)))
        # Correct header lines in raw xmfa output
        xmfa_corrector = correctXMFA.CorrectXMFA(self.dist, path.join(out_dir, self.prefix, "parsnpAligner.xmfa"),
                                        self.min_si,  self.max_si, self.input_contigs_dict, self.write_log)
        corrected_headers = xmfa_corrector.correct_xmfa()
        # Check for error code
        try:
            int(corrected_headers)
            self.write_log("Error while correcting xmfa. Bummer :(\n")
            return -2
        except TypeError:
            pass
        # Rename corrected xmfa output file and move from tmp folder to output directory
        os.rename(path.join(out_dir, self.prefix, "parsnpAligner.xmfa.corr"),
                  path.join(out_dir, "{0}.xmfa".format(self.prefix)))
        # Generate unaligned sequence files
        useq_generator = generateUseq.GenerateUseq(self.seq_db, corrected_headers, self.input_contigs_dict,
                                                     self.input_contigs_org_names, out_dir, self.prefix, self.write_log)
        try:
            useq_generator.generate_useqs(generate_useq, generate_icstats)
        except Exception:
            self.write_log("Error while creating unaligned file xmfa. Bummer :(\n")
            return -3
        # Clean up
        os.remove(path.join(out_dir, self.prefix, "parsnpAligner.log"))
        os.remove(path.join(out_dir, self.prefix, "parsnpAligner.xmfa"))
        os.rmdir(path.join(out_dir, self.prefix))
        os.remove(path.join(out_dir, "{0}parsnp_config.ini".format(self.prefix)))
        time_elapsed = (dt.now() - start).total_seconds()
        self.write_log("Finished aligning {0}+1 sequences on {1} CPUs in {2} seconds.\n".
                       format(len(self.files_p), self.threads, time_elapsed))
        return 0





if __name__ == '__main__':
    # Test run
    # Input format: dist, ref, core, out_dir, prefix, sample_file(s)
    if len(sys.argv[1:]) != 6:
        print("Input format: dist, ref, core, out_dir, prefix, sample_file(s)")
        exit(3)
    simple_snp = SimpleParSNP()
    simple_snp.set_dist(sys.argv[1])
    simple_snp.set_reference(sys.argv[2])
    simple_snp.set_threads(sys.argv[3])
    simple_snp.add_files(sys.argv[6:])
    simple_snp.set_prefix(sys.argv[5])
    simple_snp.run_parsnp(sys.argv[4])
