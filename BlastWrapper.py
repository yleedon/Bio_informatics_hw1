import os
import Bio as bo
import numpy as np


class blastWrapper:

    def __init__(self):
        self.blast_db_path = "/Users/yanivleedon/Desktop/BioInformaticsAssignment1/blastDB"
        self.local_directory = os.getcwd()
        self.initial_set = "/usr/local/ncbi/blast/bin/"


    def make_blast_db(self, path_to_sequences = "/Users/yanivleedon/Desktop/BioInformaticsAssignment1/RawData/celegance_mrna_unspliced.fasta", title="celegansMRNA", output_directory: str = "/Users/yanivleedon/Desktop/BioInformaticsAssignment1/blastDB/"):
        cmd = "{0}makeblastdb -in {1} -title {2} -dbtype nucl -out {3}/{2}".format(self.initial_set, path_to_sequences, title, output_directory)
        os.system(cmd)

    def run_blast_search(self, blast_db_name, mi_rna_file_path = "/Users/yanivleedon/Desktop/BioInformaticsAssignment1/RawData/celegans-precursor.fasta", e_value_threshold: str = "1e-10", out_format:  str = "6", outputfile_path: str = "/Users/yanivleedon/Desktop/BioInformaticsAssignment1/res_blastn_compact"):
        cmd = "{0}blastn -db {1} -query {2} -out {3} -outfmt {4} -evalue {5}  ".format(self.initial_set, blast_db_name, mi_rna_file_path, outputfile_path, out_format, e_value_threshold)
        os.system(cmd)
