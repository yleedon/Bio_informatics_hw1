from BlastWrapper import blastWrapper as bw


bw = bw()

# bw.make_blast_db(path_to_sequences="/Users/yanivleedon/Desktop/BioInformaticsAssignment1/RawData/celegance_mrna_unspliced.fasta", title="CElegans", output_directory="/Users/yanivleedon/Desktop/BioInformaticsAssignment1/blastDB/C_ElegansBlastDB")

mi_rna_file_path = "/Users/yanivleedon/Desktop/BioInformaticsAssignment1/RawData/celegans-precursor.fasta"
e_value_threshold = "1e-10"
out_format = "10"
blast_db_name = "/Users/yanivleedon/Desktop/BioInformaticsAssignment1/blastDB/C_ElegansBlastDB/CElegans"
output_file_path = "/Users/yanivleedon/Desktop/BioInformaticsAssignment1/csv_ready_results"
bw.run_blast_search(blast_db_name, mi_rna_file_path, e_value_threshold, out_format, output_file_path)
# # #
# # output_file_path = "/Users/yanivleedon/Desktop/BioInformaticsAssignment1/res_blastn_compact"
# # blast_db_name = "/Users/yanivleedon/Desktop/BioInformaticsAssignment1/blastDB/celegansMRNA"
# # bw.run_blast_search(blast_db_name, mi_rna_file_path, e_value_threshold, out_format)
#
#
# from Bio import SearchIO
# uncommented = output_file_path
# # qresult = SearchIO.read(uncommented, 'blast-tab')
# print(qresult)