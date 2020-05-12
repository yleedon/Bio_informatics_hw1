import pandas as pd
from Bio import SearchIO
from BlastWrapper import blastWrapper as bw


bw = bw()

# bw.make_blast_db(path_to_sequences="/Users/yanivleedon/Desktop/BioInformaticsAssignment1/RawData/briggsae-pre.fasta", title="cbriggsae", output_directory="/Users/yanivleedon/Desktop/BioInformaticsAssignment1/blastDB/C_BriggsaeBlastDB")

# mi_rna_file_path = "/Users/yanivleedon/Desktop/BioInformaticsAssignment1/RawData/Mature_mirnas_paired_with_precursurs.fasta"
# e_value_threshold = "1e-3"
# out_format = "6"
# blast_db_name = "/Users/yanivleedon/Desktop/BioInformaticsAssignment1/blastDB/C_ElegansBlastDB/CElegans"
output_file_path = "/Users/yanivleedon/Desktop/BioInformaticsAssignment1/csv_ready_results"
# bw.run_blast_search(blast_db_name, mi_rna_file_path, e_value_threshold, out_format, output_file_path, additional_arguments=["-task blastn"])
# #
# output_file_path = "/Users/yanivleedon/Desktop/BioInformaticsAssignment1/res_blastn_compact"
# blast_db_name = "/Users/yanivleedon/Desktop/BioInformaticsAssignment1/blastDB/celegansMRNA"
# bw.run_blast_search(blast_db_name, mi_rna_file_path, e_value_threshold, out_format)


# uncommented = output_file_path
# qresult = SearchIO.read(uncommented, 'blast-tab')
# print(qresult)

# def pair_mirna_with_pre_mirna(path1: str, path2: str):
#
#     all_mirnas = {}
#     mirnas_paired_with_precursur = {}
#     with open(path1, "r") as file:
#         line = file.readline()
#         while line != "\n" and not line == "" and line is not None:
#             if line.startswith(">"):
#                 mirnaID = line.replace("\n", "")
#                 mirnaSeq = file.readline().replace("\n", "")
#                 all_mirnas[mirnaID] = mirnaSeq
#                 line = file.readline()
#                 continue
#
#     with open("Mature_mirnas_paired_with_precursurs.fasta", "w") as new_file:
#         with open(path2, "r") as file:
#             line = file.readline()
#             mirnas_ids = list(all_mirnas.keys())
#             while line != "\n" and not line == "" and line is not None:
#                 if line.startswith(">"):
#                     precursur_id = line.replace("_pre", "").replace("\n", "")
#                     for key in mirnas_ids:
#                         if key.__contains__(precursur_id):
#                             new_file.write("{}\n{}\n".format(key, all_mirnas[key]))
#                             break
#
#                 line = file.readline()
#
#
# pair_mirna_with_pre_mirna("/Users/yanivleedon/Desktop/BioInformaticsAssignment1/RawData/Celeganse_Mature_miRNA.fasta","/Users/yanivleedon/Desktop/BioInformaticsAssignment1/RawData/celegans-precursor.fasta")

def get_seq_dict(path: str):
    ans = dict()
    with open(path, "r") as file:
        line = file.readline().replace("\n", "")
        while line is not None:
            if line == "":
                line = file.readline().replace("\n", "")
                if line == "":
                    break
                continue
            if line.startswith(">"):
                new_line = ""
                mir_id = line.replace(">", "")
                seq_line = file.readline().replace("\n", "")
                while seq_line != "\n" and not seq_line == "" and seq_line is not None and not seq_line.startswith(">"):
                    new_line += seq_line
                    seq_line = file.readline().replace("\n", "")

                line = seq_line
                ans[mir_id] = new_line

    return ans


def compare_seq(utr_seq, mirna_seq):
    if utr_seq == 'Sequence unavailable' or mirna_seq == 'Sequence unavailable':
        return '-'

    mirna_seq = mirna_seq[::-1].replace('U', 'T') # reverse mirnae
    mirna_seq = mirna_seq[1:8]
    utr_seq = utr_seq[1:8]

    return 'yes' if mirna_seq == utr_seq else 'no'






def check_if_targets_utr(gene_id: str, mirna_id, utr3_ids_to_seq, mature_mir_ids_to_seq_dic):

    if gene_id not in utr3_ids_to_seq or mirna_id not in mature_mir_ids_to_seq_dic:
        return '-'
    return compare_seq(utr3_ids_to_seq[gene_id], mature_mir_ids_to_seq_dic[mirna_id] )




def parser(file_path, out_path):
    utr3_ids_to_seq = get_seq_dict('/Users/yanivleedon/Desktop/BioInformaticsAssignment1/RawData/celegance_3utr')
    mature_mir_ids_to_seq_dic = get_seq_dict("/Users/yanivleedon/Desktop/BioInformaticsAssignment1/RawData/Celeganse_Mature_miRNA.fasta")
    unspliced_gene_id_to_seq_dic = get_seq_dict("/Users/yanivleedon/Desktop/BioInformaticsAssignment1/RawData/celegance_mrna_unspliced.fasta")

    results_dic = dict({"C.elegans Mature Name":[], "C.elegans Mature Sequence":[], "C.elegans Pre-miRNA Name":[], "C.elegans Pre-mRNA Sequence":[], "Host Gene Name":[], "Targets the Host Gene":[], "Conserved in C.Briggsae":[]})
    with open(file_path, 'r') as file:
        line = file.readline()
        while line != "\n" and not line == "" and line is not None:
            split_line = line.split('\t')
            if split_line[2] != '100.000':
                line = file.readline()
                continue
            mirna_id = split_line[0]
            gene_id = split_line[1]
            results_dic["C.elegans Mature Name"].append(mirna_id)
            results_dic["C.elegans Mature Sequence"].append(mature_mir_ids_to_seq_dic[mirna_id])
            results_dic["C.elegans Pre-miRNA Name"].append("{}_pre".format(mirna_id.split("_")[0]))
            results_dic["C.elegans Pre-mRNA Sequence"].append(unspliced_gene_id_to_seq_dic[gene_id])
            results_dic["Host Gene Name"].append(gene_id)
            results_dic["Targets the Host Gene"].append(check_if_targets_utr(gene_id, mirna_id, utr3_ids_to_seq, mature_mir_ids_to_seq_dic))
            results_dic["Conserved in C.Briggsae"].append(None)
            line = file.readline()

    # todo: add all mirna Ids that have no host gene with appropriate data
    return pd.DataFrame(results_dic)



parser('/Users/yanivleedon/Desktop/BioInformaticsAssignment1/csv_ready_results', '')
