import os



file_path = "/Users/yanivleedon/Desktop/BioInformaticsAssignment1/RawData/Celeganse_Mature_miRNA.fasta"
new_file_lines = []
with open(file_path, "r") as file:
    line = file.readline()
    while line != "\n" and not line == "" and line is not None:
        if line.startswith(">"):
            new_file_lines.append(line)
            line = file.readline()
            continue
        else:
            new_line = line.replace("U", "T")
            new_file_lines.append(new_line)
            line = file.readline()

with open(file_path, "w") as file:
    for line in new_file_lines:
        file.write("{}".format(line))

