with open("02_protein_order.txt", "r") as file1:
    protein_order = [line.replace("\n", "") for line in file1.readlines()]

num_dic = {}
with open("08_drug_number_per_prot.txt", "r") as file2:
    for line in file2.readlines():
        data = line.replace("\n", "").split("\t")
        num_dic[data[0]] = int(data[1])


with open("new_prot2prot_proportionality.txt", "w") as o1:
    print("Prot_1", "Drug_1#", "Prot_2", "Drug_2#", "Proportionality", sep="\t", file=o1)
    with open("prot2prot_proportionality.txt", "r") as file3:
        for line in file3.readlines():
            data = line.replace("\n", "").split("\t")
            print(protein_order[int(data[0])].replace("_HUMAN", ""), num_dic[protein_order[int(data[0])]], protein_order[int(data[1])].replace("_HUMAN", ""), num_dic[protein_order[int(data[1])]], data[2], sep="\t", file=o1)
