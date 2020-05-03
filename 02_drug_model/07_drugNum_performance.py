drug_num_dic = {}
with open("..//01_raw_data//08_drug_number_per_prot.txt", "r") as file1:
    for line in file1.readlines():
        data = line.replace("\n", "").split("\t")
        
        drug_num_dic[data[0]] = int(data[1])


prot_zscore_dic = {}
with open("zscore_all_model.txt", "r") as file2:
    for line in file2.readlines()[1:]:
        data = line.replace("\n", "").split("\t")
        prot_zscore_dic[data[0]] = float(data[1])

with open("drug_num_zscore.txt", "w") as o1:
    print("protein", "drug_number", "model_zscore", sep="\t", file=o1)
    for key in prot_zscore_dic:
        print(key, drug_num_dic[key], prot_zscore_dic[key], sep="\t", file=o1)