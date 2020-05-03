import os

prot_dic = {}

file_list = os.listdir(".//all_model_scc_zscore//")


for files in file_list:
    prot = files.replace("_SVR_0.7_0.3_1000loop.txt", "")
    with open(".//all_model_scc_zscore//" + files, "r") as file1:
        data = file1.readlines()
        drug_num = int(data[0].replace("\n", "").split(" ")[-1])
        median_zscore = float(data[7].replace("\n", "").split(" ")[-1])
        prot_dic[prot] = [drug_num, median_zscore]
        
with open("protNum_n_medianZ.txt", "w") as o1:
    print("protein", "drug_number", "median_zscore", sep="\t", file=o1)
    for key in prot_dic:
        print(key, prot_dic[key][0], prot_dic[key][1], sep="\t", file=o1)
