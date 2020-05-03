from statistics import mean
import os
from scipy.stats import spearmanr


# load GO term
GO_dic = {}
file_list = os.listdir(".//prot_raw_BP_GO")
for files in file_list:
    with open(".//prot_raw_BP_GO//" + files, "r") as gofile:
        go_terms = set(line.replace("\n", "") for line in gofile)
        GO_dic[files.replace(".txt", "")] = go_terms

# jaccard
def jaccard(path_list):
    GO_path_list = [GO_dic[key] for key in path_list]
    
    intersections = len(set.intersection(*GO_path_list))
    union_list = []
    
    if len(GO_path_list) > 2:
        for i in range(0, len(GO_path_list)-1):
            for j in range(i+1, len(GO_path_list)):
                union_list.append(set.intersection(set(GO_path_list[i]), set(GO_path_list[j])))
        unions = len(set().union(*union_list))
    else:
        unions = len(set().union(set(GO_path_list[0]), set(GO_path_list[1])))
    return 1. * intersections / unions


prot_pair_list = []
coef_list = []
GOsim_list = []

# load feature coef
with open(".//model_to_model_coeff_scc.txt", "r") as file1:
    f1 = file1.readlines()[1:]
    for line in f1:
        data = line.replace("\n", "").split("\t")
        coef = data[-1]

        prot_pair_list.append(data[:2])
        coef_list.append(coef)
        GOsim_list.append(jaccard(data[:2]))
        
# write file
with open("feature_coef_GOsim.txt", "w") as o1:
    print("SCC of feature coef & GOsim", spearmanr(coef_list, GOsim_list), sep="\n", file=o1)
    print("model_1", "model_2", "SCC_of_coef", "GO_sim", sep="\t", file=o1)
    for i in range(len(prot_pair_list)):
        print(prot_pair_list[i][0], prot_pair_list[i][1], coef_list[i], GOsim_list[i], sep="\t", file=o1)