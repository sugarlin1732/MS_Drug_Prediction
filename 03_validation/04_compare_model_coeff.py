from sklearn.svm import SVR
from joblib import dump, load
from scipy.stats import spearmanr, zscore
import numpy as np
import os

# load model
model_list = os.listdir("..//02_drug_model//prot_model//")

# load prot name
with open("..//02_drug_model//02_protein_order.txt", "r") as file1:
    prot_name_lsit = [line.replace("\n", "") for line in file1]

# make coeff dic
coeff_matrix = []
for model in model_list:
    prot_name = model.replace("_model_MACCS_properties", "")
    svr = load("..//02_drug_model//prot_model//" + model)
    
    coeff_matrix.append(svr.coef_[0])

#np.savetxt("coef_of_model.txt", coeff_matrix, fmt = "%.5e", delimiter=",")

# load model drug number
drug_num_dic = {}
with open("..//01_raw_data//08_drug_number_per_prot.txt", "r") as file2:
    f2 = file2.readlines()
    for line in f2:
        data = line.replace("\n", "").split("\t")
        drug_num_dic[data[0]] = int(data[1])


counter = 1
with open("model_to_model_coeff_scc.txt", "w") as o1:
    print("model_1", "model_2", "drug_1#", "drug_2#", "mean_drug#", "SCC_of_coeff", "pvalue", sep="\t", file=o1)
    for i in range(0, len(coeff_matrix)-1):
        print(counter)
        for j in range(i+1, len(coeff_matrix)):
            scc, pval = spearmanr(coeff_matrix[i], coeff_matrix[j])
            num_1 = drug_num_dic[prot_name_lsit[i]]
            num_2 = drug_num_dic[prot_name_lsit[j]]
            print(prot_name_lsit[i], prot_name_lsit[j], num_1, num_2, str((num_1+num_2)/2), scc, pval, sep="\t", file=o1)
        counter += 1