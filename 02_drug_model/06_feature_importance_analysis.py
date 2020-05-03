from sklearn.svm import SVR
from joblib import dump, load
from scipy.stats import spearmanr, zscore
import numpy as np
import os

prot = "DRD3_HUMAN"
properties = ["acd_logd", "acd_logp", "alogp", "aromatic_rings", "full_mwt", "hba", "hba_lipinski", "hbd", "hbd_lipinski", "heavy_atoms", "mw_freebase", "mw_monoisotopic", "num_lipinski_ro5_violations", "num_ro5_violations", "psa", "qed_weighted", "rtb"]

# load prot name
with open(".//02_protein_order.txt", "r") as file1:
    prot_name_lsit = [line.replace("\n", "") for line in file1]

# make coeff dic
model = load(".//prot_model//" + prot + "_model_MACCS_properties")

with open("properties_coef_" + prot + ".txt", "w") as o1:
    print("weights of property feature", file=o1)
    for i in range(-17, 0):
        print(float(model.coef_[0][i]), properties[i+17], sep="\t", file=o1)
    



