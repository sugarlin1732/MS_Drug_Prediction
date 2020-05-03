# use a model to predict its own drugs and all drugs (#86475)
# compare the distribution of predicted -logKi
# check difference by cohen's d
import numpy as np
from scipy import stats
import statistics
import os
from sklearn import preprocessing
import pandas as pd
from joblib import dump, load

# cohen's d
def cohen_d(x,y):
    nx = len(x)
    ny = len(y)
    dof = nx + ny - 2
    return (np.mean(x) - np.mean(y)) / np.sqrt(((nx-1)*np.std(x, ddof=1) ** 2 + (ny-1)*np.std(y, ddof=1) ** 2) / dof)


# check single model======================================================================
"""
protein ="AL5AP"
with open(protein + "_HUMAN_model_" + protein + "_HUMAN_drug.txt", "r") as file1:
    own_drug = [float(i.replace("\n", "").split("\t")[-1]) for i in file1.readlines()[2:]]


with open(protein + "_HUMAN_model_all_database_drug.txt", "r") as file2:
    all_drug = [float(i.replace("\n", "")) for i in file2.readlines()[1:]]


with open(protein + "cohens_d_of_model_drug_N_all_drug.txt", "w") as o1:
    print("model", "drug_number", "all_drug_number", "cohen's_d", sep="\t", file=o1)
    print(protein, len(own_drug), len(all_drug), cohen_d(own_drug, all_drug), sep="\t", file=o1)

"""

# check all models======================================================================
model_list = os.listdir("..//02_drug_model//prot_all_model//")

# load drug info data (for standardization)
all_df = pd.read_csv("..//01_raw_data//05_all_drug_info.txt", sep="\t", header=None, lineterminator="\n")
all_properties = np.array(list(all_df[all_df.columns[-1]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
scaler = preprocessing.StandardScaler().fit(all_properties)
all_MACCS = np.array(list(all_df[all_df.columns[-2]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
all_properties = scaler.transform(all_properties)
all_feature = np.array(np.append(all_MACCS, all_properties, axis=1))


# write output
counter = 1
with open("cohens_d_of_model_drug_N_all_drug.txt", "w") as o1:
    print("model", "drug_number", "all_drug_number", "cohen's_d", sep="\t", file=o1)

    for prot in model_list:
        print(counter)
        protein = prot.replace("_HUMAN_model_MACCS_properties", "")
        svr = load("..//02_drug_model//prot_all_model//" + prot)
        df = pd.read_csv("..//02_drug_model//drug_per_prot//" + protein + "_HUMAN.txt", sep="\t", header=None, lineterminator="\n")
        MACCS = np.array(list(df[df.columns[-2]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
        properties = np.array(list(df[df.columns[-1]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
        # set standardization scale
        properties = scaler.transform(properties)
        feature = np.array(np.append(MACCS, properties, axis=1))
        own_drug = svr.predict(feature)
        all_drug = svr.predict(all_feature)

        print(protein, len(own_drug), len(all_drug), cohen_d(own_drug, all_drug), sep="\t", file=o1)
        counter += 1