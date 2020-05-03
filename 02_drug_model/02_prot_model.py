import sys
from sklearn.svm import SVR
from scipy.stats import spearmanr
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
from sklearn.externals import joblib
import numpy as np
import pandas as pd
import math
import os
np.set_printoptions(threshold=sys.maxsize)

# get protin list
file_list = os.listdir(".//drug_per_prot")


# load drug info data (for standardization)
all_df = pd.read_csv("..//01_raw_data//05_all_drug_info.txt", sep="\t", header=None, lineterminator="\n")
all_properties = np.array(list(all_df[all_df.columns[-1]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
scaler = preprocessing.StandardScaler().fit(all_properties)

counter = 1
# # loop
# for protein_file in file_list:
#     print(counter)
#     # loading protein_drug data
#     df = pd.read_csv(".//drug_per_prot//" + protein_file, sep="\t", header=None, lineterminator="\n")
#     if len(df) > 10:  
#         MACCS = np.array(list(df[df.columns[-2]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
#         properties = np.array(list(df[df.columns[-1]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
#         Ki = np.array(list(df[df.columns[3]]))
    
#         # data standardization
#         properties = scaler.transform(properties)
#         MACCS_properties = np.array(np.append(MACCS, properties, axis=1))

#         # train test
#         x_train = MACCS_properties
#         y_train = Ki
#         x_test = MACCS_properties
#         y_test = Ki

#         # setting model parameters
#         svr_linear = SVR(kernel="linear")

#         # training
#         svr_linear.fit(x_train, y_train)

#         # save model
#         protein_name = protein_file.split(".")[0]
#         joblib.dump(svr_linear, ".//prot_model//" + protein_name + "_model_MACCS_properties")
#     counter += 1




for protein_file in file_list:
    print(counter)
    # loading protein_drug data
    df = pd.read_csv(".//drug_per_prot//" + protein_file, sep="\t", header=None, lineterminator="\n")

    MACCS = np.array(list(df[df.columns[-2]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
    properties = np.array(list(df[df.columns[-1]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
    Ki = np.array(list(df[df.columns[3]]))
    
    # data standardization
    properties = scaler.transform(properties)
    MACCS_properties = np.array(np.append(MACCS, properties, axis=1))

    # train test
    x_train = MACCS_properties
    y_train = Ki
    x_test = MACCS_properties
    y_test = Ki

    # setting model parameters
    svr_linear = SVR(kernel="linear")

    # training
    svr_linear.fit(x_train, y_train)

    # save model
    protein_name = protein_file.split(".")[0]
    joblib.dump(svr_linear, ".//prot_all_model//" + protein_name + "_model_MACCS_properties")
    counter += 1
