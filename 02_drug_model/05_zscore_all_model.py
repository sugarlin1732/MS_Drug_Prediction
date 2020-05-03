import numpy as np
import pandas as pd
import os
from sklearn.svm import SVR
from sklearn import preprocessing
from joblib import dump, load
from scipy.stats import spearmanr
import math
import statistics
import sys
np.set_printoptions(threshold=sys.maxsize)

# get model list
model_list = os.listdir(".//prot_model")

# load drug info data (for standardization)
all_df = pd.read_csv("..//01_raw_data//05_all_drug_info.txt", sep="\t", header=None, lineterminator="\n")
all_properties = np.array(list(all_df[all_df.columns[-1]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
scaler = preprocessing.StandardScaler().fit(all_properties)


with open("zscore_all_model.txt", "w") as o1:
    print("model", "zscore", sep="\t", file=o1)
    counter = 1
    for model in model_list:
        
        print(counter)
        prot = model.replace("_model_MACCS_properties", "")
        
        # loading protein_drug data
        df = pd.read_csv(".//drug_per_prot//" + prot + ".txt", sep="\t", header=None, lineterminator="\n")
        MACCS = np.array(list(df[df.columns[-2]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
        properties = np.array(list(df[df.columns[-1]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
        Ki = np.array(list(df[df.columns[3]]))

        # data standardization
        properties = scaler.transform(properties)
        MACCS_properties = np.array(np.append(MACCS, properties, axis=1))
        x_train = MACCS_properties
        y_train = Ki
        x_test = MACCS_properties
        y_test = Ki

        # setting model parameters
        svr_linear = SVR(kernel="linear")

        # load model
        SVR_linear = load(".//prot_model//" + model)
        y_pred = SVR_linear.predict(x_test)

        def fisher_z(scc, data_number):
            arct = np.arctanh(scc)
            return arct*(math.pow(((data_number-3)/1.06), 0.5))

        # SCC
        testing_scc = spearmanr(y_test, y_pred)

        # Fisher Zscore
        testing_arctanh = np.arctanh(testing_scc)
        testing_z = fisher_z(testing_scc[0], len(y_pred))
        print(prot, testing_z, sep="\t", file=o1)

        counter += 1