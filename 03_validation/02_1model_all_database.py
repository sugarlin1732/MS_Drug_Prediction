import sys
from sklearn.svm import SVR
from scipy.stats import spearmanr
from sklearn import preprocessing
from sklearn.externals import joblib
import numpy as np
import pandas as pd
import math
np.set_printoptions(threshold=sys.maxsize)

model_prot = "AL5AP_HUMAN"


# load drug info data (for standardization)
all_df = pd.read_csv("..//01_raw_data//05_all_drug_info.txt", sep="\t", header=None, lineterminator="\n")
all_properties = np.array(list(all_df[all_df.columns[-1]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
scaler = preprocessing.StandardScaler().fit(all_properties)


# loading protein_drug data
df = pd.read_csv("..//01_raw_data//05_all_drug_info.txt", sep="\t", header=None, lineterminator="\n")
MACCS = np.array(list(df[df.columns[-2]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
properties = np.array(list(df[df.columns[-1]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))


# data standardization
properties = scaler.transform(properties)
MACCS_properties = np.array(np.append(MACCS, properties, axis=1))


# load model
SVR_linear = joblib.load("..//02_drug_model//prot_model//" + model_prot + "_model_MACCS_properties")

x_test = MACCS_properties

# predict
result = SVR_linear.predict(x_test)
y_pred = list(result)


with open(model_prot + "_model_all_database_drug.txt", "w") as o1:
    print("predicted_value", file=o1)
    for s in range(len(y_pred)):
        print(y_pred[s], file=o1)