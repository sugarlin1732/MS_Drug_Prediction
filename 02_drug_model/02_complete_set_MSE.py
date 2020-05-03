import sys
from sklearn.svm import SVR
from sklearn.model_selection import train_test_split
from joblib import dump, load
from scipy.stats import spearmanr
from sklearn import preprocessing
from sklearn.externals import joblib
from sklearn.metrics import mean_squared_error
import random
import numpy as np
import pandas as pd
import math
import statistics
np.set_printoptions(threshold=sys.maxsize)

# set model protein
prot = "AA2AR_HUMAN"

# load drug info data (for standardization)
all_df = pd.read_csv("..//..//01_raw_data//05_all_drug_info.txt", sep="\t", header=None, lineterminator="\n")
all_properties = np.array(list(all_df[all_df.columns[-1]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
scaler = preprocessing.StandardScaler().fit(all_properties)


# loading data
df = pd.read_csv("..//drug_per_prot//" + prot + ".txt", sep="\t", header=None, lineterminator="\n")
MACCS = np.array(list(df[df.columns[-2]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
properties = np.array(list(df[df.columns[-1]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
Ki = np.array(list(df[df.columns[3]]))
# set standardization scale
properties = scaler.transform(properties)
feature = np.array(np.append(MACCS, properties, axis=1))


svr_linear = load("..//prot_all_model//" + prot + "_model_MACCS_properties")

result = svr_linear.predict(feature)
with open("complete_dataset_MSE.txt", "w") as o1:
    print("MSE of complete dataset", file=o1)
    print(mean_squared_error(Ki, result), file=o1)
