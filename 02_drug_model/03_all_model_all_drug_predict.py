import numpy as np
import pandas as pd
import os
from sklearn.svm import SVR
from sklearn import preprocessing
from sklearn.externals import joblib
import math
import statistics
import sys
np.set_printoptions(threshold=sys.maxsize)

# get model list
model_list = os.listdir(".//prot_model")


# load drug info data
all_df = pd.read_csv("..//01_raw_data//05_all_drug_info.txt", sep="\t", header=None, lineterminator="\n")
all_properties = np.array(list(all_df[all_df.columns[-1]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
MACCS = np.array(list(all_df[all_df.columns[-2]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
properties = np.array(list(all_df[all_df.columns[-1]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))

# data standardization
scaler = preprocessing.StandardScaler().fit(all_properties)

properties = scaler.transform(properties)
MACCS_properties = np.array(np.append(MACCS, properties, axis=1))

# set testing data
x_test = MACCS_properties


# generate predicted -logKi matrix

# run each model
counter = 1
with open("protein_order.txt", "w") as o1:
    with open("03_predicted_ki_matrix.txt", "w") as o2:
        print(counter)
        for model in model_list:
            print(model.split("_model")[0], file=o1)

            # load model
            SVR_linear = joblib.load(".//prot_model//" + model)

            result = SVR_linear.predict(x_test)
            for i in range(len(result)):
                if result[i] == result[-1]:
                    print(result[i], file=o2)
                else:
                    print(result[i], end=",", file=o2)

            counter += 1
print("END")
