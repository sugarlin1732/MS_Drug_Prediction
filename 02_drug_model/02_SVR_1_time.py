import sys
from sklearn.svm import SVR
from scipy.stats import spearmanr
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
from sklearn.externals import joblib
import numpy as np
import pandas as pd
import math
np.set_printoptions(threshold=sys.maxsize)

# set model protein
prot = "AA2AR_HUMAN"

# load drug info data (for standardization)
all_df = pd.read_csv("..//..//01_raw_data//05_all_drug_info.txt", sep="\t", header=None, lineterminator="\n")
all_properties = np.array(list(all_df[all_df.columns[-1]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
scaler = preprocessing.StandardScaler().fit(all_properties)



# loading protein_drug data
df = pd.read_csv("..//drug_per_prot//" + prot + ".txt", sep="\t", header=None, lineterminator="\n")
MACCS = np.array(list(df[df.columns[-2]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
properties = np.array(list(df[df.columns[-1]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
Ki = np.array(list(df[df.columns[3]]))


# data standardization
properties = scaler.transform(properties)
MACCS_properties = np.array(np.append(MACCS, properties, axis=1))


######################    MACCS    ############################
# self train or data split

x_train, x_test, y_train, y_test = train_test_split(
    MACCS_properties, Ki, test_size=.3)

"""
x_train = MACCS_properties
y_train = Ki
x_test = MACCS_properties
y_test = Ki
"""

# setting model parameters
svr_linear = SVR(kernel="linear")

# training
svr_linear.fit(x_train, y_train)

train_result = svr_linear.predict(x_train)
train_result = list(train_result)

test_result = svr_linear.predict(x_test)
y_pred = list(test_result)

def fisher_z(scc, data_number):
    arct = np.arctanh(scc)
    return arct*(math.pow(((data_number-3)/1.06), 0.5))

# SCC
training_scc = spearmanr(y_train, train_result)
testing_scc = spearmanr(y_test, y_pred)

# Fisher Zscore
training_arctanh = np.arctanh(training_scc)
training_z = fisher_z(training_scc[0], len(train_result))

testing_arctanh = np.arctanh(testing_scc)
testing_z = fisher_z(testing_scc[0], len(y_pred))

# write file
with open(prot + "_linearSVR_0.7_0.3_1time.txt", "w") as o1:
    print("==================================70% training set==================================", file=o1)
    print(training_scc, file=o1)
    print("Fisher zscore = ", training_z, file=o1)
    print("true_value", "pred_value", sep="\t", file=o1)
    for s in range(len(train_result)):
        print(y_train[s], train_result[s], sep="\t", file=o1)

    print(file=o1)

    print("==================================30% testing set==================================", file=o1)
    print(testing_scc, file=o1)
    print("Fisher zscore = ", testing_z, file=o1)
    print("true_value", "pred_value", sep="\t", file=o1)
    for s in range(len(y_pred)):
        print(y_test[s], y_pred[s], sep="\t", file=o1)





"""
print("MACCS, MACCS_properties, train/test_7/3", file=o1)
print("R_square =", r2_score(y_test, y_pred), file=o1)
print(spearmanr(y_test, y_pred), file=o1)
for s in range(len(y_test)):
    print(y_test[s], y_pred[s], file=o1)

"""


# save model
#joblib.dump(svr_rbf, "A2A_newest_model_MACCS_properties")

#o1.close
# o2.close
