import sys
from sklearn.svm import SVR
from sklearn.model_selection import train_test_split
from scipy.stats import spearmanr
from sklearn import preprocessing
from sklearn.externals import joblib
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

train_rho_result = list()
train_pvalue = list()

rho_result = list()
train_zscore = list()
test_zscore = list()
pvalue = list()

for i in range(0, 1000):
    print(i)
    # data split
    portion = np.random.rand(len(df)) < 0.7
    train = df[portion]
    test = df[~portion]

    MACCS_train = np.array(list(train[df.columns[-2]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
    MACCS_test = np.array(list(test[df.columns[-2]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
    properties_train = np.array(list(train[df.columns[-1]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
    properties_test = np.array(list(test[df.columns[-1]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))

    Ki_train = np.array(list(train[train.columns[3]]))
    Ki_test = np.array(list(test[test.columns[3]]))

    # set standardization scale
    properties_train = scaler.transform(properties_train)
    properties_test = scaler.transform(properties_test)

    #fingerprint + properties
    x_train = np.array(np.append(MACCS_train, properties_train, axis=1))
    x_test = np.array(np.append(MACCS_test, properties_test, axis=1))

    y_train = Ki_train
    y_test = Ki_test

    # setting model parameters
    svr_linear = SVR(kernel="linear")

    # training
    svr_linear.fit(x_train, y_train)


    def fisher_z(scc, data_number):
        arct = np.arctanh(scc)
        return arct*(math.pow(((data_number-3)/1.06), 0.5))

    #testing result
    result = svr_linear.predict(x_test)
    y_pred = list(result)

    #training result
    training_result = svr_linear.predict(x_train)
    train_predict = list(training_result)

    # record training result
    train_rho, train_pval = spearmanr(y_train, train_predict)
    train_z = fisher_z(float(train_rho), len(train_predict))

    train_rho_result.append(float(train_rho))
    train_pvalue.append(float(train_pval))
    train_zscore.append(train_z)

    # record testing result
    rho, pval = spearmanr(y_test, y_pred)
    test_z = fisher_z(float(rho), len(y_pred))

    test_zscore.append(test_z)
    rho_result.append(float(rho))
    pvalue.append(float(pval))




with open(prot + "_linearSVR_0.7_0.3_1000loop.txt", "w") as o1:

    print("min/max training spearman rho results", min(train_rho_result), max(train_rho_result), file=o1)
    print("min/max training spearman p value", min(train_pvalue), max(train_pvalue), file=o1)
    print("average training spearman results", statistics.mean(train_rho_result), file=o1)
    print("mdeian training spearman results", statistics.median(train_rho_result), file=o1)

    print("min/max training zscore", min(train_zscore), max(train_zscore), file=o1)
    print("average training zscore", statistics.mean(train_zscore), file=o1)
    print("mdeian training zscore", statistics.median(train_zscore), file=o1)

    print("SCC", "zscore", sep="\t", file=o1)
    for i in range(len(train_rho_result)):
        print(train_rho_result[i], train_zscore[i], file=o1)

    print("min/max testing spearman rho results", min(rho_result), max(rho_result), file=o1)
    print("min/max testing spearman p value", min(pvalue), max(pvalue), file=o1)
    print("average testing spearman results", statistics.mean(rho_result), file=o1)
    print("mdeian testing spearman results", statistics.median(rho_result), file=o1)

    print("min/max testing zscore", min(test_zscore), max(test_zscore), file=o1)
    print("average testing zscore", statistics.mean(test_zscore), file=o1)
    print("mdeian testing zscore", statistics.median(test_zscore), file=o1)

    print("SCC", "zscore", sep="\t", file=o1)
    for j in range(len(rho_result)):
        print(rho_result[j], test_zscore[j], file=o1)
