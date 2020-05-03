import sys
import os
from sklearn.svm import SVR
from sklearn.model_selection import train_test_split
from scipy.stats import spearmanr
from sklearn import preprocessing
import random
import numpy as np
import pandas as pd
import math
import statistics
np.set_printoptions(threshold=sys.maxsize)


# load drug info data (for standardization)
all_df = pd.read_csv("..//01_raw_data//05_all_drug_info.txt", sep="\t", header=None, lineterminator="\n")
all_properties = np.array(list(all_df[all_df.columns[-1]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
scaler = preprocessing.StandardScaler().fit(all_properties)

model_list = os.listdir(".//prot_all_model//")

counter = 1
for model in model_list:
    print("model " + str(counter))
    print(model)
    prot = model.replace("_model_MACCS_properties", "")
    
    # loading data
    df = pd.read_csv("drug_per_prot//" + prot + ".txt", sep="\t", header=None, lineterminator="\n")
    drug_num = len(df.index)
    if not drug_num > 13:
        continue

    # init result list
    train_rho_result = list()
    train_pvalue = list()
    rho_result = list()
    train_zscore = list()
    test_zscore = list()
    pvalue = list()

    for i in range(0, 1000):

        MACCS = np.array(list(df[df.columns[-2]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))
        properties = np.array(list(df[df.columns[-1]].apply(lambda x: x[1:-1].replace("[", "").replace("]", "").split(","))))

        Ki = np.array(list(df[df.columns[3]]))

        # set standardization scale
        properties = scaler.transform(properties)
        feature = np.array(np.append(MACCS, properties, axis=1))

        # shuffle data
        temp = list(zip(feature, Ki))
        random.shuffle(temp)
        feature, Ki = zip(*temp)
        feature = np.asarray(feature)
        Ki = np.asarray(Ki)

        # data split
        x_train, x_test, y_train, y_test = train_test_split(feature, Ki, test_size=.3)

        # setting model parameters
        svr_linear = SVR(kernel="linear")

        # training
        svr_linear.fit(x_train, y_train)

        # data number must > 3
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


    with open(".//all_model_scc_zscore//" + prot.replace("_HUMAN", "") + "_SVR_0.7_0.3_1000loop.txt", "w") as o1:
        print("drug_number = ", drug_num, file=o1)
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

    counter += 1

print("END")